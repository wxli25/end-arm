/*
 * m4cycles -- QEMU TCG plugin: Cortex-M4 cycle estimate for a Thumb-2 guest (qemu-system-arm).
 *
 * Every executed instruction gets a base cost from the Cortex-M4 TRM (r0p1, Table 3-1 and
 * section 3.3 "Load/store timings"); taken branches add a pipeline refill detected at run time
 * (the next translation block does not start where the previous one ended).
 *
 *   single LDR/LDRB/LDRH/LDRS*   2, or 1 when it directly follows a single load whose
 *                                 destination it does not use for its address (pipelined)
 *   STR with immediate offset    1;  STR with register offset  2 (1 after a single load)
 *   LDRD/STRD                    3   LDM/STM/PUSH/POP   1+N   VLDM/VSTM/VPUSH/VPOP  1+N
 *   VLDR/VSTR                    2   VMOV two core regs 2
 *   MLA/MLS                      mla (2)    SDIV/UDIV  div (7; the TRM gives 2..12)
 *   TBB/TBH                      3   LDR pc  2          IT  it (1)
 *   everything else              1   (incl. MUL, (U|S)MULL, (U|S)MLAL, SMLA*, SMUAD, ...)
 *   taken direct branch          +pdir (2)   taken indirect branch (bx, pop/ldm/ldr pc, tbb)  +pind (3)
 *
 * The model assumes zero wait-state memory (pqm4 clocks the NUCLEO-L4R5ZI so that flash has no
 * wait states).  It is an estimate: validate it against measured boards (run_m4model.sh, then m4model_report.py --board).
 *
 * Arguments (key=value): pdir, pind, mla, div, it, profile=on|off.
 * Output (qemu -d plugin -D file):  "insns N", "cycles N", class counters, and with profile=on
 * one line per function "sym <name> <insns> <cycles>" (a translation block is attributed to the
 * symbol of its first instruction; the refill after a taken branch to the branching block).
 */
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <qemu-plugin.h>

QEMU_PLUGIN_EXPORT int qemu_plugin_version = QEMU_PLUGIN_VERSION;

static unsigned p_dir = 2, p_ind = 3, c_mla = 2, c_div = 7, c_it = 1;
static int profile;

enum { EXIT_NONE, EXIT_DIR, EXIT_IND };
enum { K_ALU, K_LD, K_ST_IMM, K_ST_REG, K_OTHER_MEM };
enum { S_LD, S_LD_PIPE, S_ST, S_MULTI_WORDS, S_LDRD, S_MLA, S_DIV, S_BR_DIR, S_BR_IND, S_TAKEN_DIR,
       S_TAKEN_IND, S_FP_MEM, S_IT, S_NSTAT };
static const char *const stat_name[S_NSTAT] = {
    "load", "load_pipelined", "store", "ldm_stm_words", "ldrd_strd", "mla_mls", "div", "branch_direct",
    "branch_indirect", "taken_direct", "taken_indirect", "fp_load_store", "it"};

struct tbinfo {
    uint64_t start, end;
    uint32_t cycles, ninsn;
    uint32_t stat[S_NSTAT];
    int exit_kind, sym;
};

/* ---- symbol table (pointer-keyed) -------------------------------------------------------- */
#define SYM_MAX 8192
static const char *sym_name[SYM_MAX];
static uint64_t sym_insn[SYM_MAX], sym_cyc[SYM_MAX];
static int nsym;

static int sym_index(const char *s)
{
    int i;
    if (!s) s = "?";
    for (i = 0; i < nsym; i++)          /* translation time only; the list stays short */
        if (sym_name[i] == s || strcmp(sym_name[i], s) == 0) return i;
    if (nsym == SYM_MAX) return SYM_MAX - 1;
    sym_name[nsym] = s;
    return nsym++;
}

/* ---- run-time state ----------------------------------------------------------------------- */
static uint64_t tot_cyc, tot_insn, tot_stat[S_NSTAT];
static uint64_t last_end;
static int last_exit = EXIT_NONE, last_sym = -1;

static void tb_exec(unsigned int vcpu, void *udata)
{
    const struct tbinfo *t = udata;
    int i;
    (void)vcpu;
    if (last_sym >= 0 && t->start != last_end) {           /* control flow changed: refill */
        unsigned p = last_exit == EXIT_DIR ? p_dir : p_ind;
        tot_cyc += p;
        sym_cyc[last_sym] += p;
        tot_stat[last_exit == EXIT_DIR ? S_TAKEN_DIR : S_TAKEN_IND]++;
    }
    tot_cyc += t->cycles;
    tot_insn += t->ninsn;
    for (i = 0; i < S_NSTAT; i++) tot_stat[i] += t->stat[i];
    sym_cyc[t->sym] += t->cycles;
    sym_insn[t->sym] += t->ninsn;
    last_end = t->end;
    last_exit = t->exit_kind;
    last_sym = t->sym;
}

/* ---- Thumb-2 classification ---------------------------------------------------------------- */
struct dec {
    unsigned cost;
    int kind;            /* K_* : what the load/store pipelining rule needs */
    int rt, rn, rm;      /* -1 when absent */
    int exit_kind;       /* EXIT_* if the instruction can change control flow */
    int stat;            /* S_* or -1 */
    unsigned words;      /* LDM/STM words */
};

static unsigned popcount16(unsigned x) { return (unsigned)__builtin_popcount(x & 0xFFFFu); }

static void decode(const uint8_t *b, size_t len, struct dec *d)
{
    unsigned hw1 = b[0] | (b[1] << 8);
    memset(d, 0, sizeof *d);
    d->cost = 1;
    d->kind = K_ALU;
    d->rt = d->rn = d->rm = -1;
    d->stat = -1;
    if (len == 2) {
        if ((hw1 & 0xF800) == 0x4800) {                              /* LDR literal */
            d->kind = K_LD; d->rt = (hw1 >> 8) & 7; d->rn = 15;
        } else if ((hw1 & 0xF000) == 0x5000) {                       /* register offset */
            unsigned opb = (hw1 >> 9) & 7;
            d->rm = (hw1 >> 6) & 7; d->rn = (hw1 >> 3) & 7; d->rt = hw1 & 7;
            d->kind = opb >= 3 ? K_LD : K_ST_REG;
        } else if ((hw1 & 0xE000) == 0x6000 || (hw1 & 0xF000) == 0x8000) {   /* imm word/byte/half */
            d->rn = (hw1 >> 3) & 7; d->rt = hw1 & 7;
            d->kind = (hw1 & 0x0800) ? K_LD : K_ST_IMM;
        } else if ((hw1 & 0xF000) == 0x9000) {                       /* SP relative */
            d->rn = 13; d->rt = (hw1 >> 8) & 7;
            d->kind = (hw1 & 0x0800) ? K_LD : K_ST_IMM;
        } else if ((hw1 & 0xF000) == 0xC000) {                       /* LDM/STM */
            d->words = popcount16(hw1 & 0xFF);
            d->cost = 1 + d->words; d->kind = K_OTHER_MEM; d->stat = S_MULTI_WORDS;
        } else if ((hw1 & 0xFE00) == 0xB400) {                       /* PUSH */
            d->words = popcount16(hw1 & 0x1FF);
            d->cost = 1 + d->words; d->kind = K_OTHER_MEM; d->stat = S_MULTI_WORDS;
        } else if ((hw1 & 0xFE00) == 0xBC00) {                       /* POP */
            d->words = popcount16(hw1 & 0x1FF);
            d->cost = 1 + d->words; d->kind = K_OTHER_MEM; d->stat = S_MULTI_WORDS;
            if (hw1 & 0x100) d->exit_kind = EXIT_IND;
        } else if ((hw1 & 0xF500) == 0xB100) {                       /* CBZ/CBNZ */
            d->exit_kind = EXIT_DIR;
        } else if ((hw1 & 0xF000) == 0xD000 && ((hw1 >> 8) & 0xF) < 0xE) {   /* B<c> */
            d->exit_kind = EXIT_DIR;
        } else if ((hw1 & 0xF800) == 0xE000) {                       /* B */
            d->exit_kind = EXIT_DIR;
        } else if ((hw1 & 0xFC00) == 0x4400) {                       /* hi-reg ADD/CMP/MOV, BX/BLX */
            unsigned op = (hw1 >> 8) & 3, rd = ((hw1 >> 4) & 8) | (hw1 & 7);
            if (op == 3 || ((op == 0 || op == 2) && rd == 15)) d->exit_kind = EXIT_IND;
        } else if ((hw1 & 0xFF00) == 0xBF00 && (hw1 & 0xF)) {        /* IT */
            d->cost = c_it; d->stat = S_IT;
        }
    } else {
        unsigned hw2 = b[2] | (b[3] << 8);
        if ((hw1 & 0xFE00) == 0xF800) {                              /* single load/store */
            int load = (hw1 >> 4) & 1;
            unsigned size = (hw1 >> 5) & 3;
            d->rn = hw1 & 0xF; d->rt = hw2 >> 12;
            if (!(hw1 & 0x80) && d->rn != 15 && (hw2 & 0x0FC0) == 0) d->rm = hw2 & 0xF;
            if (load && d->rt == 15) {
                if (size == 2) { d->cost = 2; d->exit_kind = EXIT_IND; d->kind = K_OTHER_MEM; }
                /* else PLD/PLI hint: 1 cycle */
            } else if (load) {
                d->kind = K_LD;
            } else {
                d->kind = d->rm >= 0 ? K_ST_REG : K_ST_IMM;
            }
        } else if ((hw1 & 0xFE40) == 0xE800) {                       /* LDM/STM/PUSH/POP (T2) */
            d->words = popcount16(hw2 & 0xDFFF);
            d->cost = 1 + d->words; d->kind = K_OTHER_MEM; d->stat = S_MULTI_WORDS;
            if ((hw1 & 0x10) && (hw2 & 0x8000)) d->exit_kind = EXIT_IND;
        } else if ((hw1 & 0xFE40) == 0xE840) {                       /* dual / exclusive / table branch */
            if ((hw1 & 0xFFF0) == 0xE8D0 && (hw2 & 0xFFE0) == 0xF000) {
                d->cost = 3; d->exit_kind = EXIT_IND; d->kind = K_OTHER_MEM;
            } else if (hw1 & 0x0120) {                                /* LDRD/STRD */
                d->cost = 3; d->kind = K_OTHER_MEM; d->stat = S_LDRD;
            } else {                                                  /* LDREX/STREX family */
                d->rn = hw1 & 0xF; d->rt = hw2 >> 12;
                d->kind = (hw1 & 0x10) ? K_LD : K_ST_REG;
            }
        } else if ((hw1 & 0xFF80) == 0xFB00) {                       /* 32-bit multiply / MAC */
            unsigned op1 = (hw1 >> 4) & 7, ra = hw2 >> 12, op2 = (hw2 >> 4) & 3;
            if (op1 == 0 && ((op2 == 0 && ra != 15) || op2 == 1)) { d->cost = c_mla; d->stat = S_MLA; }
        } else if ((hw1 & 0xFF80) == 0xFB80) {                       /* long multiply / divide */
            unsigned op1 = (hw1 >> 4) & 7;
            if (op1 == 1 || op1 == 3) { d->cost = c_div; d->stat = S_DIV; }
        } else if ((hw1 & 0xF800) == 0xF000 && (hw2 & 0x8000)) {     /* branches, misc control */
            unsigned t = hw2 & 0x5000;
            if (t != 0 || ((hw1 >> 6) & 0xE) != 0xE) d->exit_kind = EXIT_DIR;
        } else if ((hw1 & 0xEE00) == 0xEC00) {                       /* coprocessor load/store */
            if ((hw1 & 0xFFE0) == 0xEC40) {
                d->cost = 2;                                          /* VMOV 2 core <-> 2 S / 1 D */
            } else if ((hw1 & 0x0120) == 0x0100) {
                d->cost = 2; d->kind = K_OTHER_MEM; d->stat = S_FP_MEM;   /* VLDR/VSTR */
            } else {
                d->words = hw2 & 0xFF;
                d->cost = 1 + d->words; d->kind = K_OTHER_MEM; d->stat = S_MULTI_WORDS;
            }
        }
    }
    if (d->exit_kind) d->stat = d->exit_kind == EXIT_DIR ? S_BR_DIR : S_BR_IND;
}

static void tb_trans(qemu_plugin_id_t id, struct qemu_plugin_tb *tb)
{
    size_t n = qemu_plugin_tb_n_insns(tb), i;
    struct tbinfo *t = calloc(1, sizeof *t);
    int prev_kind = K_ALU, prev_rt = -1;
    (void)id;
    t->start = qemu_plugin_tb_vaddr(tb);
    t->ninsn = (uint32_t)n;
    t->sym = sym_index(n ? qemu_plugin_insn_symbol(qemu_plugin_tb_get_insn(tb, 0)) : NULL);
    for (i = 0; i < n; i++) {
        struct qemu_plugin_insn *insn = qemu_plugin_tb_get_insn(tb, i);
        uint8_t buf[4] = {0};
        size_t len = qemu_plugin_insn_data(insn, buf, sizeof buf);
        struct dec d;
        unsigned cost;
        decode(buf, len, &d);
        cost = d.cost;
        if (d.kind == K_LD || d.kind == K_ST_REG) {
            int dep = prev_rt >= 0 && (d.rn == prev_rt || d.rm == prev_rt);
            if (prev_kind == K_LD && !dep) {
                cost = 1;
                if (d.kind == K_LD) t->stat[S_LD_PIPE]++;
            } else {
                cost = 2;
            }
        }
        if (d.kind == K_LD) t->stat[S_LD]++;
        if (d.kind == K_ST_IMM || d.kind == K_ST_REG) t->stat[S_ST]++;
        if (d.stat == S_MULTI_WORDS) t->stat[S_MULTI_WORDS] += d.words;
        else if (d.stat >= 0) t->stat[d.stat]++;
        t->cycles += cost;
        prev_kind = d.kind;
        prev_rt = d.kind == K_LD ? d.rt : -1;
        if (i + 1 == n) {
            t->end = qemu_plugin_insn_vaddr(insn) + qemu_plugin_insn_size(insn);
            t->exit_kind = d.exit_kind;
        }
    }
    qemu_plugin_register_vcpu_tb_exec_cb(tb, tb_exec, QEMU_PLUGIN_CB_NO_REGS, t);
}

static void at_exit(qemu_plugin_id_t id, void *p)
{
    char line[512];
    int i;
    (void)id; (void)p;
    snprintf(line, sizeof line, "insns %" PRIu64 "\ncycles %" PRIu64 "\n", tot_insn, tot_cyc);
    qemu_plugin_outs(line);
    for (i = 0; i < S_NSTAT; i++) {
        snprintf(line, sizeof line, "stat %s %" PRIu64 "\n", stat_name[i], tot_stat[i]);
        qemu_plugin_outs(line);
    }
    if (profile)
        for (i = 0; i < nsym; i++) {
            if (!sym_insn[i] && !sym_cyc[i]) continue;
            snprintf(line, sizeof line, "sym %s %" PRIu64 " %" PRIu64 "\n", sym_name[i], sym_insn[i], sym_cyc[i]);
            qemu_plugin_outs(line);
        }
}

QEMU_PLUGIN_EXPORT int qemu_plugin_install(qemu_plugin_id_t id, const qemu_info_t *info, int argc, char **argv)
{
    int i;
    (void)info;
    for (i = 0; i < argc; i++) {
        char *eq = strchr(argv[i], '=');
        if (!eq) { fprintf(stderr, "m4cycles: bad argument %s\n", argv[i]); return -1; }
        *eq = 0;
        const char *k = argv[i], *v = eq + 1;
        if (!strcmp(k, "pdir")) p_dir = (unsigned)strtoul(v, NULL, 0);
        else if (!strcmp(k, "pind")) p_ind = (unsigned)strtoul(v, NULL, 0);
        else if (!strcmp(k, "mla")) c_mla = (unsigned)strtoul(v, NULL, 0);
        else if (!strcmp(k, "div")) c_div = (unsigned)strtoul(v, NULL, 0);
        else if (!strcmp(k, "it")) c_it = (unsigned)strtoul(v, NULL, 0);
        else if (!strcmp(k, "profile")) {
            bool on;
            if (!qemu_plugin_bool_parse(k, v, &on)) { fprintf(stderr, "m4cycles: bad profile=%s\n", v); return -1; }
            profile = on;
        } else { fprintf(stderr, "m4cycles: unknown option %s\n", k); return -1; }
    }
    qemu_plugin_register_vcpu_tb_trans_cb(id, tb_trans);
    qemu_plugin_register_atexit_cb(id, at_exit, NULL);
    return 0;
}
