/*
 * Minimal Cortex-M startup: vector table, .data/.bss initialisation, FPU enable, semihosting init,
 * stack painting (for stack-usage measurement), main(), exit via semihosting (terminates QEMU).
 */
#include <stdint.h>
#include <string.h>
#include <stdlib.h>

extern uint32_t _sidata, _sdata, _edata, _sbss, _ebss, _estack, _sstack, _end, _heap_limit;
extern int main(int argc, char **argv);
extern void initialise_monitor_handles(void);

/* semihosting SYS_GET_CMDLINE (0x15): fetch the command line passed by QEMU "-semihosting-config ...,arg=..." */
static char cmdline[512];
static char *argv_buf[16];
static int get_args(void)
{
    struct { char *buf; uint32_t len; } blk = { cmdline, sizeof cmdline - 1 };
    register uint32_t r0 __asm__("r0") = 0x15;
    register void *r1 __asm__("r1") = &blk;
    int argc = 0;
    char *p;
    __asm__ volatile("bkpt 0xAB" : "+r"(r0) : "r"(r1) : "memory");
    if (r0 != 0) return 0;
    cmdline[blk.len < sizeof cmdline ? blk.len : sizeof cmdline - 1] = 0;
    for (p = cmdline; *p && argc < 15; ) {
        while (*p == ' ') p++;
        if (!*p) break;
        argv_buf[argc++] = p;
        while (*p && *p != ' ') p++;
        if (*p) *p++ = 0;
    }
    argv_buf[argc] = 0;
    return argc;
}

#define STACK_PAINT 0xDEADBEEFu

void Reset_Handler(void);
void Default_Handler(void) { for (;;) { } }

__attribute__((section(".vectors"), used))
const void *vector_table[16] = {
    (void *)&_estack, (void *)Reset_Handler, Default_Handler, Default_Handler, Default_Handler, Default_Handler,
    Default_Handler, 0, 0, 0, 0, Default_Handler, Default_Handler, 0, Default_Handler, Default_Handler
};

/* stack painting: fill the stack region below the current SP with a pattern; stack_used() counts overwritten words */
static void paint_stack(void)
{
    uint32_t *p = &_sstack;
    uint32_t sp;
    __asm__ volatile("mov %0, sp" : "=r"(sp));
    while ((uint32_t)p < sp - 64) *p++ = STACK_PAINT;
}
uint32_t platform_stack_used(void)
{
    uint32_t *p = &_sstack;
    while ((uint32_t)p < (uint32_t)&_estack && *p == STACK_PAINT) p++;
    return (uint32_t)&_estack - (uint32_t)p;
}

void Reset_Handler(void)
{
    uint32_t *src = &_sidata, *dst = &_sdata;
    int r;
    while (dst < &_edata) *dst++ = *src++;
    for (dst = &_sbss; dst < &_ebss; dst++) *dst = 0;
#if defined(__ARM_FP)
    /* enable CP10/CP11 (FPU) */
    *(volatile uint32_t *)0xE000ED88 |= (0xFu << 20);
    __asm__ volatile("dsb; isb");
#endif
    paint_stack();
    initialise_monitor_handles();
    r = main(get_args(), argv_buf);
    exit(r);
}

/* newlib hooks that rdimon does not provide for a custom startup */
void *_sbrk(int incr)
{
    static char *heap = 0;
    char *prev;
    if (!heap) heap = (char *)&_end;
    if (heap + incr > (char *)&_heap_limit) return (void *)-1;
    prev = heap; heap += incr;
    return prev;
}
