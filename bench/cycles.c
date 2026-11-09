/*
 * Cycle / time counters for the benchmarks.
 *   x86-64            : rdtsc
 *   AArch64 (Linux)   : PMCCNTR_EL0 if user access was enabled (-DEND_AARCH64_PMU), else CNTVCT_EL0 (generic timer)
 *   ARMv7E-M / v8-M   : DWT->CYCCNT (0 under QEMU, which does not model cycles)
 *   other             : clock()
 */
#include <stdint.h>
#include <time.h>
#include "cycles.h"

#if defined(__x86_64__) || defined(_M_X64)
uint64_t cycles_now(void)
{
    uint32_t lo, hi;
    __asm__ __volatile__("rdtsc" : "=a"(lo), "=d"(hi));
    return ((uint64_t)hi << 32) | lo;
}
const char *cycles_unit(void) { return "rdtsc cycles"; }
#elif defined(__aarch64__)
uint64_t cycles_now(void)
{
    uint64_t v;
#ifdef END_AARCH64_PMU
    __asm__ __volatile__("mrs %0, pmccntr_el0" : "=r"(v));
#else
    __asm__ __volatile__("isb; mrs %0, cntvct_el0" : "=r"(v));
#endif
    return v;
}
#ifdef END_AARCH64_PMU
const char *cycles_unit(void) { return "PMU cycles"; }
#else
const char *cycles_unit(void) { return "cntvct ticks"; }
#endif
#elif defined(__arm__) && defined(END_NO_DWT)
/* QEMU does not model the DWT: report 0 cycles (instruction counts come from the TCG plugin instead) */
uint64_t cycles_now(void) { return 0; }
const char *cycles_unit(void) { return "cycles (DWT disabled)"; }
#elif defined(__arm__) && (defined(__ARM_ARCH_7EM__) || defined(__ARM_ARCH_8M_MAIN__))
#define DWT_CTRL   (*(volatile uint32_t *)0xE0001000)
#define DWT_CYCCNT (*(volatile uint32_t *)0xE0001004)
#define DEMCR      (*(volatile uint32_t *)0xE000EDFC)
static int dwt_init = 0;
uint64_t cycles_now(void)
{
    if (!dwt_init) { DEMCR |= 1u << 24; DWT_CYCCNT = 0; DWT_CTRL |= 1u; dwt_init = 1; }
    return DWT_CYCCNT;
}
const char *cycles_unit(void) { return "DWT cycles"; }
#else
uint64_t cycles_now(void) { return (uint64_t)clock(); }
const char *cycles_unit(void) { return "clock() ticks"; }
#endif

#if !defined(__arm__) || defined(__linux__)
uint64_t wall_ns(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (uint64_t)ts.tv_sec * 1000000000ull + (uint64_t)ts.tv_nsec;
}
#else
uint64_t wall_ns(void) { return 0; }
#endif
