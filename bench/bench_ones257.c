/* Standalone Cortex-M prefix-product benchmark. Its reference is the same
 * private function used by expkey_from_sk. Including that translation unit
 * here gives this test access without changing the scheme's interface. */
#include <stdio.h>
#include <string.h>
#include "end_params.h"
static void __attribute__((noinline)) ones_product_mod257(int16_t out[END_N], const int16_t in[END_N]);
#include "../src/common/end_kem.c"
#include "cycles.h"

void END_NS(bench_ones257_dsp)(int16_t out[END_N], const int16_t in[END_N]);
typedef void (*ones_fn)(int16_t *, const int16_t *);

static int16_t input[END_N] END_ALIGN, reference[END_N] END_ALIGN, output[END_N] END_ALIGN;

static void fill_input(unsigned pattern)
{
    uint32_t state = 1u + pattern;
    for (unsigned i = 0; i < END_N; i++) {
        state = state * 1664525u + 1013904223u;
        switch (pattern) {
        case 0: input[i] = 0; break;
        case 1: input[i] = 256; break;
        case 2: input[i] = (i & 1u) ? 256 : 0; break;
        case 3: input[i] = (i % 4 == 3 || i % 4 == 0) ? 256 : 1; break;
        default: input[i] = (int16_t)(state % 257u); break;
        }
    }
}

#ifndef END_NO_DWT
static void measure(const char *label, ones_fn function)
{
    /* Each call remains observable even when the input is unchanged. */
    ones_fn volatile run = function;
    uint64_t start = cycles_now();
    for (unsigned i = 0; i < 64; i++) run(output, input);
    uint64_t elapsed = cycles_now() - start;
    printf("%s: %lu %s per call\n", label, (unsigned long)(elapsed / 64), cycles_unit());
}
#endif

int main(int argc, char **argv)
{
    if (argc == 2) {
        /* One operation for the model plugin, with identical input/setup.
         * 'none' is the subtraction baseline used by m4model_kernels.sh. */
        fill_input(4);
        if (!strcmp(argv[1], "ref")) ones_product_mod257(output, input);
        else if (!strcmp(argv[1], "dsp")) END_NS(bench_ones257_dsp)(output, input);
        else if (strcmp(argv[1], "none")) return 2;
        __asm__ volatile("" : : "m"(output) : "memory");
        return 0;
    }
    if (argc != 1) return 2;
    for (unsigned pattern = 0; pattern < 20; pattern++) {
        fill_input(pattern);
        ones_product_mod257(reference, input);
        END_NS(bench_ones257_dsp)(output, input);
        for (unsigned i = 0; i < END_N; i++) {
            if (reference[i] != output[i]) {
                printf("prefix mismatch at pattern %u, coefficient %u\n", pattern, i);
                return 1;
            }
        }
    }
    printf("END-%d prefix product: all coefficients match\n", END_N);
#ifndef END_NO_DWT
    measure("reference", ones_product_mod257);
    measure("DSP", END_NS(bench_ones257_dsp));
#endif
    return 0;
}
