# END KEM for ARM -- top-level Makefile (host builds, AArch64 Linux cross builds)
#
#   make test_kernels                       # portable C backend, END-512, host compiler
#   make END_N=1024 test_all                # END-1024
#   make BACKEND=aarch64 CC=aarch64-linux-gnu-gcc CFLAGS="-O3 -march=armv8-a+simd" test_all
#   make BACKEND=ref ARCH=aarch64 ...
# Cortex-M bare-metal builds live in platform/cortex-m/Makefile.

END_N    ?= 512
BACKEND  ?= ref
ARCH     ?= host
CC       ?= gcc
CFLAGS   ?= -O2
WARN     ?= -Wall -Wextra -Wno-unused-parameter -Wno-unused-function
STD      ?= -std=gnu99
EXTRA_CPPFLAGS ?=
CPPFLAGS  = -DEND_N=$(END_N) -DEND_BACKEND_NAME=\"$(BACKEND)\" -Isrc/common -Isrc/common/ntrugen -Isrc/$(BACKEND) -Isrc/ref $(EXTRA_CPPFLAGS)
LDFLAGS  ?=
LDLIBS   ?=
RUN      ?=

BUILD    := build/$(ARCH)-$(BACKEND)-$(END_N)

KERNELS      := ntt257 ntt3329 ntt12289 polyvec baseinv257
KERNEL_SRCS  := $(foreach k,$(KERNELS),$(if $(wildcard src/$(BACKEND)/$(k).c),src/$(BACKEND)/$(k).c,src/ref/$(k).c))
KERNEL_ASM   := $(wildcard src/$(BACKEND)/*.S)
COMMON_SRCS  := src/common/end_codec.c src/common/end_decode.c src/common/end_kem.c src/common/end_keygen.c \
                src/common/fips202.c src/common/randombytes.c
NTRUGEN_SRCS := src/common/ntrugen/ntru.c src/common/ntrugen/module_p.c src/common/ntrugen/big_integer.c \
                src/common/ntrugen/fixed_point.c src/common/ntrugen/poly_big.c
LIB_SRCS     := $(KERNEL_SRCS) $(COMMON_SRCS) $(NTRUGEN_SRCS)
LIB_OBJS     := $(patsubst %.c,$(BUILD)/%.o,$(LIB_SRCS)) $(patsubst %.S,$(BUILD)/%.o,$(KERNEL_ASM))
KERNEL_OBJS  := $(patsubst %.c,$(BUILD)/%.o,$(KERNEL_SRCS)) $(patsubst %.S,$(BUILD)/%.o,$(KERNEL_ASM))

DEPFLAGS = -MMD -MP

$(BUILD)/%.o: %.c
	@mkdir -p $(dir $@)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(WARN) $(STD) $(DEPFLAGS) -c $< -o $@

$(BUILD)/%.o: %.S
	@mkdir -p $(dir $@)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(BUILD)/libend.a: $(LIB_OBJS)
	$(AR) rcs $@ $^

lib: $(BUILD)/libend.a

$(BUILD)/test_kernels: test/test_kernels.c $(KERNEL_OBJS)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(WARN) $(STD) $^ -o $@ $(LDFLAGS) $(LDLIBS)

$(BUILD)/test_kat: test/test_kat.c $(BUILD)/libend.a
	$(CC) $(CPPFLAGS) $(CFLAGS) $(WARN) $(STD) $< -o $@ $(BUILD)/libend.a $(LDFLAGS) $(LDLIBS)

$(BUILD)/test_roundtrip: test/test_roundtrip.c $(BUILD)/libend.a
	$(CC) $(CPPFLAGS) $(CFLAGS) $(WARN) $(STD) $< -o $@ $(BUILD)/libend.a $(LDFLAGS) $(LDLIBS)

$(BUILD)/test_entropy: test/test_entropy.c platform/cortex-m/randombytes_m.c $(BUILD)/libend.a
	$(CC) $(CPPFLAGS) $(CFLAGS) $(WARN) $(STD) $^ -o $@ $(LDFLAGS) $(LDLIBS)

$(BUILD)/test_wraps: test/test_wraps.c $(BUILD)/libend.a
	$(CC) $(CPPFLAGS) $(CFLAGS) $(WARN) $(STD) $< -o $@ $(BUILD)/libend.a $(LDFLAGS) $(LDLIBS)

$(BUILD)/bench: bench/bench.c bench/cycles.c $(BUILD)/libend.a
	$(CC) $(CPPFLAGS) $(CFLAGS) $(WARN) $(STD) bench/bench.c bench/cycles.c -o $@ $(BUILD)/libend.a $(LDFLAGS) $(LDLIBS)

test_kernels: $(BUILD)/test_kernels
	$(RUN) $(BUILD)/test_kernels

test_kat: $(BUILD)/test_kat
	$(RUN) $(BUILD)/test_kat

test_roundtrip: $(BUILD)/test_roundtrip
	$(RUN) $(BUILD)/test_roundtrip

test_wraps: $(BUILD)/test_wraps
	$(RUN) $(BUILD)/test_wraps kat/kat_end$(END_N)_ref.txt

test_entropy: $(BUILD)/test_entropy
	$(RUN) $(BUILD)/test_entropy

test_all: test_kernels test_kat test_roundtrip test_wraps test_entropy

bench: $(BUILD)/bench
	$(RUN) $(BUILD)/bench

clean:
	rm -rf build

-include $(wildcard $(BUILD)/src/*/*.d $(BUILD)/src/*/*/*.d)

.PHONY: lib test_kernels test_kat test_roundtrip test_wraps test_entropy test_all bench clean
