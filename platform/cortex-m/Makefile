# END KEM -- Cortex-M bare-metal build (arm-none-eabi-gcc, newlib + semihosting)
#
#   make                                   # END-512, Cortex-M4, portable kernels, all test/bench ELFs
#   make END_N=1024 BACKEND=armv7em        # END-1024 with the ARMv7E-M optimised kernels
#   make CPU=cortex-m33 QEMU_MACHINE=mps2-an505
#   make run-test_kat                      # run under QEMU (semihosting; cwd = project root for the KAT files)
#
# Real boards: pass FLASH_ORIGIN/FLASH_LENGTH/RAM_ORIGIN/RAM_LENGTH/STACK_SIZE via LDDEFS and replace
# randombytes_m.c / cycles with the board support.

END_N        ?= 512
BACKEND      ?= ref
CPU          ?= cortex-m4
ifeq ($(CPU),cortex-m4)
FPUFLAGS     ?= -mfloat-abi=hard -mfpu=fpv4-sp-d16
QEMU_MACHINE ?= mps2-an386
else
FPUFLAGS     ?= -mfloat-abi=hard -mfpu=fpv5-sp-d16
QEMU_MACHINE ?= mps2-an505
# QEMU mps2-an505 boots the M33 in Secure state: use the secure aliases of code SRAM / SSRAM2-3
LDDEFS       ?= -Wl,--defsym=FLASH_ORIGIN=0x10000000 -Wl,--defsym=RAM_ORIGIN=0x38000000
endif
OPT          ?= -O3
ROOT         := ../..
CC           := arm-none-eabi-gcc
OBJCOPY      := arm-none-eabi-objcopy
SIZE         := arm-none-eabi-size
QEMU         ?= qemu-system-arm
LDDEFS       ?=

CFLAGS   = -mcpu=$(CPU) -mthumb $(FPUFLAGS) $(OPT) -g -ffunction-sections -fdata-sections -fno-common \
           -Wall -Wextra -Wno-unused-parameter -Wno-unused-function -std=gnu99
# REAL_BOARD=1 enables the DWT cycle counter (QEMU does not model it -> END_NO_DWT by default)
ifeq ($(REAL_BOARD),1)
DWTFLAG =
else
DWTFLAG = -DEND_NO_DWT
endif
CPPFLAGS = -DEND_N=$(END_N) -DEND_BACKEND_NAME=\"$(BACKEND)\" -DEND_RANDOMBYTES_EXTERNAL -DEND_STATIC_EXPKEY $(DWTFLAG) \
           -DEND_STATIC_KEYGEN_TMP -DEND_PLATFORM_STACK -I$(ROOT)/src/common -I$(ROOT)/src/common/ntrugen -I$(ROOT)/src/$(BACKEND) \
           -I$(ROOT)/src/ref -I$(ROOT)/bench -I.
LDFLAGS  = -T link.ld --specs=rdimon.specs --specs=nano.specs -Wl,--gc-sections -nostartfiles $(LDDEFS) \
           -Wl,-Map=$@.map
LDLIBS   = -lc -lrdimon -lm

BUILD := build/$(BACKEND)-$(END_N)-$(CPU)

KERNELS      := ntt257 ntt3329 ntt12289 polyvec baseinv257
KERNEL_SRCS  := $(foreach k,$(KERNELS),$(if $(wildcard $(ROOT)/src/$(BACKEND)/$(k).c),$(ROOT)/src/$(BACKEND)/$(k).c,$(ROOT)/src/ref/$(k).c))
KERNEL_ASM   := $(wildcard $(ROOT)/src/$(BACKEND)/*.S)
COMMON_SRCS  := $(ROOT)/src/common/end_codec.c $(ROOT)/src/common/end_decode.c $(ROOT)/src/common/end_kem.c \
                $(ROOT)/src/common/end_keygen.c $(ROOT)/src/common/fips202.c
NTRUGEN_SRCS := $(ROOT)/src/common/ntrugen/ntru.c $(ROOT)/src/common/ntrugen/module_p.c $(ROOT)/src/common/ntrugen/big_integer.c \
                $(ROOT)/src/common/ntrugen/fixed_point.c $(ROOT)/src/common/ntrugen/poly_big.c
PLAT_SRCS    := startup.c randombytes_m.c $(ROOT)/bench/cycles.c
LIB_SRCS     := $(KERNEL_SRCS) $(COMMON_SRCS) $(NTRUGEN_SRCS) $(PLAT_SRCS)
LIB_OBJS     := $(patsubst %.c,$(BUILD)/%.o,$(subst $(ROOT)/,,$(LIB_SRCS))) $(patsubst %.S,$(BUILD)/%.o,$(subst $(ROOT)/,,$(KERNEL_ASM)))

TESTS := test_kernels test_kat test_roundtrip test_wraps bench
ELFS  := $(addprefix $(BUILD)/,$(addsuffix .elf,$(TESTS)))

all: $(ELFS)
	$(SIZE) $(ELFS)

$(BUILD)/%.o: $(ROOT)/%.c
	@mkdir -p $(dir $@)
	$(CC) $(CPPFLAGS) $(CFLAGS) -MMD -MP -c $< -o $@

$(BUILD)/%.o: $(ROOT)/%.S
	@mkdir -p $(dir $@)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(BUILD)/%.o: %.c
	@mkdir -p $(dir $@)
	$(CC) $(CPPFLAGS) $(CFLAGS) -MMD -MP -c $< -o $@

$(BUILD)/libend.a: $(LIB_OBJS)
	arm-none-eabi-ar rcs $@ $^

$(BUILD)/test_kernels.elf: $(ROOT)/test/test_kernels.c $(BUILD)/libend.a
	$(CC) $(CPPFLAGS) $(CFLAGS) $< $(BUILD)/libend.a -o $@ $(LDFLAGS) $(LDLIBS)
$(BUILD)/test_kat.elf: $(ROOT)/test/test_kat.c $(BUILD)/libend.a
	$(CC) $(CPPFLAGS) $(CFLAGS) $< $(BUILD)/libend.a -o $@ $(LDFLAGS) $(LDLIBS)
$(BUILD)/test_roundtrip.elf: $(ROOT)/test/test_roundtrip.c $(BUILD)/libend.a
	$(CC) $(CPPFLAGS) $(CFLAGS) $< $(BUILD)/libend.a -o $@ $(LDFLAGS) $(LDLIBS)
$(BUILD)/test_wraps.elf: $(ROOT)/test/test_wraps.c $(BUILD)/libend.a
	$(CC) $(CPPFLAGS) $(CFLAGS) $< $(BUILD)/libend.a -o $@ $(LDFLAGS) $(LDLIBS)
$(BUILD)/bench.elf: $(ROOT)/bench/bench.c $(BUILD)/libend.a
	$(CC) $(CPPFLAGS) $(CFLAGS) $< $(BUILD)/libend.a -o $@ $(LDFLAGS) $(LDLIBS)

# run under QEMU with semihosting; file paths in the tests are relative to the project root
QEMU_ARGS = -M $(QEMU_MACHINE) -cpu $(CPU) -nographic -semihosting-config enable=on,target=native -monitor none -serial none
run-%: $(BUILD)/%.elf
	cd $(ROOT) && $(QEMU) $(QEMU_ARGS) -kernel platform/cortex-m/$< ; true

clean:
	rm -rf build

-include $(wildcard $(BUILD)/*.d $(BUILD)/*/*.d $(BUILD)/*/*/*.d $(BUILD)/*/*/*/*.d)
.PHONY: all clean
