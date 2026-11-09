/*
 * Linker script for Cortex-M4/M33 bare-metal builds.
 * Default memory map = QEMU mps2-an386 / mps2-an505 (4 MB "flash" at 0, 4 MB RAM at 0x20000000).
 * For a real board (e.g. NUCLEO-L4R5ZI: 2 MB flash @ 0x08000000, 640 KB SRAM @ 0x20000000) override
 * FLASH_ORIGIN/FLASH_LENGTH/RAM_ORIGIN/RAM_LENGTH with -Wl,--defsym or edit the MEMORY block.
 */
FLASH_ORIGIN = DEFINED(FLASH_ORIGIN) ? FLASH_ORIGIN : 0x00000000;
FLASH_LENGTH = DEFINED(FLASH_LENGTH) ? FLASH_LENGTH : 4M;
RAM_ORIGIN   = DEFINED(RAM_ORIGIN)   ? RAM_ORIGIN   : 0x20000000;
RAM_LENGTH   = DEFINED(RAM_LENGTH)   ? RAM_LENGTH   : 4M;
STACK_SIZE   = DEFINED(STACK_SIZE)   ? STACK_SIZE   : 128K;

MEMORY
{
    FLASH (rx)  : ORIGIN = FLASH_ORIGIN, LENGTH = FLASH_LENGTH
    RAM   (rwx) : ORIGIN = RAM_ORIGIN,   LENGTH = RAM_LENGTH
}

ENTRY(Reset_Handler)

SECTIONS
{
    .text :
    {
        KEEP(*(.vectors))
        *(.text*)
        *(.rodata*)
        . = ALIGN(4);
    } > FLASH

    .ARM.extab : { *(.ARM.extab* .gnu.linkonce.armextab.*) } > FLASH
    .ARM.exidx : { __exidx_start = .; *(.ARM.exidx* .gnu.linkonce.armexidx.*) __exidx_end = .; } > FLASH
    . = ALIGN(4);
    _etext = .;

    .data : AT(_etext)
    {
        _sdata = .;
        *(.data*)
        . = ALIGN(4);
        _edata = .;
    } > RAM
    _sidata = LOADADDR(.data);

    .bss :
    {
        _sbss = .;
        *(.bss*)
        *(COMMON)
        . = ALIGN(8);
        _ebss = .;
    } > RAM

    . = ALIGN(8);
    _end = .;
    end = _end;                     /* heap start for newlib's sbrk */
    _heap_limit = ORIGIN(RAM) + LENGTH(RAM) - STACK_SIZE;
    _estack = ORIGIN(RAM) + LENGTH(RAM);
    _sstack = _estack - STACK_SIZE;
}
