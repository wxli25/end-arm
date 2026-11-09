#ifndef END_CYCLES_H
#define END_CYCLES_H
#include <stdint.h>
uint64_t cycles_now(void);
const char *cycles_unit(void);
uint64_t wall_ns(void);
#endif
