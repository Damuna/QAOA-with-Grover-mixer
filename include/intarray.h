#ifndef INTARRAY_H
#define INTARRAY_H

#include <stdio.h>
#include <stdlib.h>

typedef struct {
    size_t n;
    size_t bits;
    uint64_t* part;
} array_t;

#define sw_init(A, B) \
    do { \
        (A).bits = (B); \
        size_t BITS = (B); \
        size_t N = 0; \
        do { \
            BITS >>= 6; \
            ++N; \
        } while(0 != BITS); \
        (A).n = N; \
        (A).part = calloc((A).n, sizeof(uint64_t)); \
    } while (0)

#define sw_setbit(A, B) \
    do { \
        size_t WHICHPART = (B) >> 6; \
        size_t WHICHBIT = (B) - (WHICHPART << 6); \
        (A).part[WHICHPART] |= (1ULL << WHICHBIT); \
    } while (0)

#define sw_clrbit(A, B) \
    do { \
        size_t WHICHPART = (B) >> 6; \
        size_t WHICHBIT = (B) - (WHICHPART << 6); \
        (A).part[WHICHPART] &= (~(1ULL << WHICHBIT)); \
    } while (0)

#define sw_copy(A, B) \
    do { \
        sw_init((A), (B).bits); \
        for (size_t LOOPINDEX = 0; LOOPINDEX < (B).bits; ++LOOPINDEX) { \
            if (sw_tstbit((B), LOOPINDEX) == 1) { \
                sw_setbit((A), LOOPINDEX); \
            } \
        } \
    } while (0)

#define sw_tstbit(A, B) \
    ({ \
        size_t which_part = (B) >> 6; \
        size_t which_bit = (B) - (which_part << 6); \
        int res = (((A).part[which_part] & (1ULL << which_bit)) >> which_bit); \
        res; \
    })

#define sw_set_ui_0(A) \
    do { \
        for (size_t LOOPINDEX = 0; LOOPINDEX < (A).n; ++LOOPINDEX) { \
            (A).part[LOOPINDEX] = 0; \
        } \
    } while (0)

#define sw_print(A) \
    do { \
        for (size_t LOOPINDEX = 0; LOOPINDEX < (A).bits; ++LOOPINDEX) { \
            printf("%d", sw_tstbit((A), LOOPINDEX)); \
        } \
    } while (0)

#define sw_println(A) \
    do { \
        for (size_t LOOPINDEX = 0; LOOPINDEX < (A).bits; ++LOOPINDEX) { \
            printf("%d", sw_tstbit((A), LOOPINDEX)); \
        } \
        putchar('\n'); \
    } while (0)

#define sw_clear(A) (free((A).part))

#define sw_cmp(A1, A2) \
    ({ \
        bool result = true; \
        for (size_t LOOPINDEX = 0; LOOPINDEX < (A1).n; ++LOOPINDEX) { \
            if ((A1).part[LOOPINDEX] != (A2).part[LOOPINDEX]) { \
                result = false; \
                break; \
            } \
        } \
        result; \
    })

#endif
