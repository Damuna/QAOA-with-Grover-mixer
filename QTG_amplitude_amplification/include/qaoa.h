#ifndef QAOA_H
#define QAOA_H

/* 
 * =============================================================================
 *                            includes
 * =============================================================================
 */

#include "stategen.h"
#include <complex.h>

/* 
 * =============================================================================
 *                            C++ check
 * =============================================================================
 */

#ifdef __cplusplus
extern "C" {
#endif

/* 
 * =============================================================================
 *                            type definitions
 * =============================================================================
 */

typedef double complex      cmplx;

typedef struct metastate_amplitude {
    choice_profit_t choice_profit;
    cmplx amplitude;
} metastate_amplitude_t;

typedef struct metastate_probability {
    choice_profit_t choice_profit;
    double probability;
} metastate_probability_t;

typedef struct qaoa_result {
    double optimal_value;
    metastate_probability_t *metastate_probability;
} qaoa_result_t;

/* 
 * =============================================================================
 *                            Quasi-Adiabatic Evolution
 * =============================================================================
 */

void phase_separation_unitary(metastate_amplitude_t*, num_t);


void mixing_unitary(metastate_amplitude_t*, double*, num_t);


void quasi_adiabatic_evolution(metastate_amplitude_t*, double*, num_t*);




#ifdef __cplusplus
}
#endif

#endif
