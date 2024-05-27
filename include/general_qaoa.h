#ifndef GENERAL_QAOA_H
#define GENERAL_QAOA_H


/*
 * =============================================================================
 *                            includes
 * =============================================================================
 */

#include "knapsack.h"
#include <complex.h>


/*
 * =============================================================================
 *                            type definitions
 * =============================================================================
 */

typedef double complex      cmplx;


/*
 * Struct:              solution_amplitude_t
 * ---------------------------
 * Description:     This struct equips a feasible solution, identified via its profit, with a complex amplitude.
 * Contents:
 *      profit:     Profit (objective function value) of the computational basis state / feasible solution.
 *      amplitude:  Complex amplitude associated with this computational basis state / feasible solution.
 */
typedef struct cbs {
    num_t profit;
    cmplx amplitude;
} cbs_t;


/*
 * =============================================================================
 *                            Quasi-Adiabatic Evolution
 * =============================================================================
 */

/*
 * Function:        quasiadiabatic_evolution
 * --------------------
 * Description:     This function emulates one execution of the circuit representing the quasi-adiabatic evolution.
 *                  First, it defines separate pointers for the gamma and the beta values from the input pointer to the
 *                  angles. Afterward, it allocates the required memory and transforms the output from the initial QTG
 *                  application to the format we need for applying the two unitaries, i.e. the real-valued probabilities
 *                  are mapped to complex numbers with vanishing imaginary part. Finally, the phase separation and the
 *                  mixing unitaries are applied alternately where the number of iterations is determined by the depth.
 * Parameters:
 *      angles:     Pointer to list of angles with length equaling twice the depth.
 * Returns:         The state with updated amplitudes after the alternating application.
 * Side Effect:     Allocates dynamically; pointer should eventually be freed.
 */
cbs_t* quasiadiabatic_evolution(
    void (*initial_state_prep)(cbs_t*),
    void (*phase_separation_unitary)(cbs_t*, double),
    void (*mixing_unitary)(cbs_t*, double),
    const num_t depth,
    const size_t num_states,
    const double* angles
);


/*
 * =============================================================================
 *                          Evaluation & Optimization
 * =============================================================================
 */

/*
 * Function:            expectation_value
 * --------------------
 * Description:         Computes the expectation value of the objective Hamiltonian in a given state based on the actual
                        amplitudes of the state.
 * Parameters:
 *      angle_state:    Pointer to the state to compute the expectation value for.
 *      num_states:     The number of states included in the angle state.
 * Returns:             The sum of all terms making the expectation value.
 */
double expectation_value(const size_t num_states, const cbs_t* angle_state);


/*
 * Function:            angles_to_value
 * --------------------
 * Description:         Wrapper function that computes the expectation value by executing the quasi-adiabatic evolution
 *                      first and processing the result afterwards.
 * Parameters:
 *      angles:         Pointer to list of angles with length equaling twice the depth.
 *      num_states:     The number of states included in the angle state.
 * Returns:             The negative expectation value corresponding to the specified angles.
 */
double angles_to_value(
    void (*initial_state_prep)(cbs_t*),
    void (*phase_separation_unitary)(cbs_t*, double),
    void (*mixing_unitary)(cbs_t*, double),
    const num_t depth,
    const size_t num_states,
    const double* angles
);


#endif //GENERAL_QAOA_H
