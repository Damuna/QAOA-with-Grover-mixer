#ifndef QAOA_H
#define QAOA_H

/*
 * =============================================================================
 *                            includes
 * =============================================================================
 */

#include "stategen.h"
#include <complex.h>
#include <nlopt.h>
#include <time.h>

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


/*
 * Struct:              metastate_amplitude_t
 * ---------------------------
 * Description:         This struct equips a computational basis state (identified via its corresponding bitstring and
 *                      its profit) with a complex amplitude.
 * Contents:
 *      choice_profit:  The meta information of a computational basis state consisting of variable assignment and profit
 *      amplitude:      Complex amplitude associated with this computational basis state.
 */
typedef struct metastate_amplitude {
    choice_profit_t choice_profit;
    cmplx amplitude;
} metastate_amplitude_t;


/*
 * Struct:              metastate_probability_t
 * ---------------------------
 * Description:         This struct equips a computational basis state (identified via its corresponding bitstring and
 *                      its profit) with a probability.
 * Contents:
 *      choice_profit:  The meta information of a computational basis state consisting of variable assignment and profit
*       probability:    Probability associated with this computational basis state.
 */
typedef struct metastate_probability {
    choice_profit_t choice_profit;
    double probability;
} metastate_probability_t;


/*
 * Struct:              approxratio_probability_t
 * ---------------------------
 * Description:         This struct is used to assign a probability to an approximation ratio.
 * Contents:
 *      approx_ratio:   The value of the approximation ratio.
 *      probability:    Probability associated with the approximation ratio.
 */
typedef struct approxratio_probability {
    double approx_ratio;
    double probability;
} approxratio_probability_t;


/*
 * Struct:              qaoa_result_t
 * ---------------------------
 * Description:         This struct equips a computational basis state (identified via its corresponding bitstring and
 *                      its profit) with a complex amplitude.
 * Contents:
 *      solution_value:     The value returned by the classical optimization routine employed by the QTG-induced QAOA.
 *      approxratio_prob:   List of approximation ratio - probability pairs.
 */
typedef struct qaoa_result {
    double solution_value;
    approxratio_probability_t *approxratio_prob;
} qaoa_result_t;

/*
 * enum:                OptimizationType
 * ------------------------------------
 * Description:         Choose the optimization that you want to use in qaoa
 *
 * Contents:            Three types of optimization considered
 */
enum OptimizationType {
    BFGS,
    NELDER_MEAD,
    POWELL,
};

/*
 * =============================================================================
 *                            Free newly defined types
 * =============================================================================
 */


/*
 * Function:    free_metastates_probability
 * ----------------------
 * Description: This function frees a dynamically allocated list of struct metastate_probability_t together with the
 *              vector contained two levels below.
 * Parameter:   Pointer to the metastate_probability_t struct to be freed.
 */
void free_metastates_probability(metastate_probability_t*, size_t);


/*
 * Function:    free_metastates_probability
 * ----------------------
 * Description: This function frees a dynamically allocated list of struct metastate_amplitude_t together with the
 *              vector contained two levels below.
 * Parameter:   Pointer to the metastate_amplitude_t struct to be freed.
 */
void free_metastates_amplitude(metastate_amplitude_t*, size_t);



/*
 * =============================================================================
 *                            Quasi-Adiabatic Evolution
 * =============================================================================
 */


/*
 * Function:        phase_separation_unitary
 * --------------------
 * Description:     Classical emulation of the application of the phase separation unitary. Its underlying formula is
 *                  based on theoretical considerations. It owes its simplicity to the objective Hamiltonian being
 *                  diagonal in the computational basis by design.
 * Parameters:
 *      parameter1: Pointer to the current state before the application; will be updated.
 *      parameter2: Value of the angle gamma that parametrizes the unitary.
 */
void phase_separation_unitary(metastate_amplitude_t*, double);


/*
 * Function:        mixing_unitary
 * --------------------
 * Description:     Classical emulation of the application of the mixing unitary. Its underlying formula is based on
 *                  theoretical considerations. First, it calculates the scalar product of the state returned by the
 *                  initial QTG application and the current state. The result is used in an expression that comes out
 *                  when cleverly re-writing the action of the mixing unitary.
 * Parameters:
 *      parameter1: Pointer to the current state before the application; will be updated.
 *      parameter2: Value of the angle beta that parametrizes the unitary.
 */
void mixing_unitary(metastate_amplitude_t*, double);


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
 *      parameter1: Pointer to list of angles with length equaling twice the depth.
 * Returns:         The state with updated amplitudes after the alternating application.
 * Side Effect:     Allocates dynamically; pointer should eventually be freed.
 */
metastate_amplitude_t* quasiadiabatic_evolution(double*);


/*
 * =============================================================================
 *                          Measurement & Expectation Value
 * =============================================================================
 */


/*
 * Function:        measurement
 * --------------------
 * Description:     Measures a given state in the sense that it draws from its associated probability distribution via
 *                  the inverse transform sampling.
 * Parameters:
 *      parameter1: Pointer to the state to measure.
 * Returns:         The measured basis state identified via its index the original state or -1 in case of error.
 */
int measurement(metastate_amplitude_t*);


/*
 * Function:        sample_for_probabilities
 * --------------------
 * Description:     This function first transfers the metadata of each computational basis state (i.e. bitstring and
 *                  associated profit) to the less memory requiring format that suffices to store probabilities.
 *                  Afterward, it measures the given state multiple times and updates - based on the measurement
 *                  outcome - the probability of each basis state (identified via its index) accordingly. The number of
 *                  samples is determined by the global variable num_smpls.
 * Parameters:
 *      parameter1: Pointer to the state to sample from.
 * Returns:         A dictionary in which each computational basis state (identified via its bitstring and profit) is
 *                  assigned a probability to get measured.
 * Side Effect:     Allocates memory dynamically for the probability dictionary; pointer should eventually be freed.
 */
metastate_probability_t* sample_for_probabilities(metastate_amplitude_t*);


/*
 * Function:        expectation_value
 * --------------------
 * Description:     Computes the expectation value of the objective Hamiltonian in a given state. To this end, the
 *                  function first samples from the state and sums up the product of obtained probability and profit
 *                  for each computational basis state.
 * Parameters:
 *      parameter1: Pointer to the state to sample from.
 * Returns:         The sum of all terms making the expectation value.
 * Side Effect:     Frees the memory allocated for the probability dict in sample_for_probabilities.
 */
double expectation_value(metastate_amplitude_t*);


/*
 * Function:        angles_to_value
 * --------------------
 * Description:     This function computes the expectation value for a given set of angles. To do so, it calls
 *                  quasiadiabatic_evolution as well as expectation_value. This is the function that shall be optimized
 *                  (minimized) by the chosen classical routine.
 * Parameters:
 *      parameter1: Pointer to list of angles with length equaling twice the depth.
 * Returns:         The negative expectation value corresponding to the specified angles.
 * Side Effect:     Frees the memory allocated for the angle state in quasiadiabatic_evolution.
 */
double angles_to_value(double*);



/*
 * =============================================================================
 *                            Combination to full QAOA
 * =============================================================================
 */


/*
 * Function:        qaoa_qtg
 * --------------------
 * Description:     This is the main function to execute the QTG-induced QAOA, wrapping up all other functions defined
 *                  here. First, it updates the global variables holding the depth and the number of samples. Then, it
 *                  sorts the input knapsack by the relative profit of each item (i.e. the profit divided by the weight)
 *                  and computes an integer greedy solution for the given knapsack instance. The resulting lower bound
 *                  is used in the following initial application of the QTG for biasing. The output is transformed to
 *                  a different shape as the remaining cost is no longer necessary to be stored. With this, the chosen
 *                  classical optimizing routine minimizes the function angles_to_value with randomly chosen initial
 *                  conditions for the angles whose length equals twice the depth. The negative optimization result is
 *                  the solution of the QAOA. For creating expressive graphics, the quasi-adiabatic evolution must be
 *                  run a last time with the optimal angle values as input. Ultimately, it is sampled again from the
 *                  resulting angle state. In order to transform the not instance-agnostic profits in the probability
 *                  dictionary to approximation ratios, combo is executed once to compute the optimal solution of the
 *                  given knapsack instance. Also, the cycles required to execute the Grover mixing unitary are counted.
 * Parameters:
 *      parameter1: Pointer to the knapsack on which the QAOA is to be applied.
 *      parameter2: The depth of the QAOA.
 *      parameter3: The bias for the QTG.
 *      parameter4: The number of samples determining the number of measurements of the angle state in each iteration.
 * Returns:         The solution value returned by the classical optimization routine together with a dictionary in
 *                  which approximation ratios are mapped to probabilities.
 * Side Effect:     Allocates memory dynamically for the more memory-efficient QTG output.
 * Side Effect:     Frees the memory allocated in the QTG after transforming it to the more memory-efficient shape.
 *                  Frees the memory allocated earlier for the more memory-efficient QTG output.
 */
qaoa_result_t qaoa_qtg(knapsack_t*, num_t, size_t, size_t, enum OptimizationType);





#ifdef __cplusplus
}
#endif


#endif //QAOA_H
