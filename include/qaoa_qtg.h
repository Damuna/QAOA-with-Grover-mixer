#ifndef QAOA_H
#define QAOA_H

/*
 * =============================================================================
 *                            includes
 * =============================================================================
 */

#include "combowrp.h"
#include "stategen.h"
#include "qtgcount.h"
#include <complex.h>
#include <nlopt.h>
#include <time.h>

/*
 * =============================================================================
 *                            C++ check
 * =============================================================================
 */

extern num_t dpth;
extern size_t num_smpls;
extern size_t num_states;
extern node_t *qtg_nodes;

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
 * Struct:              solution_amplitude_t
 * ---------------------------
 * Description:     This struct equips a feasible solution, identified via its profit, with a complex amplitude.
 * Contents:
 *      profit:     Profit (objective function value) of the computational basis state / feasible solution.
 *      amplitude:  Complex amplitude associated with this computational basis state / feasible solution.
 */
typedef struct feassol_ampl {
    num_t profit;
    cmplx amplitude;
} feassol_ampl_t;


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
typedef enum optimization_type {
    BFGS,
    NELDER_MEAD,
    POWELL,
} optimization_type_t;


/*
 * =============================================================================
 *                                     Utils
 * =============================================================================
 */

/*
 * Function:    random_value_on_windows_or_linux
 * ----------------------
 * Description: This function generates a random number in the interval [0,1) after seeding the random generator via the
 *              time. This is done via built-in functions according to the used operating system.
 * Returns:     The generated random number.
 */
double random_value_on_windows_or_linux();


/*
 * Function:                        write_plot_data_to_file
 * ----------------------
 * Description:                     This function is used to extract pairs of approximation ratio and associated
 *                                  probability from the final QAOA state and write them to an external file. The
 *                                  approximation ratios are calculated via the specified optimal solution value and the
 *                                  probabilities are obtained from the amplitudes of the feasible solutions.
 * Parameters:
 *      angle_state:                Pointer to the angle state whose data is to be exported.
 *      solution_value:             Final solution value returned by the QAOA.
 *      optimal_solution_value:     Optimal solution of the knapsack instance at hand.
 */
void write_plot_data_to_file(feassol_ampl_t*, double, num_t);


/*
 * =============================================================================
 *                            Quasi-Adiabatic Evolution
 * =============================================================================
 */

/*
 * Function:            phase_separation_unitary
 * --------------------
 * Description:         Classical emulation of the application of the phase separation unitary. Its underlying formula
 *                      is based on theoretical considerations. It owes its simplicity to the objective Hamiltonian
 *                      being diagonal in the computational basis by design.
 * Parameters:
 *      angle_state:    Pointer to the current state before the application; will be updated.
 *      gamma:          Value of the angle gamma that parametrizes the unitary.
 */
void phase_separation_unitary(feassol_ampl_t*, double);


/*
 * Function:            mixing_unitary
 * --------------------
 * Description:         Classical emulation of the application of the mixing unitary. Its underlying formula is based on
 *                      theoretical considerations. First, it calculates the scalar product of the state returned by the
 *                      initial QTG application and the current state. The result is used in an expression that comes
 *                      out when cleverly re-writing the action of the mixing unitary.
 * Parameters:
 *      angle_state:    Pointer to the current state before the application; will be updated.
 *      beta:           Value of the angle beta that parametrizes the unitary.
 */
void mixing_unitary(feassol_ampl_t*, double);


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
feassol_ampl_t* quasiadiabatic_evolution(const double*);


/*
 * =============================================================================
 *                          Evaluation & Optimization
 * =============================================================================
 */

/*
 * Function:            expectation_value
 * --------------------
 * Description:         Computes the expectation value of the objective Hamiltonian in a given state. To this end, the
 *                      function samples from the state and iteratively sums up the relative objective function value of
 *                      the measured feasible solutions (i.e. their associated profit divided by the number of samples).
 * Parameters:
 *      angle_state:    Pointer to the state to compute the expectation value for.
 * Returns:             The sum of all terms making the expectation value.
 */
double expectation_value(feassol_ampl_t*);


/*
 * Function:            angles_to_value
 * --------------------
 * Description:         This function computes the expectation value for a given set of angles. To do so, it calls
 *                      quasiadiabatic_evolution as well as expectation_value. This is the function that shall be
 *                      optimized (minimized) by the chosen classical routine. To this end, it needs the present format.
 * Parameters:
 *      n:              The number of parameters to optimize (not used).
 *      angles:         Pointer to list of angles with length equaling twice the depth.
 *      grad:           Gradient (not used).
 *      my_func_data:   Function data, e.g. parameters that shall not be optimized over (not used).
 * Returns:             The negative expectation value corresponding to the specified angles.
 * Side Effect:         Frees the memory allocated for the angle state in quasiadiabatic_evolution.
 */
double angles_to_value(unsigned, const double*, double*, void*);



/*
 * =============================================================================
 *                            Combination to full QAOA
 * =============================================================================
 */

/*
 * Function:                qaoa_qtg
 * --------------------
 * Description:             This is the main function to execute the QTG-induced QAOA, wrapping up all other functions
 *                          defined here. First, it updates the global variables holding the depth and the number of
 *                          samples. Then, it sorts the input knapsack by the relative profit of each item (i.e. the
 *                          profit divided by the weight) and computes an integer greedy solution for the given knapsack
 *                          instance. The resulting lower bound is used in the following initial application of the QTG
 *                          for biasing. Then, the chosen classical optimizing routine minimizes the function
 *                          angles_to_value. The negative optimization result is the solution of the QAOA. For creating
 *                          expressive graphics, the quasi-adiabatic evolution must be run a last time with the optimal
 *                          angle values as input. In order to transform the not instance-agnostic profits in the
 *                          probability dictionary to approximation ratios, combo is executed once to compute the
 *                          optimal solution of the given knapsack instance. Both the final angle state and the exact
 *                          solution are passed to the function write_plot_data_to_file which exports the needed data to
 *                          an external file. Ultimately, the cycles required to execute the mixing unitary are counted.
 * Parameters:
 *      k:                  Pointer to the knapsack on which the QAOA is to be applied.
 *      depth:              The depth of the QAOA.
 *      bias:               The bias for the QTG.
 *      num_samples:        The number of samples determining the number of measurements in each iteration.
 *      optimizationType:   The classical optimizer that shall be used for the optimization.
 *      m:                  Precision of the fine grid search
 * Returns:                 The negative solution value obtained from inserting the optimized angles returned by the
 *                          classical optimization routine.
 * Side Effect:             Frees the memory allocated in path_rep for the integer greedy solution.
 *                          Frees the memory allocated in qtg for its output nodes.
 *                          Frees the memory allocated in quasiadiabatic_evolution for the final QAOA state obtained
 *                          from inserting the optimized angles.
 */
double qaoa_qtg(knapsack_t* k, num_t depth, size_t bias, size_t num_samples, optimization_type_t optimization_type, int m);





#ifdef __cplusplus
}
#endif


#endif //QAOA_H
