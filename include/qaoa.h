#ifndef QAOA_H
#define QAOA_H


/*
 * =============================================================================
 *                                includes
 * =============================================================================
 */

#include <complex.h>
#include <nlopt.h>
#include <time.h>
#include "general_count.h"
#include "qtg_count.h"
#include "copula_count.h"
#include "stategen.h"
#include "combowrp.h"


/*
 * =============================================================================
 *                                C++ check
 * =============================================================================
 */

#ifdef __cplusplus
extern "C" {
#endif


/*
 * =============================================================================
 *                              Type definitions
 * =============================================================================
 */

typedef double complex      cmplx;


/*
 * Struct:          cbs_t
 * ---------------------------
 * Description:     This struct equips a solution, identified via its profit, with a complex amplitude.
 * Contents:
 *      profit:     Profit (objective function value) of the computational basis state / feasible solution.
 *      amplitude:  Complex amplitude associated with this computational basis state / feasible solution.
 */
typedef struct cbs {
    num_t profit;
    cmplx amplitude;
} cbs_t;


/*
 * enum:            qaoa_type_t
 * ------------------------------------
 * Description:     Choose the QAOA method that you want to use, i.e. the mixing strategy.
 *
 * Contents:        Two types of QAOA mixers, QTG and Copula via van Dam.
 */
typedef enum qaoa_type {
    QTG,
    COPULA
} qaoa_type_t;


/*
 * enum:                opt_t
 * ------------------------------------
 * Description:         Choose the optimization that you want to use in qaoa
 *
 * Contents:            Three types of optimization considered
 */
typedef enum opt {
    BFGS,
    NELDER_MEAD,
    POWELL,
} opt_t;


/*
 * =============================================================================
 *                              Global variables
 * =============================================================================
 */

extern knapsack_t* kp;
extern qaoa_type_t qaoa_type;
extern size_t bias;
extern num_t depth;
extern double k;
extern double theta;
extern size_t num_states;
extern node_t *qtg_nodes;


/*
 * =============================================================================
 *                              QTG-specific
 * =============================================================================
 */

/*
* Function:            qtg_initial_state_prep
* --------------------
* Description:         Classical emulation of the application of the QTG to prepare the initial state of the QTG QAOA.
* Parameters:
*      angle_state:    Pointer to the current state before the application; will be updated.
*/

void qtg_initial_state_prep(cbs_t*);


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
void qtg_grover_mixer(cbs_t*, double);


/*
 * =============================================================================
 *                                 Copula-specific
 * =============================================================================
 */

/*
 * Function:            prob_dist
 * --------------------
 * Description:         Computes the value of the probability distribution that the initial state is to be biased
 *                      towards in the Copula-QAOA for a certain index.
 * Parameters:
 *      index:          Item in the knapsack whose value is wanted.
 *      kp:             Pointer to the underlying knapsack.
 *      k:              Hyperparameter of the probability distribution.
 * Returns:             The value of the probability distribution.
 */
double prob_dist(bit_t index);

/*
 * Function:            modified_objective_func
 * --------------------
 * Description:         Computes the value of the modified objective function, which returns zero in case the given
 *                      solution is infeasible and its standard objective function value otherwise.
 * Parameters:
 *      kp:             Pointer to the underlying knapsack.
 *      solution:       Variable assignment whose modified objective function value shall be calculated.
 * Returns:             The corresponding value of the modified objective function.
 */
int modified_objective_func(int solution);


/*
 * Function:            copula_initial_state_prep
 * --------------------
 * Description:         Prepares the initial state of the Copula-QAOA based on the underlying probability distribution.
 * Parameters:
 *      angle_state:    Pointer to the current state before the application; will be updated.
 */
void copula_initial_state_prep(cbs_t* angle_state);


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
void phase_separation_unitary(cbs_t*, double);

/*
 * Function:        quasiadiabatic_evolution
 * --------------------
 * Description:     Depending on the type of the QAOA (either QTG- or Copula-based), this function performs the
 *                  corresponding quasi-adiabatic evolution. In case of the QTG-approach, it emulates one execution of
 *                  the full circuit, whereas it actually simulates the gates for the Copula-QAOA. The basic structure
 *                  in both cases is as follows: After initializing the angle state in the first place, the initial
 *                  state is prepared. Afterwards, the depth specifies the number of alternating repitions of calling
 *                  the phase separation and mixing unitaries, respectively.
 * Parameters:
 *      angles:     Pointer to list of angles with length equaling twice the depth.
 * Returns:         The state with updated amplitudes after the alternating application.
 * Side Effect:     Allocates dynamically; pointer should eventually be freed.
 */
cbs_t* quasiadiabatic_evolution(const double* angles);


/*
 * =============================================================================
 *                                 Evaluation
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
double expectation_value(const cbs_t* angle_state);


/*
 * Function:            angles_to_value
 * --------------------
 * Description:         Computes the expectation value, given a set of angles, by executing the quasi-adiabatic
 *                      evolution first and processing the result afterwards.
 * Parameters:
 *      angles:         Pointer to list of angles with length equaling twice the depth.
 * Returns:             The negative expectation value corresponding to the specified angles.
 */
double angles_to_value(const double* angles);


/*
 * =============================================================================
 *                                 Optimization
 * =============================================================================
 */

/*
 * Function:            angles_to_value_nlopt
 * --------------------
 * Description:         Wrapper function to compute the expectation value, given a set of angles, in a shape that can be
 *                      processed to the NLOPT optimizer.
 * Parameters:
 *      n:              The number of parameters to optimize (not used).
 *      angles:         Pointer to list of angles with length equaling twice the depth.
 *      grad:           Gradient (not used).
 *      my_func_data:   Function data, e.g. parameters that shall not be optimized over (not used).
 * Returns:             The negative expectation value corresponding to the specified angles.
 */
double angles_to_value_nlopt(unsigned n, const double* angles, double* grad, void* my_func_data);


/*
 * =============================================================================
 *                               Export results
 * =============================================================================
 */

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
void write_plot_data_to_file(const cbs_t* angle_state, double solution_value, num_t optimal_solution_value);


/*
 * =============================================================================
 *                            Combination to full QAOA
 * =============================================================================
 */

/*
 * Function:                qaoa
 * --------------------
 * Description:             This is the main function for executing the QTG-induced or Copula-based QAOA, depending on
 *                          the input QAOA type. It wraps up all other functions defined in this file. First, it updates
 *                          all the global variables, i.e. the underlying knapsack, the QAOA type, the depth, the bias
 *                          and the hyperparameter k of the probability distribution neeeded in the Copula approach.
 *                          Then, it sorts the input knapsack by the relative profit of each item (i.e. the
 *                          profit divided by the weight). Depending on the context, it computes an integer greedy
 *                          solution for the given knapsack instance, which is needed for the application of the QTG.
 *                          The resulting lower bound is used in the following initial application of the QTG
 *                          for biasing. Then, the chosen classical optimizing routine minimizes the function
 *                          angles_to_value. The negative optimization result is the solution of the QAOA. For creating
 *                          expressive graphics, the quasi-adiabatic evolution must be run a last time with the optimal
 *                          angle values as input. In order to transform the not instance-agnostic profits in the
 *                          probability dictionary to approximation ratios, combo is executed once to compute the
 *                          optimal solution of the given knapsack instance. Both the final angle state and the exact
 *                          solution are passed to the function write_plot_data_to_file which exports the needed data to
 *                          an external file. Ultimately, the resources required to run the routine are counted.
 * Parameters:
 *      input_kp:           Pointer to the knapsack on which the QAOA is to be applied.
 *      input_qaoa_type:    The type of the QAOA, i.e. QTG or Copula.
 *      input_depth:        The depth of the QAOA.
 *      opt_type:           The classical method that shall be used for the optimization.
 *      input_bias:         The bias for the QTG.
 *      copula_k:           The hyperparameter k for the probability distribution in the Copula ansatz.
 * Returns:                 The negative solution value obtained from inserting the optimized angles returned by the
 *                          classical optimization routine, i.e. a positive result to the given knapsack problem.
 * Side Effect:             Frees the memory allocated in path_rep for the integer greedy solution.
 *                          Frees the memory allocated in qtg for its output nodes.
 *                          Frees the memory allocated in quasiadiabatic_evolution for the final QAOA state obtained
 *                          from inserting the optimized angles.
 */
double qaoa(
    knapsack_t* input_kp,
    qaoa_type_t input_qaoa_type,
    num_t input_depth,
    size_t input_bias,
    opt_t opt_type,
    double copula_k,
    double copula_theta
);



#ifdef __cplusplus
}
#endif

#endif //QAOA_H
