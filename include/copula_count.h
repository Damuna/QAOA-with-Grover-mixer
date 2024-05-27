#ifndef COPULA_COUNT_H
#define COPULA_COUNT_H

/*
 * =============================================================================
 *                            includes
 * =============================================================================
 */

#include "knapsack.h"



/*
 * Function: 		qubit_count_copula_mixer
 * --------------------------------
 * Description:		This function counts the number of qubits that is
 *                  needed for an application of the copula mixer.
 * Parameters:
 *      k:          Underlying knapsack.
 * Returns:			Number of qubits.
 */
bit_t qubit_count_copula_mixer(const knapsack_t* k);


/*
 * Function: 		cycle_count_copula_mixer
 * --------------------------------
 * Description:		This function counts the number of cycles used
 *                  for an application of the copula mixer.
 * Parameters:
 *      k:          Underlying knapsack.
 * Returns:			Number of cycles.
 */
count_t cycle_count_copula_mixer(const knapsack_t* k);


/*
 * Function: 		gate_count_copula_mixer
 * --------------------------------
 * Description:		This function counts the number of gates used
 *                  for an application of the copula mixer.
 * Parameters:
 *      k:          Underlying knapsack.
 * Returns:			Number of gates.
 */
count_t gate_count_copula_mixer(const knapsack_t* k);


/*
 * Function: 		qubit_count_copula_qaoa
 * --------------------------------
 * Description:		This function counts the number of qubits that is
 *                  needed for an application of the full copula QAOA.
 * Parameters:
 *      k:          Underlying knapsack.
 * Returns:			Number of qubits.
 */
bit_t qubit_count_copula_qaoa(const knapsack_t* k);


/*
 * Function: 		cycle_count_copula_qaoa
 * --------------------------------
 * Description:		This function counts the number of cycles used for
 *                  for an application of the full copula QAOA.
 * Parameters:
 *      k:          Underlying knapsack.
 *      depth:      The desired depth.
 * Returns:			Number of cycles.
 */
count_t cycle_count_copula_qaoa(const knapsack_t* k, num_t depth);


/*
 * Function: 		gate_count_copula_qaoa
 * --------------------------------
 * Description:		This function counts the number of gates used for
 *                  for an application of the full copula QAOA.
 * Parameters:
 *      k:          Underlying knapsack.
 *      depth:      Desired depth.
 * Returns:			Number of gates.
 */
count_t gate_count_copula_qaoa(const knapsack_t* k, num_t depth);



#endif //COPULA_COUNT_H
