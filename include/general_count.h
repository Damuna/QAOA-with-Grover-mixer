#ifndef GENERAL_COUNT_H
#define GENERAL_COUNT_H

/*
 * =============================================================================
 *                            includes
 * =============================================================================
 */

#include "knapsack.h"


/*
 * =============================================================================
 *                            type definitions
 * =============================================================================
 */

/*
 * Struct:      resource_t
 * -----------------------
 * Description: This struct represents a resource counter for a complete run of
 *              the QMaxSearch algorithm.
 * Contents:
 *      qubit_count:        Number of qubits.
 *      cycle_count:        Number of quantum cycles.
 *      gate_count:         Number of quantum gates.
 *      cycle_count_decomp: Number of quantum cycles (decomposed Toffolis).
 *      gate_count_decomp:  Number of quantum gates (decomposed Toffolis).
 */
typedef struct resource {
    bit_t qubit_count;
    count_t cycle_count;
    count_t gate_count;
    count_t cycle_count_decomp;
    count_t gate_count_decomp;
} resource_t;


/*
 * Function: 		cycle_count_phase_separator
 * --------------------------------
 * Description:		This function counts the number of cycles used for an
 *                  application of the standard phase separation unitary.
 * Returns:			Number of cycles.
 */
count_t cycle_count_phase_separator();


/*
 * Function: 		gate_count_phase_separator
 * --------------------------------
 * Description:		This function counts the number of gates used for an
 *                  application of the standard phase separation unitary.
 * Parameters:
 *      k:          The underlying knapsack.
 * Returns:			Number of gates.
 */
count_t gate_count_phase_separator(const knapsack_t* k);


/*
 * Function: 		cycle_count_qaoa
 * --------------------------------
 * Description:		       This function provides an interface for the
 *                         number of cycles used for one execution of a
 *                         QAOA circuit.
 * Parameters:
 *      depth:             Desired depth of the QAOA circuit.
 *      cycle_count_mixer: Number of cycles needed for one application of
 *                         the mixing unitary.
 * Returns:			       Number of cycles.
 */
count_t cycle_count_qaoa(num_t depth, count_t cycle_count_mixer);


/*
 * Function: 		gate_count_qaoa
 * --------------------------------
 * Description:		       This function provides an interface for the
 *                         number of gates used for one execution of a
 *                         QAOA circuit.
 * Parameters:
 *      k:                 Underlying knapsack.
 *      depth:             Desired depth of the QAOA circuit.
 *      gate_count_mixer:  Number of gates needed for one application of
 *                         the mixing unitary.
 * Returns:			       Number of gates.
 */
count_t gate_count_qaoa(const knapsack_t* k, num_t depth, count_t gate_count_mixer);





#endif //GENERAL_COUNT_H
