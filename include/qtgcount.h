#ifndef QTGCOUNT_H
#define QTGCOUNT_H

/* 
 * =============================================================================
 *                            includes
 * =============================================================================
 */

#include "knapsack.h"

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

typedef enum qft {
    COPPERSMITH,
} qft_t;

typedef enum mc {
    TOFFOLI,
} mc_t;

/* 
 * =============================================================================
 *                            enum names
 * =============================================================================
 */

/*
 * Function:    get_qft_name
 * -------------------------
 * Description: This function returns a QFT method's name as a string.
 * Parameter:   Implementation method for the QFT.
 * Returns:     QFT method's name as a string.
 */
const char* get_qft_name(qft_t);

/*
 * Function:    get_mc_name
 * ------------------------
 * Description: This function returns a decomposition method's name for
 *              multi-controls as a string.
 * Parameter:   Decomposition method for the multi-control.
 * Returns:     Decomposition method's name as a string.
 */
const char* get_mc_name(mc_t);

/* 
 * =============================================================================
 *                                  utils
 * =============================================================================
 */

/*
 * Function: 	num_bits
 * ---------------------
 * Description:	This function returns the number of bits of the given input.
 * Parameter:	Number whose amount of bits should be evaluated.
 * Returns:		Number of bits of the given input.
 */
bit_t num_bits(num_t);

/*
 * Function: 	num_set_bits
 * ---------------------
 * Description:	This function returns the number of set bits of the given input, i.e. those with value 1.
 * Parameter:	Number whose amount of set bits should be evaluated.
 * Returns:		Number of set bits of the given input.
 */
bit_t num_set_bits(num_t);

/*
 * Function: 	lso
 * ----------------
 * Description:	This function returns the index of the least significant one.
 * Parameter:	Number whose least significant one should be evaluated.
 * Returns:		Least significant one of the input number.
 */
bit_t lso(num_t);

/* 
 * =============================================================================
 *                            counts
 * =============================================================================
 */

/*
 * Function:    path_reg_size
 * --------------------------
 * Description: This function calculates the size of the path register used for
 *              the QTG on a knapsack instance. Its size is given by the number
 *              of items assigned to the knapsack.
 * Parameters:  Pointer to knapsack that should be considered.
 * Returns:     Size (number of qubits) of the path register.
 */
bit_t path_reg_size(const knapsack_t*);

/*
 * Function:    cost_reg_size
 * --------------------------
 * Description: This function calculates the size of the cost register used for
 *              the QTG on a knapsack instance. Its size is given by the number
 *              of bits necessary to express the capacity.
 * Parameters:  Pointer to knapsack that should be considered.
 * Returns:     Size (number of qubits) of the cost register.
 */
bit_t cost_reg_size(const knapsack_t*);

/*
 * Function: 		anc_count_qft
 * ------------------------------
 * Description:		This function counts the ancilla qubits of the given QFT
 *                  implementation on a given register. 
 * Parameters:
 *		parameter1:	Size of the register the QFT is performed on.
 *		parameter2:	Implementation method for the QFT.
 * Returns:			Number of ancilla qubits of the QFT implementation.
 */
bit_t anc_count_qft(bit_t, qft_t);

/*
 * Function: 		cycle_count_qft
 * --------------------------------
 * Description:		This function counts the cycles of the given QFT
 *                  implementation on a given register. 
 * Parameters:
 *		parameter1:	Size of the register the QFT is performed on.
 *		parameter2:	Implementation method for the QFT.
 *      parameter3: Whether toffoli gates should be decomposed or not.
 * Returns:			Number of cycles of the QFT implementation.
 */
count_t cycle_count_qft(bit_t, qft_t, bool_t);

/*
 * Function: 		gate_count_qft
 * -------------------------------
 * Description:		This function counts the gates of the given QFT
 *                  implementation on a given register. 
 * Parameters:
 *		parameter1:	Size of the register the QFT is performed on.
 *		parameter2:	Implementation method for the QFT.
 *      parameter3: Whether toffoli gates should be decomposed or not.
 * Returns:			Number of gates of the QFT implementation.
 */
count_t gate_count_qft(bit_t, qft_t, bool_t);

/*
 * Function: 		anc_count_mc
 * -----------------------------
 * Description:		This function counts the ancilla qubits for decomposing a
 *                  multi-controlled single-qubit-gate. 
 * Parameters:
 *		parameter1:	Number of control qubits.
 *		parameter2:	Decomposition method for the multi-control.
 * Returns:			Number of ancilla qubits of decomposed multi-control.
 */
bit_t anc_count_mc(bit_t, mc_t);

/*
 * Function: 		cycle_count_mc
 * -------------------------------
 * Description:		This function counts the cycles for decomposing a
 *                  multi-controlled single-qubit-gate.  
 * Parameters:
 *		parameter1:	Number of control qubits.
 *		parameter2:	Decomposition method for the multi-control.
 *      parameter3: Whether toffoli gates should be decomposed or not.
 * Returns:			Number of cycles of the multi-control decomposition.
 */
count_t cycle_count_mc(bit_t, mc_t, bool_t);

/*
 * Function: 		gate_count_mc
 * ------------------------------
 * Description:		This function counts the gates for decomposing a
 *                  multi-controlled single-qubit-gate. 
 * Parameters:
 *		parameter1:	Number of control qubits.
 *		parameter2:	Decomposition method for the multi-control.
 *      parameter3: Whether toffoli gates should be decomposed or not.
 * Returns:			Number of gates of the multi-control decomposition.
 */
count_t gate_count_mc(bit_t, mc_t, bool_t);

/*
 * Function:        anc_count_comp
 * -------------------------------
 * Description:     This function counts the ancilla qubits of an integer
 *                  comparison on a given register. 
 * Parameters:
 *      parameter1: Size of the compared register.
 *      parameter2: Integer that is being compared.
 *      parameter3: Decomposition method for the arising multi-control.
 *      parameter4: Whether the comparison or its negation should be tested.
 * Returns:         Number of ancilla of the integer comparison method.
 */
bit_t anc_count_comp(bit_t, num_t, mc_t, bool_t);

/*
 * Function:        cycle_count_comp
 * ---------------------------------
 * Description:     This function counts the cycles of an integer comparison on
 *                  a given register. 
 * Parameters:
 *      parameter1: Size of the compared register.
 *      parameter2: Integer that is being compared.
 *      parameter3: Decomposition method for the arising multi-control.
 *      parameter4: Whether the comparison or its negation should be tested.
 *      parameter5: Whether toffoli gates should be decomposed or not.
 * Returns:         Number of cycles of the integer comparison method.
 */
count_t cycle_count_comp(bit_t, num_t, mc_t, bool_t, bool_t);

/*
 * Function:        gate_count_comp
 * --------------------------------
 * Description:     This function counts the gates of an integer comparison on a
 *                  given register. 
 * Parameters:
 *      parameter1: Size of the compared register.
 *      parameter2: Integer that is being compared.
 *      parameter3: Decomposition method for the arising multi-control.
 *      parameter4: Whether the comparison or its negation should be tested.
 *      parameter5: Whether toffoli gates should be decomposed or not.
 * Returns:         Number of gates of the integer comparison method.
 */
count_t gate_count_comp(bit_t, num_t, mc_t, bool_t, bool_t);

/*
 * Function: 		qubit_count_qtg
 * --------------------------------
 * Description:		This function counts the qubits used for an application of
 *                  the QTG. The implementation of the QTG's subroutines alsos
 *                  have to be specified. The main registers as well as any
 *                  ancilla register contribute to the total number.
 * Parameters:
 *		parameter:	Pointer to knapsack that should be considered.
 * Returns:			Number of qubits of one application of the QTG.
 */
bit_t qubit_count_qtg(const knapsack_t*);

/*
 * Function: 		cycle_count_qtg
 * --------------------------------
 * Description:		This function counts the number of cycles used for an
 *                  application of the QTG. The implementation of the QTG's
 *                  subroutines also have to be specified. Subroutines on
 *                  distinct registers are parallelized whenever possible.
 * Parameters:
 *      parameter1: Pointer to knapsack that should be considered.
 *      parameter2: Implementation method for the QFT.
 *      parameter3: Decomposition method for any multi-controls.
 *      parameter4: Whether toffoli gates should be decomposed or not.
 * Returns:			Number of cycles of one application of the QTG.
 */
count_t cycle_count_qtg(const knapsack_t*, qft_t, mc_t, bool_t);

/*
 * Function: 		gate_count_qtg
 * -------------------------------
 * Description:		This function counts the number of gates used for an
 *                  application of the QTG. The implementation of the QTG's
 *                  subroutines also have to be specified.
 * Parameters:
 *      parameter1: Pointer to knapsack that should be considered.
 *      parameter2: Implementation method for the QFT.
 *      parameter3: Decomposition method for any multi-controls.
 *      parameter4: Whether toffoli gates should be decomposed or not.
 * Returns:			Number of gates of one application of the QTG.
 */
count_t gate_count_qtg(const knapsack_t*, qft_t, mc_t, bool_t);

/*
 * Function: 		qubit_count_qtg_mixer
 * --------------------------------
 * Description:		This function counts the qubits used for an application of
 *                  the Grover mixer induced by the QTG. The main registers as
 *                  well as any ancilla register contribute to the total number.
 * Parameters:
 *		k:	        Pointer to knapsack that should be considered.
 * Returns:			Number of qubits of one application of the QTG mixer.
 */
bit_t qubit_count_qtg_mixer(const knapsack_t*);

/*
 * Function: 		cycle_count_qtg_mixer
 * --------------------------------
 * Description:		This function counts the number of cycles used for an
 *                  application of the Grover mixer induced by the QTG. The
 *                  implementation of the QTG's subroutines also have to be
 *                  specified. Subroutines on distinct registers are parallelized
 *                  whenever possible.
 * Parameters:
 *      k:          Pointer to knapsack that should be considered.
 *      method_qft: Implementation method for the QFT.
 *      method_mc:  Decomposition method for any multi-controls.
 *      tof_decomp: Whether toffoli gates should be decomposed or not.
 * Returns:			Number of cycles of one application of the QTG mixer.
 */
count_t cycle_count_qtg_mixer(const knapsack_t*, qft_t, mc_t, bool_t);

/*
 * Function: 		gate_count_qtg_mixer
 * -------------------------------
 * Description:		This function counts the number of gates used for an
 *                  application of the Grover mixer induced by the QTG. The
 *                  implementation of the QTG's subroutines also have to be
 *                  specified.
 * Parameters:
 *      k:          Pointer to knapsack that should be considered.
 *      method_qft: Implementation method for the QFT.
 *      method_mc:  Decomposition method for any multi-controls.
 *      tof_decomp: Whether toffoli gates should be decomposed or not.
 * Returns:			Number of gates of one application of the QTG mixer.
 */
count_t gate_count_qtg_mixer(const knapsack_t*, qft_t, mc_t, bool_t);

/*
 * Function:        print-qtg_counts
 * ---------------------------------
 * Description:     This function prints the entire QTG statistics for a given
 *                  Knapsack instance with the specified circuit construction
 *                  methods.
 * Parameters:
 *      parameter1: Pointer to knapsack that should be considered.
 *      parameter2: Implementation method for the QFT.
 *      parameter3: Decomposition method for any multi-controls.
 *      parameter4: Whether toffoli gates should be decomposed or not.
 */
void print_qtg_counts(const knapsack_t*, qft_t, mc_t, bool_t);

#ifdef __cplusplus
}
#endif

#endif
