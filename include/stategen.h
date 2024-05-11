#ifndef STATEGEN_H
#define STATEGEN_H

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
 * Struct:      node_t
 * -------------------
 * Description: This struct represents a node in a knapsack's decision tree.
 * Contents:
 *      path:   Partial path associated with the node. Items of unexplored
 *              layers are conventionally not set.
 *      prob:   Probability of traversing this node.
 */
typedef struct node {
    path_t path;
    double prob;
} node_t;


/* 
 * =============================================================================
 *                            free nodes
 * =============================================================================
 */

void free_nodes(node_t*, size_t);

/* 
 * =============================================================================
 *                            branch probability
 * =============================================================================
 */

/*
 * Function:        branch_prob
 * ----------------------------
 * Description:     This function returns the factor by which a child node's
 *                  probability differs from its parent node's probability while
 *                  traversing the decision tree of a knapsack instance.
 * Parameters:
 *      parameter1: Current item/layer at which the branching happens.
 *      parameter2: Bias towards certain branch.
 *      parameter3: Whether left or right subtree is considered.
 *      parameter4: Bit string representation of current optimal path.
 * Returns:         Branching probability.
 */
double branch_prob(bit_t, size_t, bool_t, array_t);

/* 
 * =============================================================================
 *                            Quantum Tree Generator
 * =============================================================================
 */

/*
 * Function:        qtg
 * --------------------
 * Description:     This function simulates applying the quantum tree generator
 *                  to the initial state |0> |Z> |0>, where Z is the capacity of
 *                  the specified knapsack. The simulation is conducted via
 *                  breadth-first search: The decision tree is traversed layer
 *                  by layer. Following the design of the QTG, it is first
 *                  checked whether the item corresponding to the current layer
 *                  can be included or not. If not, no branching occurs.
 *                  Otherwise, the current node is expanded into two child nodes
 *                  which carry the measurement probability of the parent node,
 *                  but updated by the specified branching rule. This rule
 *                  directly corresponds to the way the Hadamard gates are
 *                  biased in the proposed design. For simulation purposes, at
 *                  each node, the local optimum of its subtree is calculated by
 *                  applying Combo to the corresponding knapsack subinstance. If
 *                  the local optimum falls below a specified threshold, the
 *                  branch is cut. However, the probability of the remaining
 *                  branch is still updated. Therefore, this function only
 *                  collects paths whose total profit lies above the specified
 *                  threshold, while still assigning them the correct
 *                  probabilities.
 * Parameters:
 *      parameter1: Pointer to knapsack whose decision tree should be traversed.
 *      parameter2: Bias towards certain branch.
 *      parameter3: Bit string representation of current solution used for biasing.
 *      parameter4: Pointer to states counter; will be updated.
 * Returns:         Array of paths, together with their sampling probabilities,
 *                  whose total profit lie above the specified threshold.  
 * Side Effect:     Allocates dynamically; pointer should eventually be freed. 
 */
node_t* qtg(const knapsack_t*, size_t, array_t, size_t*);

#ifdef __cplusplus
}
#endif

#endif
