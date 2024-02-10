/* 
 * =============================================================================
 *                            includes
 * =============================================================================
 */

#include "stategen.h"

/* 
 * =============================================================================
 *                            macros
 * =============================================================================
 */

#define TRUE                    1
#define FALSE                   0

#define SWAP(a, b, T)           do { register T q; q = *(a); *(a) = *(b); \
                                *(b) = q; } while(0)


/* 
 * =============================================================================
 *                            free nodes
 * =============================================================================
 */

 void
 free_nodes(node_t nodes[], size_t num_nodes) {
    for (size_t i = 0; i < num_nodes; ++i) {
        sw_clear(nodes[i].path.vector);
    }
    free(nodes);
 }


/* 
 * =============================================================================
 *                            branch probability
 * =============================================================================
 */

double
//remove single and default
branch_prob(const knapsack_t* k, bit_t i, size_t bias, bool_t left, \
            array_t cur_sol) {

    if (left) {
        return (1. + (1 - sw_tstbit(cur_sol, i)) * bias) \
                        / (bias + 2);
    } else {
        return (1. + sw_tstbit(cur_sol, i) * bias) \
                        / (bias + 2);
        }

}

/* 
 * =============================================================================
 *                            Quantum Tree Generator
 * =============================================================================
 */

node_t*
qtg(const knapsack_t* k, size_t bias, \
    array_t cur_sol, size_t* num_states) {
    
	*num_states = 1; /* start from the root */
    size_t a = 0; /* start from leftmost node */
    /* initialize root node as single-element node_t array */
    node_t* parent = malloc(sizeof(node_t));
    parent->path.remain_cost = k->capacity;
    parent->path.tot_profit = 0;
    sw_init(parent->path.vector, k->size);
    parent->prob = 1.;
    
    for (bit_t i = 0; i < k->size; a = 0, ++i) { /* start from leftmost node */
        /*
         * The size of the child layer is upper bounded by twice the parent
         * layer's size which is given as the current number of states.
         */
        //a is counting the states in the child (current) layer (we have to figure out a)
        //j is counting the states in the parent layer
        //KEEP
        node_t* child = malloc(2 * (*num_states) * sizeof(node_t));
        for (size_t j = 0; j < (*num_states); ++j) {
            if (parent[j].path.remain_cost < k->items[i].cost) {
                /* item cannot be included, thus no branching */
                child[a].path.remain_cost = parent[j].path.remain_cost;
                child[a].path.tot_profit = parent[j].path.tot_profit;
                sw_init(child[a].path.vector, k->size);
                sw_set(child[a].path.vector, parent[j].path.vector);
                child[a].prob = parent[j].prob;
                // printf("----------------\n");
                // printf("Node info:\n");
                // printf("Remaining cost: %ld\n", child[a].path.remain_cost);
                // printf("Total profit: %ld\n", child[a].path.tot_profit);
                // gmp_printf("Vector: %Zd\n", child[a].path.vector);
                // printf("Probability: %.16f\n", child[a].prob);
                // printf("----------------\n");
                ++a;
                continue;
            }
            // The item can be included. 
            
            //LEFT BRANCH
            /* remaining cost, total profit, and vector do not change */
            child[a].path.remain_cost = parent[j].path.remain_cost;
            child[a].path.tot_profit = parent[j].path.tot_profit;
            sw_init(child[a].path.vector, k->size);
            sw_set(child[a].path.vector, parent[j].path.vector);
            /* update probability, then increase child index */
            child[a].prob = parent[j].prob * branch_prob(k, i, bias, \
                                                TRUE, cur_sol);
            // printf("----------------\n");
            // printf("Node info:\n");
            // printf("Remaining cost: %ld\n", child[a].path.remain_cost);
            // printf("Total profit: %ld\n", child[a].path.tot_profit);
            // gmp_printf("Vector: %Zd\n", child[a].path.vector);
            // printf("Probability: %.16f\n", child[a].prob);
            // printf("----------------\n");
            ++a;

            //RIGHT BRANCH
            /* update remaining cost */
            child[a].path.remain_cost = parent[j].path.remain_cost \
                                        - k->items[i].cost;
            /* update total profit */
            child[a].path.tot_profit = k->items[i].profit \
                                       + parent[j].path.tot_profit;
            /* include item: set the corresponding bit to 1 */
            sw_init(child[a].path.vector, k->size);
            sw_set(child[a].path.vector, parent[j].path.vector);
            sw_setbit(child[a].path.vector, i);
            
            /* update probability, then increase child index */
            child[a].prob = parent[j].prob * branch_prob(k, i, \
                                      bias, FALSE, cur_sol);
            // printf("----------------\n");
            // printf("Node info:\n");
            // printf("Remaining cost: %ld\n", child[a].path.remain_cost);
            // printf("Total profit: %ld\n", child[a].path.tot_profit);
            // gmp_printf("Vector: %Zd\n", child[a].path.vector);
            // printf("Probability: %.16f\n", child[a].prob);
            // printf("----------------\n");
            ++a;

        }
        /* swap pointer to parent and child layer */
        SWAP(&parent, &child, node_t*);
        /* resize new parent layer and free former parent layer */
        parent = realloc(parent, a * sizeof(node_t));
        free_nodes(child, *num_states);
        *num_states = a;
        // printf("---------------------------------\n");
        // printf("DONE WITH LAYER\n");
        // printf("---------------------------------\n");
    }
    /* final layer comprises all feasible paths above threshold */
    // printf("Number of states after QTG: %zu\n", *num_states);
    return parent;
}
