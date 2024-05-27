/*
 * =============================================================================
 *                            includes
 * =============================================================================
 */

#include "copula_count.h"
#include "general_count.c"



bit_t
qubit_count_copula_mixer(const knapsack_t* k) {
    return k->size;
}

count_t
cycle_count_copula_mixer(const knapsack_t* k) {
    const count_t cycle_count_two_init = 3;
    const count_t cycle_count_two_copula = 2 * cycle_count_two_init + 1;
    if (k->size % 2 == 0) {
        return 2 * cycle_count_two_copula;
    } else { // TODO Why else redundant?
        return 3 * cycle_count_two_copula;
    }
}

count_t
gate_count_copula_mixer(const knapsack_t* k) {
    const count_t gate_count_two_init = 3;
    const count_t gate_count_two_copula = 2 * gate_count_two_init + 2;
    return k->size * gate_count_two_copula;
}

bit_t
qubit_count_copula_qaoa(const knapsack_t* k) {
    return qubit_count_qaoa(k, qubit_count_copula_mixer(k));
}

count_t
cycle_count_copula_qaoa(const knapsack_t* k, const num_t depth) {
    return cycle_count_qaoa(depth, cycle_count_copula_mixer(k));
}

count_t
gate_count_copula_qaoa(const knapsack_t* k, const num_t depth) {
    return gate_count_qaoa(k, depth, gate_count_copula_mixer(k));
}