/*
 * =============================================================================
 *                            includes
 * =============================================================================
 */

#include "general_count.h"
#include "knapsack.c"


bit_t qubit_count_phase_separator(const knapsack_t* k) {
    return k->size;
}

count_t
cycle_count_phase_separator() {
    return 1;
}

count_t
gate_count_phase_separator(const knapsack_t* k) {
    return k->size;
}

bit_t
qubit_count_qaoa(const knapsack_t* k, const bit_t qubit_count_mixer) {
    return MAX(qubit_count_phase_separator(k), qubit_count_mixer);
}

count_t
cycle_count_qaoa(const num_t depth, const count_t cycle_count_mixer) {
    return depth * (cycle_count_phase_separator() + cycle_count_mixer);
}

count_t
gate_count_qaoa(const knapsack_t* k, const num_t depth, const count_t gate_count_mixer) {
    return depth * (gate_count_phase_separator(k) + gate_count_mixer);
}
