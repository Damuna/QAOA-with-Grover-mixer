/* 
 * =============================================================================
 *                            includes
 * =============================================================================
 */

#if defined(_WIN32) || defined(_WIN64)
    #include "..\include\qaoa.h"
#else
    #include "../include/qaoa.h"
#endif

#include <time.h>

/* 
 * =============================================================================
 *                            Quasi-Adiabatic Evolution
 * =============================================================================
 */

void
phase_separation_unitary(metastate_amplitude_t *metastate, num_t gamma) {
    for (int idx = 0; idx < sizeof(metastate); ++idx) {
        metastate[idx].amplitude *= cexp(-I * gamma * metastate[idx].choice_profit.tot_profit);
    }
}


void 
mixing_unitary(metastate_amplitude_t *metastate, double *qtg_amplitudes, num_t beta) {
    cmplx scalar_product = 0.0;
    for (int idx = 0; idx < sizeof(qtg_amplitudes); ++idx) {
        scalar_product += qtg_amplitudes[idx] * metastate[idx].amplitude;
    }

    for (int idx = 0; idx < sizeof(metastate); ++idx) {
        metastate[idx].amplitude += (cexp(-I * beta) - 1.0) * scalar_product * qtg_amplitudes[idx];
    }
}


void 
quasi_adiabatic_evolution(metastate_amplitude_t *metastate, double *qtg_amplitudes, num_t *angles) {
    double *gamma_values = angles;
    double *beta_values = angles + 1;

    for (int j = 0; j < sizeof(angles) / 2; ++j) {
        phase_separation_unitary(metastate, gamma_values[j]);
        mixing_unitary(metastate, qtg_amplitudes, beta_values[j]);
    }
}


/* 
 * =============================================================================
 *                            Evaluation & Optimization
 * =============================================================================
 */

int measurement(metastate_amplitude_t *metastate_amplitude) {
    // Seed the random generator and generate a random value in the range [0,1)
    srand((unsigned int)time(NULL));
    double random_number = (double)rand() / RAND_MAX;

    // Use the inverse transform sampling to draw a basis state index
    double cumulative_probability = 0.0;
    for (int idx = 0; idx < sizeof(metastate_amplitude); ++idx) {
        cumulative_probability += cabs(metastate_amplitude[idx].amplitude) * cabs(metastate_amplitude[idx].amplitude);

        if (cumulative_probability > random_number) {
            return idx;
        }
    }
}

metastate_probability_t* sample_for_probabilities(metastate_amplitude_t *metastate_amplitude, num_t num_samples) {
    metastate_probability_t* metastate_probability = malloc(sizeof(metastate_amplitude) * sizeof(metastate_probability_t));
    for (int idx = 0; idx < sizeof(metastate_amplitude); ++idx) {
        metastate_probability[idx].choice_profit = metastate_amplitude->choice_profit;
    }

    for (int sample = 0; sample < num_samples; ++sample) {
        int measured_basis_state = measurement(metastate_amplitude);
        metastate_probability[measured_basis_state].probability += 1.0 / num_samples;
    }

    return metastate_probability;
}

void free_metastates_probabilities() {
    // TODO
}

double expectation_value(metastate_amplitude_t *metastate_amplitude, num_t num_samples) {
    metastate_probability_t* metastate_probability = sample_for_probabilities(metastate_amplitude, num_samples);
    double expectation_value = 0.0;
    for (int idx = 0; idx < metastate_amplitude; ++idx) {
        expectation_value += metastate_probability[idx].choice_profit.tot_profit * metastate_probability[idx].probability;
    }
    return expectation_value;
}

double angles_to_value(metastate_amplitude_t *metastate_amplitude, double *qtg_amplitudes, double *angles, num_t num_samples) {
    quasi_adiabatic_evolution(metastate_amplitude, qtg_amplitudes, angles);
    return - expectation_value(metastate_amplitude, num_samples); // TODO: The function to optimize should only have the parameters to optimize as input variables, right?
}



/* 
 * =============================================================================
 *                            Combination to full QAOA
 * =============================================================================
 */

qaoa_result_t qaoa(knapsack_t* k, int depth) {

    node_t* qtg_nodes;

    sort_knapsack(k, RATIO);

    qtg_nodes = qtg(k); // Continue after merging changes from main, i.e. threshold not existing anymore

    if (qtg_nodes != NULL) {
        free_nodes(qtg_nodes);
    }
}