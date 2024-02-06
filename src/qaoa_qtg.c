/*
 * =============================================================================
 *                            includes
 * =============================================================================
 */


#include "qaoa_qtg.h"

/*
 * =============================================================================
 *                            Global variables
 * =============================================================================
 */

num_t dpth;
size_t num_smpls;
size_t num_states;
metastate_probability_t* qtg_output;


/*
 * =============================================================================
 *                            Free newly defined types
 * =============================================================================
 */

void
free_metastates_probability(metastate_probability_t metastates_probability[], size_t num_metastates) {
    for (size_t idx = 0; idx < num_metastates; ++idx) {
        sw_clear(metastates_probability[idx].choice_profit.vector);
    }
    free(metastates_probability);
}

void
free_metastates_amplitude(metastate_amplitude_t metastates_amplitude[], size_t num_metastates) {
    for (size_t idx = 0; idx < num_metastates; ++idx) {
        sw_clear(metastates_amplitude[idx].choice_profit.vector);
    }
    free(metastates_amplitude);
}


/*
 * =============================================================================
 *                            Quasi-Adiabatic Evolution
 * =============================================================================
 */

void
phase_separation_unitary(metastate_amplitude_t *angle_state, double gamma) {
    for (size_t idx = 0; idx < num_states; ++idx) {
        angle_state[idx].amplitude *= cexp(-I * gamma * angle_state[idx].choice_profit.tot_profit);
    }
}


void
mixing_unitary(metastate_amplitude_t *angle_state, double beta) {
    cmplx scalar_product = 0.0;
    for (size_t idx = 0; idx < num_states; ++idx) {
        scalar_product += qtg_output[idx].probability * angle_state[idx].amplitude;
    }

    for (int idx = 0; idx < num_states; ++idx) {
        angle_state[idx].amplitude += (cexp(-I * beta) - 1.0) * scalar_product * qtg_output[idx].probability;
    }
}


metastate_amplitude_t*
quasiadiabatic_evolution(double *angles) {
    //Put the odd values of angles as gammas and the even ones as betas

    double gamma_values[dpth];
    size_t gamma_count = 0;
    double beta_values[dpth];
    size_t beta_count = 0;

    for (size_t i = 0; i < 2*dpth; ++i) {
        if ((i + 1) % 2) {  // if (i + 1) % 2 is nonzero, i + 1 is odd
            gamma_values[gamma_count++] = angles[i];
        } else {
            beta_values[beta_count++] = angles[i];
        }
    }
    metastate_amplitude_t* angle_state = malloc(num_states * sizeof(metastate_amplitude_t));
    for (size_t idx = 0; idx < num_states; ++idx) {
        angle_state[idx].choice_profit = qtg_output[idx].choice_profit;
        angle_state[idx].amplitude = qtg_output[idx].probability + 0.0 * I;
    }

    for (int j = 0; j < dpth; ++j) {
        phase_separation_unitary(angle_state, gamma_values[j]);
        mixing_unitary(angle_state, beta_values[j]);
    }

    return angle_state;
}


/*
 * =============================================================================
 *                          Measurement & Expectation Value
 * =============================================================================
 */

int
measurement(metastate_amplitude_t *angle_state) {
    // Seed the random generator and generate a random value in the range [0,1)
    srand((unsigned int)time(NULL));
    double random_number = (double)rand() / RAND_MAX;

    // Use the inverse transform sampling to draw a basis state index
    double cumulative_probability = 0.0;
    for (int idx = 0; idx < num_states; ++idx) {
        cumulative_probability += cabs(angle_state[idx].amplitude) * cabs(angle_state[idx].amplitude);

        if (cumulative_probability > random_number) {
            return idx;
        }
    }

    return -1; // Should not happen, used for error handling
}

metastate_probability_t*
sample_for_probabilities(metastate_amplitude_t *angle_state) {
    metastate_probability_t* probs_dict = malloc(num_states * sizeof(metastate_probability_t));
    for (int idx = 0; idx < num_states; ++idx) {
        probs_dict[idx].choice_profit = angle_state->choice_profit;
    }

    for (size_t sample = 0; sample < num_smpls; ++sample) {
        int measured_basis_state = measurement(angle_state);
        if (measured_basis_state != -1) {
            probs_dict[measured_basis_state].probability += 1.0 / num_smpls;
        }
    }

    return probs_dict;
}

double
expectation_value(metastate_amplitude_t *angle_state) {
    metastate_probability_t* probs_dict = sample_for_probabilities(angle_state);
    double expectation_value = 0.0;
    for (int idx = 0; idx < num_states; ++idx) {
        expectation_value += probs_dict[idx].choice_profit.tot_profit * probs_dict[idx].probability;
    }
    if (probs_dict != NULL) {
        free_metastates_probability(probs_dict, num_states);
    }
    return expectation_value;
}

double
angles_to_value(double *angles) {
    metastate_amplitude_t *angle_state = quasiadiabatic_evolution(angles);
    double exp_value = expectation_value(angle_state);
    if (angle_state != NULL) {
        free_metastates_amplitude(angle_state, num_states);
    }
    return - exp_value;
}

double objective(unsigned n, const double *angles, double *grad, void *my_func_data) {
    n = num_states;
    metastate_amplitude_t *angle_state = quasiadiabatic_evolution(angles);
    double exp_value = expectation_value(angle_state);
    if (angle_state != NULL) {
        free_metastates_amplitude(angle_state, num_states);
    }
    // grad is NULL bcs both Nelder Mead and Powell are derivative-free algorithms
    return - exp_value;
}

double * nelder_mead(){
    nlopt_opt opt = nlopt_create(NLOPT_LN_NELDERMEAD, 2 * dpth);

    // Set your optimization parameters
    nlopt_set_xtol_rel(opt, 1e-6);
    void *f_data = NULL;
    // Set the objective function
    nlopt_set_min_objective(opt, objective, f_data);

    // Set initial guess
    // Set all gammas to 0 and betas to sum up to pi/2
    double x[2 * dpth] ;
    for (size_t i = 0; i < dpth; i++)
        x[i] = (i + 1) % 2 ? 0 : M_PI / (2 * dpth);

    // Run the optimization
    nlopt_result result = nlopt_optimize(opt, x, NULL);

    // Check the result and print the optimized values
    if (result < 0) {
        printf("NLOpt failed with code %d\n", result);
    } else {
        printf("Found minimum at f(%g, %g) = %g\n", x[0], x[1], objective(num_states, x, NULL, NULL));
    }

    // Clean up
    nlopt_destroy(opt);
    return x;
}

/*
 * =============================================================================
 *                            Combination to full QAOA
 * =============================================================================
 */

qaoa_result_t
qaoa_qtg(knapsack_t* k, num_t depth, size_t bias, size_t num_samples, enum OptimizationType optimizationType) {

    path_t* int_greedy_sol;
    node_t* qtg_nodes;
    double * results;

    dpth = depth;
    num_smpls = num_samples;

    sort_knapsack(k, RATIO);
    apply_int_greedy(k);
    int_greedy_sol = path_rep(k);
    remove_all_items(k);

    qtg_nodes = qtg(k, bias, int_greedy_sol->choice_profit.vector, &num_states);

    qtg_output = malloc(num_states * sizeof(metastate_probability_t));
    for (int idx = 0; idx < num_states; ++idx) {
        qtg_output[idx].choice_profit = qtg_nodes[idx].path.choice_profit;
        qtg_output[idx].probability = qtg_nodes[idx].prob;
    }

    if (qtg_nodes != NULL) {
        free_nodes(qtg_nodes, num_states);
    }

    // TODO Optimize angles via function angles_to_value here -> Negative optimization result is final qaoa result

    switch (optimizationType) {
        case BFGS:
            //
            break;
        case NELDER_MEAD:
            results = nelder_mead();
            break;
        case POWELL:
            break;
    }
    // TODO Insert optimal angles into quasiadiabatic_evolution to obtain the final state based on optimal angles
    // TODO Call sample_for_probabilities with the final angle state to get a probability dictionary (mapping: profit <-> prob)
    // TODO -> Ask Lennart, Sören about better way instead of re-running the circuit-emulating function

    if (qtg_output != NULL) {
        free_metastates_probability(qtg_output, num_states);
    }

    // TODO Gate count

    // TODO Execute combo to transform probabilities obtained in last quasiadiabatic_evolution call to approximation ratios
    // TODO Collect and combine probability dictionary entries with the same approximation ratios
    //  -> This could also be done earlier in sample_for_probabilities (for equal profits)

    // TODO Wrap these in the shape of qaoa_result_t and return the result
}
