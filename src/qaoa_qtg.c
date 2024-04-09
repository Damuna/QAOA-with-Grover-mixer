/*
 * =============================================================================
 *                                  includes
 * =============================================================================
 */

#include "qaoa_qtg.h"


/*
 * =============================================================================
 *                         Macros & global variables
 * =============================================================================
 */

#define TRUE        1
#define FALSE       0


/*
 * =============================================================================
 *                                   Utils
 * =============================================================================
 */

num_t dpth;
size_t num_smpls;
size_t num_states;
node_t *qtg_nodes;

double
random_value_on_windows_or_linux() {
#if defined(_WIN32) || defined(_WIN64)
    srand(time(NULL));
    return rand() / RAND_MAX;
#else
    srand48(time(NULL));
    return drand48();
#endif
}


void
write_plot_data_to_file(feassol_ampl_t *angle_state, double solution_value, num_t optimal_solution_value) {
    FILE *file;

    create_dir(""); // TODO Create meaningful directories according to generated instances
    file = fopen("testfile", "w"); // TODO Create meaningful filenames according to generated instances

    fprintf(file, "%llu\n", num_states); // Save number of states for easier Python access
    fprintf(file, "%ld\n", optimal_solution_value); // Save optimal solution value for documentation
    double tot_approx_ratio = solution_value / optimal_solution_value;
    fprintf(file, "%f\n", tot_approx_ratio); // Save final approximation ratio for documentation

    for (size_t idx = 0; idx < num_states; ++idx) {
        double approx_ratio = (double) angle_state[idx].profit / optimal_solution_value;
        double prob = cabs(angle_state[idx].amplitude) * cabs(angle_state[idx].amplitude);
        fprintf(file, "%f %f\n", approx_ratio, prob);
    }

    fclose(file);
}


/*
 * =============================================================================
 *                            Quasi-Adiabatic Evolution
 * =============================================================================
 */

void
phase_separation_unitary(feassol_ampl_t *angle_state, double gamma) {
    for (size_t idx = 0; idx < num_states; ++idx) {
        angle_state[idx].amplitude *= cexp(-I * gamma * angle_state[idx].profit);
    }
}


void
mixing_unitary(feassol_ampl_t *angle_state, double beta) {
    cmplx scalar_product = 0.0;
    for (size_t idx = 0; idx < num_states; ++idx) {
        scalar_product += sqrt(qtg_nodes[idx].prob) * angle_state[idx].amplitude;
    }

    for (int idx = 0; idx < num_states; ++idx) {
        angle_state[idx].amplitude += (cexp(-I * beta) - 1.0) * scalar_product * sqrt(qtg_nodes[idx].prob);
    }
}


feassol_ampl_t *
quasiadiabatic_evolution(const double *angles) {

    feassol_ampl_t *angle_state = malloc(num_states * sizeof(feassol_ampl_t));
    for (size_t idx = 0; idx < num_states; ++idx) {
        angle_state[idx].profit = qtg_nodes[idx].path.tot_profit;
        angle_state[idx].amplitude = sqrt(qtg_nodes[idx].prob);
    }

    for (int j = 0; j < dpth; ++j) {
        phase_separation_unitary(angle_state,
                                 angles[2 * j]); // gamma values are even positions in angles since starting at index 0
        mixing_unitary(angle_state,
                       angles[2 * j + 1]); // beta values are odd positions in angles since starting at index 0
    }

    return angle_state;
}


/*
 * =============================================================================
 *                            Evaluation & Optimization
 * =============================================================================
 */

double
expectation_value(feassol_ampl_t *angle_state) {
    double exp_val = 0.0;

    for (size_t sample = 0; sample < num_smpls; ++sample) {
        double rand_val = random_value_on_windows_or_linux();

        // Draw from the probability distribution (= measure) via the inverse transform sampling
        double cum_prob = 0.0;
        for (size_t idx = 0; idx < num_states; ++idx) {
            cum_prob += cabs(angle_state[idx].amplitude) * cabs(angle_state[idx].amplitude);

            if (cum_prob > rand_val) {
                // Update the expectation value by the relative objective function value (profit) of the measured state
                exp_val += (double) angle_state[idx].profit / num_smpls;
                break;
            }
        }
    }

    return exp_val;
}


double
angles_to_value(unsigned n, const double *angles, double *grad, void *my_func_data) {
    feassol_ampl_t *angle_state = quasiadiabatic_evolution(angles);
    double exp_value = expectation_value(angle_state);
    if (angle_state != NULL) {
        free(angle_state);
    }
    // grad is NULL bcs both Nelder Mead and Powell are derivative-free algorithms
    return -exp_value;
}


nlopt_algorithm map_enum_to_nlopt_algorithm(optimization_type_t opt_type) {
    switch (opt_type) {
        case BFGS:
            return 0;
        case NELDER_MEAD:
            return NLOPT_LN_NELDERMEAD;
        case POWELL:
            return NLOPT_LN_BOBYQA; // or another suitable algorithm for POWELL
        default:
            // Handle invalid enum value
            return NLOPT_LN_BOBYQA; // Default to BOBYQA or any other suitable default
    }
}


double *
nlopt_optimizer(optimization_type_t optimization_type) {
    double *angles = malloc(2 * dpth);

    nlopt_algorithm nlopt_optimization_algorithm = map_enum_to_nlopt_algorithm(optimization_type);
    nlopt_opt opt = nlopt_create(nlopt_optimization_algorithm, 2 * dpth);

    // Set your optimization parameters
    nlopt_set_xtol_rel(opt, 1e-6);

    // Set the objective function
    nlopt_set_min_objective(opt, angles_to_value, NULL);

    // Set initial guess
    // Set all gammas to random and betas to 0
    for (size_t j = 0; j < dpth; j++) {
        if ((j + 1) % 2) {
            angles[j] = 0;
//            angles[j] = random_value_on_windows_or_linux() * 2.0 * M_PI;
        } else {
            angles[j] = 0;
//            angles[j] = random_value_on_windows_or_linux() * 2.0 * M_PI;
        }
    }
    // Run the optimization
    // opt_f must not be NULL
    double obj = 0;
    nlopt_result result = nlopt_optimize(opt, angles, &obj);

    // Check the result and print the optimized values
    if (result < 0) {
        printf("NLOpt failed with code %d\n", result);
    } else {
        printf("Found minimum at f(%g, %g) = %g\n", angles[0], angles[1], angles_to_value(num_states, angles, NULL,
                                                                                          NULL)); // TODO There will be more than 2 angles, so printing x[0] and x[1] is meaningless
    }

    // Clean up
    nlopt_destroy(opt);
    return angles;
}


/*
 * =============================================================================
 *                            Combination to full QAOA
 * =============================================================================
 */

double
qaoa_qtg(knapsack_t *k, num_t depth, size_t bias, size_t num_samples, optimization_type_t optimization_type) {

    double *opt_angles;

    dpth = depth;
    num_smpls = num_samples;

    sort_knapsack(k, RATIO);
    apply_int_greedy(k);
    path_t *int_greedy_sol = path_rep(k);
    printf("greedy = %ld\n", int_greedy_sol->tot_profit);
    remove_all_items(k);

    printf("stategen\n");
    qtg_nodes = qtg(k, bias, int_greedy_sol->vector, &num_states);
    printf("states = %zu\n", num_states);

    free_path(int_greedy_sol);

    printf("angle opt\n");
    opt_angles = nlopt_optimizer(optimization_type);

    printf("evolution\n");
    fflush(stdout);
    feassol_ampl_t *opt_angle_state = quasiadiabatic_evolution(opt_angles);

    printf("exp value\n");
    if (qtg_nodes != NULL) {
        free_nodes(qtg_nodes, num_states);
    }

    double sol_val = -expectation_value(opt_angle_state);
    num_t optimal_sol_val = combo_wrap(k, 0, k->capacity, FALSE, FALSE, TRUE, FALSE);
    printf("optimal = %ld\n", optimal_sol_val);
    write_plot_data_to_file(opt_angle_state, sol_val, optimal_sol_val);

    if (opt_angle_state != NULL) {
        free(opt_angle_state);
    }

    resource_t res = {
            .qubit_count = qubit_count_qtg_mixer(k),
            .cycle_count = cycle_count_qtg_mixer(k, COPPERSMITH, TOFFOLI, FALSE),
            .gate_count = gate_count_qtg_mixer(k, COPPERSMITH, TOFFOLI, FALSE),
            .cycle_count_decomp = cycle_count_qtg_mixer(k, COPPERSMITH, TOFFOLI, TRUE),
            .gate_count_decomp = gate_count_qtg_mixer(k, COPPERSMITH, TOFFOLI, TRUE)
    };

    return sol_val;
}
