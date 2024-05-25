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
#define POW2(X) (1 << (X))

/*
 * =============================================================================
 *                                   Utils
 * =============================================================================
 */

num_t depth;
size_t numStates;
mixer_t mixerStrategy;
node_t *qtgNodes;

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
write_plot_data_to_file(cbs_t *angle_state, double solution_value, num_t optimal_solution_value) {
    FILE *file;

    create_dir(""); // TODO Create meaningful directories according to generated instances
    file = fopen("testfile", "w"); // TODO Create meaningful filenames according to generated instances

    fprintf(file, "%zu\n", numStates); // Save number of states for easier Python access
    fprintf(file, "%ld\n", optimal_solution_value); // Save optimal solution value for documentation
    double tot_approx_ratio =  solution_value / (double) optimal_solution_value;
    fprintf(file, "%f\n", tot_approx_ratio); // Save final approximation ratio for documentation

    for (size_t idx = 0; idx < numStates; ++idx) {
        double approx_ratio = (double) angle_state[idx].profit / (double) optimal_solution_value;
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

void apply_ry(cbs_t *angle_state, int qubit, double prob) {
    size_t blockDistance = POW2(qubit + 1);
    size_t flipDistance = POW2(qubit);
    cmplx tmp;
    for (size_t i = 0; i < numStates; i += blockDistance) {
        for (size_t j = i; j < i + flipDistance; ++j) {
            tmp = (angle_state + j)->amplitude;
            (angle_state + j)->amplitude = sqrt(1 - prob * prob) * tmp \
                                            - prob * (angle_state + j + flipDistance)->amplitude;
            (angle_state + j + flipDistance)->amplitude = prob * tmp + sqrt(1 - prob * prob) \
                                                           * (angle_state + j + flipDistance)->amplitude;
        }
    }
}

void apply_cry(cbs_t *angle_state, int control, int target, bool_t condition, double prob) {
    size_t blockDistance = POW2(target + 1);
    size_t flipDistance = POW2(target);
    cmplx tmp;
    for (size_t i = 0; i < numStates; i += blockDistance) {
        for (size_t j = i; j < i + flipDistance; ++j) {
            if ((condition && (j & POW2(control))) || (!condition && !(j & POW2(control)))) {
                tmp = (angle_state + j)->amplitude;
                (angle_state + j)->amplitude = sqrt(1 - prob * prob) * tmp \
                                        - prob * (angle_state + j + flipDistance)->amplitude;
                (angle_state + j + flipDistance)->amplitude = prob * tmp + sqrt(1 - prob * prob) \
                                                       * (angle_state + j + flipDistance)->amplitude;
            }
        }
    }
}

void apply_rz(cbs_t *angle_state, int qubit, double angle) {
    size_t blockDistance = POW2(qubit + 1);
    size_t flipDistance = POW2(qubit);
    for (size_t i = 0; i < numStates; i += blockDistance) {
        for (size_t j = i; j < i + flipDistance; ++j) {
            (angle_state + j)->amplitude *= cos(angle) - I * sin(angle);
            (angle_state + j + flipDistance)->amplitude *= cos(angle) + I * sin(angle);
        }
    }
}

void
phase_separation_unitary(cbs_t *angle_state, double gamma) {
    for (size_t idx = 0; idx < numStates; ++idx) {
        angle_state[idx].amplitude *= cexp(-I * gamma * angle_state[idx].profit);
    }
}


void
grover_mixer(cbs_t *angle_state, double beta) {
    cmplx scalar_product = 0.0;
    for (size_t idx = 0; idx < numStates; ++idx) {
        scalar_product += sqrt(qtgNodes[idx].prob) * angle_state[idx].amplitude;
    }

    for (int idx = 0; idx < numStates; ++idx) {
        angle_state[idx].amplitude += (cexp(-I * beta) - 1.0) * scalar_product * sqrt(qtgNodes[idx].prob);
    }
}

void
vandam_mixer(cbs_t *angle_state, double beta) {

}

cbs_t *
quasiadiabatic_evolution(const double *angles) {

    cbs_t *angle_state = malloc(numStates * sizeof(cbs_t));
    switch (mixerStrategy) {
        case QTG:
            for (size_t idx = 0; idx < numStates; ++idx) {
                angle_state[idx].profit = qtgNodes[idx].path.tot_profit;
                angle_state[idx].amplitude = sqrt(qtgNodes[idx].prob);
            }
            break;
        case VANDAM:
            for (size_t idx = 0; idx < numStates; ++idx) {
                angle_state[idx].profit = 0; // TODO: invoke classical objective function for idx as bit string
                angle_state[idx].amplitude = 0; // TODO: invoke initial state method for getting amplitudes
            }
    }


    for (size_t j = 0; j < depth; ++j) {
        phase_separation_unitary(angle_state,
                                 angles[2 * j]); // gamma values are even positions in angles since starting at index 0
        switch (mixerStrategy) {
            case QTG:
                grover_mixer(angle_state, angles[2 * j + 1]);
                break;
            case VANDAM:
                vandam_mixer(angle_state, angles[2 * j + 1]);
                break;
        }
    }

    return angle_state;
}


/*
 * =============================================================================
 *                            Evaluation & Optimization
 * =============================================================================
 */

double
expectation_value(cbs_t *angle_state) {
    double exp_val = 0.0;
    for (size_t idx = 0; idx < numStates; ++idx) {
        exp_val += (double) angle_state[idx].profit \
                    * cabs(angle_state[idx].amplitude) \
                    * cabs(angle_state[idx].amplitude);
    }
    return exp_val;
}


double
angles_to_value(unsigned n, const double *angles, double *grad, void *my_func_data) {
    cbs_t *angle_state = quasiadiabatic_evolution(angles);
    double exp_value = expectation_value(angle_state);
    if (angle_state != NULL) {
        free(angle_state);
    }
    // grad is NULL bcs both Nelder Mead and Powell are derivative-free algorithms
    return -exp_value;
}


nlopt_algorithm map_enum_to_nlopt_algorithm(opt_t opt_type) {
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
nlopt_optimizer(opt_t optimizationType) {
    double *angles = malloc(2 * depth * sizeof(double));

    nlopt_algorithm nlopt_optimization_algorithm = map_enum_to_nlopt_algorithm(optimizationType);
    nlopt_opt opt = nlopt_create(nlopt_optimization_algorithm, 2 * depth);

    // Set your optimization parameters
    nlopt_set_xtol_rel(opt, 1e-6);

    // Set the objective function
    nlopt_set_max_objective(opt, angles_to_value, NULL);

    // Set initial guess
    // Set all gammas to random and betas to 0
    for (size_t j = 0; j < depth; j++) {
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
        printf("Found minimum at f(%g, %g) = %g\n", angles[0], angles[1], angles_to_value(numStates, angles, NULL,
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
qaoa(knapsack_t *k, num_t inputDepth, size_t bias, mixer_t inputMixerStrategy, opt_t optimizationType) {

    double *opt_angles;

    depth = inputDepth;
    mixerStrategy = inputMixerStrategy;

    sort_knapsack(k, RATIO);

    switch (mixerStrategy) {
        case QTG:
            apply_int_greedy(k);
            path_t *int_greedy_sol = path_rep(k);
            printf("greedy = %ld\n", int_greedy_sol->tot_profit);
            remove_all_items(k);

            printf("stategen\n");
            qtgNodes = qtg(k, bias, int_greedy_sol->vector, &numStates);
            printf("states = %zu\n", numStates);

            free_path(int_greedy_sol);
            break;
        case VANDAM:
            numStates = POW2(k->size);
            break;
    }

    printf("angle opt\n");
    opt_angles = nlopt_optimizer(optimizationType);

    printf("evolution\n");
    fflush(stdout);
    cbs_t *opt_angle_state = quasiadiabatic_evolution(opt_angles);

    printf("exp value\n");
    if (qtgNodes != NULL) {
        free_nodes(qtgNodes, numStates);
    }

    double sol_val = expectation_value(opt_angle_state);
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
