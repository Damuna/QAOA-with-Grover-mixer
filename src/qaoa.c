/*
 * =============================================================================
 *                            includes
 * =============================================================================
 */

#include "qaoa.h"


/*
 * =============================================================================
 *                                  Macros
 * =============================================================================
 */

#define TRUE        1
#define FALSE       0

#define POW2(X) (1 << (X))


/*
 * =============================================================================
 *                              Global variables
 * =============================================================================
 */

// Variables to be set from the outset
knapsack_t* kp;
qaoa_type_t qaoa_type;
size_t bias;
num_t depth;
double k;
double theta;

// Variables that are initialized later
size_t num_states;
node_t* qtg_nodes;
num_t* sol_profits;
double* prob_dist_vals;
bool_t* sol_feasibilities;


/*
 * =============================================================================
 *                                  Utils
 * =============================================================================
 */

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


/*
 * =============================================================================
 *                              Gate application
 * =============================================================================
 */

void
apply_ry(cbs_t* angle_state, const int qubit, const double prob) {
    const size_t blockDistance = POW2(qubit + 1);
    const size_t flipDistance = POW2(qubit);
    for (size_t i = 0; i < num_states; i += blockDistance) {
        for (size_t j = i; j < i + flipDistance; ++j) {
            const cmplx tmp = (angle_state + j)->amplitude;
            (angle_state + j)->amplitude = sqrt(1 - prob) * tmp \
                                            - sqrt(prob) * (angle_state + j + flipDistance)->amplitude;
            (angle_state + j + flipDistance)->amplitude = sqrt(prob) * tmp + sqrt(1 - prob) \
                                                           * (angle_state + j + flipDistance)->amplitude;
        }
    }
}


void
apply_ry_inv(cbs_t* angle_state, const int qubit, const double prob) {
    const size_t blockDistance = POW2(qubit + 1);
    const size_t flipDistance = POW2(qubit);
    for (size_t i = 0; i < num_states; i += blockDistance) {
        for (size_t j = i; j < i + flipDistance; ++j) {
            const cmplx tmp = (angle_state + j)->amplitude;
            (angle_state + j)->amplitude = sqrt(1 - prob) * tmp \
                                            + sqrt(prob) * (angle_state + j + flipDistance)->amplitude;
            (angle_state + j + flipDistance)->amplitude = - sqrt(prob) * tmp + sqrt(1 - prob) \
                                                           * (angle_state + j + flipDistance)->amplitude;
        }
    }
}


void
apply_cry(cbs_t* angle_state, const int control, const int target, const bool_t condition, const double prob) {
    const size_t blockDistance = POW2(target + 1);
    const size_t flipDistance = POW2(target);
    for (size_t i = 0; i < num_states; i += blockDistance) {
        for (size_t j = i; j < i + flipDistance; ++j) {
            if ((condition && (j & POW2(control))) || (!condition && !(j & POW2(control)))) {
                const cmplx tmp = (angle_state + j)->amplitude;
                (angle_state + j)->amplitude = sqrt(1 - prob) * tmp \
                                        - sqrt(prob) * (angle_state + j + flipDistance)->amplitude;
                (angle_state + j + flipDistance)->amplitude = sqrt(prob) * tmp + sqrt(1 - prob) \
                                                       * (angle_state + j + flipDistance)->amplitude;
            }
        }
    }
}


void
apply_cry_inv(cbs_t* angle_state, const int control, const int target, const bool_t condition, const double prob) {
    const size_t blockDistance = POW2(target + 1);
    const size_t flipDistance = POW2(target);
    for (size_t i = 0; i < num_states; i += blockDistance) {
        for (size_t j = i; j < i + flipDistance; ++j) {
            if ((condition && (j & POW2(control))) || (!condition && !(j & POW2(control)))) {
                const cmplx tmp = (angle_state + j)->amplitude;
                (angle_state + j)->amplitude = sqrt(1 - prob) * tmp \
                                        + sqrt(prob) * (angle_state + j + flipDistance)->amplitude;
                (angle_state + j + flipDistance)->amplitude = - sqrt(prob) * tmp + sqrt(1 - prob) \
                                                       * (angle_state + j + flipDistance)->amplitude;
            }
        }
    }
}


void
apply_rz(cbs_t* angle_state, const int qubit, const double angle) {
    const size_t blockDistance = POW2(qubit + 1);
    const size_t flipDistance = POW2(qubit);
    for (size_t i = 0; i < num_states; i += blockDistance) {
        for (size_t j = i; j < i + flipDistance; ++j) {
            (angle_state + j)->amplitude *= cos(angle) - I * sin(angle);
            (angle_state + j + flipDistance)->amplitude *= cos(angle) + I * sin(angle);
        }
    }
}


/*
 * =============================================================================
 *                                 QTG-specific
 * =============================================================================
 */

void
qtg_initial_state_prep(cbs_t* angle_state) {
    for (size_t idx = 0; idx < num_states; ++idx) {
        angle_state[idx].profit = qtg_nodes[idx].path.tot_profit;
        angle_state[idx].amplitude = sqrt(qtg_nodes[idx].prob);
    }
}


void
qtg_grover_mixer(cbs_t *angle_state, double beta) {
    cmplx scalar_product = 0.0;
    for (size_t idx = 0; idx < num_states; ++idx) {
        scalar_product += sqrt(qtg_nodes[idx].prob) * angle_state[idx].amplitude;
    }

    for (int idx = 0; idx < num_states; ++idx) {
        angle_state[idx].amplitude += (cexp(-I * beta) - 1.0) * scalar_product * sqrt(qtg_nodes[idx].prob);
    }
}


/*
 * =============================================================================
 *                                 Copula-specific
 * =============================================================================
 */

double
prob_dist(const bit_t index) {
    const double c = cost_sum(kp) / kp->capacity - 1;
    const double r_st = stop_item_ratio(kp);
    const double r = kp->items[index].profit / kp->items[index].cost;
    return 1 / (1 + c * exp(-k * (r - r_st)));
}


void
copula_initial_state_prep(cbs_t* angle_state) {
    for (bit_t bit = 0; bit < kp->size; ++bit) {
        angle_state->profit = sol_profits[bit];

        const double d = prob_dist_vals[bit];
        apply_ry(angle_state, bit, d);
    }
}


void
apply_r_dist(
    cbs_t* angle_state,
    const num_t qubit1,
    const num_t qubit2,
    const double d1,
    const double d2given1,
    const double d2givennot1
) {
    apply_ry(angle_state, qubit1, d1);
    apply_cry(angle_state, qubit1, qubit2, 1, d2given1);
    apply_cry(angle_state, qubit1, qubit2, 0, d2givennot1);
}


void
apply_r_dist_inv(
    cbs_t* angle_state,
    const num_t qubit1,
    const num_t qubit2,
    const double d1,
    const double d2given1,
    const double d2givennot1
) {
    apply_cry_inv(angle_state, qubit1, qubit2, 0, d2givennot1);
    apply_cry_inv(angle_state, qubit1, qubit2, 1, d2given1);
    apply_ry_inv(angle_state, qubit1, d1);
}


void
apply_two_copula(cbs_t* angle_state, const int qubit1, const int qubit2, const double beta) {
    const double d1 = prob_dist_vals[qubit1];
    const double d2 = prob_dist_vals[qubit2];

    const double d2given1 = d2 + theta * d2 * (1 - d1) * (1 - d2);
    const double d2givennot1 = d2 - theta * d1 * d2 * (1 - d2);

    apply_r_dist_inv(angle_state, qubit1, qubit2, d1, d2given1, d2givennot1);
    apply_rz(angle_state, qubit1, 2 * beta);
    apply_rz(angle_state, qubit2, 2 * beta);
    apply_r_dist(angle_state, qubit1, qubit2, d1, d2given1, d2givennot1);
}


void
copula_mixer(cbs_t* angle_state, double beta) {
    const bool_t kp_size_even = {kp->size % 2 == 0};
    for (bit_t qubit = 0; qubit < kp->size - 1; ++qubit) {
        const bit_t odd_qubit = 2 * qubit + 1;
        apply_two_copula(angle_state, odd_qubit, odd_qubit + 1, beta);
    }
    if (!kp_size_even) {
        apply_two_copula(angle_state, kp->size, 1, beta);
    }
    for (bit_t qubit = 0; qubit < kp->size - 1; ++qubit) {
        const bit_t even_qubit = 2 * qubit;
        apply_two_copula(angle_state, even_qubit, even_qubit + 1, beta);
    }
    if (kp_size_even) {
        apply_two_copula(angle_state, kp->size, 1, beta);
    }
}


/*
 * =============================================================================
 *                            Quasi-Adiabatic Evolution
 * =============================================================================
 */

void
phase_separation_unitary(cbs_t *angle_state, double gamma) {
    for (size_t idx = 0; idx < num_states; ++idx) {
        angle_state[idx].amplitude *= cexp(-I * gamma * angle_state[idx].profit);
    }
}


cbs_t *
quasiadiabatic_evolution(const double *angles) {
    void (*initial_state_prep)(cbs_t*);
    void (*mixing_unitary)(cbs_t*, double);
    switch (qaoa_type) {
        case QTG:
            initial_state_prep = qtg_initial_state_prep;
            mixing_unitary = qtg_grover_mixer;
            break;
        case COPULA:
            initial_state_prep = copula_initial_state_prep;
            mixing_unitary = copula_mixer;
            break;
    }

    cbs_t* angle_state = malloc(num_states * sizeof(cbs_t));
    initial_state_prep(angle_state);

    for (int j = 0; j < depth; ++j) {
        // gamma values are even positions in angles since starting at index 0
        phase_separation_unitary(angle_state, angles[2 * j]);
        // beta values are odd positions in angles since starting at index 0
        mixing_unitary(angle_state, angles[2 * j + 1]);

    }
    return angle_state;
}


/*
 * =============================================================================
 *                                  Evaluation
 * =============================================================================
 */

double
expectation_value(const cbs_t* angle_state) {
    double exp_val = 0;

    for (size_t idx = 0; idx < num_states; ++idx) {
        if (qaoa_type == COPULA) {
            if (!sol_feasibilities[idx])
                continue; // Add 0 in case that solution is infeasible (modified objective function)
        }

        const double prob = cabs(angle_state[idx].amplitude) * cabs(angle_state[idx].amplitude);
        exp_val += prob * angle_state[idx].profit;
    }
    return exp_val;
}


double
angles_to_value(const double* angles) {
    cbs_t *angle_state = quasiadiabatic_evolution(angles);

    const double exp_value = expectation_value(angle_state);
    if (angle_state != NULL) {
        free(angle_state);
    }
    return -exp_value;
}


/*
 * =============================================================================
 *                                Optimization
 * =============================================================================
 */

double
angles_to_value_nlopt(unsigned n, const double *angles, double *grad, void *my_func_data) {
    // grad is NULL bcs both Nelder Mead and Powell are derivative-free algorithms
    return angles_to_value(angles);
}


nlopt_algorithm
map_enum_to_nlopt_algorithm(const opt_t opt_type) {
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


// Helper function to perform a fine grid search
void fine_grid_search(const int m, double* best_angles, double* best_value) {
    const int total_angles = 2 * depth;
    const double step_size = 2 * M_PI / m;
    double angles[total_angles];
    *best_value = -INFINITY;

    for (int i = 0; i < (int)pow(m, total_angles); ++i) {
        int index = i;
        for (int j = 0; j < total_angles; ++j) {
            angles[j] = (index % m) * step_size;
            index /= m;
        }
        const double value = angles_to_value(angles);
        if (value > *best_value) {
            *best_value = value;
            for (int k = 0; k < total_angles; ++k) {
                best_angles[k] = angles[k];
            }
        }
    }
}



double*
nlopt_optimizer(const opt_t optimization_type, const int m) {
    double *angles = malloc(2 * depth * sizeof(double));
    double best_value;

    nlopt_algorithm nlopt_optimization_algorithm = map_enum_to_nlopt_algorithm(optimization_type);
    nlopt_opt opt = nlopt_create(nlopt_optimization_algorithm, 2 * depth);

    // Set your optimization parameters
    nlopt_set_xtol_rel(opt, 1e-6);

    // Set the objective function
    nlopt_set_min_objective(opt, angles_to_value_nlopt, NULL);

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
    printf("Performing fine-grid search...\n");
    // Perform fine grid search before optimizing
    fine_grid_search(m, angles, &best_value);

    printf("Starting NLOpt...\n");
    // Run the optimization
    // opt_f must not be NULL
    double obj = 0;
    const nlopt_result result = nlopt_optimize(opt, angles, &obj);

    // Check the result and print the optimized values
    if (result < 0) {
        printf("NLOpt failed with code %d\n", result);
    } else {
        char* string_to_adjust = "value";
        if (depth > 1) {
            string_to_adjust = strcat(string_to_adjust, "s");
        }
        printf("Done! NLOPT returned optimal gamma ");
        printf("%s (", string_to_adjust);
        for (size_t j = 0; j < depth; j++) {
            printf("%g", angles[2 * j]);
        }
        printf(") and optimal beta ");
        printf("%s (", string_to_adjust);
        for (size_t j = 0; j < depth; j++) {
            printf("%g", angles[2 * j + 1]);
        }
        printf(") with objective function value %g\n", angles_to_value(angles));
    }

    // Clean up
    nlopt_destroy(opt);
    return angles;
}


/*
 * =============================================================================
 *                                  Export data
 * =============================================================================
 */

char*
path_for_instance(const char* instance) {
    char qaoa_type_str[16];
    switch (qaoa_type) {
        case QTG:
            strcpy(qaoa_type_str, "qtg");
            break;
        case COPULA:
            strcpy(qaoa_type_str, "copula");
            break;
    }
    char *path = calloc(1024, sizeof(char));
    sprintf(path, "../instances/%s/%s/", instance, qaoa_type_str);
    return path;
}


void
export_results(
    const char* instance, const cbs_t* angle_state, const double solution_value, const num_t optimal_solution_val
) {
    FILE* file = fopen(strcat(path_for_instance(instance), "results"), "w");

    fprintf(file, "%llu\n", num_states); // Save number of states for easier Python access
    fprintf(file, "%ld\n", optimal_solution_val); // Save optimal solution value for documentation
    const double tot_approx_ratio = solution_value / optimal_solution_val;
    fprintf(file, "%f\n", tot_approx_ratio); // Save final approximation ratio for documentation

    for (size_t idx = 0; idx < num_states; ++idx) {
        double const approx_ratio = (double) angle_state[idx].profit / optimal_solution_val;
        double const prob = cabs(angle_state[idx].amplitude) * cabs(angle_state[idx].amplitude);
        fprintf(file, "%f %f\n", approx_ratio, prob);
    }

    fclose(file);
}

void
export_resources(const char* instance, const resource_t res) {
    FILE* file = fopen(strcat(path_for_instance(instance), "resources"), "w");

    fprintf(file, "%d\n", res.qubit_count);
    fprintf(file, "%u\n", res.cycle_count);
    fprintf(file, "%u\n", res.gate_count);
    fprintf(file, "%u\n", res.gate_count_decomp);

    fclose(file);
}


/*
 * =============================================================================
 *                            Combination to full QAOA
 * =============================================================================
 */

double
qaoa(
    const char* instance,
    knapsack_t* input_kp,
    const qaoa_type_t input_qaoa_type,
    const num_t input_depth,
    const opt_t opt_type,
    const num_t m,
    const size_t input_bias,
    const double copula_k,
    const double copula_theta
) {
    kp = input_kp;
    qaoa_type = input_qaoa_type;
    depth = input_depth;
    bias = input_bias;
    k = copula_k;
    theta = copula_theta;

    sort_knapsack(kp, RATIO);

    // Initialize algorithms
    printf("\n===== Preparation ======\n");
    switch (qaoa_type) {
        case QTG:
            apply_int_greedy(kp);
            path_t *int_greedy_sol = path_rep(kp);
            printf("Integer greedy solution = %ld\n", int_greedy_sol->tot_profit);
            remove_all_items(kp);

            printf("Generating states via QTG...\n");
            qtg_nodes = qtg(kp, bias, int_greedy_sol->vector, &num_states);
            printf("Done! Number of states = %zu\n", num_states);

            free_path(int_greedy_sol);
            break;

        case COPULA:
            num_states = POW2(kp->size);

            prob_dist_vals = malloc(kp->size * sizeof(double));
            for (bit_t bit = 0; bit < kp->size; ++bit) {
                prob_dist_vals[bit] = prob_dist(bit);
            }

            sol_profits = malloc(num_states * sizeof(num_t));
            for (size_t idx = 0; idx < num_states; ++idx) {
                sol_profits[idx] = objective_func(kp, idx);
            }

            sol_feasibilities = malloc(num_states * sizeof(bool_t));
            for (size_t idx = 0; idx < num_states; ++idx) {
                sol_feasibilities[idx] = sol_cost(kp, idx) <= kp->capacity;
            }
            break;
    }

    // Optimize
    printf("\n===== Optimization =====\n");
    printf("Optimize angles...\n");
    const double* opt_angles = nlopt_optimizer(opt_type, m);

    printf("\n===== Final phase of evaluation =====\n");

    printf("Quasi-adiabatic evolution of optimal angles...\n");
    fflush(stdout);
    cbs_t *opt_angle_state = quasiadiabatic_evolution(opt_angles);

    printf("Compute final expectation value...\n");
    if (qtg_nodes != NULL) {
        free_nodes(qtg_nodes, num_states);
    }

    const double sol_val = -expectation_value(opt_angle_state);

    // Export results data
    const num_t optimal_sol_val = combo_wrap(kp, 0, kp->capacity, FALSE, FALSE, TRUE, FALSE);
    printf("Optimal solution value (COMBO) = %ld\n", optimal_sol_val);
    export_results(instance, opt_angle_state, sol_val, optimal_sol_val);

    // Free optimal-angles solution
    if (opt_angle_state != NULL) {
        free(opt_angle_state);
    }

    // Export resources data
    resource_t res;
    switch (qaoa_type) {
        case QTG:
            res.qubit_count = qubit_count_qtg_qaoa(kp);
            res.cycle_count = cycle_count_qtg_qaoa(kp, depth, COPPERSMITH, TOFFOLI, FALSE);
            res.gate_count = gate_count_qtg_qaoa(kp, depth, COPPERSMITH, TOFFOLI, FALSE);
            res.cycle_count_decomp = cycle_count_qtg_qaoa(kp, depth, COPPERSMITH, TOFFOLI, TRUE);
            res.gate_count_decomp = gate_count_qtg_qaoa(kp, depth, COPPERSMITH, TOFFOLI, TRUE);

        case COPULA:
            res.qubit_count = qubit_count_copula_qaoa(kp);
            res.cycle_count = cycle_count_copula_qaoa(kp, depth);
            res.gate_count = gate_count_copula_qaoa(kp, depth);
            res.cycle_count_decomp = res.cycle_count;
            res.gate_count_decomp = res.gate_count;
    }
    export_resources(instance, res);

    return sol_val;
}
