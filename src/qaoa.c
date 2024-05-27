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

knapsack_t* kp;
qaoa_type_t qaoa_type;
size_t bias;
num_t depth;
double k;
size_t num_states;
node_t* qtg_nodes;


/*
 * =============================================================================
 *                              Utils
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


int
modified_objective_func(const int solution) {
    int tot_profit = 0;
    int tot_cost = 0;
    for (bit_t idx = 0; idx < kp->size; ++idx) {
        if (solution & (1 << idx) == 1) {
            tot_profit += kp->items[idx].profit;
            tot_cost += kp->items[idx].cost;
            if (tot_cost > kp->capacity) {
                return 0;
            }
        }
    }
    return tot_profit;
}


void
copula_initial_state_prep(cbs_t* angle_state) {
    for (size_t idx = 0; idx < num_states; ++idx) {
        angle_state->profit = modified_objective_func(idx);

        const double d = prob_dist(idx);
        bit_t num_set_bits = num_set_bits(idx);
        angle_state->amplitude = pow(sqrt(d), num_set_bits) * pow(sqrt(1-d), kp->size - num_set_bits);
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
        case COPULA:
            initial_state_prep = copula_initial_state_prep;
            // TODO assign Copula mixer once ready
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


double*
nlopt_optimizer(const opt_t optimization_type) {
    double *angles = malloc(2 * depth * sizeof(double));

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
    // Run the optimization
    // opt_f must not be NULL
    double obj = 0;
    nlopt_result result = nlopt_optimize(opt, angles, &obj);

    // Check the result and print the optimized values
    if (result < 0) {
        printf("NLOpt failed with code %d\n", result);
    } else {
        printf("Found minimum at f(%g, %g) = %g\n", angles[0], angles[1], angles_to_value_nlopt(num_states, angles, NULL,
                                                                                          NULL)); // TODO There will be more than 2 angles, so printing x[0] and x[1] is meaningless
    }

    // Clean up
    nlopt_destroy(opt);
    return angles;
}


/*
 * =============================================================================
 *                                  Export results
 * =============================================================================
 */

void
write_plot_data_to_file(
    const cbs_t* angle_state, const double solution_value, const num_t optimal_solution_value
) {
    create_dir(""); // TODO Create meaningful directories according to generated instances
    FILE* file = fopen("testfile", "w"); // TODO Create meaningful filenames according to generated instances

    fprintf(file, "%llu\n", num_states); // Save number of states for easier Python access
    fprintf(file, "%ld\n", optimal_solution_value); // Save optimal solution value for documentation
    const double tot_approx_ratio = solution_value / optimal_solution_value;
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
 *                            Combination to full QAOA
 * =============================================================================
 */

double
qaoa(
    knapsack_t* input_kp,
    const qaoa_type_t input_qaoa_type,
    const num_t input_depth,
    const opt_t opt_type,
    const size_t input_bias,
    const double copula_k
) {
    kp = input_kp;
    qaoa_type = input_qaoa_type;
    depth = input_depth;
    bias = input_bias;
    k = copula_k;

    sort_knapsack(kp, RATIO);

    switch (qaoa_type) {
        case QTG:
            apply_int_greedy(kp);
            path_t *int_greedy_sol = path_rep(kp);
            printf("greedy = %ld\n", int_greedy_sol->tot_profit);
            remove_all_items(kp);

            printf("stategen\n");
            qtg_nodes = qtg(kp, bias, int_greedy_sol->vector, &num_states);
            printf("states = %zu\n", num_states);

            free_path(int_greedy_sol);

        case COPULA:
            num_states = POW2(kp->size);
    }

    printf("angle opt\n");
    const double* opt_angles = nlopt_optimizer(opt_type);

    printf("evolution\n");
    fflush(stdout);
    cbs_t *opt_angle_state = quasiadiabatic_evolution(opt_angles);

    printf("exp value\n");
    if (qtg_nodes != NULL) {
        free_nodes(qtg_nodes, num_states);
    }

    const double sol_val = -expectation_value(opt_angle_state);
    const num_t optimal_sol_val = combo_wrap(kp, 0, kp->capacity, FALSE, FALSE, TRUE, FALSE);
    printf("optimal = %ld\n", optimal_sol_val);
    write_plot_data_to_file(opt_angle_state, sol_val, optimal_sol_val);

    if (opt_angle_state != NULL) {
        free(opt_angle_state);
    }

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

    return sol_val;
}
