/*
 * =============================================================================
 *                            includes
 * =============================================================================
 */


#include "qaoa_qtg.h"

/*
 * =============================================================================
*                           Macros & global variables
 * =============================================================================
 */

#define TRUE        1
#define FALSE       0

num_t dpth;
size_t num_smpls;
size_t num_states;
node_t* qtg_nodes;


/*
 * =============================================================================
 *                              Free metastates
 * =============================================================================
 */

void
free_metastates(metastate_t metastates[], size_t num_metastates) {
    for (size_t idx = 0; idx < num_metastates; ++idx) {
        sw_clear(metastates[idx].choice_profit.vector);
    }
    free(metastates);
}


/*
 * =============================================================================
*                                      Utils
 * =============================================================================
 */

double
randomValueBetween0And1OnWindowsOrLinux() {
#if defined(_WIN32) || defined(_WIN64)
    srand(time(NULL));
    return rand() / RAND_MAX;
#else
    srand48(time(NULL));
        return drand48();
#endif
}

void
writePlotDataToFile(metastate_t *angle_state, num_t optimal_solution_value) {
    FILE* stream;

    create_dir(""); // TODO Create meaningful directories according to generated instances
    stream = fopen("testfile", "w"); // TODO Create meaningful filenames according to generated instances

    for (size_t idx = 0; idx < num_states; ++idx) {
        double approx_ratio = (double)angle_state[idx].choice_profit.tot_profit / optimal_solution_value;
        fprintf(stream, "%f %f\n", approx_ratio, angle_state[idx].probability);
    }
    fclose(stream);
}


/*
 * =============================================================================
 *                            Quasi-Adiabatic Evolution
 * =============================================================================
 */

void
phase_separation_unitary(metastate_t *angle_state, double gamma) {
    for (size_t idx = 0; idx < num_states; ++idx) {
        angle_state[idx].amplitude *= cexp(-I * gamma * angle_state[idx].choice_profit.tot_profit);
    }
}


void
mixing_unitary(metastate_t *angle_state, double beta) {
    cmplx scalar_product = 0.0;
    for (size_t idx = 0; idx < num_states; ++idx) {
        scalar_product += qtg_nodes[idx].prob * angle_state[idx].amplitude;
    }

    for (int idx = 0; idx < num_states; ++idx) {
        angle_state[idx].amplitude += (cexp(-I * beta) - 1.0) * scalar_product * qtg_nodes[idx].prob;
    }
}


metastate_t*
quasiadiabatic_evolution(const double *angles) {

    metastate_t* angle_state = malloc(num_states * sizeof(metastate_t));
    for (size_t idx = 0; idx < num_states; ++idx) {
        angle_state[idx].choice_profit = qtg_nodes[idx].path.choice_profit;
        angle_state[idx].amplitude = qtg_nodes[idx].prob + I * 0.0;
    }

    for (int j = 0; j < dpth; ++j) {
        phase_separation_unitary(angle_state, angles[2 * j]); // gamma values are even positions in angles since starting at index 0
        mixing_unitary(angle_state, angles[2 * j + 1]); // beta values are odd positions in angles since starting at index 0
    }

    return angle_state;
}


/*
 * =============================================================================
 *                          Measurement & Expectation Value
 * =============================================================================
 */

int
measurement(metastate_t *angle_state) {
    double random_number = randomValueBetween0And1OnWindowsOrLinux();

    // Use the inverse transform sampling to draw a basis state index
    double cumulative_probability = 0.0;
    for (int idx = 0; idx < num_states; ++idx) {
        cumulative_probability += cabs(angle_state[idx].amplitude) * cabs(angle_state[idx].amplitude);

        if (cumulative_probability > random_number) {
            return idx;
        }
    }

    return -1; // Should not happen, indicates an error
}

double
expectation_value(metastate_t *angle_state) {
    // Sample from the probability distribution induced by angle_state to get probabilities for the feasible states
    for (size_t sample = 0; sample < num_smpls; ++sample) {
        int measured_basis_state = measurement(angle_state);
        if (measured_basis_state != -1) {
            angle_state[measured_basis_state].probability += 1.0 / num_smpls;
        }
    }

    // Calculate expectation value according to obtained probabilities
    double expectation_value = 0.0;
    for (int idx = 0; idx < num_states; ++idx) {
        expectation_value += angle_state[idx].choice_profit.tot_profit * angle_state[idx].probability;
    }

    return expectation_value;
}

/*
 * =============================================================================
 *                            Optimization Functions
 * =============================================================================
 */


double objective(unsigned n, const double *angles, double *grad, void *my_func_data) {
    metastate_t *angle_state = quasiadiabatic_evolution(angles);
    double exp_value = expectation_value(angle_state);
    if (angle_state != NULL) {
        free_metastates(angle_state, num_states);
    }
    // grad is NULL bcs both Nelder Mead and Powell are derivative-free algorithms
    return - exp_value;
}

double * nelder_mead(){
    double randomValue;
    double x[2 * dpth] ;
    void *f_data = NULL;
    nlopt_opt opt = nlopt_create(NLOPT_LN_NELDERMEAD, 2 * dpth);

    // Set your optimization parameters
    nlopt_set_xtol_rel(opt, 1e-6);

    // Set the objective function
    nlopt_set_min_objective(opt, objective, f_data);

    // Set initial guess
    // Set all gammas to random and betas to 0

    for (size_t i = 0; i < dpth; i++)
        if((i + 1) % 2){
            randomValue = randomValueBetween0And1OnWindowsOrLinux() * 2.0 * M_PI;
            x[i] = randomValue;
            }
        else{
            x[i] = 0;
        }

    // Run the optimization
    nlopt_result result = nlopt_optimize(opt, x, NULL);

    // Check the result and print the optimized values
    if (result < 0) {
        printf("NLOpt failed with code %d\n", result);
    } else {
        printf("Found minimum at f(%g, %g) = %g\n", x[0], x[1], objective(num_states, x, NULL, NULL)); // TODO There will be more than 2 angles, so printing x[0] and x[1] is meaningless
    }

    // Clean up
    nlopt_destroy(opt);
    return x;
}

double * powell(){
    double randomValue;
    double x[2 * dpth] ;
    void *f_data = NULL;

    nlopt_opt opt = nlopt_create(NLOPT_LN_BOBYQA, 2 * dpth);

    // Set your optimization parameters
    nlopt_set_xtol_rel(opt, 1e-6);

    // Set the objective function
    nlopt_set_min_objective(opt, objective, f_data);

    // Set initial guess
    // Set all gammas to random and betas to 0

    for (size_t i = 0; i < dpth; i++)
        if((i + 1) % 2){
            randomValue = randomValueBetween0And1OnWindowsOrLinux() * 2.0 * M_PI;
            x[i] = randomValue;
        }
        else{
            x[i] = 0;
        }

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

double
qaoa_qtg(knapsack_t* k, num_t depth, size_t bias, size_t num_samples, enum OptimizationType optimizationType) {

    double* opt_angles;

    dpth = depth;
    num_smpls = num_samples;

    sort_knapsack(k, RATIO);
    apply_int_greedy(k);
    path_t* int_greedy_sol = path_rep(k);
    remove_all_items(k);

    qtg_nodes = qtg(k, bias, int_greedy_sol->choice_profit.vector, &num_states);

    free_path(int_greedy_sol);

    switch (optimizationType) {
        case BFGS:
            //
            break;
        case NELDER_MEAD:
            opt_angles = nelder_mead();
            break;
        case POWELL:
            opt_angles = powell();
            break;
    }

    metastate_t* opt_angle_state = quasiadiabatic_evolution(opt_angles);

    if (qtg_nodes != NULL) {
        free_nodes(qtg_nodes, num_states);
    }

    num_t optimal_sol_val = combo_wrap(k, 0, k->capacity, FALSE, FALSE, TRUE, FALSE);
    writePlotDataToFile(opt_angle_state, optimal_sol_val);

    double sol_val = - expectation_value(opt_angle_state);

    if (opt_angle_state != NULL) {
        free_metastates(opt_angle_state, num_states);
    }

    // TODO Gate count

    // TODO Collect and combine probability dictionary entries with the same approximation ratios -> Do this in python to avoid the need of looping twice through the array?

    // TODO Wrap these in the shape of qaoa_result_t and return the result -> If writing the probabilities to an external file, we could return only the solution value itself

    return sol_val;
}
