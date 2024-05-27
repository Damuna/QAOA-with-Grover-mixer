/*
 * =============================================================================
 *                            includes
 * =============================================================================
 */

#include "general_qaoa.h"



/*
 * =============================================================================
 *                          Evaluation & Optimization
 * =============================================================================
 */

cbs_t *
quasiadiabatic_evolution(
    void (*initial_state_prep)(cbs_t*),
    void (*phase_separation_unitary)(cbs_t*, double),
    void (*mixing_unitary)(cbs_t*, double),
    const num_t depth,
    const size_t num_states,
    const double *angles
) {

    cbs_t *angle_state = malloc(num_states * sizeof(cbs_t));
    initial_state_prep(angle_state);

    for (int j = 0; j < depth; ++j) {
        phase_separation_unitary(angle_state,
                                 angles[2 * j]); // gamma values are even positions in angles since starting at index 0
        mixing_unitary(angle_state,
                       angles[2 * j + 1]); // beta values are odd positions in angles since starting at index 0
    }

    return angle_state;
}


double
expectation_value(const size_t num_states, const cbs_t* angle_state) {
    double exp_val = 0;

    for (size_t idx = 0; idx < num_states; ++idx) {
       const double prob = cabs(angle_state[idx].amplitude) * cabs(angle_state[idx].amplitude);
       exp_val += prob * angle_state[idx].profit;
    }
    return exp_val;
}


double
angles_to_value(
    void (*initial_state_prep)(cbs_t*),
    void (*phase_separation_unitary)(cbs_t*, double),
    void (*mixing_unitary)(cbs_t*, double),
    const num_t depth,
    const size_t num_states,
    const double *angles
) {
    cbs_t *angle_state = quasiadiabatic_evolution(
        initial_state_prep, phase_separation_unitary, mixing_unitary, depth, num_states, angles
    );

    const double exp_value = expectation_value(num_states, angle_state);
    if (angle_state != NULL) {
        free(angle_state);
    }
    return -exp_value;
}
