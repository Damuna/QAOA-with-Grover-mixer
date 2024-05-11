//
// Created by Sören Wilkening on 09.04.24.
//

#include <stdio.h>
#include "stategen.h"
#include "knapsack.h"
#include "include/qaoa_qtg.h"

int main() {
    knapsack_t *k = create_empty_knapsack(4, 10);
    k->items[0].profit = 5;
    k->items[0].cost = 2;

    k->items[1].profit = 12;
    k->items[1].cost = 9;

    k->items[2].profit = 3;
    k->items[2].cost = 3;

    k->items[3].profit = 2;
    k->items[3].cost = 6;

    // Define the states, that should be generated
    node_t *should_be = calloc(8, sizeof(node_t));
    should_be[0].path.tot_profit = 0;
    should_be[0].prob = 4. / 81;
    sw_init(should_be[0].path.vector, 4);
    should_be[1].path.tot_profit = 2;
    should_be[1].prob = 2. / 81;
    sw_init(should_be[1].path.vector, 4);
    sw_setbit(should_be[1].path.vector, 3);
    should_be[2].path.tot_profit = 3;
    should_be[2].prob = 8. / 81;
    sw_init(should_be[2].path.vector, 4);
    sw_setbit(should_be[2].path.vector, 2);
    should_be[3].path.tot_profit = 5;
    should_be[3].prob = 4. / 81;
    sw_init(should_be[3].path.vector, 4);
    sw_setbit(should_be[3].path.vector, 2);
    sw_setbit(should_be[3].path.vector, 3);
    should_be[4].path.tot_profit = 12;
    should_be[4].prob = 1. / 9;
    sw_init(should_be[4].path.vector, 4);
    sw_setbit(should_be[4].path.vector, 1);
    should_be[5].path.tot_profit = 5;
    should_be[5].prob = 4. / 27;
    sw_init(should_be[5].path.vector, 4);
    sw_setbit(should_be[5].path.vector, 0);
    should_be[6].path.tot_profit = 7;
    should_be[6].prob = 2. / 27;
    sw_init(should_be[6].path.vector, 4);
    sw_setbit(should_be[6].path.vector, 0);
    sw_setbit(should_be[6].path.vector, 3);
    should_be[7].path.tot_profit = 8;
    should_be[7].prob = 4. / 9;
    sw_init(should_be[7].path.vector, 4);
    sw_setbit(should_be[7].path.vector, 0);
    sw_setbit(should_be[7].path.vector, 2);

    apply_int_greedy(k);
    array_t cur;
    sw_init(cur, 4);
    for (int i = 0; i < 4; ++i) { if (k->items[i].included == 1) sw_setbit(cur, i); }
    qtgNodes = qtg(k, 1, cur, &numStates);

    // Check, if the routine "qtg" worked properly
    // bias = 1
    if (numStates == 8) printf("Correct number of states!\n");
    int correct_prob = 1, correct_profit = 1, correct_vector = 1;
    for (int i = 0; i < numStates; ++i) {
        if (qtgNodes[i].path.tot_profit != should_be[i].path.tot_profit) correct_profit = 0;
        if (qtgNodes[i].prob != should_be[i].prob) correct_prob = 0;
        if (!sw_cmp(should_be[i].path.vector, qtgNodes[i].path.vector)) correct_vector = 0;
    }
    if (correct_prob) printf("Correct probabilities!\n");
    else printf("Incorrect Probabilities!\n");
    if (correct_profit) printf("Correct profits!\n");
    else printf("Incorrect Profits!\n");
    if (correct_vector) printf("Correct vectors!\n");
    else printf("Incorrect Vectors!\n");

    // Check, if we get the same expectation value, if we use the angles=(0,0)
    depth = 1;
    double opt_angles[2];
    opt_angles[0] = 0;
    opt_angles[1] = 0;
    cbs_t *opt_angle_state = quasiadiabatic_evolution(opt_angles);
    double exp = 0;
    for (int l = 0; l < numStates; ++l) {
        exp += pow(cabs(opt_angle_state[l].amplitude), 2) * (double) opt_angle_state[l].profit;
    }
    if (fabs(exp - 6.74074) < pow(10, -5)) printf("Correct Expectation for p=1 angles=(0,0)!\n");
    else printf("Incorrect Expectation for p=1 angles=(0,0)!\n");

    // Check, if we get the same expectation value, if we use the angles=(0.11,0.22)
    opt_angles[0] = 0.11;
    opt_angles[1] = 0.22;
    opt_angle_state = quasiadiabatic_evolution(opt_angles);
    exp = 0;
    for (int l = 0; l < numStates; ++l) {
        exp += pow(cabs(opt_angle_state[l].amplitude), 2) * (double) opt_angle_state[l].profit;
    }
    printf("%f\n", exp);
//    if (fabs(exp - 6.74074) < pow(10, -5)) printf("Correct Expectation for p=1 angles=(0,0)!\n");
}