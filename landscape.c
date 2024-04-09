//
// Created by Sören Wilkening on 12.03.24.
//

#include "include/qaoa_qtg.h"


int main() {
    knapsack_t *k;
    int p, samples, bias;
    char instance[1023];
    char opt_type[32];
    char line[1023];
    FILE *file = fopen("../benchmark_instances.txt", "r");
    while (fgets(line, sizeof(line), file)) { // all the instances will be considered
        if (line[0] != '#') { // lines startin with '#' are ignored
            sscanf(line, "%s %d %s %d %d\n", instance, &p, opt_type, &bias, &samples);
            printf("p = %d\n", p);
            printf("bias = %d\n", bias);
            printf("samples = %d\n", samples);

            printf("%s\n", instance);
            k = create_jooken_knapsack(instance);
        }
    }

    dpth = 5;
    num_smpls = samples;

    sort_knapsack(k, RATIO);
    apply_int_greedy(k);
    path_t *int_greedy_sol = path_rep(k);
    printf("greedy = %ld\n", int_greedy_sol->tot_profit);
//    remove_all_items(k);

    qtg_nodes = qtg(k, bias, int_greedy_sol->vector, &num_states);

    int max = 0;
    int index = 0;
    for (int i = 0; i < num_states; ++i) {
        if (max < qtg_nodes[i].path.tot_profit) {
            max = qtg_nodes[i].path.tot_profit;
            index = i;
        }
//        printf("%ld\n", qtg_nodes[i].path.tot_profit);
    }
//    printf("%ld\n", max);
//    printf("%f\n", qtg_nodes[index].prob);
//    printf("%ld\n", qtg_nodes[num_states - 3].path.tot_profit);
//    printf("%.20f\n", qtg_nodes[num_states - 3].prob);

    double exp = 0;
    for (int l = 0; l < num_states; ++l) {
//                        exp += pow(cabs(opt_angle_state[l].amplitude), 2) * (double) opt_angle_state[l].profit;
//                        if (int_greedy_sol->tot_profit < opt_angle_state[l].profit)
        exp += qtg_nodes[l].prob;
    }
    printf("%f\n", exp);

    double opt_angles[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    double best = 890.608767;
    double init = 0;

    int index1 = 0;
    int index2 = 0;

    int steps = 50;
    for (int i = 0; i < steps; ++i) {
        for (int j = 0; j < steps; ++j) {
            opt_angles[0] = 2 * 3.14159 / steps * i;
            opt_angles[1] = 2 * 3.14159 / steps * j;
            feassol_ampl_t *opt_angle_state = quasiadiabatic_evolution(opt_angles);
            double exp = 0;
            for (int l = 0; l < num_states; ++l) {
                exp += pow(cabs(opt_angle_state[l].amplitude), 2) * (double) opt_angle_state[l].profit;
            }
            if (i == 0 && j == 0) {
                best = exp;
                init = exp;
            }
            if (exp > best) {
                index1 = i;
                index2 = j;
                best = exp;
            }
            printf("\r%3d %3d %3d %3d %f", i, 1, j, 1, best);
            fflush(stdout);
        }
    }
    printf("\nbest = %f", best);
    printf("\ninit = %f\n", init);
    printf("\n%d %d\n", index1, index2);
}