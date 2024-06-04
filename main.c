#include <stdio.h>
#include "knapsack.h"
#include "qaoa.h"
#include "stategen.h"
#include "qaoa.h"

int main(int argc, const char **argv) {

    int p, m, bias;
    double k, theta;
    char instance[1023];
    char input_qaoa_type[16];
    char input_opt_type[32];
    char line[1023];
    FILE *file = fopen("../benchmark_instances.txt", "r");
    while (fgets(line, sizeof(line), file)) { // all the instances will be considered
        if (line[0] != '#') { // lines startin with '#' are ignored
            sscanf(
                line,
                "%s %s %d %s %d %d %.0f %.0f\n",
                instance, input_qaoa_type, &p, input_opt_type, &m, &bias, &k,  &theta
            );
            printf("p = %d\n", p);
            printf("bias = %d\n", bias);
            printf("k = %.0f\n", k);
            printf("theta = %.0f\n", theta);

            printf("%s\n", instance);
            knapsack_t* kp = create_jooken_knapsack(instance);
//            print_knapsack(k);

            printf("QAOA type = %s", input_qaoa_type);
            qaoa_type_t qaoa_type;
            if (strcmp(input_qaoa_type, "qtg") == 0) {
                qaoa_type = QTG;
            } else if (strcmp(input_qaoa_type, "copula") == 0) {
                qaoa_type = COPULA;
            } else {
                printf("Error: Input for QAOA type does not match any of the permitted values.");
                return -1;
            }

            printf("Optimization type = %s", input_opt_type);
            opt_t opt_type;
            if (strcmp(input_opt_type, "bfgs") == 0) {
                opt_type = BFGS;
            } else if (strcmp(input_opt_type, "nelder_mead") == 0) {
                opt_type = NELDER_MEAD;
            } else if (strcmp(input_opt_type, "powell") == 0) {
                opt_type = POWELL;
            } else {
                printf("Error: Input for optimization type does not match any of the permitted values.");
                return -1;
            }

            char path[1023] = "../instances/";
            strcat(path, instance);
            create_dir(strcat(path, input_qaoa_type));

            qaoa(instance, kp, qaoa_type, depth, opt_type, m, bias, k, theta);
        }
    }
    fclose(file);

    return 0;
}

