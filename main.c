#include <stdio.h>
#include "knapsack.h"
#include "qaoa.h"
#include "stategen.h"
#include "qaoa.h"

int main(int argc, const char **argv) {

    int p, bias;
    double k;
    char instance[1023];
    char input_qaoa_type[16];
    char input_opt_type[32];
    char line[1023];
    FILE *file = fopen("../benchmark_instances.txt", "r");
    while (fgets(line, sizeof(line), file)) { // all the instances will be considered
        if (line[0] != '#') { // lines startin with '#' are ignored
            sscanf(line, "%s %s %d %s %d %.0f\n", instance, input_qaoa_type, &p, input_opt_type, &bias, &k);
            printf("p = %d\n", p);
            printf("bias = %d\n", bias);
            printf("k = %.0f\n", k);

            printf("%s\n", instance);
            knapsack_t* kp = create_jooken_knapsack(instance);
//            print_knapsack(k);

            printf("QAOA type = %s", input_qaoa_type);
            qaoa_type_t qaoa_type;
            switch (input_qaoa_type) {
                case "qtg":
                    qaoa_type = QTG;
                case "copula":
                    qaoa_type = COPULA;
                default:
                    printf("Error: Input for QAOA type does not match any of the permitted values.");
                    return -1;
            }

            printf("Optimization type = %s", input_opt_type);
            opt_t opt_type;
            switch (input_opt_type) {
                case "bfgs":
                    opt_type = BFGS;
                case "nelder_mead":
                    opt_type = NELDER_MEAD;
                case "powell":
                    opt_type = POWELL;
                default:
                    printf("Error: Input for optimization type does not match any of the permitted values.");
                    return -1;
            }

            qaoa(kp, qaoa_type, depth, opt_type, bias, k);
        }
    }
    fclose(file);

    return 0;
}

