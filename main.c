#include <stdio.h>
#include "knapsack.h"
#include "qaoa.h"
#include "stategen.h"
#include "qaoa.h"

int main(int argc, const char **argv) {

    int p, m, bias, memory_size;
    double k, theta;
    char instance[1023];
    char input_qaoa_type[16], input_opt_type[16];
    char line[1023];

    const char *benchmark_instance = argv[1];
    char *path_to_benchmark = calloc(1024, sizeof(char));
    sprintf(path_to_benchmark, "../benchmark_instances/%s.txt", benchmark_instance);
    FILE *file = fopen(path_to_benchmark, "r");

    while (fgets(line, sizeof(line), file)) { // all the instances will be considered
        if (line[0] != '#') { // lines startin with '#' are ignored
            sscanf(
                line,
                "%s %s %d %s %d %d %lf %lf %d\n",
                instance, input_qaoa_type, &p, input_opt_type, &m, &bias, &k,  &theta, &memory_size
            );
            printf("\n===== Input parameters =====\n");

            printf("Instance = %s\n", instance);
            char path_to_instance[1024];
            sprintf(path_to_instance, "..%cinstances%c%s%c", path_sep(), path_sep(), instance, path_sep());
            strcat(path_to_instance, "test.in");
            knapsack_t* kp = create_jooken_knapsack(path_to_instance);

            printf("QAOA type = %s\n", input_qaoa_type);
            qaoa_type_t qaoa_type;
            if (strcmp(input_qaoa_type, "qtg") == 0) {
                qaoa_type = QTG;
            } else if (strcmp(input_qaoa_type, "copula") == 0) {
                qaoa_type = COPULA;
            } else {
                printf("Error: Input for QAOA type does not match any of the permitted values.");
                return -1;
            }

            printf("p = %d\n", p);

            printf("Optimization type = %s\n", input_opt_type);
            opt_t opt_type;
            if (strcmp(input_opt_type, "powell") == 0) {
                opt_type = POWELL;
            } else if (strcmp(input_opt_type, "nelder-mead") == 0) {
                opt_type = NELDER_MEAD;
            } else if (strcmp(input_opt_type, "bfgs") == 0) {
                opt_type = BFGS;
            } else {
                printf("Error: Input for optimization type does not match any of the permitted values.");
                return -1;
            }

            printf("m = %d\n", m);

            switch (qaoa_type) {
                case QTG:
                    printf("bias = %d\n", bias);
                case COPULA:
                    printf("k = %.1f\n", k);
                    printf("theta = %.1f\n", theta);
            }

            if (opt_type == BFGS) {
                printf("Memory size for BFGS = %d\n", memory_size);
            }

            strcat(path_to_instance, input_qaoa_type);
            create_dir(path_to_instance);
            char depth_string[16];
            sprintf(depth_string, "%cp_%d%c", path_sep(), p, path_sep());
            strcat(path_to_instance, depth_string);
            create_dir(path_to_instance);
            strcat(path_to_instance, input_opt_type);
            create_dir(path_to_instance);

            qaoa(instance, kp, qaoa_type, p, opt_type, m, bias, k, theta, memory_size);
        }
    }
    fclose(file);

    return 0;
}

