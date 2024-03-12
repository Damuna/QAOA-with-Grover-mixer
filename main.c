#include <stdio.h>
#include "knapsack.h"
#include "stategen.h"
#include "qaoa_qtg.h"

int main(int argc, const char **argv) {

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
//            print_knapsack(k);

            qaoa_qtg(k, p, bias, samples, POWELL);
        }
    }
    fclose(file);

    return 0;
}

