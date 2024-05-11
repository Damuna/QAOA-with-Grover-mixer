#include <stdio.h>
#include "knapsack.h"
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
            sscanf(line, "%s %d %s %d\n", instance, &p, opt_type, &bias);
            printf("p = %d\n", p);
            printf("bias = %d\n", bias);

            printf("%s\n", instance);
            k = create_jooken_knapsack(instance);
//            print_knapsack(k);

            qaoa(k, p, bias, QTG, POWELL);
        }
    }
    fclose(file);

    return 0;
}

