#include <stdio.h>
#include "knapsack.h"
#include "stategen.h"
#include "qaoa_qtg.h"

int main(int argc, const char **argv) {

    knapsack_t *k;

    char instance[1023];
    char opt_type[32];
    int p, samples, bias;
    char line[1023];
    FILE *file = fopen("../benchmark_instances.txt", "r");
    while (fgets(line, sizeof(line), file)) { // all the instances will be considered
        if (line[0] != '#') { // lines startin with '#' are ignored
            sscanf(line, "%s %d %s %d %d\n", instance, &p, opt_type, &bias, &samples);
        }
    }
    fclose(file);

    return 0;
}

