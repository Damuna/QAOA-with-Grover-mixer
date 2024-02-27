/*
 * =============================================================================
 *                            includes
 * =============================================================================
 */

#include "combowrp.h"

/* 
 * =============================================================================
 *                            macros
 * =============================================================================
 */

#define TRUE                    1
#define FALSE                   0

/* 
 * =============================================================================
 *                            Combo wrapper
 * =============================================================================
 */

num_t
combo_wrap(const knapsack_t *k, bit_t first_item, num_t capacity, bool_t def, \
           bool_t relx, bool_t exe_combo, bool_t save) {

    num_t opt_sol;

    // knapsack_t* k_copy = create_empty_knapsack(k->size - first_item, capacity);
    // memcpy(k_copy->items, k->items, k_copy->size * sizeof(item_t));
    // strcpy(k_copy->name, k->name);

    knapsack_t k_copy = {.size = k->size - first_item, .capacity = capacity, \
                         .remain_cost = capacity, .tot_profit = 0, \
                         .items = k->items + first_item, .name = k->name};
    /* check whether instance is trivial */
    if (is_trivial(&k_copy, &opt_sol)) {
        return opt_sol;
    }

    item *f;
    item *l;

    /* Set lower and upper bound */
    // num_t lbi = int_greedy(&k_new, RATIO);
    // num_t ubi = frac_greedy(k, RATIO);

    /* conversation of item_t structure to Combo's item structure */
    item items[k->size - first_item];
    for (size_t i = 0; i < k_copy.size; ++i) {
        items[i].p = k_copy.items[i].profit;
        items[i].w = k_copy.items[i].cost;
        items[i].x = k_copy.items[i].included;
    }
    f = items;
    l = items + k_copy.size - 1;
    /* either start combo or return 0 */
    if (exe_combo) {
        char pathname[256];
        char filename_ndef[256];
        snprintf(pathname, sizeof(pathname), "instances%c%s%ccombo%c", path_sep(), \
             k->name, path_sep(), path_sep());
        create_dir(pathname);
        if (def == 0) {
            snprintf(filename_ndef, sizeof(filename_ndef), \
             "%scombo_counts_def=false.csv", pathname);
        } else {
            snprintf(filename_ndef, sizeof(filename_ndef), \
             "%scombo_counts_def=true.csv", pathname);
        }
        opt_sol = combo(f, l, k_copy.capacity, 0, 0, def, relx);
    } else {
        opt_sol = 0;
    }
    // free(k_copy);
    return opt_sol;
}

/* 
 * =============================================================================
 *                            Combo data
 * =============================================================================
 */

num_t
combo_data(const knapsack_t *k, bit_t first_item, num_t capacity, bool_t def, \
           bool_t relx, bool_t read) {

    num_t opt_sol;

    /* Check whether instance's solution was already calculated by Combo */
    FILE *stream;
    char instancename[256];
    char pathname[256];
    char filename[256];
    char line[128];
    snprintf(instancename, sizeof(instancename), "instances%c%s", path_sep(), \
             k->name);
    snprintf(pathname, sizeof(pathname), "%s%ccombo", instancename, path_sep());
    snprintf(filename, sizeof(filename), "%s%csize=%"PRIu64"_capacity=" \
             "%"PRIu64".txt", pathname, path_sep(), \
             (uint64_t) (k->size - first_item), (uint64_t) capacity);
    if (file_exists(filename) && read) {
        stream = fopen(filename, "r");
        fscanf(stream, "%ld", &opt_sol);
        fclose(stream);
        return opt_sol;
    }

    opt_sol = combo_wrap(k, first_item, capacity, def, relx, TRUE, 0);

    /* save the result */
    if (!file_exists(filename)) {
        create_dir(instancename);
        create_dir(pathname);
        FILE *file = fopen(filename, "w");
        fprintf(file, "%"PRIu64"\n", (uint64_t) opt_sol);
        fclose(file);
    }
    return opt_sol;
}
