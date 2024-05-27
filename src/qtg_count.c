/* 
 * =============================================================================
 *                            includes
 * =============================================================================
 */

#include "qtg_count.h"
#include "general_count.c"

/* 
 * =============================================================================
 *                            macros
 * =============================================================================
 */

#define TRUE                    1
#define FALSE                   0

#define SWAP(a, b, T)           do { register T q; q = *(a); *(a) = *(b); \
                                *(b) = q; } while(0)

/* 
 * =============================================================================
 *                            enum names
 * =============================================================================
 */

const char*
get_qft_name(qft_t method) {
    switch (method) {
        case COPPERSMITH: {
            return "Coppersmith";
        }

        default: {
            return "Unspecified QFT method!";
        }
    }
}

const char*
get_mc_name(mc_t method) {
    switch (method) {
        case TOFFOLI: {
            return "Toffoli";
        }

        default: {
            return "Unspecified multi-control decomposition!";
        }
    }
}

/* 
 * =============================================================================
 *                                  utils
 * =============================================================================
 */

bit_t
num_bits(num_t number) {
    bit_t result = 0;
    do {
        number >>= 1; /* removes last bit */
        ++result; /* increases counter */
    } while (0 != number);

    return result;
}

bit_t
num_set_bits(num_t number) {
    bit_t result = 0;
    while (0 != number) {
        if (number & 1) {
            ++result; // least significant bit is set, i.e. has a value of 1
        }
        number >>= 1; // move to the next bit by right shifting number by 1
    }

    return result;
}

bit_t
lso(num_t number) {
    bit_t result = 0;
    while (!(number & 1)) {
        number >>= 1;
        ++result;
    }

    return result;
}

/* 
 * =============================================================================
 *                            counts
 * =============================================================================
 */

bit_t
path_reg_size(const knapsack_t* k) {
    return k->size;
}

bit_t
cost_reg_size(const knapsack_t* k) {
    return num_bits(k->capacity);
}

bit_t
anc_count_qft(bit_t reg_size, qft_t method) {
    switch (method) {
        case COPPERSMITH: { /* original QFT version without swaps */
            return 0;
        }

        default: {
            printf("Unspecified QFT method!");
            return -1;
        }
    }
}

count_t
cycle_count_qft(bit_t reg_size, qft_t method, bool_t tof_decomp) {
    switch (method) {
        case COPPERSMITH: { /* original QFT version without swaps */
            return 2 * reg_size - 1; 
        }

        default: {
            printf("Unspecified QFT method!");
            return -1;
        }
    }
}

count_t
gate_count_qft(bit_t reg_size, qft_t method, bool_t tof_decomp) {
    switch (method) {
        case COPPERSMITH: { /* original QFT version without swaps */
            return reg_size * (reg_size + 1) / 2; 
        }

        default: {
            printf("Unspecified QFT method!");
            return -1;
        }
    }
}

bit_t
anc_count_mc(bit_t control_qubits, mc_t method) {
    switch (method) {
        case TOFFOLI: { /* decompose into cascade of Toffoli gates */
            return control_qubits - 1;
        }

        default: {
            printf("Unspecified multi-control decomposition!");
            return -1;
        }
    }
}

count_t
cycle_count_mc(bit_t control_qubits, mc_t method, bool_t tof_decomp) {
    switch (method) {
        case TOFFOLI: { /* decompose into cascade of Toffoli gates */
            count_t toffoli_cycles = 2 * num_bits(control_qubits - 1);
            return (tof_decomp) ? 5 * toffoli_cycles + 1 : toffoli_cycles + 1;
            /* decomposing Toffolis yields 5 addtional gates per Toffoli gate*/
        }

        default: {
            printf("Unspecified multi-control decomposition!");;
            return -1;
        }
    }
}

count_t
gate_count_mc(bit_t control_qubits, mc_t method, bool_t tof_decomp) {
    switch (method) {
        case TOFFOLI: { /* decompose into cascade of Toffoli gates */
            count_t toffolis = 2 * (control_qubits - 1);
            return (tof_decomp) ? 5 * toffolis + 1 : toffolis + 1;
            /* decomposing Toffolis yields 5 addtional gates per Toffoli gate*/
        }

        default: {
            printf("Unspecified multi-control decomposition!");
            return -1;
        }
    }
}

bit_t
anc_count_comp(bit_t reg_size, num_t to_compare, mc_t method, \
               bool_t unnegated) {
    bit_t i = 0;
    if (unnegated) {
        /* check for rightmost zero in binary represented (to_compare - 1) */
        while (((to_compare - 1) >> i) & 1) {
            ++i;
        }
    } else {
        /* check for rightmost one in binary represented to_compare */
        while (!((to_compare >> i) & 1)) {
            ++i;
        }
    }
    return anc_count_mc(reg_size - i, method);
}

count_t
cycle_count_comp(bit_t reg_size, num_t to_compare, mc_t method, \
                 bool_t unnegated, bool_t tof_decomp) {
    count_t cycle_count;

    if (unnegated) {
        /* check whether being strictly larger than (to_compare - 1) */
        cycle_count = 0;
        for (bit_t i = 0; i < num_bits(to_compare); ++i) {
            if (!(((to_compare - 1)>> i) & 1)) {
                /* invoke multi-controlled single-qubit gate if bit is a zero */
                cycle_count += cycle_count_mc(reg_size - i, method, tof_decomp);
            }
        }
    } else {
        /* check whether being strictly smaller than to_compare */
        cycle_count = 1; /* adding one single-qubit gate in the beginning */
        for (bit_t i = 0; i < num_bits(to_compare); ++i) {
            if ((to_compare >> i) & 1) {
                /* invoke multi-controlled single-qubit gate if bit is a one */
                cycle_count += cycle_count_mc(reg_size - i, method, tof_decomp);
            }
        }
    }
    return cycle_count;
}

count_t
gate_count_comp(bit_t reg_size, num_t to_compare, mc_t method, \
                bool_t unnegated, bool_t tof_decomp) {
    count_t gate_count;

    if (unnegated) {
        /* check whether being strictly larger than (to_compare - 1) */
        gate_count = 0;
        for (bit_t i = 0; i < num_bits(to_compare); ++i) {
            if (!(((to_compare - 1)>> i) & 1)) {
            /* invoke multi-controlled single-qubit gate if bit is a zero */
                gate_count += gate_count_mc(reg_size - i, method, tof_decomp);
            }
        }
    } else {
        /* check whether being strictly smaller than to_compare */
        gate_count = 1; /* adding one single-qubit gate in the beginning */
        for (bit_t i = 0; i < num_bits(to_compare); ++i) {
            if ((to_compare >> i) & 1) {
                /* invoke multi-controlled single-qubit gate if bit is a one */
                gate_count += gate_count_mc(reg_size - i, method, tof_decomp);
            }
        }
    }
    return gate_count;
}

bit_t
qubit_count_qtg(const knapsack_t* k) {
    bit_t reg_a = path_reg_size(k);
    bit_t reg_b = cost_reg_size(k);

    bit_t reg_anc = MAX(reg_a, reg_b);

    return reg_a + reg_b + reg_anc - 1; // TODO Is this -1 still correct?
}

count_t
cycle_count_qtg(const knapsack_t* k, qft_t method_qft, mc_t method_mc, bool_t tof_decomp) {
    bit_t reg_b = cost_reg_size(k);
    bit_t break__item = break_item(k);
    count_t cost_qft = cycle_count_qft(reg_b, method_qft, tof_decomp);

    /*
     * first block consisting of:
     * - first comparison on the cost register for controlled operation on the
     *   path register
     * - QFT on the cost register
     * - first subtraction on the cost register
     * - inverse QFT on the cost register
     */
    count_t first_block;
    /* check whether unnegated or negated first comparison is more efficient */

    /* first item can always be included => no comparison needed */
    if (break__item == 0) {
        first_block = 2 * cost_qft + 1;
    } else {
        first_block = cost_qft + 1;
    }

    /* 
     * second block consisting of k->size - 2 similar blocks, each consisting of:
     * - comparison on the cost register for controlled operation on the path
     *   register (only after break item)
     * - QFT on the cost register (only after break item)
     * - subtraction on the cost register
     * - inverse QFT on the cost register (starting at break item)
     */
    count_t second_block = 0;
    count_t second_comp;

    for (bit_t i = 1; i < break__item; ++i) {
        second_block += 2 * num_bits(reg_b - lso(k->items[i].cost) - 1) + 2;
        /* 
         * for each item: Since comparison and subtraction on the cost register
         * cannot be parallelized, both their lengths have to be added.
         */
    }

    second_block += num_bits(reg_b - lso(k->items[break__item].cost) - 1) + 1 + cost_qft;
    

    for (bit_t i = break__item + 1; i < k->size - 1; ++i) {
        /*
         * for each item, check whether unnegated or negated comparison is more
         * efficient
         */
        second_comp = MIN(cycle_count_comp(reg_b, ((k->items)[i]).cost, \
                                           method_mc, TRUE, tof_decomp), \
                          cycle_count_comp(reg_b, ((k->items)[i]).cost, \
                                           method_mc, FALSE, tof_decomp));

        second_block += second_comp + 2 * cost_qft + 1;
        /* 
         * for each item: Since comparison and subtraction on the cost register
         * cannot be parallelized, both their lengths have to be added.
         */
    }

    /*
     * third block consisting of:
     * - last comparison on the cost register for controlled operation on the
     *   path register
     * Note: the last subtraction together with the QFT and inverse QFT on the
     *       cost register can be dropped, since no further comparisons have to
     *       be executed and the cost register is not of independent interest
     *       for further applications.
     */

    /* check whether unnegated or negated last comparison is more efficient */
    count_t third_block = MIN(cycle_count_comp(reg_b, \
                                             k->items[k->size - 1].cost, \
                                             method_mc, TRUE, tof_decomp), \
                            cycle_count_comp(reg_b, \
                                             k->items[k->size - 1].cost, \
                                             method_mc, FALSE, tof_decomp));

    return first_block + second_block + third_block;
}

count_t
gate_count_qtg(const knapsack_t* k, qft_t method_qft, mc_t method_mc, bool_t tof_decomp) {
    bit_t reg_b = cost_reg_size(k);
    count_t sub_gates = 0;
    count_t comp_gates = 0;
    bit_t break__item = break_item(k);
    for (bit_t i = break__item + 1; i < k->size; ++i) {
        comp_gates += MIN(gate_count_comp(reg_b, k->items[i].cost, method_mc, \
                                          TRUE, tof_decomp), \
                          gate_count_comp(reg_b, k->items[i].cost, method_mc, \
                                          FALSE, tof_decomp));
    }

    for (bit_t i = 0; i < k->size - 1; ++i) {
        sub_gates += reg_b - lso(k->items[i].cost);
    }

    return comp_gates /* all comparison gates */ \
           + (2 * ((count_t)(k->size) - 1 - break__item - 1) + 1)
             * gate_count_qft(reg_b, method_qft, tof_decomp) /* (un)QFT-ing cost register */ \
           + sub_gates; /* all subtraction gates */
}

bit_t
qubit_count_qtg_mixer(const knapsack_t* k) {
    return qubit_count_qtg(k);
}

count_t
cycle_count_qtg_mixer(const knapsack_t* k, qft_t method_qft, mc_t method_mc, bool_t tof_decomp) {
    count_t cost_qtg = cycle_count_qtg(k, method_qft, method_mc, tof_decomp);
    count_t cost_mc = cycle_count_mc(k->size - 1, method_mc, tof_decomp);
    return 2 * (cost_qtg + 1) + cost_mc; // One additional cycle per QTG application needed to prepare the cost register
}

count_t
gate_count_qtg_mixer(const knapsack_t* k, qft_t method_qft, mc_t method_mc, bool_t tof_decomp) {
    count_t cost_qtg = gate_count_qtg(k, method_qft, method_mc, tof_decomp);
    count_t cost_prep = num_set_bits(k->capacity); // Preparation of 2nd register to encode full capacity
    count_t cost_mc = gate_count_mc(k->size - 1, method_mc, tof_decomp);
    return 2 * (cost_qtg + cost_prep) + cost_mc;
}

bit_t
qubit_count_qtg_qaoa(const knapsack_t* k) {
    return qubit_count_qaoa(k, qubit_count_qtg_mixer(k));
}

count_t
cycle_count_qtg_qaoa(
    const knapsack_t* k,
    const num_t depth,
    const qft_t method_qft,
    const mc_t method_mc,
    const bool_t tof_decomp
) {
    return cycle_count_qaoa(depth, cycle_count_qtg_mixer(k, method_qft, method_mc, tof_decomp));
}

count_t
gate_count_qtg_qaoa(
    const knapsack_t* k,
    const num_t depth,
    const qft_t method_qft,
    const mc_t method_mc,
    const bool_t tof_decomp
) {
    return gate_count_qaoa(k, depth, gate_count_qtg_mixer(k, method_qft, method_mc, tof_decomp));
}

void
print_qtg_counts(const knapsack_t* k, qft_t method_qft,
                 mc_t method_mc, bool_t tof_decomp) {
    printf("The QTG resources are estimated with the following methods:\n");
    printf("Implementation of the quantum Fourier transform (QFT): %s\n", \
           get_qft_name(method_qft));
    printf("Decomposition of multi-controlled single-qubit-gates: %s\n", \
           get_mc_name(method_mc));
    printf("Toffoli gates have%s been decomposed.\n--------\n", \
           (tof_decomp) ? "" : " not");
    printf("Number of qubits (main + ancilla registers): %u\n", \
           qubit_count_qtg(k));
    printf("Number of QTG cycles: %u\n", cycle_count_qtg(k, method_qft, method_mc, tof_decomp));
    printf("Number of QTG gates: %u\n", gate_count_qtg(k, method_qft, method_mc, tof_decomp));
    printf("Number of mixer cycles: %u\n", cycle_count_qtg_mixer(k, method_qft, method_mc, tof_decomp));
    printf("Number of mixer gates: %u\n\n", gate_count_qtg_mixer(k, method_qft, method_mc, tof_decomp));
}
