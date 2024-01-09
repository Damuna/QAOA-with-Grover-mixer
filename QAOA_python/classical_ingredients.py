import numpy as np
import random

from knapsack_problem import GenerateKnapsackProblemInstances, KnapsackProblem


class AuxiliaryFunctions:
    
    def find_selected_items(bitstring: str):
        return [idx + 1 for idx, bit in enumerate(bitstring) if bit == "1"]



class SortingProfitsAndWeights:
    """
    Class for sorting profits and weights according to the ratio of profit/weight in
    ascending order.

    Attributes:
    profits (list): the profits (values) of the items 
    weights (list): the costs (weights) of the items

    Methods:
    sorting_profits_weights: sorting according to ratio of profit over weight
    """

    def __init__(self, profits: list, weights: list):
        self.profits = profits
        self.weights = weights

    def sorting_profits_weights(self):
        if len(self.profits) == 0 or len(self.weights) == 0:
            if not len(self.profits) == 0 == len(self.weights):
                raise ValueError("The generation of subproblems may be broken - both profits and weights should be empty in case of no remaining items being affordable!")
            return {"profits": [], "weights": []}
        sorted_profits, sorted_weights = zip(*sorted(zip(self.profits, self.weights), reverse = True, key = lambda k: k[0]/k[1]))
        return {"profits": sorted_profits, "weights": sorted_weights}
    
    def sorting_permutation(self, old_profits: list, old_weights: list):
        old_profit_weight_tuples = [(old_profits[idx], old_weights[idx]) for idx in range(len(old_profits))]
        new_profit_weight_tuples = [(self.profits[idx], self.weights[idx]) for idx in range(len(self.profits))]
        sorting_permutation = []
        for new_tuple in new_profit_weight_tuples:
            first_occurrence = old_profit_weight_tuples.index(new_tuple)
            if first_occurrence in sorting_permutation:
                remaining_occurrences_of_tuple = [idx for idx in range(first_occurrence + 1, len(old_profits)) if old_profit_weight_tuples[idx] == new_tuple]
                for occurrence_idx in remaining_occurrences_of_tuple:
                    if occurrence_idx not in sorting_permutation:
                        sorting_permutation.append(occurrence_idx)
            else: 
                sorting_permutation.append(first_occurrence)
        return sorting_permutation
    


class EvaluatingProfitsAndWeights:
    """
    Class for calculating the profit and weight associated to a given (partial)
    represented as bitstring.

    Attributes:
    problem_instance (KnapsackProblem): the instance of the knapsack problem at hand
    profits_sorted (list): profits sorted according to ratio of profit/weight
    weights_sorted (list): weights sorted according to ratio of profit/weight

    Methods:
    calculate_profit: determining the profit of a given bitstring
    calculate_weight: determining the weight of a given bitstring
    """

    def __init__(self, problem_instance: KnapsackProblem):
        self.profits = SortingProfitsAndWeights(profits = problem_instance.profits, weights = problem_instance.weights).sorting_profits_weights()["profits"]
        self.weights = SortingProfitsAndWeights(profits = problem_instance.profits, weights = problem_instance.weights).sorting_profits_weights()["weights"]

    def calculate_profit(self, bitstring: str):
        as_int_list = list(map(int, list(bitstring)))
        profit = np.dot(np.array(as_int_list), np.array(self.profits[:len(as_int_list)]))
        return profit
    
    def calculate_weight(self, bitstring: str):
        as_int_list = list(map(int, list(bitstring)))
        weight = np.dot(np.array(as_int_list), np.array(self.weights[:len(as_int_list)]))
        return weight


class DynamicalSubproblems(EvaluatingProfitsAndWeights):

    def __init__(self, problem_instance: KnapsackProblem):
        EvaluatingProfitsAndWeights.__init__(self, problem_instance)
        self.capacity = problem_instance.capacity
        self.number_items = problem_instance.number_items

    def calculate_residual_capacity(self, partial_choice: str):
        return self.capacity - self.calculate_weight(partial_choice)
    
    def partial_choice_to_subproblem(self, partial_choice: str):
        residual_capacity = self.calculate_residual_capacity(partial_choice)
        if len(partial_choice) > self.number_items:
            raise ValueError("There cannot be more items specified than existing.")
        residual_profits = []
        residual_weights = []
        for idx in range(len(partial_choice), self.number_items):
            if self.weights[idx] <= residual_capacity: # No need to consider items not even fitting in the residual knapsack alone
                residual_profits.append(self.profits[idx])
                residual_weights.append(self.weights[idx])
        return KnapsackProblem(residual_profits, residual_weights, residual_capacity)

    def complete_partial_choice(self, partial_choice: str):
        number_unspecified_items = self.number_items - len(partial_choice)
        return partial_choice + "0" * number_unspecified_items


    



def main():
    random_kp_instance = GenerateKnapsackProblemInstances.generate_random_kp_instance_for_capacity_ratio_and_maximum_value(
        size = 60,
        desired_capacity_ratio = 0.023,
        maximum_value = 1e3
    )
    profits, weights = random_kp_instance.profits, random_kp_instance.weights
    sorted_profits_weights = SortingProfitsAndWeights(profits, weights).sorting_profits_weights()
    sorted_profits = sorted_profits_weights["profits"]
    sorted_weights = sorted_profits_weights["weights"]
    sorting_permutation = SortingProfitsAndWeights(sorted_profits, sorted_weights).sorting_permutation(profits, weights)
    print("Sorting permutation = ", sorting_permutation)
    


if __name__ == "__main__":
    main()
