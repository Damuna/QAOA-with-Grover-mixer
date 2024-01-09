from dataclasses import dataclass
import random
import numpy as np
from typing import List


@dataclass
class KnapsackProblem:
    """
    Class for construction and representation of instances of the knapsack problem.

    Only intended for 0-1 knapsack problems, i.e. only one knapsack, one constraint,
    integer-valued weights/costs and values/profits as well as the decision "take it 
    vs. don't take it".

    Furthermore only such instances are permitted for which all profits and all costs
    are non-negative integers, the capacity of the knapsack is a positive integer, 
    the sum of all weights exceeds the capacity and no single weight is larger than the 
    capacity (otherwise the optimal solution would trivially be given be including
    all items or a certain item could never be selected, respectively). 

    Attributes:
    profits (list): the profits (values) of the items 
    weights (list): the costs (weights) of the items
    capacity (int): the capacity of the knapsack
    total_weight (int): the sum of all weights
    number_items (int): the number of items
    """

    profits: List[int]
    weights: List[int]
    capacity: int

    def __post_init__(self):
        if len(self.profits) != len(self.weights):
            raise ValueError("Profits and weights have different length.")

        if any([p for p in self.profits if p < 0]):
            raise ValueError("Profits must not be negative.")

        if any([w for w in self.weights if w < 0]):
            raise ValueError("Weights must not be negative.")

        if self.capacity < 0:
            raise ValueError("Capacity of knapsack must not be negative.")
        
        self.profits = [np.int64(profit) for profit in self.profits]
        self.weights = [np.int64(weight) for weight in self.weights]
        self.capacity = np.int64(self.capacity)
        self.total_weight = sum(self.weights)
        self.number_items = len(self.weights)




exemplary_kp_instances = {
    "A": KnapsackProblem(profits = [1, 2], weights = [1, 1], capacity = 2),
    "B": KnapsackProblem(profits = [1, 4, 2], weights = [1, 3, 2], capacity = 3),
    "C": KnapsackProblem(profits = [3, 1, 2, 1], weights = [1, 1, 2, 2], capacity = 4),
    "D": KnapsackProblem(profits = [6, 5, 8, 9, 6, 7, 3], weights = [2, 3, 6, 7, 5, 9, 4], capacity = 9)
}




class GenerateKnapsackProblemInstances:

    def generate_random_kp_instance_for_capacity_ratio_and_maximum_value(size: int, desired_capacity_ratio: float, maximum_value: int):
        
        def adjust_weights(weights: List[int], capacity: int):
            for idx in range(len(weights)):
                if weights[idx] > capacity:
                    weights[idx] = random.randint(1, capacity)
            missing_difference = int(1/desired_capacity_ratio * capacity) - sum(weights[:-1])
            weights = [weight + int(np.ceil(missing_difference / size)) for weight in weights]
            if len([weight for weight in weights if weight > capacity]) > 0:
                weights = adjust_weights(weights, capacity)
            return weights
        
        profits = [random.randint(1, maximum_value) for _ in range(size)]
        weights = [random.randint(1, maximum_value) for _ in range(size)]
        capacity = int(np.ceil(desired_capacity_ratio * sum(weights)))
        weights = adjust_weights(weights, capacity)
        
        return KnapsackProblem(profits, weights, capacity)
    



def main():
    generated_random_kp_instance = GenerateKnapsackProblemInstances.generate_random_kp_instance_for_capacity_ratio_and_maximum_value(size = 1000, desired_capacity_ratio = 0.01, maximum_value = 1000)
    print(generated_random_kp_instance)
    print(generated_random_kp_instance.capacity / sum(generated_random_kp_instance.weights))


if __name__ == "__main__":
    main()
        

