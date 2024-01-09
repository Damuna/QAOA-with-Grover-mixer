import numpy as np
from numpy.typing import NDArray
from typing import List, Union, Dict
import itertools
import random
import scipy.optimize as opt
from copy import deepcopy
import time

from knapsack_problem import KnapsackProblem, GenerateKnapsackProblemInstances
from classical_ingredients import SortingProfitsAndWeights, EvaluatingProfitsAndWeights



class AuxiliaryFunctions:

    def sort_kp_instance(kp_instance: KnapsackProblem):
        sorting = SortingProfitsAndWeights(kp_instance.profits, kp_instance.weights).sorting_profits_weights()
        return KnapsackProblem(
            profits = sorting["profits"],
            weights = sorting["weights"],
            capacity = kp_instance.capacity
        )
    
    def is_feasible(kp_instance: KnapsackProblem, current_node: str):
        weight = EvaluatingProfitsAndWeights(kp_instance).calculate_weight(current_node)
        return False if weight > kp_instance.capacity else True
    
    def calculate_profits_of_feasible_solutions(kp_instance: KnapsackProblem, binary_values_of_feasible_solutions: List[int]):
        evaluation = EvaluatingProfitsAndWeights(kp_instance)
        return [evaluation.calculate_profit(bin(binary_value)[2:].zfill(kp_instance.number_items)) for binary_value in binary_values_of_feasible_solutions]



"""
class QTG:

    def __init__(self, kp_instance: KnapsackProblem):
        self.kp_instance = AuxiliaryFunctions.sort_kp_instance(kp_instance)


    def quantum_tree_generator(self):

        feasible_solutions = ["0" * self.kp_instance.number_items]
        feasible_solutions_amplitudes: List[Union[float, int]] = [1]
        binary_values_of_feasible_solutions = [0]

        for item_idx in range(self.kp_instance.number_items):
            feasible_solutions_at_item_idx = []
            for bitstring in feasible_solutions:
                candidate = bitstring[:item_idx] + "1" + bitstring[(item_idx + 1):]
                if AuxiliaryFunctions.is_feasible(self.kp_instance, candidate):
                    feasible_solutions_at_item_idx.append(candidate)
                    binary_values_of_feasible_solutions.append(int(candidate, 2))
                    feasible_solutions_amplitudes[feasible_solutions.index(bitstring)] *= 1/np.sqrt(2)
                    feasible_solutions_amplitudes.append(feasible_solutions_amplitudes[feasible_solutions.index(bitstring)])
            feasible_solutions += feasible_solutions_at_item_idx
        print("QTG done")
        return {"binary values of feasible solutions": binary_values_of_feasible_solutions, "feasible solutions amplitudes": feasible_solutions_amplitudes}
"""


class QuasiAdiabaticEvolution:

    def __init__(self, kp_instance: KnapsackProblem, qtg_output: Dict[str, list], depth: int):
        self.kp_instance = AuxiliaryFunctions.sort_kp_instance(kp_instance)
        self.depth = depth
        self.feasible_solutions_amplitudes = np.array(qtg_output["feasible solutions amplitudes"], dtype = "complex_")
        self.binary_values_of_feasible_solutions = qtg_output["binary values of feasible solutions"]
        self.feasible_solutions_profits = AuxiliaryFunctions.calculate_profits_of_feasible_solutions(self.kp_instance, self.binary_values_of_feasible_solutions)


    def apply_phase_separation_unitary(self, state: NDArray, gamma: Union[int, float]):
        for idx in range(len(state)):
            state[idx] = np.exp(-1j * gamma * self.feasible_solutions_profits[idx]) * state[idx]
        return state
    

    def apply_mixing_unitary(self, state: NDArray, beta: Union[int, float]):
        scalar_product = sum([self.feasible_solutions_amplitudes[idx].conjugate() * state[idx] for idx in range(len(state))])
        for idx in range(len(state)):
            state[idx] = state[idx] + (np.exp(-1j * beta) - 1) * scalar_product * self.feasible_solutions_amplitudes[idx]
        return state

    
    def apply_quasiadiabatic_evolution(self, angles: List[Union[float, int]]):
        
        if len(angles) != 2 * self.depth:
            raise ValueError("Number of provided values for gamma and beta parameters need to be consistent with specified circuit depth.")
        
        gamma_values = angles[0::2]
        beta_values = angles[1::2]
        
        state = deepcopy(self.feasible_solutions_amplitudes)
        for j in range(self.depth):
            state = self.apply_phase_separation_unitary(state, gamma_values[j])
            state = self.apply_mixing_unitary(state, beta_values[j])

        return state
    


class QAOA(QuasiAdiabaticEvolution):

    def __init__(self, kp_instance: KnapsackProblem, qtg_output: Dict[str, list], depth: int):
        self.unsorted_profits = kp_instance.profits
        self.unsorted_weights = kp_instance.weights
        QuasiAdiabaticEvolution.__init__(self, kp_instance, qtg_output, depth) 

    def measure_state(self, state: NDArray):
        probs_dict = {}
        for idx in range(len(state)):
            item_choice = bin(self.binary_values_of_feasible_solutions[idx])[2:].zfill(self.kp_instance.number_items)
            probability = np.abs(state[idx])**2
            probs_dict[f"{item_choice}"] = probability
        if not np.isclose(sum(probs_dict.values()), 1, rtol = 1e-05):
            raise ValueError("Even with tolerance all single probabilities don't sum to one.")
        return probs_dict
    
    def calculate_expectation_value(self, state: NDArray):
        probs_dict = self.measure_state(state)
        return sum([self.feasible_solutions_profits[idx] * list(probs_dict.values())[idx] for idx in range(len(probs_dict))])
    
    def angles_to_value(self, angles: List[Union[float, int]]):
        angle_state = self.apply_quasiadiabatic_evolution(angles)
        return - self.calculate_expectation_value(angle_state)
    
    def optimize(self):
        print("QAOA running...")
        gamma_range = (0, 2 * np.pi)
        beta_range = (0, np.pi)
        bounds_neldermead_format = (gamma_range, beta_range) * self.depth
        initial_guess = list(itertools.chain.from_iterable([[random.uniform(gamma_range[0], gamma_range[1])] + [random.uniform(beta_range[0], beta_range[1])] for _ in range(self.depth)]))
        optimization_result = opt.minimize(fun = self.angles_to_value, x0 = initial_guess, method = "Nelder-Mead", bounds = bounds_neldermead_format)
        optimal_angles = optimization_result.x
        optimal_value = - optimization_result.fun # Minus needed for optimization via Nelder-Mead minimum finding
        optimal_angle_state = self.apply_quasiadiabatic_evolution(optimal_angles)
        probs_dict = self.measure_state(optimal_angle_state)
        sorting_permutation = SortingProfitsAndWeights(self.kp_instance.profits, self.kp_instance.weights).sorting_permutation(self.unsorted_profits, self.unsorted_weights)
        sorted_probs_dict = {
            f"{''.join([list(bitstring)[sorting_permutation.index(idx)] for idx in range(len(sorting_permutation))])}": probs_dict[bitstring] for bitstring in list(probs_dict.keys())
        }
        return {"qaoa result": optimal_value, "solution probabilities": sorted_probs_dict}





def main():
    random_kp_instance = GenerateKnapsackProblemInstances.generate_random_kp_instance_for_capacity_ratio_and_maximum_value(
        size = 60,
        desired_capacity_ratio = 0.023,
        maximum_value = 1e3
    )
    """
    print("Random KP instance = ", random_kp_instance)
    print("Qubits required = ", random_kp_instance.number_items + int(np.floor(np.log2(random_kp_instance.capacity)) + 1))
    qtg_output = QTG(random_kp_instance).quantum_tree_generator()
    start_time = time.time()
    print("QAOA result = ", QAOA(random_kp_instance, qtg_output, depth = 1).optimize()["qaoa result"])
    print("Elapsed time = ", time.time() - start_time)
    """


if __name__ == "__main__":
    main()
