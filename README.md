# QTG-based AAM-QAOA

This repository contains the source code for the QTG-based Amplitude Amplification-Mixer QAOA (AAM-QAOA) as described in the paper "[Quantum tree generator improves QAOA state-of-the-art for the knapsack problem](https://arxiv.org/abs/2411.00518)". Please cite this paper as:

Paul Christiansen, Lennart Binkowski, Debora Ramacciotti, Sören Wilkening (2024), arXiv:2411.00518 [quant-ph]

The corresponding BibTex entry is:

@misc{Christiansen2024QuantumTreeGeneratorImprovesQAOAStateOfTheArtForTheKnapsackProblem,  
&nbsp;&nbsp; author={Christiansen, Paul and Binkowski, Lennart and Ramacciotti, Debora and Wilkening, Sören},  
&nbsp;&nbsp; year={2024},  
&nbsp;&nbsp; title={{Quantum tree generator improves QAOA state-of-the-art for the knapsack problem}},  
&nbsp;&nbsp; eprint={2411.00518},  
&nbsp;&nbsp; archivePrefix={arXiv},  
&nbsp;&nbsp; primaryClass={quant-ph}  
}

## Background 
The AAM-QAOA is a straight-forward extension of the Grover-mixer QAOA by 
[Baertschi and Eidenbenz](https://arxiv.org/abs/2006.00354), who suggested to shift the complexity in a feasibility 
preserving QAOA from designing a suitable mixer to finding a state preparation routine that generates a uniform 
superposition of all feasible states. The Quantum Tree Generator (QTG) is a technique proposed by 
[Wilkening et al.](https://arxiv.org/abs/2310.06623) as part of a promising new quantum algorithm for the 0-1 Knapsack
(KP) problem. Here, we use the QTG as a state preparation unitary. However, the superposition created by the QTG 
features unequal amplitudes that are influenced by the convention according to which the items of the given KP instance 
are sorted. In analogy to how amplitude amplification generalizes Grover's algorithm, this approach gives rise to a 
more general QAOA framework, which we name AAM-QAOA. Also, we take over their high-level simulation technique, which 
rather turns out as medium-level in this case, as we may not discard states whose profit is below a certain threshold. 
Nevertheless, we are able to classically emulate the QTG-based AAM-QAOA simulation up to 35 variables (qubits).

We identified the Copula-QAOA by [van Dam et al.](https://arxiv.org/abs/2108.08805) as the state-of-the-art competitor
against which we wish to compare the performance of our algorithm. As various attempts of making contact with the 
authors were not successful, we implemented their QAOA from scratch. In contrast to our algorithm, the simulation of
the Copula-QAOA is based on gate applications.


## Structure

### `main.c`

Main file for executing the desired QAOA on a certain KP instance. Command-line argument is of the form 
`benchmark_instance_n[value]_g[value]` where the value of $n$ specifies the number of items and $g$ is an indicator for
the complexity of the instance. The actual values of $n$ and $g$ have to correspond to an existing instance (more on 
instance creation below).

### `benchmark_instances`

Contains one instruction file for every instance that has been created via `generator.cpp` in the `source` directory.
The creation of benchmark instances is based on the KP instance generator by 
[Jooken, Leyman and de Causmaecker](https://www.sciencedirect.com/science/article/abs/pii/S037722172101016X). The 
created instances are stored in the `instances` folder where each subfolder carries the name of an instance and contains
a `test.in` file that comprises the number of items in that instance, its profit-weight pairs as well as the capacity. 
Each file in the `benchmark_instances` directory needs to specify the instance via its name, the QAOA type (QTG or 
Copula) that shall be run, the desired depth $p$, the optimization type (BFGS, Powell or Nelder-Mead) one wants to use 
and a number $m$ of discretization steps for the fine-grid search. Additionally, the QTG-QAOA needs a bias for the QTG 
application. On the other hand, the Copula-QAOA requires parameter values for $k$ and $\theta$. In case that BFGS is 
chosen as classical optimization method, a memory size may be provided as final hyper-parameter.

### `instances`

Contains one folder for every instance that has been created. Next to the defining `test.in` file that was mentioned 
above, each folder holds the results of the simulations conducted. On the highest level, we collect results based on
the QAOA type; next criterion is the depth, and finally the classical optimizer. Running the main file leads to
(over-)writing three files: `resources`, stored on the same level as the classical optimizer subdirectories, holds 
information about the qubit count as well as gate and cycle counts (with and without parallelization). On the deepest
level, `results` contains the number of states in the simulation, the solution value of integer Greedy, the total 
approximation ratios of Greedy and QAOA, and the probability of measuring a (feasible) state whose profit is larger than
the value returned by Greedy. Next to this file, `raw_data` stores the pairs of approximation ratio and probability for 
all involved states in order to not lose information from the simulation.

### `src`

#### `qaoa.c`

Most important file that contains major part of the logic for simulating both QTG-QAOA or Copula-QAOA.

#### `qtg_count.c`

Counts the resources required by the QTG-QAOA for a given KP instance.

#### `copula_count.c`

Counts the resources required by the Copula-QAOA for a given KP instance.

#### `general_count.c`

Wrapper functionality for resource counts that are independent of the QAOA type.

#### `stategen.c`

Applies the QTG to a given KP instance.

#### `generator.cpp`

Creates benchmark instances (see above).

#### `knapsack.c`

Holds functionality that is related to KP instances.

#### `combo.c`

COMBO algorithm described by [Martello, Pisinger and Toth](https://www.jstor.org/stable/2634886).

#### `combowrp.c`

Wrapper functionality for the COMBO algorithm. COMBO is used to determine the optimal solution of the respective KP 
instance at hand, which is needed, e.g., for calculating approximation ratios.

#### `syslinks.c`

Contains simple functionality to make the code OS-agnostic.


## Conclusion

In case you are interested in taking up on our code for your own project, we are always open for suggestions, 
discussions or cooperations. Please don't hesitate to contact us if you have any questions, e.g., via 
[this email](mailto:p.chr.1999@gmail.com).



  
