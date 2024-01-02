from time import time
import gurobipy as gp
from gurobipy import GRB
import numpy as np

def read_instance(instance):
    file2 = open(instance, "r").read().split("\n")
    file2.pop()

    n = int(file2[0])
    Z = int(file2[-1])

    file2.pop(0)
    file2.pop()

    p, w = [0] * n, [0] * n

    a = 0
    for i in file2:
        j = i.split()
        p[a] = int(j[1])
        w[a] = int(j[2])
        a += 1

    return n, Z, p, w

def solve_knapsack(instance):

    n, capacity, values, weights = read_instance(instance)

    num_items = len(weights)

    # Create a new Gurobi model
    model = gp.Model("knapsack")

    # mute model
    model.setParam('OutputFlag', 0)

    # Create binary decision variables for each item
    x = model.addVars(num_items, vtype=GRB.BINARY, name="x")

    # Set the objective to maximize the total value
    model.setObjective(sum(values[i] * x[i] for i in range(num_items)), sense=GRB.MAXIMIZE)

    # Add the constraint: total weight <= capacity
    model.addConstr(sum(weights[i] * x[i] for i in range(num_items)) <= capacity, name="capacity")


    # Optimize the model
    t1 = time()
    model.optimize()
    t2 = time() - t1
    #
    # if model.status == GRB.OPTIMAL:
    #     print("Objective value:", model.objVal)
    #     print(*[int(x[i].X) for i in range(num_items)])
    # else:
    #     print("No feasible solution found.")

    # Dispose of the model
    model.dispose()
    return t2 * 2.6 * 10 ** 9

# file = open("benchmark_instances.txt", "r").read().split("\n")
# file.pop()
file = ['instances/problemInstances/n_30_c_10000000000_g_2_f_0.3_eps_1e-05_s_300/test.in',
 'instances/problemInstances/n_30_c_10000000000_g_2_f_0.3_eps_0_s_200/test.in',
 'instances/problemInstances/n_30_c_10000000000_g_2_f_0.3_eps_0_s_300/test.in',
 'instances/problemInstances/n_30_c_10000000000_g_2_f_0.3_eps_0_s_100/test.in',
 'instances/problemInstances/n_30_c_10000000000_g_2_f_0.3_eps_1e-05_s_200/test.in',
 'instances/problemInstances/n_30_c_10000000000_g_2_f_0.3_eps_1e-05_s_100/test.in',
 'instances/problemInstances/n_30_c_10000000000_g_4_f_0.3_eps_0_s_100/test.in',
 'instances/problemInstances/n_30_c_10000000000_g_4_f_0.3_eps_0_s_200/test.in',
 'instances/problemInstances/n_30_c_10000000000_g_4_f_0.3_eps_0_s_300/test.in',
 'instances/problemInstances/n_30_c_10000000000_g_4_f_0.3_eps_1e-05_s_100/test.in',
 'instances/problemInstances/n_30_c_10000000000_g_4_f_0.3_eps_1e-05_s_200/test.in',
 'instances/problemInstances/n_30_c_10000000000_g_4_f_0.3_eps_1e-05_s_300/test.in',
 'instances/problemInstances/n_30_c_10000000000_g_6_f_0.3_eps_0_s_200/test.in',
 'instances/problemInstances/n_30_c_10000000000_g_6_f_0.3_eps_0_s_100/test.in',
 'instances/problemInstances/n_30_c_10000000000_g_6_f_0.3_eps_0_s_300/test.in',
 'instances/problemInstances/n_30_c_10000000000_g_6_f_0.3_eps_1e-05_s_300/test.in',
 'instances/problemInstances/n_30_c_10000000000_g_6_f_0.3_eps_1e-05_s_100/test.in',
 'instances/problemInstances/n_30_c_10000000000_g_6_f_0.3_eps_1e-05_s_200/test.in',
 'instances/problemInstances/n_30_c_10000000000_g_8_f_0.3_eps_0_s_100/test.in',
 'instances/problemInstances/n_30_c_10000000000_g_8_f_0.3_eps_0_s_300/test.in',
 'instances/problemInstances/n_30_c_10000000000_g_8_f_0.3_eps_0_s_200/test.in',
 'instances/problemInstances/n_30_c_10000000000_g_8_f_0.3_eps_1e-05_s_100/test.in',
 'instances/problemInstances/n_30_c_10000000000_g_8_f_0.3_eps_1e-05_s_200/test.in',
 'instances/problemInstances/n_30_c_10000000000_g_8_f_0.3_eps_1e-05_s_300/test.in',
 'instances/problemInstances/n_30_c_10000000000_g_10_f_0.3_eps_0_s_100/test.in',
 'instances/problemInstances/n_30_c_10000000000_g_10_f_0.3_eps_1e-05_s_200/test.in',
 'instances/problemInstances/n_30_c_10000000000_g_10_f_0.3_eps_0_s_300/test.in',
 'instances/problemInstances/n_30_c_10000000000_g_10_f_0.3_eps_0_s_200/test.in',
 'instances/problemInstances/n_30_c_10000000000_g_10_f_0.3_eps_1e-05_s_300/test.in',
 'instances/problemInstances/n_30_c_10000000000_g_10_f_0.3_eps_1e-05_s_100/test.in']

reps = 50
times = np.array([[0] * len(file) for i in range(reps)])

for i in range(reps):
    for j in range(len(file)):
        times[i, j] = solve_knapsack(file[j])

for i in range(len(file)):
	print(int(np.mean(times[:, i])))

print("\n\n\n")
for i in range(len(file)):
	print(int(np.std(times[:, i])))










