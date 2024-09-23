import matplotlib.pyplot as plt
import os
import numpy as np

n_qtg = []
n_copula = []
rat = []
p = "instances"
for i in os.listdir(p):
    dirs = os.listdir(f"{p}/{i}")
    if "qtg" in dirs:
        f = open(f"{p}/{i}/qtg/results/global_results.txt").read().split("\n")[2:]
        for j in f:
            sp = j.split()
            if "Powell" in sp and "50" in sp:
                n_qtg.append((int(i.split("_")[1]), float(sp[-1])))
    if "copula" in dirs:
        f = open(f"{p}/{i}/copula/results/global_results.txt").read().split("\n")[2:]
        for j in f:
            sp = j.split()
            if "Powell" in sp and "50" in sp:
                n_copula.append((int(i.split("_")[1]), float(sp[-1])))


n_qtg.sort()
n_copula.sort()


plt.plot([i for i, _ in n_qtg], [i for _, i in n_qtg], ".--", label="QTG")
plt.plot([i for i, _ in n_copula], [i for _, i in n_copula], ".--", label="Copula $p=1$")
plt.ylabel("Approximation Ratio")
plt.xlabel("$n$")
plt.legend()
plt.tight_layout()
plt.savefig("approx_ratio.pdf")
plt.show()