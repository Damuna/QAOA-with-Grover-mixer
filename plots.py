import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd
import seaborn as sns

sns.set_theme()

df = pd.DataFrame({
    "n": [],
    "impl": [],
    "p": [],
    "approx": [],
    "succ": [],
    "triv": []
})


def extrect_data(df, n, i, p, implementation):
    d = [x for x in os.listdir(f"{p}/{i}/{implementation}/") if x not in ['resources', 'backup', "results"]]
    if not d: return df

    p_list = [int(x.split("_")[1]) for x in d]
    for ps in p_list:
        try:
            f = open(f"{p}/{i}/{implementation}/p_{int(ps)}/bfgs/results.txt").read().split()

            if float(f[-3]) == 1.:
                triv = True
            else:
                triv = False

            new = {"n": n, "impl": implementation, "p": ps, "approx": float(f[-2]), "succ": float(f[-1]), "triv": triv}
            df = pd.concat([df, pd.DataFrame([new])], ignore_index=True)
        except:
            df = df
    return df


n_qtg = []
n_copula = []
rat = []
p = "instances"
for i in os.listdir(p):
    dirs = os.listdir(f"{p}/{i}")
    n = int(i.split("_")[1])
    if "qtg" in dirs:
        df = extrect_data(df, n, i, p, "qtg")
    if "copula" in dirs:
        df = extrect_data(df, n, i, p, "copula")

# print(df)
# print(df[df["n"] == 15][["n", "p", "impl"]])


sns.set_palette("tab20", 11)
def plot_approx(p, implementation):
    sns.lineplot(df[(df["p"] == p) & (df["impl"] == implementation)], x="n", y="approx", label=f"{implementation} $p={p}$")

def plot_succ(p, implementation):
    sns.lineplot(df[(df["p"] == p) & (df["impl"] == implementation)], x="n", y="succ", label=f"{implementation} $p={p}$")

plot_approx(1, "qtg")
for i in range(1, 11):
    plot_approx(i, "copula")
plt.legend()
plt.ylabel("$\left< P_{QAOA}\\right> / P_{OPT}$")
plt.xlabel("$n$")
plt.tight_layout()
plt.savefig("approximation_ratio.pdf")
plt.show()

# plot_succ(1, "qtg")
# for i in range(1, 11):
#     plot_succ(i, "copula")
# plt.yscale("log")
# plt.legend()
# plt.tight_layout()
# plt.show()
