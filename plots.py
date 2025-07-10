import matplotlib.pyplot as plt
import os
import pandas as pd
import seaborn as sns
import matplotlib

sns.set_theme()

sns.set_theme()
font = {'family': 'serif', 'size': 16}
matplotlib.rc('font', **font)

plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = "Times New Roman"
plt.rcParams['mathtext.fontset'] = 'stix'

df = pd.DataFrame({
    "n": [],
    "impl": [],
    "p": [],
    "approx": [],
    "succ": [],
    "triv": []
})


def extrect_data(df, n, i, p, implementation, optimizer):
    d = [x for x in os.listdir(f"{p}/{i}/{implementation}/") if x not in ['resources', 'backup', "results"]]
    if not d: return df

    p_list = [int(x.split("_")[1]) for x in d]
    for ps in p_list:
        try:
            f = open(f"{p}/{i}/{implementation}/p_{int(ps)}/{optimizer}/results.txt").read().split()

            if float(f[-3]) == 1.:
                triv = True
            else:
                triv = False

            new = {"n": n, "impl": implementation, "p": ps, "approx": float(f[-2]), "succ": float(f[-1]), "triv": triv}
            df = pd.concat([df, pd.DataFrame([new])], ignore_index=True)
        except:
            df = df
    return df


optim = "powell"
n_qtg = []
n_copula = []
rat = []
p = "instances"
for i in [i for i in os.listdir(p) if i != ".DS_Store"]:
    dirs = os.listdir(f"{p}/{i}")
    n = int(i.split("_")[1])
    if "qtg" in dirs:
        df = extrect_data(df, n, i, p, "qtg", optim)
    if "copula" in dirs:
        df = extrect_data(df, n, i, p, "copula", optim)

df = df[(df["impl"] == "copula") | (df["impl"] == "qtg")]


sns.set_palette("tab20", 11)
def plot_apprrox(p, implementation):
    sns.lineplot(df[(df["p"] == p) & (df["impl"] == implementation)], x="n", y="approx", label=f"{implementation} $p={p}$")

def plot_succ(p, implementation):
    sns.lineplot(df[(df["p"] == p) & (df["impl"] == implementation)], x="n", y="succ", label=f"{implementation} $p={p}$")

plot_apprrox(1, "qtg")
for i in range(1, 11):
    plot_apprrox(i, "copula")
plt.legend()
plt.ylabel("$\left< f\\left(x\\right)\\right> / f_{OPT}$")
plt.xlabel("$n$")
plt.tight_layout()
plt.savefig(f"approximation_ratio_{optim}.pdf")
plt.show()

plot_succ(1, "qtg")
for i in range(1, 11):
    plot_succ(i, "copula")
plt.legend()
plt.ylabel("$P\\left(f\\left(x\\right) > f_{VG}\\right)$")
plt.tight_layout()
plt.savefig(f"better_than_greedy_{optim}.pdf")
plt.show()
