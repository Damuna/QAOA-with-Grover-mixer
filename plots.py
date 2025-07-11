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
    l = "QTG"
    if implementation == "copula":
        l = "Copula"
    sns.lineplot(df[(df["p"] == p) & (df["impl"] == implementation) & (df["succ"] != 0)], x="n", y="approx", label=f"{l} $q={p}$")

def plot_succ(p, implementation):
    l = "QTG"
    if implementation == "copula":
        l = "Copula"
    else: print(len(df[(df["p"] == p) & (df["impl"] == implementation) & (df["succ"] == 0)]))
    
    sns.lineplot(df[(df["p"] == p) & (df["impl"] == implementation) & (df["succ"] != 0)], x="n", y="succ", label=f"{l} $q={p}$")

# plot_apprrox(1, "qtg")
# for i in range(1, 11):
#     plot_apprrox(i, "copula")
# plt.legend()
# plt.ylabel("$\left< f\\right> / f_{OPT}$")
# plt.xlabel("$n$")
# plt.tight_layout()
# plt.savefig(f"approximation_ratio_{optim}.pdf")
# plt.show()
#
# plot_succ(1, "qtg")
# for i in range(1, 11):
#     plot_succ(i, "copula")
# plt.legend()
# plt.ylabel("$P\\left(f > f_{VG}\\right)$")
# plt.xlabel("$n$")
# # plt.yscale("log")
# plt.tight_layout()
# plt.savefig(f"better_than_greedy_{optim}.pdf")
# plt.show()



# plot of the resources
for impl in ["qtg", "copula"]:
    qtg_res = [[] for _ in range(11)]
    # ns = [[] for _ in range(11)]
    i_ = os.listdir("instances")
    try:
        i_.remove(".DS_Store")
    except: pass
    i_.sort()
    for i in i_:
        try:
            d_ = [x for x in os.listdir(f"{p}/{i}/{impl}/") if x not in ['resources', 'backup', "results"]]
            # print(d_)
            if len(d_) > 0:
                for d in d_:
                    p_ = int(d.split("_")[1])
                    if impl == "qtg" and p_ == 1 or impl == "copula" and int(i.split("_")[1]) <= 20:
                        f = int(open(f"{p}/{i}/{impl}/{d}/resources").read().split()[1])
                        # ns[p_].append(int(i.split("_")[1]))
                        qtg_res[p_].append((int(i.split("_")[1]), f))
        except: pass

    lab = "QTG"
    if impl == "copula": lab = "Copula"
    for i in range(11):
        qtg_res[i].sort(key = lambda x: x[0])
        if len(qtg_res[i]) != 0:
            a, b = [i for i, j in qtg_res[i]], [j for i, j in qtg_res[i]]
            plt.plot(a, b, ".--", label=f"{lab} $q={i}$")

plt.yscale("log")
plt.xlabel("$n$")
plt.ylabel("Cycle count")
plt.legend()
plt.tight_layout()
plt.savefig("cycle_counts.pdf")
plt.show()