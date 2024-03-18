import matplotlib.pyplot as plt
import numpy as np
from thermopack.cubic import SoaveRedlichKwong
from plottools.cyclers import NormedCmap

eos = SoaveRedlichKwong('C1,C3')

zc_lst = np.linspace(0, 1, 100, endpoint=True)
Tc_lst = np.empty_like(zc_lst)
pc_lst = np.empty_like(zc_lst)
for i, z in enumerate(zc_lst):
    Tc_lst[i], _, pc_lst[i] = eos.critical([z, 1 - z])

T_cuts = (150, 180, 200, 250, 300, 350, 365)
T_lst = np.linspace(150, 400, 100)
tcmap = NormedCmap('cool', T_lst)

def plot_cutouts():
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection='3d')

    ax.scatter(Tc_lst, pc_lst / 1e5, zc_lst, color='r')
    tcmap = NormedCmap('cool', T_lst)
    has_cuts = []
    for T in T_lst:
        alpha = 0.2
        for Ti in T_cuts:
            if abs(T - Ti) < 10:
                if Ti in has_cuts:
                    break
                alpha = 1
                has_cuts.append(Ti)
                break

        _, lve, _ = eos.get_binary_pxy(T)
        if lve.x is not None:
            if len(lve.x) == 0:
                continue
            ax.plot([T for _ in lve.x], lve.p / 1e5, lve.x, color=tcmap(T), alpha=alpha)
            ax.plot([T for _ in lve.x], lve.p / 1e5, lve.y, color=tcmap(T), alpha=alpha)

    ax.set_xlabel(r'$T$ [K]')
    ax.set_ylabel(r'$p$ [bar]')
    ax.set_zlabel(r'$z, x, y$ [-]')
    plt.tight_layout()
    plt.show()

def plot_cuts():
    plt.plot(zc_lst, pc_lst / 1e5, color='r')
    for T in T_cuts[::-1]:
        _, lve, _ = eos.get_binary_pxy(T, minimum_pressure=1e3)
        if lve.x is not None:
            plt.plot(lve.x, lve.p / 1e5, color=tcmap(T), label=int(T))
            plt.plot(lve.y, lve.p / 1e5, color=tcmap(T))

        #
        plt.legend(title=r'$T$ [K]')

    plt.xlabel(r'$x, y$ [-]')
    plt.ylabel(r'$p$ [bar]')
    plt.xlim(0, 1)
    plt.savefig(f'T_cut_all.pdf')
    plt.show()

plot_cutouts()
# plot_cuts()