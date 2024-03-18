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

p_cuts = (1e5, 20e5, 40e5, 60e5, 90e5, 100e5)
p_lst = np.linspace(1e5, max(pc_lst), 100)

pmax = min(pc_lst) - 1e5

class TxyEquilibrium:
    def __init__(self, x, y, T):
        self.x = np.array(x)
        self.y = np.array(y)
        self.T = np.array(T)

def get_binary_txy(p):
    if p < pmax:
        return eos.get_binary_txy(p)

    T_range = np.linspace(min(Tc_lst) - 10, max(Tc_lst) + 10, 200)
    x_lst = []
    y_lst = []
    T_lst = []
    z_lst = np.linspace(0, 1, 50, endpoint=True)
    for T in T_range:
        for z in z_lst:
            flsh = eos.two_phase_tpflash(T, p, [z, 1 - z])
            if flsh.phase == eos.TWOPH:
                x_lst.append(flsh.x[0])
                y_lst.append(flsh.y[0])
                T_lst.append(T)
                break
    if len(T_lst) > 0:
        y_lst = [x_lst[0]] + y_lst
        x_lst = [y_lst[1]] + x_lst
        y_lst.append(x_lst[-1])
        x_lst.append(y_lst[-2])
        T_lst.append(T_lst[-1])
        T_lst = [T_lst[0]] + T_lst
    return None, TxyEquilibrium(x_lst, y_lst, T_lst), None

pcmap = NormedCmap('viridis', p_lst)

def plot_cutouts():
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection='3d')

    ax.scatter(Tc_lst, pc_lst / 1e5, zc_lst, color='r')
    have_cuts = []
    for p in p_lst:
        alpha = 0.2
        for pi in p_cuts:
            if pi in have_cuts:
                continue
            if abs((p - pi) / 1e5) < 5:
                alpha = 1
                have_cuts.append(pi)
                break

        _, lve, _ = get_binary_txy(p)
        if lve.x is not None:
            plt.plot(lve.T, [p / 1e5 for _ in lve.T], lve.x, color=pcmap(p), alpha=alpha)
            plt.plot(lve.T, [p / 1e5 for _ in lve.T], lve.y, color=pcmap(p), alpha=alpha)

    ax.set_xlabel(r'$T$ [K]')
    ax.set_ylabel(r'$p$ [bar]')
    ax.set_zlabel(r'$z, x, y$ [-]')
    plt.tight_layout()
    plt.show()

def plot_cuts():
    plt.plot(zc_lst, Tc_lst, color='r')
    for p in list(p_cuts)[::-1]:
        _, lve, _ = get_binary_txy(p)
        if lve.x is not None:
            plt.plot(lve.x, lve.T, color=pcmap(p), label=int(p / 1e5))
            plt.plot(lve.y, lve.T, color=pcmap(p))


        plt.legend(title=r'$p$ [bar]')

    plt.xlabel(r'$x, y$ [-]')
    plt.ylabel(r'$T$ [K]')
    plt.tight_layout()
    plt.xlim(0, 1)
    plt.savefig(f'p_cut_all.pdf')
    plt.show()

plot_cutouts()
# plot_cuts()