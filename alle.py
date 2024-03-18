import matplotlib.pyplot as plt
import numpy as np
from thermopack.cubic import SoaveRedlichKwong
from plottools.cyclers import NormedCmap

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, projection='3d')
# ax2 = fig.add_subplot(1, 2, 2)

eos = SoaveRedlichKwong('C1,C3')

zc_lst = np.linspace(0, 1, 100, endpoint=True)
Tc_lst = np.empty_like(zc_lst)
pc_lst = np.empty_like(zc_lst)
for i, z in enumerate(zc_lst):
    Tc_lst[i], _, pc_lst[i] = eos.critical([z, 1 - z])

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

def plot_z_ends():
    z_lst = (0, 1)
    zcmap = NormedCmap('winter', z_lst)
    for z in z_lst:
        T, p = eos.get_envelope_twophase(1e5, [z, 1 - z], minimum_temperature=120, maximum_pressure=1e9)
        ax.plot(T, p / 1e5, [z for _ in T], color=zcmap(z))

def plot_z():
    # ax.scatter(Tc_lst, pc_lst / 1e5, zc_lst, color='r', marker='.')
    z_lst = (0, 0.25, 0.5, 0.75, 0.95, 1)
    zcmap = NormedCmap('winter', z_lst)
    for z in z_lst:
        T, p = eos.get_envelope_twophase(1e5, [z, 1 - z], minimum_temperature=150)
        ax.plot(T, p / 1e5, [z for _ in T], color=zcmap(z))

def plot_T():
    T_lst = (150, 180, 200, 250, 300, 350, 365)
    tcmap = NormedCmap('cool', T_lst)
    for T in T_lst:
        _, lve, _ = eos.get_binary_pxy(T, minimum_pressure=1e2)
        if lve.x is not None:
            ax.plot([T for _ in lve.x], lve.p / 1e5, lve.x, color=tcmap(T))
            ax.plot([T for _ in lve.x], lve.p / 1e5, lve.y, color=tcmap(T))

    plot_z_ends()

def plot_p():
    # ax.scatter(Tc_lst, pc_lst / 1e5, zc_lst, color='r')
    p_lst = (1e5, 20e5, 40e5, 60e5, 90e5, 100e5)
    pcmap = NormedCmap('viridis', p_lst)
    for p in p_lst:
        _, lve, _ = get_binary_txy(p)
        if lve.x is not None:
            plt.plot(lve.T, [p / 1e5 for _ in lve.T], lve.x, color=pcmap(p))
            plt.plot(lve.T, [p / 1e5 for _ in lve.T], lve.y, color=pcmap(p))

    plot_z_ends()

# plot_p()
plot_T()
# plot_z()
ax.set_xlabel(r'$T$ [K]')
ax.set_ylabel(r'$p$ [bar]')
ax.set_zlabel(r'$z, x, y$ [-]')
plt.tight_layout()
plt.show()