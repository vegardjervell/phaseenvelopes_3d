import matplotlib.pyplot as plt
import numpy as np
from thermopack.cubic import SoaveRedlichKwong
from plottools.cyclers import NormedCmap

eos = SoaveRedlichKwong('C1,C3')

z_cuts = (0, 0.25, 0.5, 0.75, 0.95, 1)
z_lst = np.linspace(0, 1, 500, endpoint=True)
print(z_lst)
zcmap = NormedCmap('winter', z_lst)

zc_lst = np.linspace(0, 1, 100, endpoint=True)
Tc_lst = np.empty_like(zc_lst)
pc_lst = np.empty_like(zc_lst)
for i, z in enumerate(zc_lst):
    Tc_lst[i], _, pc_lst[i] = eos.critical([z, 1 - z])

def plot_cutouts():
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection='3d')

    zc_lst = np.linspace(0, 1, 100, endpoint=True)
    Tc_lst = np.empty_like(zc_lst)
    pc_lst = np.empty_like(zc_lst)
    for i, z in enumerate(zc_lst):
        Tc_lst[i], _, pc_lst[i] = eos.critical([z, 1 - z])

    ax.scatter(Tc_lst, pc_lst / 1e5, zc_lst, color='r')
    has_cuts = []
    for z in z_lst:
        T, p = eos.get_envelope_twophase(1e5, [z, 1 - z])
        alpha = 0.05
        for zc in z_cuts:
            if zc in has_cuts:
                continue
            if abs(z - zc) < 0.01:
                alpha = 1
                print(f'cut {z}, {zc}, {has_cuts}')
                has_cuts.append(zc)
                break

        ax.plot(T, p / 1e5, [z for _ in T], color=zcmap(z), alpha=alpha)

    ax.set_xlabel(r'$T$ [K]')
    ax.set_ylabel(r'$p$ [bar]')
    ax.set_zlabel(r'$z, x, y$ [-]')
    plt.show()

def plot_cuts():
    for z in z_cuts:
        T, p = eos.get_envelope_twophase(1e5, [z, 1 - z])
        plt.plot(T, p / 1e5, color=zcmap(z), label=f'{z:.2f}')
        plt.plot(Tc_lst, pc_lst / 1e5, color='r')
        plt.legend(title=r'$z$ [-]')

    plt.xlabel(r'$T$ [K]')
    plt.ylabel(r'$p$ [bar]')
    plt.savefig(f'z_cut_all.pdf')
    plt.show()

plot_cutouts()
plot_cuts()