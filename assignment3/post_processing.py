#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  2 18:19:55 2025

@author: utku
"""

import beadspring as bsa
import matplotlib.pyplot as plt
import MDAnalysis as mda
import numpy as np

from cycler import cycler


plt.style.use('../scripts/default.mplstyle')

plt.rcParams['axes.prop_cycle'] = plt.cycler(cycler(color = ['#332288', 
                                    '#CC6677',
                                    '#88CCEE',
                                    '#DDCC77', 
                                    '#117733', 
                                    '#882255', 
                                    '#44AA99', 
                                    '#999933', 
                                    '#AA4499',
                                    '#DDDDDD'
                                ]))



def main():
    u, u2 = setup(topology="final.data", trajectory="productionLog.dat")
    shells, selection_indices, L = define_shells(u2)
    time, positions = fetch_positions(u, shells, selection_indices)
    msds = compute_shell_msds(positions)
    TAU_P = 1.0
    dwfs, _ = compute_shell_dwf(time, msds, tau_p=TAU_P)
    
    plot_msds(time, msds, tau_p=TAU_P)
    plt.clf()
    print(shells)

    plot_dwfs(dwfs)
    plt.clf()

    r_mid, density = compute_density(u2, L)

    plot_density(r_mid, density)


def setup(topology, trajectory):
    u = bsa.setup_universe(topology, trajectory)
    u2 = mda.Universe(topology)
    return u, u2

def rg_filler(u2):
    filler = u2.select_atoms("type 1")
    _, eigvals = bsa.compute_gyration_tensor(filler.positions.copy())
    rg = np.sum(eigvals) ** 0.5
    print("filler R_g:", rg)

def define_shells(u2, num_shells=10):
    filler = u2.select_atoms("type 1")
    rbound, _ = filler.bsphere()
    L = u2.dimensions[0]
    dr = np.linspace(rbound, L / 2, num_shells)
    shell_bounds = list(zip(dr[:-1], dr[1:]))
    shells = [f"{r0:.2f}-{r1:.2f}" for r0, r1 in shell_bounds]
    indices = [
        u2.select_atoms(f"sphlayer {r0} {r1} type 1").indices for r0, r1 in shell_bounds
    ]
    return shells, indices, L


def fetch_positions(u, shells, selection_indices):
    positions = {shell: [] for shell in shells}
    time = np.zeros(len(u.trajectory))

    for ts, frame in enumerate(u.trajectory):
        time[ts] = frame.data["time"]
        for shell, indices in zip(shells, selection_indices):
            positions[shell].append(u.atoms[indices].positions.copy())

    return time, positions


def compute_shell_msds(positions):
    return {shell: bsa.compute_msd(pos) for shell, pos in positions.items()}

def compute_shell_dwf(time, msds, tau_p=1.0):
    tau_p = np.where(time[1:] == tau_p)[0][0]
    return {shell: msd[tau_p] for shell, msd in msds.items()}, np.array([msd[tau_p] for msd in msds.values()])


def plot_msds(time, msds, tau_p = 1.0):
    for i, (shell, msd) in enumerate(msds.items()):
        plt.loglog(time[1:], msd, "o-", label=f"Shell {i+1}")
    plt.vlines(tau_p, *plt.ylim(), linestyles=":" , colors = 'k')
    plt.legend(ncols = 2)
    plt.xlabel("Time")
    plt.ylabel("MSD")
    #plt.title("Mean Square Displacement by Shell")
    plt.tight_layout()
    plt.savefig("msd.png", dpi=300)
    #plt.show()
def compute_dwf_boundary(shells, dwfs, start_idx = 0):
    bulk_dwf = np.mean(dwfs[start_idx:])
    boundary_idx = np.where(dwfs >= 0.64 * bulk_dwf)[0][0]
    boundary = shells[boundary_idx]
    return boundary

def plot_dwfs(dwfs):
    shells = np.array([np.mean(np.array(shell.split('-'), dtype=float)) for shell in dwfs.keys()])
    dwf = [val for val in dwfs.values()]
    

    boundary = compute_dwf_boundary(shells, dwf, start_idx = 1)
    print(boundary)    
    plt.plot(shells, dwf, "o-")
    y_lims = plt.ylim()
    plt.vlines(boundary, *y_lims, linestyles=":", colors='k')
    plt.ylim(y_lims)
    xticks = np.arange(5.0, 21, 2.5)
    plt.xlabel("Radius")
    # plt.ylabel(r"$\langle u^2\rangle$")
    plt.ylabel("Debye-Waller Factor")
    plt.xticks(xticks, [f"{x:.1f}" if int(x) != x else int(x) for x in xticks])
    plt.tight_layout()
    plt.savefig("dwf.png", dpi = 300)
    #plt.show()

def compute_density(u2, L, num_points=50):
    dr = np.linspace(0, L / 2, num_points)
    deltar = np.diff(dr)[0]
    r_pairs = zip(dr[:-1], dr[1:])
    r_mid = 0.5 * (dr[:-1] + dr[1:])

    density = [
        u2.select_atoms(f"sphlayer {r1} {r2} type 1").n_atoms
        / ((4 * np.pi / 3) * ((r1 + deltar) ** 3 - r1**3))
        for r1, r2 in r_pairs
    ]
    return r_mid, density

def compute_struct_boundary(r_mid, density):
    # determining bound layer 
    sorted_densities = np.sort(density)
    print(sorted_densities)
    boundary_density = sorted_densities[sorted_densities <= 0.6 * max(density)][-1]

    boundary_idx = np.where(density == boundary_density)[0][0]
    return r_mid[boundary_idx]

def plot_density(r_mid, density):

    struct_boundary = compute_struct_boundary(r_mid, density)
    print(struct_boundary)
    
    plt.plot(r_mid, density, "o-")
    ylim = plt.ylim()
    plt.xlabel("Radius")
    plt.ylabel("Radial Density")
    plt.vlines(struct_boundary, *ylim, linestyles=":", colors='k' )
    plt.ylim(ylim)
    #plt.title("Radial Density Profile")
    plt.tight_layout()
    plt.savefig("radial_density.png", dpi=300)
    #plt.show()


if __name__ == "__main__":
    main()
