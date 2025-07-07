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


def main():
    u, u2 = setup(topology="final.data", trajectory="productionLog.dat")
    shells, selection_indices, L = define_shells(u2)
    time, positions = fetch_positions(u, shells, selection_indices)
    msds = compute_shell_msds(positions)

    plot_msds(time, msds)

    r_mid, density = compute_density(u2, L)
    plot_density(r_mid, density)


def setup(topology, trajectory):
    u = bsa.setup_universe(topology, trajectory)
    u2 = mda.Universe(topology)
    return u, u2


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


def plot_msds(time, msds):
    for i, (shell, msd) in enumerate(msds.items()):
        plt.loglog(time[1:], msd, "o-", label=f"Shell {i+1}")
    plt.legend()
    plt.xlabel("Time")
    plt.ylabel("MSD")
    plt.title("Mean Square Displacement by Shell")
    plt.show()


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


def plot_density(r_mid, density):
    plt.plot(r_mid, density, "o-")
    plt.xlabel("Radius")
    plt.ylabel("Density")
    plt.title("Radial Density Profile")
    plt.show()


if __name__ == "__main__":
    main()
