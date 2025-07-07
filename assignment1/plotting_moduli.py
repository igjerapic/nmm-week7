import numpy as np
import matplotlib.pyplot as plt
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
    file = "elastic_moduli.txt"
    ORIGINAL_EPSILON = 2.1683

    epsilons, moduli = np.loadtxt(file).T

    
    mask = epsilons != ORIGINAL_EPSILON 
    og_idx = np.where(epsilons == ORIGINAL_EPSILON)[0][0]
    print(moduli[og_idx])
    plt.plot(epsilons[mask], moduli[mask], "o", label = "Additional points")
    plt.scatter(ORIGINAL_EPSILON, moduli[og_idx], marker ="o", label = r"Original $\epsilon$", color = "C1")
    plt.xlabel("Interaction Engery (eV)")
    plt.ylabel("Bulk Modulus (GPa)")
    plt.legend()
    plt.tight_layout()
    plt.savefig("bulk-modulus_vs_epsilons.png", dpi=300)
    #plt.show()

if __name__ == "__main__":
    main()