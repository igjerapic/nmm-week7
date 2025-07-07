from cycler import cycler

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from sklearn.metrics import auc 

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


def calculate_strain_stress_props(strain, stress, tolerence = 0.01, min_strain_elongation = 0.1):
    max_mask = stress > max(stress) - 2 * tolerence
    max_stress = np.mean(stress[max_mask])

    # elongation
    max_stress_idx = np.where(stress == max(stress))
    mask = (stress <= tolerence) & (stress >= -tolerence) & (strain > strain[max_stress_idx])
    elongation = strain[mask][0]

    # fracture energy
    fract_eng = auc(strain, stress)

    return max_stress, elongation, fract_eng



def main():
    file_name = "stressCurve.dat"
    data = np.zeros((3, 4))
    for i in range(1,4):
        #plt.clf()
        fig, ax = plt.subplots()
        for j in range(1,4):
            dir_name = f"deform_film{j}/"

            strain, stress = np.loadtxt(dir_name + file_name).T
            
            max_stress, elongation, fract_eng = calculate_strain_stress_props(strain, stress)
            if i == 1:
                data[j-1] = np.array([j, max_stress, elongation, fract_eng ])

            
            print(f"Film{i}")
            print(f"\t max stress: {max_stress:.2f},  elongation: {elongation:.2f},  Fraction energy: {fract_eng:.2f}")

            alpha = 1.0 if j == i else 0.5
            zorder = 10 if j == i else 2
            ax.plot(strain, stress, label = f"film{j}", alpha = alpha, zorder = zorder)
        if i == 1:
            data = data.T
            df = pd.DataFrame({"film": data[0], "max stress": data[1], "elongation": data[2], "fraction energy": data[3]} )
            df.to_csv("stress-strain_films.txt", index=False, sep='\t')
            df.to_latex("stress-strain_films.tex", 
                        index = False, formatters={"name": str.upper}, 
                        float_format="{:.2f}".format,)
        
        ax.set_xlabel("Strain")
        ax.set_ylabel("Stress ($\epsilon / \sigma ^3$)")
        # ax.legend()
        plt.tight_layout()
        plt.savefig(f"stress-strain_film{i}.png", dpi=300)
        #plt.show()


if __name__ == "__main__":
    main()