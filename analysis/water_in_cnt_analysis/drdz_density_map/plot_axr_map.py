import numpy as np
import matplotlib.pyplot as plt
import argparse

def plot_from_csv(csv_name, figure_title):
    density = np.loadtxt(csv_name, delimiter=",")

    z_length = 14 
    r_max = 6.0

    plt.imshow(density.T, extent=(-z_length/2, +z_length/2, 0, r_max), aspect='auto', origin='lower', cmap='viridis')
    plt.colorbar(label='Density (N/Å³)')
    plt.xlabel('Axial Distance (Å)')
    plt.ylabel('Radial Distance (Å)')
    plt.title(figure_title)
    plt.tight_layout()
    plt.savefig(f"{figure_title}.svg", transparent=True)
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot water density from a CSV file.")
    parser.add_argument("csv_name", help="Name of the CSV file containing the density data.")
    parser.add_argument("figure_title", help="Title for the generated figure.")
    
    args = parser.parse_args()

    plot_from_csv(args.csv_name, args.figure_title)

