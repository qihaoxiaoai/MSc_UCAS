import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
import argparse

def plot_water_density_around_CNT_dp(u, figure_title, csv_name):
    o_atoms = u.select_atoms('type OT')
    cnt_c_atoms = u.select_atoms('resname CNT and name C')

    z_mean = cnt_c_atoms.positions[:, 2].mean()
    cnt_up_h_atoms = u.select_atoms(f'resname CNT and name H and prop z > {z_mean}')
    cnt_down_h_atoms = u.select_atoms(f'resname CNT and name H and prop z <= {z_mean}')

    z_length = cnt_up_h_atoms.positions[:, 2].mean() - cnt_down_h_atoms.positions[:, 2].mean()
    z_center_z = (cnt_up_h_atoms.positions[:, 2].mean() + cnt_down_h_atoms.positions[:, 2].mean()) / 2
    z_up = z_center_z + z_length / 2
    z_down = z_center_z - z_length / 2

    r_max = 6.0
    dr = 0.1
    dz = 0.1

    r_bins = np.arange(0, r_max + dr, dr)
    z_bins = np.arange(z_down, z_up + dz, dz)

    density = np.zeros((len(z_bins) - 1, len(r_bins) - 1))
    cnt_center_xy = cnt_c_atoms.positions[:, :2].mean(axis=0)

    for ts in u.trajectory[5:30]:
        o_positions = o_atoms.positions.copy()
        o_positions[:, :2] -= cnt_center_xy

        radial_dists = np.sqrt(o_positions[:, 0]**2 + o_positions[:, 1]**2)
        axial_dists = o_positions[:, 2]

        for i, z in enumerate(z_bins[:-1]):
            z_next = z_bins[i + 1]
            axial_mask = (axial_dists >= z) & (axial_dists < z_next)

            for j, r in enumerate(r_bins[:-1]):
                r_next = r_bins[j + 1]
                radial_mask = (radial_dists >= r) & (radial_dists < r_next)

                mask = axial_mask & radial_mask
                ring_volume = np.pi * (r_next**2 - r**2) * dz
                density[i, j] += (mask.sum() / ring_volume)

    density /= len(u.trajectory)
    np.savetxt(csv_name, density, delimiter=",")
    
    plt.imshow(density.T, extent=(-z_length/2, +z_length/2, 0, r_max), aspect='auto', origin='lower', cmap='viridis')
    plt.colorbar(label='Density (N/Å³)')
    plt.xlabel('Axial Distance (Å)')
    plt.ylabel('Radial Distance (Å)')
    plt.title(figure_title)
    plt.tight_layout()
    plt.savefig(f"{figure_title}.svg", transparent=True)
    plt.show()

    return density

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot water density around CNT.")
    parser.add_argument("trjtpr", help="Path to the TPR file.")
    parser.add_argument("trjxtc", help="Path to the XTC file.")
    parser.add_argument("csv_name", help="Name of the CSV file to save the density data.")
    parser.add_argument("figure_title", help="Title for the generated figure.")
    
    args = parser.parse_args()

    u = mda.Universe(args.trjtpr, args.trjxtc)
    plot_water_density_around_CNT_dp(u, args.figure_title, args.csv_name)

