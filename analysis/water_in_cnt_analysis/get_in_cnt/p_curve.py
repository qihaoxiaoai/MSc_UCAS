import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde

def plot_smooth_density(data, ax, xlabel, title, color):
    density = gaussian_kde(data)
    xs = np.linspace(min(data), max(data), 1000)
    ax.plot(xs, density(xs), color=color)
    ax.set_xlabel(xlabel, fontsize=14)
    ax.set_title(title, fontsize=16)
    ax.grid(True)

def plot_water_count_vs_time(df_output):
    plt.figure(figsize=(10, 6))
    plt.plot(df_output["Frame"] * 0.5, df_output["Number_in_CNT"], marker='o', linestyle='-')
    plt.xlabel('Time (ps)',fontsize=16)
    plt.ylabel('Number of Water Molecules in CNT',fontsize=14)
    plt.title('Number of Water Molecules in CNT vs Time',fontsize=16)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig('P1_water_count_vs_time.svg',transparent=True)

def plot_angle_bisector_vs_time(df_output):
    time_interval = 0.5  # ps
    times = df_output["Frame"] * time_interval
    plt.figure(figsize=(10, 6))
    plt.plot(times, df_output["bsecz_mean"], label="bsecz_mean", color="blue")
    plt.xlabel("Time (ps)", fontsize=14)
    plt.ylabel("Angle (°)", fontsize=14)
    plt.title("Angle Bisector of Water Molecule with Z-Axis", fontsize=16)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig('P2_angle_bisecz_vs_time.svg',transparent=True)

def time_plot(df_output):
    fig, axs = plt.subplots(2, 1, figsize=(10, 12), sharex=True)

    # Plot Number of Water Molecules in CNT vs Time
    axs[0].plot(df_output["Frame"] * 0.5, df_output["Number_in_CNT"], marker='o', linestyle='-')
    axs[0].set_ylabel('Number of Water Molecules in CNT',fontsize=14)
    axs[0].set_title('Number of Water Molecules in CNT vs Time',fontsize=16)
    axs[0].grid(True)

    # Plot Angle Bisector of Water Molecule with Z-Axis vs Time
    time_interval = 0.5  # ps
    times = df_output["Frame"] * time_interval
    axs[1].plot(times, df_output["bsecz_mean"], label="bsecz_mean", color="blue")
    axs[1].set_xlabel("Time (ps)", fontsize=14)
    axs[1].set_ylabel("Angle (°)", fontsize=14)
    axs[1].set_title("Angle Bisector of Water Molecule with Z-Axis", fontsize=16)
    axs[1].legend()
    axs[1].grid(True)

    plt.tight_layout()
    plt.savefig('P3_all_time_related_plots.svg',transparent=True)

def plot_orientation_consistency(df_output):
    fig, ax = plt.subplots(figsize=(10, 6))
    plot_smooth_density(df_output["Orient_Prob"], ax, "Water Orientation Consistency in CNT", "Orientation Consistency Distribution", "green")
    plt.tight_layout()
    plt.savefig('P4_orientation_consistency.svg',transparent=True)

def plot_bond_length_distribution(df_output):
    fig, ax = plt.subplots(figsize=(10, 6))
    plot_smooth_density(df_output["bond_mean"], ax, "Bond length of water in CNT", "Bond Length Distribution", "blue")
    plt.tight_layout()
    plt.savefig('P5_incnt_water_bond_length_distribution.svg',transparent=True)

def plot_angle_distribution(df_output):
    fig, ax = plt.subplots(figsize=(10, 6))
    plot_smooth_density(df_output["hohang_mean"], ax, "Angle of water in CNT", "Angle Distribution", "red")
    plt.tight_layout()
    plt.savefig('P6_incnt_water_angle_distribution.svg',transparent=True)

def plot_number_in_CNT_distribution(df_output):
    fig, ax = plt.subplots(figsize=(10, 6))
    plot_smooth_density(df_output["Number_in_CNT"], ax, "Number of Water Molecules in CNT", "Number in CNT Distribution", "purple")
    plt.tight_layout()
    plt.savefig('P7_number_in_CNT_distribution.svg',transparent=True)

def plot_all_distributions_together(df_output):
    fig, axs = plt.subplots(2, 2, figsize=(12, 8))
    
    plot_smooth_density(df_output["Orient_Prob"], axs[0, 0], "Orientation Consistency", "Orientation Consistency Distribution", "green")
    plot_smooth_density(df_output["bond_mean"], axs[0, 1], "Bond Length", "Bond Length Distribution", "blue")
    plot_smooth_density(df_output["hohang_mean"], axs[1, 0], "Angle", "Angle Distribution", "red")
    plot_smooth_density(df_output["Number_in_CNT"], axs[1, 1], "Number in CNT", "Number in CNT Distribution", "purple")
    
    plt.tight_layout()
    plt.savefig('P8_all_distributions_together.svg',transparent=True)

if __name__ == "__main__":
    df_output = pd.read_csv('02_water_in_cnt_info.csv')
    plot_water_count_vs_time(df_output)
    plot_angle_bisector_vs_time(df_output)
    time_plot(df_output)
    plot_orientation_consistency(df_output)
    plot_bond_length_distribution(df_output)
    plot_angle_distribution(df_output)
    plot_number_in_CNT_distribution(df_output)
    plot_all_distributions_together(df_output)
    print("Successful!")

