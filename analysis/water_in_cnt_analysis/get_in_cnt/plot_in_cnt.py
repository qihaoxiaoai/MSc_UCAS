import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde
from scipy.signal import find_peaks

def plot_histogram(data, ax, xlabel, title, color):
    ax.hist(data, bins=40, density=True, color=color, alpha=0.7, label='Histogram')
    ax.set_xlabel(xlabel, fontsize=14)
    ax.set_ylabel("Probability", fontsize=14)
    ax.set_title(title, fontsize=16)
    ax.grid(True)
    ax.legend()

def save_data_to_csv(data, filename="data.csv"):
    data.to_csv(filename, index=False)

def plot_histogram_and_density(data, ax, xlabel, title, color):
    ax.hist(data, bins='auto', density=True, color=color, alpha=0.7, label='Histogram')
    density = gaussian_kde(data)
    xs = np.linspace(min(data), max(data), 1000)
    ys = density(xs)
    ax.plot(xs, ys, color=color, label='Density')
    ax.set_xlabel(xlabel, fontsize=14)
    ax.set_ylabel("Probability", fontsize=14)
    ax.set_title(title, fontsize=16)
    ax.grid(True)
    ax.legend()
    peaks, _ = find_peaks(ys)
    peak_values = [(xs[p], ys[p]) for p in peaks]
    for x, y in peak_values:
        ax.annotate(f'{x:.2f}', (x, y), textcoords="offset points", xytext=(0,10), ha='center')
    return peak_values

def plot_histogram_and_density_for_CNT(data, ax, xlabel, title, color):
    counts, bins, patches = ax.hist(data, bins=np.arange(data.min(), data.max() + 2) - 0.5, density=True, color=color, alpha=0.7, label='Histogram')
    bin_centers = 0.5 * (bins[:-1] + bins[1:])
    ax.plot(bin_centers, counts, color=color, label='Density', marker='o', linestyle='-')
    ax.set_xlabel(xlabel, fontsize=14)
    ax.set_ylabel("Probability", fontsize=14)
    ax.set_title(title, fontsize=16)
    ax.grid(True)
    ax.legend()
    for x, y in zip(bin_centers, counts):
        ax.annotate(f'{x:.0f}', (x, y), textcoords="offset points", xytext=(0,10), ha='center')
    return list(zip(bin_centers, counts))

def save_peak_values_to_file(peak_values_dict, filename="peak_values.txt"):
    with open(filename, "w") as f:
        for key, values in peak_values_dict.items():
            f.write(f"{key}:\n")
            for x, y in values:
                f.write(f"  Peak at x={x:.2f} with y={y:.2f}\n")
            f.write("\n")

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
    axs[0].plot(df_output["Frame"] * 0.5, df_output["Number_in_CNT"], marker='o', linestyle='-')
    axs[0].set_ylabel('Number of Water Molecules in CNT',fontsize=14)
    axs[0].set_title('Number of Water Molecules in CNT vs Time',fontsize=16)
    axs[0].grid(True)
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

def plot_all_distributions_together(df_output):
    fig, axs = plt.subplots(2, 2, figsize=(12, 8))
    plot_histogram(df_output["Orient_Prob"], axs[0, 0], "Orientation Consistency", "Orientation Consistency Distribution", "green")
    plot_histogram_and_density(df_output["bond_mean"], axs[0, 1], "Bond Length", "Bond Length Distribution", "blue")
    plot_histogram_and_density(df_output["hohang_mean"], axs[1, 0], "Angle", "Angle Distribution", "red")
    data_for_csv = plot_histogram_and_density_for_CNT(df_output["Number_in_CNT"], axs[1, 1], "Number in CNT", "Number in CNT Distribution", "purple")
    plt.tight_layout()
    plt.savefig('P8_all_distributions_together.svg',transparent=True)
    df_data_for_csv = pd.DataFrame(data_for_csv, columns=["Bin_Center", "Count"])
    save_data_to_csv(df_data_for_csv, "Number_in_CNT_data.csv")

if __name__ == "__main__":
    df_output = pd.read_csv('02_water_in_cnt_info.csv')
    plot_water_count_vs_time(df_output)
    plot_angle_bisector_vs_time(df_output)
    time_plot(df_output)
    plot_all_distributions_together(df_output)
    fig, ax = plt.subplots(figsize=(6, 4))
    plot_histogram(df_output["Orient_Prob"], ax, "Orientation Consistency", "Orientation Consistency Distribution", "green")
    plt.tight_layout()
    plt.savefig('P4_Orientation_Consistency.svg', transparent=True)
    print("Successful!")

