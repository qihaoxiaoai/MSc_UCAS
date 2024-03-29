import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

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
    plt.figure(figsize=(10, 6))
    sns.histplot(df_output["Orient_Prob"], kde=True, bins=30, color="green", stat="probability")
    plt.xlabel("Water Orientation Consistency in CNT", fontsize=14)
    plt.ylabel("Probability Density (%)", fontsize=14)
    plt.grid(True)
    plt.tight_layout()
    y_ticks = plt.gca().get_yticks()
    plt.gca().set_yticklabels(['{:.0f}%'.format(y*100) for y in y_ticks])
    plt.savefig('P4_orientation_consistency.svg',transparent=True)

def plot_bond_length_distribution(df_output):
    plt.figure(figsize=(10, 6))
    sns.histplot(df_output["bond_mean"], kde=True, bins=100, color="blue", stat="probability")
    plt.xlabel("Bond length of water in CNT", fontsize=14)
    plt.ylabel("Probability Density (%)", fontsize=14)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig('P5_incnt_water_bond_length_distribution.svg',transparent=True)

def plot_angle_distribution(df_output):
    plt.figure(figsize=(10, 6))
    sns.histplot(df_output["hohang_mean"], kde=True, bins=100, color="red", stat="probability")
    plt.xlabel("Angle of water in CNT", fontsize=14)
    plt.ylabel("Probability Density (%)", fontsize=14)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig('P6_incnt_water_angle_distribution.svg',transparent=True)

if __name__ == "__main__":
    df_output = pd.read_csv('02_water_in_cnt_info.csv')
    plot_water_count_vs_time(df_output)
    plot_angle_bisector_vs_time(df_output)
    time_plot(df_output)
    plot_orientation_consistency(df_output)
    plot_bond_length_distribution(df_output)
    plot_angle_distribution(df_output)
    print("Successful!")
