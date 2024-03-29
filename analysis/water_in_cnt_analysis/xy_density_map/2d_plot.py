import MDAnalysis as mda
from MDAnalysis.analysis import density
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from scipy.optimize import leastsq
from scipy.ndimage import gaussian_filter
from matplotlib.lines import Line2D
import os
import sys

data_path=os.getcwd()
outdir =os.getcwd()

if not os.path.exists(outdir):
    os.makedirs(outdir)

os.chdir(outdir)

trjtpr=data_path+ "/cnt.tpr"
trjxtc=data_path+ "/nopbc.xtc"
u = mda.Universe(trjtpr,trjxtc)

# 获取盒子的尺寸
box_dims = u.dimensions
x_box_length = box_dims[0]
y_box_length = box_dims[1]

water_OT = u.select_atoms('type OT')
cnt_c_atoms = u.select_atoms('resname CNT and name C')
z_mean = cnt_c_atoms.positions[:,2].mean()
cnt_up_h_atoms = u.select_atoms(f'resname CNT and name H and prop z > {z_mean}')
cnt_down_h_atoms= u.select_atoms(f'resname CNT and name H and prop z <= {z_mean}')
z_mean_up = cnt_up_h_atoms.positions[:,2].mean()
z_mean_down = cnt_down_h_atoms.positions[:,2].mean()

selected_OT = water_OT.select_atoms(f'prop z > {z_mean_down} and prop z < {z_mean_up}')
D = mda.analysis.density.DensityAnalysis(selected_OT, delta=0.4).run()
average_density = D.results.density.grid.mean(axis=2)

cnt_c_atoms = u.select_atoms('resname CNT and name C')
num_frames = len(u.trajectory)
selected_frames = np.random.choice(num_frames, 50, replace=False)

all_c_xy_coords = []
for frame in selected_frames:
    u.trajectory[frame]
    all_c_xy_coords.extend(cnt_c_atoms.positions[:, :2])
all_c_xy_coords = np.array(all_c_xy_coords)

def calc_R(xc, yc):
    return np.sqrt((all_c_xy_coords[:,0]-xc)**2 + (all_c_xy_coords[:,1]-yc)**2)

def f(c):
    Ri = calc_R(*c)
    return Ri - Ri.mean()

center_estimate = np.mean(all_c_xy_coords[:,0]), np.mean(all_c_xy_coords[:,1])
center, _ = leastsq(f, center_estimate)
xc, yc = center
Ri = calc_R(*center)
R = Ri.mean()

smoothed_density = gaussian_filter(average_density, sigma=0.4)

plt.figure(figsize=(9.8,8))
img = plt.imshow(smoothed_density.T, extent=[D.results.density.edges[0][0], D.results.density.edges[0][-1], D.results.density.edges[1][0], D.results.density.edges[1][-1]], origin='lower', cmap='Blues', aspect='auto')
cbar = plt.colorbar(img, format='%.0e')
cbar.set_label(r'Density ($N/\mathring{A}^3$)')
plt.xlabel(r'$X \ (\mathring{A})$')
plt.ylabel(r'$Y \ (\mathring{A})$')
plt.title('Density Distribution of Water OW in XY-Plane', pad=16)
circle = plt.Circle((xc, yc), R, color='green', fill=False)
plt.gca().add_patch(circle)
legend_elements = [Line2D([0], [0], color='green', lw=2, label='Fitted Circle for CNT')]
plt.legend(handles=legend_elements)
plt.xlim(0, x_box_length)
plt.ylim(0, y_box_length)
plt.tight_layout()
plt.savefig(sys.argv[1],transparent=True)

probability_density = smoothed_density / smoothed_density.sum()

plt.figure(figsize=(9.8,8))
img = plt.imshow(probability_density.T, extent=[D.results.density.edges[0][0], D.results.density.edges[0][-1], D.results.density.edges[1][0], D.results.density.edges[1][-1]], origin='lower', cmap='Blues', aspect='auto')
cbar = plt.colorbar(img, format='%.0e')
cbar.set_label('Probability Density')
plt.xlabel(r'$X \ (\mathring{A})$')
plt.ylabel(r'$Y \ (\mathring{A})$')
plt.title('Probability Density Distribution of Water OW in XY-Plane', pad=16)
circle = plt.Circle((xc, yc), R, color='green', fill=False)
plt.gca().add_patch(circle)
legend_elements = [Line2D([0], [0], color='green', lw=2, label='Fitted Circle for CNT')]
plt.legend(handles=legend_elements)
plt.xlim(0, x_box_length)
plt.ylim(0, y_box_length)
plt.savefig(sys.argv[2],transparent=True)

