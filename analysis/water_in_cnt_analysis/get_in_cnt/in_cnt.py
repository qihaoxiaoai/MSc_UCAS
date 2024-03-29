import os
import sys
import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.lib import distances

# 获取CNT内的水分子的函数
def get_water_inside_CNT(u, df_water, cylinder_radius, cnt_length):
    """Retrieve water molecules inside the CNT."""
    ow_inside_cylinder = u.select_atoms(f"type OT and cyzone {cylinder_radius} {cnt_length/2} {-cnt_length/2} resname CNT")
    water_molecules = []
    for oxygen in ow_inside_cylinder:
        if oxygen.index in df_water["OT"].values:
            water_info = df_water[df_water["OT"] == oxygen.index].iloc[0].tolist()
            ot_coords = oxygen.position.tolist()
            ht1_coords = u.atoms[int(water_info[1])].position.tolist()
            ht2_coords = u.atoms[int(water_info[2])].position.tolist()
            water_info.extend([ot_coords, ht1_coords, ht2_coords])
            water_molecules.append(water_info)
    return water_molecules

# 计算角平分线的函数
def calculate_bisector_angle(oxygen, h1, h2):
    """Calculate the bisector angle between water molecule and z-axis."""
    bisector = (h1.position + h2.position) / 2 - oxygen.position
    bisector /= np.linalg.norm(bisector)  # Normalize
    z_axis = np.array([0, 0, 1])
    cos_angle = np.dot(bisector, z_axis)
    angle_rad = np.arccos(np.clip(cos_angle, -1.0, 1.0))
    angle_deg = np.degrees(angle_rad)
    return angle_deg

def main(cnt_radius):
    # 设置路径
    data_path = '/home/jxhe/works/dp-cnt/cnt-dpmd-data/6-6-6-300'
    outdir = '/home/jxhe/works/dp-cnt/cnt-dpmd-data/work_plat'
    
    # 创建输出目录
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    print("Output directory:", outdir)
    print("Data directory:", data_path)
    os.chdir(outdir)
    print("Current working directory:", os.getcwd())
    
    # 加载轨迹
    trjtpr = os.path.join(data_path, "cnt.tpr")
    trjxtc = os.path.join(data_path, "nopbc.xtc")
    u = mda.Universe(trjtpr, trjxtc)
    total_frames = len(u.trajectory)
    print(f"Total frames in the trajectory: {total_frames}")

    # CNT键统计
    cnt_c_atoms = u.select_atoms("resname CNT and type CA")
    cnt_h_atoms = u.select_atoms("resname CNT and type HP")
    oxygens = u.select_atoms("type OT")
    cc_bond_count = 0
    ch_bond_count = 0
    ht_ot_bond_count = 0
    cc_bond_indices = set()
    ch_bond_indices = set()
    ht_ot_bond_indices = set()

    for c_atom in cnt_c_atoms:
        near_c_atoms = u.select_atoms(f"(type CA and around 1.6 index {c_atom.index}) and not index {c_atom.index}")
        for near_c in near_c_atoms:
            pair = tuple(sorted([c_atom.index, near_c.index]))
            if pair not in cc_bond_indices:
                cc_bond_indices.add(pair)
                cc_bond_count += 1
        near_h_atoms = u.select_atoms(f"type HP and around 1.3 index {c_atom.index}")
        for near_h in near_h_atoms:
            pair = (c_atom.index, near_h.index)
            if pair not in ch_bond_indices:
                ch_bond_indices.add(pair)
                ch_bond_count += 1

    # 水分子键统计和键长、键角计算
    water_data = []
    for oxygen in oxygens:
        near_h_sel = u.select_atoms(f"type HT and around 1.1 index {oxygen.index}")
        if len(near_h_sel) == 2:
            for near_h in near_h_sel:
                pair = (oxygen.index, near_h.index)
                if pair not in ht_ot_bond_indices:
                    ht_ot_bond_indices.add(pair)
                    ht_ot_bond_count += 1
            OH1_length = distances.calc_bonds(oxygen.position, near_h_sel[0].position)
            OH2_length = distances.calc_bonds(oxygen.position, near_h_sel[1].position)
            angle_rad = distances.calc_angles(near_h_sel[0].position, oxygen.position, near_h_sel[1].position)
            angle_deg = np.degrees(angle_rad)
            ow_hw_data = [oxygen.index] + [near_h_sel[0].index, near_h_sel[1].index] + [angle_deg, OH1_length, OH2_length]
            water_data.append(ow_hw_data)

    # 创建DataFrame
    index = [f"SOL{i+1}" for i in range(len(water_data))]
    columns = ["OT", "HT1", "HT2", "Angle", "OH1", "OH2"]
    df_water = pd.DataFrame(water_data, index=index, columns=columns)

    # 保存到CSV文件
    df_water.to_csv('00.all_water_information.csv', index_label="Water Molecule")
    bond_data = list(cc_bond_indices) + list(ch_bond_indices) + list(ht_ot_bond_indices)
    df_bonds = pd.DataFrame(bond_data, columns=["Atom1", "Atom2"])
    df_bonds.to_csv('01.all_bond_indices.csv', index=False)

    print(f"CNT have C-C Bond: {cc_bond_count}")
    print(f"CNT have C-H Bond: {ch_bond_count}")
    print(f"Water's HT-OT Bond: {ht_ot_bond_count}")
    print(f"Sum {len(water_data)} Water Molecules。")
    print("00.all_water_information.csv: Index OT HT1 HT2 Angle lenOH1 lenOH2")
    print("01.all_bond_indices.csv: Atom1 Atom2")

    cnt_c_atoms = u.select_atoms('resname CNT and name C')
    z_mean = cnt_c_atoms.positions[:,2].mean()
    cnt_up_h_atoms = u.select_atoms(f'resname CNT and name H and prop z > {z_mean}')
    cnt_down_h_atoms = u.select_atoms(f'resname CNT and name H and prop z <= {z_mean}')
    z_mean_up = cnt_up_h_atoms.positions[:,2].mean()
    z_mean_down = cnt_down_h_atoms.positions[:,2].mean()
    cnt_length = z_mean_up - z_mean_down - 1.2

    output_data = []
    for ts in u.trajectory:
        frame_data = get_water_inside_CNT(u, df_water, cnt_radius, cnt_length)
        bisector_angles = []
        hoh_angles = []
        bond_lengths = []
        same_orientation_count = 0
        for water_info in frame_data:
            oxygen = u.atoms[int(water_info[0])]
            h1 = u.atoms[int(water_info[1])]
            h2 = u.atoms[int(water_info[2])]
            
            # Calculate bisector angle
            bisector_angle = calculate_bisector_angle(oxygen, h1, h2)
            bisector_angles.append(bisector_angle)
            
            # Calculate HOH angle
            angle_rad = distances.calc_angles(h1.position, oxygen.position, h2.position)
            hoh_angles.append(np.degrees(angle_rad))
            
            # Calculate bond lengths
            OH1_length = distances.calc_bonds(oxygen.position, h1.position)
            OH2_length = distances.calc_bonds(oxygen.position, h2.position)
            bond_lengths.extend([OH1_length, OH2_length])
            
            if bisector_angle < 90:
                same_orientation_count += 1
    
        bisector_mean = np.mean(bisector_angles) if bisector_angles else 0
        hoh_angle_mean = np.mean(hoh_angles) if hoh_angles else 0
        bond_mean = np.mean(bond_lengths) if bond_lengths else 0
        orientation_prob = same_orientation_count / len(frame_data) if frame_data else 0
        if orientation_prob < 0.5:
            orientation_prob = 1 - orientation_prob
        
        # Ensure we have 9 bisector angles, fill with NaN if not
        while len(bisector_angles) < 9:
            bisector_angles.append(np.nan)
        
        output_data.append([ts.frame, len(frame_data)] + bisector_angles + [bisector_mean, orientation_prob, hoh_angle_mean, bond_mean])
    
    columns = ["Frame", "Number_in_CNT"] + [f"bsecz{i+1}" for i in range(9)] + ["bsecz_mean", "Orient_Prob", "hohang_mean", "bond_mean"]
    df_output = pd.DataFrame(output_data, columns=columns)
    df_output.to_csv('03_water_in_cnt_info.csv', index=False)
    print("02_water_in_cnt_info.csv: Frame In_Count 9-bsecz BsecZ_mean Oren_Prob Angle Bond")
    print("Successful!")
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python in_cnt.py <cnt_radius>")
        sys.exit(1)
    cnt_radius = float(sys.argv[1])
    main(cnt_radius)

