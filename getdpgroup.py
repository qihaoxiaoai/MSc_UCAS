from ase.io import read,write
from ase.visualize import view
import numpy as np

def group_atoms(vasp_file):
    atoms = read(vasp_file, format='vasp')
    z_coords = atoms.positions[:, 2]
    c_indices = [i for i, atom in enumerate(atoms) if atom.symbol == 'C']
    c_z_coords = z_coords[c_indices]
    c_zmax_index = c_indices[np.argmax(c_z_coords)]
    c_zmin_index = c_indices[np.argmin(c_z_coords)]
    max_z_diff_indices_topc = np.where(np.abs(c_z_coords - c_z_coords[c_zmax_index]) < 0.4)[0]
    min_z_diff_indices_downc = np.where(np.abs(c_z_coords - c_z_coords[c_zmin_index]) < 0.4)[0]
    c_ids=np.array(c_indices)
    max_z_H_indices = []
    for ca_index in max_z_diff_indices_topc:
        ca_atom = atoms[ca_index]
        # Find all H atoms within 1.2 angstroms of the C atom
        t_indices = [t for t, atom in enumerate(atoms) if atom.symbol == 'H' and np.linalg.norm(atom.position - ca_atom.position) <= 1.2]
        # Append the H indices to the list
        max_z_H_indices.append(t_indices[0])
    max_z_H_ids=np.array(max_z_H_indices)

    min_z_H_indices = []
    for ci_index in min_z_diff_indices_downc:
        ci_atom = atoms[ci_index]
        d_indices = [d for d, atom in enumerate(atoms) if atom.symbol == 'H' and np.linalg.norm(atom.position - ci_atom.position) <= 1.2]
        min_z_H_indices.append(d_indices[0])
    min_z_H_ids=np.array(min_z_H_indices)
    combined_arr = np.concatenate((c_ids, min_z_H_ids, max_z_H_ids))
    new_atoms = atoms[combined_arr]
    write('new_POSCAR', new_atoms, format='vasp')
    lmp_min_z_H_ids=min_z_H_ids+1
    lmp_max_z_H_ids=max_z_H_ids+1
    lmp_c_ids=c_ids+1
    str_aromaticH1 = ' '.join([str(x) for x in lmp_min_z_H_ids])
    str_aromaticH2 = ' '.join([str(x) for x in lmp_max_z_H_ids])
    str_tan=' '.join([str(x) for x in lmp_c_ids])
    aromaticH='group cntH id '+ str_aromaticH1 + ' '+ str_aromaticH2 +'\n' 
    str_cnt= ' '.join([str(x) for x in lmp_c_ids])
    cnt='group cnt id ' + str_aromaticH1 + ' '+ str_aromaticH2 +' '+ str_cnt +'\n'
    with open('group.index', 'w') as f:
        f.write(cnt)
        f.write(aromaticH)

if __name__ == '__main__':
    import sys
    vasp_file = sys.argv[1]
    print(vasp_file)
    group_atoms(vasp_file)

