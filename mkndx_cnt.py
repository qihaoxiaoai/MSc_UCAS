#!/usr/bin/env python
import sys
import MDAnalysis as mda

def get_cnt_h_atom_indices(pdb_file):
    # Load PDB file with MDAnalysis
    u = mda.Universe(pdb_file)

    # Create empty lists to store up and down H atoms in CNT residues
    cnt_up_h_atoms = []
    cnt_down_h_atoms = []

    # Get carbon atoms in CNT residues
    cnt_c_atoms = u.select_atoms('resname CNT and name C')

    # Iterate over all hydrogen atoms in the universe
    for atom in u.select_atoms('resname CNT and name H'):
        # Check if the hydrogen atom is above or below the CNT carbon atoms
        if atom.position[2] > cnt_c_atoms.positions[:,2].mean():
            cnt_up_h_atoms.append(atom.index+1)
        else:
            cnt_down_h_atoms.append(atom.index+1)

    # Create a list of the hydrogen atom indices in CNT residues
    cnt_h_atom_indices = [cnt_up_h_atoms, cnt_down_h_atoms]

    # Create a third list by concatenating cnt_up_h_atoms and cnt_down_h_atoms and incrementing each element by 1
    all_h_atoms = [index for index in cnt_up_h_atoms + cnt_down_h_atoms]

    # Save the three lists to a text file
    with open('h.ndx', 'w') as f:
        for indices,i in zip(cnt_h_atom_indices + [all_h_atoms],['down','up','all']):
            f.write(f'[ H_{i} ]  \n')
            f.write(' '.join(map(str, indices)) + '\n')


def get_cnt_c_atom_indices(pdb_file):
    # Load PDB file with MDAnalysis
    u = mda.Universe(pdb_file)

    # Create an empty list to store C atom indices in CNT residues
    cnt_c_atom_indices = []

    # Iterate over all atoms in the universe
    for atom in u.atoms:
        # Check if atom is a carbon atom and is in a residue with name "CNT"
        if atom.name == 'C' and atom.resname == 'CNT':
            # Add atom index to the list
            cnt_c_atom_indices.append(atom.index+1)

    # Save the list of C atom indices to a text file
    with open('cnt.ndx', 'w') as f:
        f.write("[ CNT_C ]  \n")
        f.write('  '.join(map(str, cnt_c_atom_indices)))


def get_water_indices(pdb_file):
    # Load the PDB file
    u = mda.Universe(pdb_file)

    # Select oxygen atoms
    oxygen_sel = u.select_atoms("name OW")

    # Create an empty list to store the OW-HW pairs
    ow_hw_pairs = []

    # Loop through each oxygen atom
    for oxygen in oxygen_sel:
        # Find the hydrogen atoms within the specified distance of the oxygen atom
        near_h_sel = u.select_atoms(f"name HW* and around 1.1 index {oxygen.index}")
        # If there are two nearby hydrogen atoms, save their indices in a list with the oxygen index
        if len(near_h_sel) == 2:
            hw_indices = [near_h_sel[0].index, near_h_sel[1].index]
            ow_hw_pair = [oxygen.index] + hw_indices
            ow_hw_pairs.append(ow_hw_pair)

    ow_hw_pairs = [[index + 1 for index in pair] for pair in ow_hw_pairs]

    # Save the OW-HW pairs in a text file
    with open("water.ndx", "w") as f:
        wt_index = 0
        for pair in ow_hw_pairs:
            wt_index +=1
            f.write(f"[ water_{wt_index} ]  \n")
            f.write("  ".join([str(index) for index in pair]) + "\n")

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print(f"Usage: python {sys.argv[0]} input.pdb output.ndx")
        sys.exit(1)
    pdb_file = sys.argv[1]
    ndx_file = sys.argv[2]
    get_cnt_h_atom_indices(pdb_file)
    get_cnt_c_atom_indices(pdb_file)
    get_water_indices(pdb_file)

    # Combine the three index files into a single file
    with open(ndx_file, 'w') as f:
        # Read and write the contents of each index file
        for index_file in ['cnt.ndx','h.ndx', 'water.ndx']:
            with open(index_file, 'r') as g:
                f.write(g.read())
            # Add a blank line between each index group
            f.write('\n')
