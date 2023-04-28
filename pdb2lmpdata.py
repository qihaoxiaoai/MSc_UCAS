'''
Author: qihaoxiaoai qihaoxiaoai@gmail.com
Date: 2023-03-22 02:46:50
LastEditors: qihaoxiaoai
LastEditTime: 2023-03-22 20:00:58
Description: Jianxin & VegeSD
'''
# resname SOL

def gen_lmp_data(pdb_file):
    import os
    pdb_name = '.'.join(pdb_file.split('.')[:-1])
    tcl_cmd = '''pbc wrap -all -compound resid;
package require topotools;
topo retypebonds;
topo -sel "resname CNT" guessbonds;
topo -sel "resname CNT SOL" guessangles;

mol reanalyze 0;
topo writelammpsdata %s-lmp-amag.data;
quit;    
''' % pdb_name

    ret = os.popen(f"echo '{tcl_cmd}' | vmd -dispdev none {pdb_file}").read()
    #print (tcl_cmd, ret) ##for debug
    return None

if __name__ == '__main__':
    import sys
    pdb_file = sys.argv[1]
    gen_lmp_data(pdb_file)
