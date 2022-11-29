import pymol, os, random, copy, glob
from pymol import cmd

dir = 'autodock_automation'

def Remover(file_path):
    files = glob.glob("%s/*" % (file_path))
    for fi in files:
        os.remove(fi)
Remover(dir + '/PDB_Files')
Remover(dir + '/PDBQT_Files')
Remover(dir + '/PDB_Files_Minimized')
Remover(dir + '/Logs')
Remover(dir + '/Outputs')

mutants = []
num = 20000
locations = ['H34', 'E66', 'W67', 'H128', 'H129', 'R254']
changes = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

for x in range(num):
    temp = copy.deepcopy(locations)
    for loc_index, location in enumerate(temp):
        random.shuffle(changes)
        temp[loc_index] = temp[loc_index] + changes[0]
    mutants.append(temp)

for mutant in mutants:
    for mutation in mutant:
        if mutation[0] == mutation[-1]:
            mutant.remove(mutation)

pymol.finish_launching(['pymol', '-cq'])
for mutation_list in mutants:
    residue_list = ['']*len(mutation_list)
    position_list = ['']*len(mutation_list)
    for x in range(len(mutation_list)):
        position_list[x] = mutation_list[x][1:-1]
        if mutation_list[x][-1] == 'A':
            residue_list[x] = 'ALA'
        elif mutation_list[x][-1] == 'R':
            residue_list[x] = 'ARG'
        elif mutation_list[x][-1] == 'N':
            residue_list[x] = 'ASN'
        elif mutation_list[x][-1] == 'D':
            residue_list[x] = 'ASP'
        elif mutation_list[x][-1] == 'C':
            residue_list[x] = 'CYS'
        elif mutation_list[x][-1] == 'Q':
            residue_list[x] = 'GLN'
        elif mutation_list[x][-1] == 'E':
            residue_list[x] = 'GLU'
        elif mutation_list[x][-1] == 'G':
            residue_list[x] = 'GLY'
        elif mutation_list[x][-1] == 'H':
            residue_list[x] = 'HIS'
        elif mutation_list[x][-1] == 'I':
            residue_list[x] = 'ILE'
        elif mutation_list[x][-1] == 'L':
            residue_list[x] = 'LEU'
        elif mutation_list[x][-1] == 'K':
            residue_list[x] = 'LYS'
        elif mutation_list[x][-1] == 'M':
            residue_list[x] = 'MET'
        elif mutation_list[x][-1] == 'F':
            residue_list[x] = 'PHE'
        elif mutation_list[x][-1] == 'P':
            residue_list[x] = 'PRO'
        elif mutation_list[x][-1] == 'S':
            residue_list[x] = 'SER'
        elif mutation_list[x][-1] == 'T':
            residue_list[x] = 'THR'
        elif mutation_list[x][-1] == 'W':
            residue_list[x] = 'TRP'
        elif mutation_list[x][-1] == 'Y':
            residue_list[x] = 'TYR'
        elif mutation_list[x][-1] == 'V':
            residue_list[x] = 'VAL'
    cmd.load(dir + '/TmD224GbFucN3-Config1.pdb') #Replace with your directory
    cmd.wizard('mutagenesis')
    for x in range(len(mutation_list)):
        cmd.get_wizard().do_select('/TmD224GbFucN3-Config1//A/%s' % (position_list[x]))
        cmd.get_wizard().set_mode("%s" % (residue_list[x]))
        cmd.get_wizard().apply()
    mutation = "+".join(mutation_list)
    cmd.save(dir + '/PDB_Files/mutant_%s.pdb' % (mutation))
    cmd.reinitialize('everything')
cmd.quit()