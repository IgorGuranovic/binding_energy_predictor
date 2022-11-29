import re, os, subprocess, copy, pickle, pymol
from pymol import cmd
import numpy as np
import pandas as pd

dir = 'autodock_automation'

for pdbFile in os.listdir(dir + '/PDB_Files_Minimized'):
    subprocess.run('cmd /c "python %s/prepare_receptor4.py -r %s/PDB_Files_Minimized/' % (dir,dir) + pdbFile + ' -o %s/PDBQT_Files/' % (dir) + pdbFile + 'qt -A checkhydrogens"')

for pdbqtFile in os.listdir(dir + '/PDBQT_Files'):
    writer = open(dir + "/config.txt", "w")
    writer.write("receptor = " + pdbqtFile + "\nligand = pnp_xylose.pdbqt\n\ncenter_x = 0\ncenter_y = 25\ncenter_z = 35\n\nsize_x = 15\nsize_y = 15\nsize_z = 15\n\nenergy_range = 4\n\nexhaustiveness = 8")
    writer.close()
    mutant_name = pdbqtFile[0:-6]
    subprocess.run('"%s/vina.exe" --receptor %s/PDBQT_Files/' % (dir,dir) + mutant_name + '.pdbqt --ligand %s/pnp_xylose.pdbqt --config %s/config.txt --log %s/Logs/log_' % (dir,dir,dir) + mutant_name + '.txt --out %s/Outputs/output_' % (dir) + mutant_name + '.pdbqt')


pymol.finish_launching(['pymol', '-cq'])
energies = []
mutants_ordered = []
for log in os.listdir(dir + '/Logs'):
    orientation = 0
    removeLog0 = log.replace("log_mutant_", "")
    removeLog = removeLog0.replace(".txt", "")
    cmd.load(dir + '/PDBQT_Files/mutant_' + removeLog + '.pdbqt', object = 'enzyme')
    cmd.load(dir + '/Outputs/output_mutant_' + removeLog + '.pdbqt', object = 'ligand')
    for i in range(9):
        cmd.set('state', i+1, 'ligand')
        verdict = 0
        for oxygenIndex in range(3):
            condition1 = cmd.distance("(/enzyme//A/FUY`1448/C1)", "(/ligand///OH%d`0/O)" % (oxygenIndex + 1)) < 4
            condition2 = cmd.distance("(/enzyme//A/GLU`266/OE1)", "(/ligand///OH%d`0/O)" % (oxygenIndex + 1)) < 6
            condition3 = cmd.distance("(/enzyme//A/GLU`266/OE2)", "(/ligand///OH%d`0/O)" % (oxygenIndex + 1)) < 6
            if condition1 and condition2 and condition3:
                verdict = 1
        if verdict:
            orientation = i+1
            break
    cmd.reinitialize('everything')
    if orientation != 0:
        reader = open(dir + "/Logs/" + log, "r")
        text = reader.read()
        lower_bound = [m.start() for m in re.finditer(' %d ' % (orientation), text)][-1]
        energy = text[lower_bound+11:lower_bound+15]
        energies.append(energy)
        mutants_ordered.append(removeLog)
        reader.close()
cmd.quit()


sorted_indices = np.argsort(energies)
sorted_indices = np.flip(sorted_indices)

matrix = []
for n in range(len(mutants_ordered)):
    matrix.append([mutants_ordered[n], energies[n]])
sorted_matrix = []
for n in sorted_indices:
    sorted_matrix.append(matrix[n])
for n in range(len(sorted_matrix)):
    sorted_matrix[n].insert(0,n+1)

df = pd.DataFrame(sorted_matrix)
df.to_csv(os.path.join(dir + '/binding_affinity_rankings.csv'), sep=',', header=None, index=None)

sorted_matrix1 = copy.deepcopy(np.array(sorted_matrix).T.tolist())

mutant_list = sorted_matrix1[1]
energy_list = sorted_matrix1[2]

for enIndex, energy in enumerate(energy_list):
    energy_list[enIndex] = float(energy)

for mutIndex, mutant in enumerate(mutant_list):
    mutant_list[mutIndex] = list(mutant.split("+"))

with open(dir + '/Pickle_files/mutants_training.pickle', 'wb') as f:
    pickle.dump(mutant_list, f)

with open(dir + '/Pickle_files/training_labels.pickle', 'wb') as f:
    pickle.dump(energy_list, f)