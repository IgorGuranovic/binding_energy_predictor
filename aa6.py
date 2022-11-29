import re, os, subprocess, pymol
from pymol import cmd
from csv import reader
import pandas as pd

dir = 'autodock_automation'

for pdbFile in os.listdir(dir + '/PDB_Files_Minimized_Validation'):
    subprocess.run('cmd /c "python %s/prepare_receptor4.py -r %s/PDB_Files_Minimized_Validation/' % (dir,dir) + pdbFile + ' -o %s/PDBQT_Files_Validation/' % (dir) + pdbFile + 'qt -A checkhydrogens"')

for pdbqtFile in os.listdir(dir + '/PDBQT_Files_Validation'):
    writer = open(dir + "/config.txt", "w")
    writer.write("receptor = " + pdbqtFile + "\nligand = pnp_xylose.pdbqt\n\ncenter_x = 0\ncenter_y = 25\ncenter_z = 35\n\nsize_x = 15\nsize_y = 15\nsize_z = 15\n\nenergy_range = 4\n\nexhaustiveness = 8")
    writer.close()
    mutant_name = pdbqtFile[0:-6]
    subprocess.run('"%s/vina.exe" --receptor %s/PDBQT_Files_Validation/' % (dir,dir) + mutant_name + '.pdbqt --ligand %s/pnp_xylose.pdbqt --config %s/config.txt --log %s/Logs_Validation/log_' % (dir,dir,dir) + mutant_name + '.txt --out %s/Outputs_Validation/output_' % (dir) + mutant_name + '.pdbqt')

pymol.finish_launching(['pymol', '-cq'])
energies = []
mutants_ordered = []
for log in os.listdir(dir + '/Logs_Validation'):
    orientation = 0
    removeLog0 = log.replace("log_mutant_", "")
    removeLog = removeLog0.replace(".txt", "")
    cmd.load(dir + '/PDBQT_Files_Validation/mutant_' + removeLog + '.pdbqt', object = 'enzyme')
    cmd.load(dir + '/Outputs_Validation/output_mutant_' + removeLog + '.pdbqt', object = 'ligand')
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
        reader_text = open(dir + "/Logs_Validation/" + log, "r")
        text = reader_text.read()
        lower_bound = [m.start() for m in re.finditer(' %d ' % (orientation), text)][-1]
        energy = text[lower_bound+11:lower_bound+15]
        energies.append(energy)
        mutants_ordered.append(removeLog)
        reader_text.close()
    else:
        reader_text = open(dir + "/Logs_Validation/" + log, "r")
        text = reader_text.read()
        energies.append(0)
        mutants_ordered.append(removeLog)
        reader_text.close()
cmd.quit()

with open (dir + '/network_predictions.csv', 'r') as networkPredictions:
    csvList = reader(networkPredictions)
    matrix = list(csvList)

mutantPredictions = []
energyPredictions = []

for index in range(25):
    mutantPredictions.append(matrix[index][1].replace(" ", ""))
    energyPredictions.append(matrix[index][2])

sorted_indices = []

for element in mutantPredictions:
    for index in range(len(mutants_ordered)):
        if mutants_ordered[index] == element:
            sorted_indices.append(index)

matrix = []
for n in range(len(mutants_ordered)):
    matrix.append([mutants_ordered[n], energies[n]])
sorted_matrix = []
for n in sorted_indices:
    sorted_matrix.append(matrix[n])
for n in range(len(sorted_matrix)):
    sorted_matrix[n].insert(0,n+1)
    sorted_matrix[n].insert(3,mutantPredictions[n])
    sorted_matrix[n].insert(4,energyPredictions[n])

df = pd.DataFrame(sorted_matrix)
df.to_csv(os.path.join(dir + '/binding_affinity_predictions.csv'), sep=',', header=None, index=None)