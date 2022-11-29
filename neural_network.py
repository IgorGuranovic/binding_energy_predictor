import tensorflow as tf
from tensorflow import keras
import numpy as np
import pandas as pd
import pickle, os, re, gc
import matplotlib.pyplot as plt

import sys
np.set_printoptions(threshold=sys.maxsize)

dir = 'D:/PyRosetta/autodock_automation/Pickle_files'

with open(dir + '/mutants_training.pickle', 'rb') as f:
    mutants = pickle.load(f)

with open(dir + '/training_labels.pickle', 'rb') as f:
    labels = np.asarray(pickle.load(f))

def Matrix (data, num):

    def ResidueOneHot (matrix, mutantInd, mutationInd, residue):
        if residue == "A":
            matrix[mutantInd,mutationInd,0] = 1
        elif residue == "R":
            matrix[mutantInd,mutationInd,1] = 1
        elif residue == "N":
            matrix[mutantInd,mutationInd,2] = 1
        elif residue == "D":
            matrix[mutantInd,mutationInd,3] = 1
        elif residue == "C":
            matrix[mutantInd,mutationInd,4] = 1
        elif residue == "Q":
            matrix[mutantInd,mutationInd,5] = 1
        elif residue == "E":
            matrix[mutantInd,mutationInd,6] = 1
        elif residue == "G":
            matrix[mutantInd,mutationInd,7] = 1
        elif residue == "H":
            matrix[mutantInd,mutationInd,8] = 1
        elif residue == "I":
            matrix[mutantInd,mutationInd,9] = 1
        elif residue == "L":
            matrix[mutantInd,mutationInd,10] = 1
        elif residue == "K":
            matrix[mutantInd,mutationInd,11] = 1
        elif residue == "M":
            matrix[mutantInd,mutationInd,12] = 1
        elif residue == "F":
            matrix[mutantInd,mutationInd,13] = 1
        elif residue == "P":
            matrix[mutantInd,mutationInd,14] = 1
        elif residue == "S":
            matrix[mutantInd,mutationInd,15] = 1
        elif residue == "T":
            matrix[mutantInd,mutationInd,16] = 1
        elif residue == "W":
            matrix[mutantInd,mutationInd,17] = 1
        elif residue == "Y":
            matrix[mutantInd,mutationInd,18] = 1
        elif residue == "V":
            matrix[mutantInd,mutationInd,19] = 1
        
    matrix = np.zeros([num,6,20])

    for mutantInd, mutant in enumerate(data):
        mutantString = ' '.join(mutant)
        ind34 = [m.start() for m in re.finditer('34', mutantString)]
        if ind34 == []:
            matrix[mutantInd,0,8] = 1
        else:
            residue = mutantString[ind34[0]+2]
            ResidueOneHot(matrix, mutantInd, 0, residue)

        ind66 = [m.start() for m in re.finditer('66', mutantString)]
        if ind66 == []:
            matrix[mutantInd,1,6] = 1
        else:
            residue = mutantString[ind66[0]+2]
            ResidueOneHot(matrix, mutantInd, 1, residue)

        ind67 = [m.start() for m in re.finditer('67', mutantString)]
        if ind67 == []:
            matrix[mutantInd,2,17] = 1
        else:
            residue = mutantString[ind67[0]+2]
            ResidueOneHot(matrix, mutantInd, 2, residue)
        
        ind128 = [m.start() for m in re.finditer('128', mutantString)]
        if ind128 == []:
            matrix[mutantInd,3,8] = 1
        else:
            residue = mutantString[ind128[0]+3]
            ResidueOneHot(matrix, mutantInd, 3, residue)
        
        ind129 = [m.start() for m in re.finditer('129', mutantString)]
        if ind129 == []:
            matrix[mutantInd,4,8] = 1
        else:
            residue = mutantString[ind129[0]+3]
            ResidueOneHot(matrix, mutantInd, 4, residue)
        
        ind254 = [m.start() for m in re.finditer('254', mutantString)]
        if ind254 == []:
            matrix[mutantInd,5,1] = 1
        else:
            residue = mutantString[ind254[0]+3]
            ResidueOneHot(matrix, mutantInd, 5, residue)
        
    return matrix
        
dataFrame = Matrix(mutants, len(labels))

n=10

test_mutants = np.asarray([x for i, x in enumerate(dataFrame) if i%n ==0])
test_labels = np.asarray([x for i, x in enumerate(labels) if i%n ==0])

training_mutants = np.asarray([x for i, x in enumerate(dataFrame) if i%n !=0])
training_labels = np.asarray([x for i, x in enumerate(labels) if i%n !=0])

model = keras.Sequential([
    keras.layers.Flatten(input_shape=(6,20)),
    keras.layers.Dense(64, activation ='relu'),
    keras.layers.Dense(16, activation ='relu'),
    keras.layers.Dense(1, activation='linear')
])

keras.optimizers.Adam(lr=.01, beta_1=0.9, beta_2=0.999, amsgrad=False)
model.compile(loss='mean_squared_error', metrics=['mse'])
fit = model.fit(training_mutants, training_labels, epochs=25, validation_split=0.2, validation_data=None, verbose=1)
test_loss, test_acc = model.evaluate(test_mutants, test_labels, verbose = 1)

print('Test error:', test_loss)

predict = model.predict(test_mutants)

predict_list = []
test_labels_list = []

for element in predict:
    predict_list.append(element[0])

for element in test_labels:
    test_labels_list.append(element)

y = predict_list
x = test_labels_list

writeFile = pd.DataFrame([x,y])
writeFile.to_csv(os.path.join('neural_network_validation.csv'), sep=',', header=None, index=None)

plt.scatter(x,y,s=10)
plt.show()

filteredEnergies = []
filteredMutants = []
filteredPlus = []

residues = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
for element1 in residues:
    predictionMutants = []
    for element2 in residues:
        for element3 in residues:
            for element4 in residues:
                for element5 in residues:
                    for element6 in residues:
                        predictionMutants.append(['H34' + element1, 'E66' + element2, 'W67' + element3, 'H128' + element4, 'H129' + element5, 'R254' + element6])

    for predictionMutant in predictionMutants:
        for mutation in predictionMutant:
            if mutation[0] == mutation[-1]:
                predictionMutant.remove(mutation)

    predictedFrame = Matrix(predictionMutants, len(predictionMutants))
    predictions = model.predict(predictedFrame)

    for energyIndex, energy in enumerate(predictions):
        if energy < -6.2:
            filteredEnergies.append(energy[0])
            filteredMutants.append(predictionMutants[energyIndex])
    del(predictionMutants)
    gc.collect()

sortedIndices = np.argsort(filteredEnergies)

filteredPlus = []
for element in filteredMutants:
    tempstring = ""
    for mutation in element:
        tempstring = tempstring + mutation + " + "
    tempstring = tempstring[0:-3]
    filteredPlus.append(tempstring)
filteredPlus = np.asarray(filteredPlus)

excel = []
for n in range(len(filteredPlus)):
    excel.append([filteredPlus[n], filteredEnergies[n]])
    
sortedExcel = []
for n in sortedIndices:
    sortedExcel.append(excel[n])
for n in range(len(sortedExcel)):
    sortedExcel[n].insert(0,n+1)

df = pd.DataFrame(sortedExcel)
df.to_csv(os.path.join('network_predictions.csv'), sep=',', header=None, index=None)