import csv

# Create dictionary to hold the data
valDic = {}
Count = 0

#years, genders, status, relation = set(), set(), set(), set()

# Read data into dictionary
with open('matrix.mtx', 'r',) as inputfile:

    reader = csv.reader(inputfile, delimiter = '\t')
    next(reader)

    for row in reader:

        key = row[0]

        cells.add(key[0])
        
        if key in valDic:
            valDic[key] +=1
            Count += 1


        
 
        


#Add missing combinations


#Prepare new CSV
newcsvfile = [["ADT", "cells"]] 

for key, val in sorted(valDic.items()):
    newcsvfile.append([valDic[key], Count])

with open('results5.csv', "w", newline='') as outputfile:
    writer = csv.writer(outputfile)
    writer.writerows(newcsvfile)  