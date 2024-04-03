import itertools
from pathlib import Path
import pandas as pd

# load dataSet and remove repeated rows
data = pd.read_csv("LIHC.csv")
data.drop_duplicates(inplace=True)
data.columns = ['Mutation_Patient']

# split Patient and Mutations by '-(' in the middle
data = data["Mutation_Patient"].str.split(r'-\(', expand=True)
data.columns = ['Mutation', 'Patient']

# replace ')' with '' at the end of patients
data["Patient"] = data["Patient"].str.replace(r'\)', '')

# write split data to file 'Splited_LIHC.csv'
# filepath = Path('Splited_LIHC.csv')
# data.to_csv(filepath)
# ------------------------------------------------------------------------------------
minSup = 2
# ------------------------------------------------------------------------------------
# frequent items with length=1
# I removed items with less than 'minsup' repeat
temp = data.groupby(['Mutation'])['Patient'].apply(set).reset_index()
temp = temp.assign(support_1=lambda x: x['Patient'].apply(len)).sort_values('support_1', ascending=False)
temp = temp[temp['support_1'] >= minSup].reset_index()
freqOneItem = temp

# write one-item frequent patterns to file 'freqOneItem.csv'
# filepath = Path('freqOneItem.csv')
# freqOneItem.to_csv(filepath)

# get one-by-one items as regard as 'freqOneItem' which keeps Mutations more equal than minsup
filtered = data[data['Mutation'].isin(freqOneItem['Mutation'])]

# write one-item frequent patterns to file '_filtered.csv'
# filepath = Path('_filtered.csv')
# filtered.to_csv(filepath)
# ------------------------------------------------------------------------------------
# construct dataframe based on patient and its mutations as well as number of mutations for each patient
patient_mutations = filtered.groupby(['Patient'])['Mutation'].apply(list).reset_index()
patient_mutations = patient_mutations.assign(length=lambda x: x['Mutation'].apply(len)).sort_values('length',
                                                                                                    ascending=True).reset_index()
# write constructed dataframe to file 'patient_mutations.csv'
# filepath = Path('patient_mutations.csv')
# patient_mutations.to_csv(filepath)
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
freqPatternDict = {}
# save frequent patterns in 'freqOneItem' to dictionary which is called 'freqPatternDict'
for index, row in freqOneItem.iterrows():
    freqPatternDict.update({row['Mutation']: row['support_1']})

# -------------------------------------------------------------------------------------
dataSet = patient_mutations.iloc[0:124, :]
maxlen = 11
lenDf = len(dataSet)
for subsetLen in range(2, maxlen + 1, 1):
    for index in range(lenDf):
        row = dataSet.iloc[index, :]
        if subsetLen <= row['length']:
            item_sets = set(itertools.combinations(row['Mutation'], subsetLen))
            for item in item_sets:
                # I check for prevent duplicate in constructing item sets
                if item not in freqPatternDict.keys():
                    cntItem = 1
                    # ----
                    # counting of item sets
                    for i in range(index + 1, len(patient_mutations)):
                        r = patient_mutations.iloc[i, :]
                        if set(item).issubset(set(r['Mutation'])):
                            cntItem += 1

                    # ----
                    if cntItem >= minSup:
                        # if an itemSet has illegal subset, it can not be as a frequent itemSet
                        # So I set its support to 0 ---> cntItem = 0
                        wrongSubsets = [dictItem for dictItem, sup in freqPatternDict.items() if
                                        set(item).issuperset({dictItem}) and sup == 0]
                        if any(wrongSubsets):
                            cntItem = 0
                        # if current itemSet is not illegal I will find closed items as regard as current item
                        if cntItem != 0:
                            # recognize itemSets witch are not closed
                            notClosedSet = [dictItem for dictItem, sup in freqPatternDict.items() if
                                            set(item).issuperset({dictItem}) and sup == cntItem]
                            dictOfItems = dict.fromkeys(notClosedSet, -1)
                            freqPatternDict.update(dictOfItems)
                        # update current itemSet's support
                        freqPatternDict.update({item: cntItem})
                    else:
                        freqPatternDict.update({item: 0})

    with open("freqPatternDict.txt", 'w') as f:
        for key, value in freqPatternDict.items():
            f.write('%s:%s\n' % (key, value))

freqPatterns = {}
tmp = [{dictItem: sup} for dictItem, sup in freqPatternDict.items() if
       sup != -1 and sup != 0]
for item in tmp:
    freqPatterns.update(item)


with open("freqPatterns.txt", 'w') as f:
    for key, value in freqPatterns.items():
        f.write('%s:%s\n' % (key, value))
