#!/usr/bin/python
import pandas as pd
import sewerannotator
import sys


def calculateFreqs(month):
    all_mutations = pd.read_csv('all_mutations_' + month + '.csv')
    knownmuts = pd.read_csv("known_muts.csv")
    b117muts = pd.read_csv("b117muts.csv")
    uniqueIDs = pd.unique(all_mutations["Sequence ID"])
    uniqueMuts = pd.unique(all_mutations["nuc name"])
    numofSequences = len(uniqueIDs)
    mutNucCount = all_mutations['nuc name'].value_counts()
    freqs = (mutNucCount / numofSequences * 100)
    freqs = freqs.T.to_dict()
    freqTable = pd.DataFrame()
    freqTable['nuc name'] = uniqueMuts
    mapping = dict(all_mutations[['nuc name', 'protein']].values)
    freqTable['protein'] = freqTable['nuc name'].map(mapping)
    mapping = dict(all_mutations[['nuc name', 'AAMutation']].values)
    freqTable['AA Mutation'] = freqTable['nuc name'].map(mapping)
    count_title = 'Count (' + str(numofSequences) + ')'
    freqTable[count_title] = freqTable['nuc name'].map(all_mutations['nuc name'].value_counts())
    freqTable['Percent of Total'] = freqTable['nuc name'].map(freqs)
    if all_mutations.groupby('nuc name')['Mut_Freq']:
        avg = all_mutations.groupby('nuc name')['Mut_Freq'].mean()
        freqTable['Average'] = freqTable['nuc name'].map(avg)
    freqTable.sort_values(by=['protein', 'Percent of Total'], ascending=False, inplace=True)
    freqTable = freqTable.loc[freqTable['Percent of Total'] >= 2]
    mapping = dict(b117muts[['nucleotide', 'lineage original']].values)
    freqTable['isUKLineage'] = freqTable['nuc name'].map(mapping)
    num_muts = len(freqTable.index)
    if num_muts <= 0:
        num_muts = 0
    geneFreq = freqTable['protein'].value_counts()
    geneFreq['total'] = freqTable['protein'].value_counts().sum()
    geneFreq = pd.DataFrame(geneFreq)
    if not geneFreq.empty:
        geneFreq['Percent of Total'] = geneFreq / num_muts * 100
    geneFreq.rename(columns={'protein': 'Count'}, inplace=True)

    writer = pd.ExcelWriter("Freq_Table_" + month + ".xlsx", engine='xlsxwriter')
    try:
        freqTable=freqTable.sort_values(by=["protein", "Average"],ascending=False)
    except:
        pass
    freqTable.to_excel(writer, sheet_name='Mutations Frequencies', index=False)
    geneFreq.to_excel(writer, sheet_name='gene count')
    writer.save()

    print("Freq_Table_" + month + ".csv is ready")


def main():
    months = ["env_feb_samples.txt", "env_mar_samples.txt", "env_apr_samples.txt", "env_may_samples.txt"]
    for month in months:
        month_str = month.split("_")[1]
        calculateFreqs(month_str)


if __name__ == '__main__':
    sewerannotator.main(sys.argv[1])
    main()
