#!/usr/bin/python
import pandas as pd
import sewerannotator
import sys


def calculateFreqs(month):
    all_mutations = pd.read_csv('all_mutations_'+month+'.csv')
    knownmuts = pd.read_csv("known_muts.csv")
    b117muts = pd.read_csv("b117muts.csv")
    uniqueIDs=pd.unique(all_mutations["Sequence ID"])
    uniqueMuts=pd.unique(all_mutations["nuc name"])
    numofSequences=len(uniqueIDs)
    mutNucCount=all_mutations['nuc name'].value_counts()
    freqs=(mutNucCount/ numofSequences*100)
    freqs=freqs.T.to_dict()
    freqTable=pd.DataFrame()
    freqTable['nuc name']=uniqueMuts
    mapping = dict(all_mutations[['nuc name', 'protein']].values)
    freqTable['protein'] = freqTable['nuc name'].map(mapping)
    mapping = dict(knownmuts[['nuc.name', 'AAMutation']].values)
    freqTable['Known mutation'] = freqTable['nuc name'].map(mapping)
    count_title='Count ('+str(numofSequences)+')'
    freqTable[count_title]=freqTable['nuc name'].map(all_mutations['nuc name'].value_counts())
    freqTable['Freq'] = freqTable['nuc name'].map(freqs)
    freqTable.sort_values(by=['protein', 'Freq'],ascending=False,inplace=True)
    freqTable=freqTable.loc[freqTable['Freq'] >= 2]
    mapping = dict(b117muts[['nucleotide', 'lineage original']].values)
    freqTable['isUKLineage'] = freqTable['nuc name'].map(mapping)
    freqTable.to_csv("Freq_Table.csv", index=False)
    print("Freq_Table_"+month+".csv is ready")


def main():
    months = ["env_feb_samples.txt", "env_mar_samples.txt", "env_apr_samples.txt", "env_may_samples.txt"]
    for month in months:
        month_str = month.split("_")[1]
        calculateFreqs(month_str)


if __name__ == '__main__':
    sewerannotator.main(sys.argv[1])
    main()

