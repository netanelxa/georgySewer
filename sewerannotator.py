import sys, os
import pandas as pd
from math import ceil
import csv
from Bio.Seq import Seq
from Bio import SeqIO
import re


def findMutPileupFiles(dir, lines):
    r = []
    subdirs = [x[0] for x in os.walk(dir)]
    for subdir in subdirs:
        if "mutationsPileups" in subdir:
            files = os.walk(subdir).__next__()[2]
            files=[file for file in files for record in lines if record in file]
            if (len(files) > 0):
                for file in files:
                    r.append(os.path.join(subdir, file))
    return r

def getRegion(i, regionsList):
    regionTitles = {}
    for region in (regionsList):
        if i + 1 in range(int(region["start"]), int(region["end"])):
            regionTitles[region["id"]] = [int(region["start"]), int(region["end"])]
        if i + 1 < int(region["start"]):
            break
    if len(regionTitles) > 1:
        try:
            regex = re.compile("^(?!ORF)")
            filteredregionTitles = {k: v for k, v in regionTitles.items() if not k.startswith('ORF')}
            region_id = str([*filteredregionTitles][0])
            return region_id, filteredregionTitles[region_id][0], filteredregionTitles[region_id][1]
        except:
            filteredregionTitles = {k: v for k, v in regionTitles.items()}
            region_id = str([*filteredregionTitles][0])
            return region_id, filteredregionTitles[region_id][0], filteredregionTitles[region_id][1]
    else:
        filteredregionTitles = {k: v for k, v in regionTitles.items()}
        if filteredregionTitles:
            region_id = str([*filteredregionTitles][0])
        else:
            filteredregionTitles["none"] = [0, 0]
            region_id = str([*filteredregionTitles][0])
        return region_id, filteredregionTitles[region_id][0], filteredregionTitles[region_id][1]


def writeToCSV(writer, record, refnuc, mutnuc, i, regionTitle):
    writer.writerow({'Sequence ID': record, 'Reference Nucleotide': refnuc,
                     'Mutation nucleotide': mutnuc, 'location': i + 1,
                     'nuc name': str(i + 1) + " " + refnuc + " -> " + mutnuc,
                     'protein': str(regionTitle)})


def main(argv):
    regiontable = "regions.csv"
    with open(regiontable, 'r') as f:
        reader = csv.reader(f)
        my_list = []
        for row in reader:
            my_list.append({'segment': row[0], 'id': row[1],
                            'region': row[2], 'start': row[3], 'end': row[4], 'function': row[5]})
        regionsList = my_list[1:]
    dirPath = argv
    with open(r'env_apr_samples.txt') as f:
        lines = f.read().splitlines()
    filesList = findMutPileupFiles(dirPath, lines)
    freq_threshold = 5
    fieldnames = ['Sequence ID', 'Reference Nucleotide', 'Mutation nucleotide', 'location', 'nuc name', 'protein']

    with open('all_mutations.csv', 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for pileupFile in filesList:
            try:
                print(str(pileupFile))
                table = pd.read_csv(pileupFile)
                for index, row in table.iterrows():
                    nucs_Freq = set(["C_freq", "A_freq", "G_freq", "T_freq","del_freq"])
                    ref = row["ref"]
                    nucs_Freq.remove(ref + "_freq")
                    for nuc in nucs_Freq:
                        if row[nuc] > freq_threshold:
                            region = getRegion(row["pos"] + 1, regionsList)
                            sampleName = pileupFile.rsplit("/", 1)[1].rsplit(".")[0]
                            writeToCSV(writer, sampleName, ref, nuc.split("_")[0], row["pos"], region[0])
            except Exception as e:
                print(e)
                print(pileupFile.rsplit("/", 1)[1].rsplit(".")[0])
                
    csvfile.close()


if __name__ == '__main__':
    main(sys.argv[1])
