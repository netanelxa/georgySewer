import sys, os
import pandas as pd
from math import ceil
import csv
from Bio.Seq import Seq
from Bio import SeqIO
import re
import os

codon_map = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"
}


def full_codon_gaps(sequence, start, end, gap='-'):
    """
    avoid partial gaps in codon and convert to whole gaps codon
    exmple: from A-- to ---
    :param sequence: Seq object (Biopython)
    :param start: start of reading frame
    :param end: end of reading frame
    :param gap: gap character. default: '-'
    :return: new sequence, with full codons
    """
    old_seq = str(sequence)
    new_seq = old_seq[:start]
    for i in range(start - 1, end, 3):
        codon = old_seq[i: i + 3]
        if '-' in codon:
            codon = '---'
        new_seq += codon
    new_seq += old_seq[end:]
    return Seq.Seq(new_seq)


# get start and end of region
# don't forget -1 in position
def translate(sequence, start, end, codon_table):
    """
    translate nucleotides sequence in given region to amino acid sequence according to codons in region start-end
    :param sequence: nucleotides sequence as str #
    :param start: position of first nucleotide
    :param end: position of last nucleotide
    :return: translated sequence (aa seq)
    """
    tranlsated = []
    for i in range(start, end, 3):
        codon = sequence[i:i + 3]
        if codon in codon_table:
            aa = codon_table[codon]
        else:
            aa = 'X'  # ignore frameshift
        tranlsated.append(aa)
    return tranlsated


def getTranslate(i, referenceSequence, record, flag, regionStart, regionEnd):
    RefAA = translate(referenceSequence, regionStart - 1, regionEnd, codon_map)
    OtherSeqAA = record
    aaMutIndexr = (i + 1 - regionStart) / 3
    if aaMutIndexr % 1 == 0:
        aaMutIndex = int(aaMutIndexr + 1)
    else:
        aaMutIndex = ceil((i + 1 - regionStart) / 3)
    refaaseq = RefAA[aaMutIndex - 2:aaMutIndex + 2]
    otheraaseq = OtherSeqAA[aaMutIndex - 2:aaMutIndex + 2]
    if refaaseq:
        try:
            AAMutToCSv = str(RefAA[aaMutIndex - 1] + str(aaMutIndex) + OtherSeqAA[aaMutIndex - 1])
        except:
            print(str(aaMutIndex) + " " + str(len(RefAA)) + " " + str(len(OtherSeqAA)))
            AAMutToCSv = "None"
        if flag == 2:
            AAMutToCSv = AAMutToCSv[:-1] + "Del"
    else:
        AAMutToCSv = "None"
    return AAMutToCSv


def findMutPileupFiles(dir, lines):
    r = []
    a = []
    subdirs = [x[0] for x in os.walk(dir)]
    for subdir in subdirs:
        if "mutationsPileups" in subdir:
            files = os.walk(subdir).__next__()[2]
            files = [file for file in files for record in lines if record in file]
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


def writeToCSV(writer, record, refnuc, mutnuc, i, regionTitle, AAMutToCSv, mutFreq):
    writer.writerow({'Sequence ID': record, 'Reference Nucleotide': refnuc,
                     'Mutation nucleotide': mutnuc, 'location': i + 1,
                     'nuc name': str(i + 1) + " " + refnuc + " -> " + mutnuc,
                     'protein': str(regionTitle), 'AAMutation': AAMutToCSv, 'Mut_Freq': mutFreq})


def findMutations(dirPath, lines, regionsList, month, refSeq):
    RefAA = ''.join(map(str, refSeq))
    OtherSeqAA = ''.join(map(str, refSeq))
    filesList = findMutPileupFiles(dirPath, lines)
    freq_threshold = 5
    fieldnames = ['Sequence ID', 'Reference Nucleotide', 'Mutation nucleotide', 'location', 'nuc name', 'protein',
                  'AAMutation', 'Mut_Freq']
    with open('all_mutations_' + month + '.csv', 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for pileupFile in filesList:
            table = pd.read_csv(pileupFile)
            for index, row in table.iterrows():
                nucs_Freq = set(["C_freq", "A_freq", "G_freq", "T_freq", "del_freq"])
                ref = row["ref"]
                nucs_Freq.remove(ref + "_freq")
                for nuc in nucs_Freq:
                    if row[nuc] > freq_threshold:
                        regionTitle, regionStart, regionEnd = getRegion(row["pos"] + 1, regionsList)
                        record = list(OtherSeqAA[int(regionStart - 1):int(regionEnd)])
                        a = row["pos"]
                        if int(regionEnd) == 0:
                            record = list(OtherSeqAA)
                        record[int(row["pos"]) - int(regionStart)+1] = nuc.split("_")[0]
                        record = translate(''.join(map(str, record)), 0, len(record), codon_map)
                        nuc_name = nuc.split("_")[0] if nuc.split("_")[0] != "del" else "-"
                        flag = 1 if nuc_name != '-' else 2
                        AAMutToCSv = getTranslate(row["pos"] + 1, RefAA, record, flag, regionStart, regionEnd)
                        sampleName = pileupFile.rsplit("/", 1)[1].rsplit(".")[0]
                        writeToCSV(writer, sampleName, ref, nuc_name, row["pos"], regionTitle, AAMutToCSv,
                                   row[nuc])
    csvfile.close()




def main(argv):
    print("Starting")
    regiontable = os.path.dirname(os.path.abspath(__file__))+"/regions.csv"
    refSeq = list(SeqIO.parse(os.path.dirname(os.path.abspath(__file__))+"/REF_NC_045512.2.fasta", "fasta"))[0].seq
    with open(regiontable, 'r') as f:
        reader = csv.reader(f)
        my_list = []
        for row in reader:
            my_list.append({'segment': row[0], 'id': row[1],
                            'region': row[2], 'start': row[3], 'end': row[4], 'function': row[5]})
        regionsList = my_list[1:]
    dirPath = str(argv)
    months = ["env_feb_samples.txt", "env_mar_samples.txt", "env_apr_samples.txt", "env_may_samples.txt"]
    for month in months:
        with open(month) as f:
            month_str = month.split("_")[1]
            lines = f.read().splitlines()
        findMutations(dirPath, lines, regionsList, month_str, refSeq)
        print("all mutations " + month_str + " is ready")


if __name__ == '__main__':
    main(sys.argv[1])
