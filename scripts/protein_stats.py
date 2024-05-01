#!/usr/bin/env python3

from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import numpy as np
import sys

def calc_protein_stats(fasta):
  stats_matrix = []
  aa_matrix = []
  for sequence in SeqIO.parse(fasta, "fasta"):
    name = str(sequence.id)
    seq = str(sequence.seq).rstrip("*") # remove trailing asterisk (stop codon)
    aa = ProteinAnalysis(seq).get_amino_acids_percent()
    mw = ProteinAnalysis(seq).molecular_weight()
    pI = ProteinAnalysis(seq).isoelectric_point()
    flex_list = ProteinAnalysis(seq).flexibility()
    if len(flex_list) == 0:
      flex = "NA"
    else:
      flex = sum(flex_list) / len(flex_list)
    gravy = ProteinAnalysis(seq).gravy()
    aroma = ProteinAnalysis(seq).aromaticity()
    aliph = (100 * aa['A']) + (2.9 * 100 * aa["V"]) + \
            (3.9 * (100 * aa["I"] + 100 * aa["L"])) # formula from Akai, 1980
    stats_matrix.append([name, mw, pI, flex, gravy, aroma, aliph])
    aa_matrix.append(list(aa.values()))
  return np.concatenate((stats_matrix, aa_matrix), axis = 1)

def main():
  files = sys.argv[1:] # list of fasta files
  for file in files:
    name = file[:file.rfind(".")]
    stats_matrix = calc_protein_stats(file)
    np.savetxt(f"{name}_stats.csv", stats_matrix, delimiter = ",",
               header = "proteinID,MolecularWeight,pI,Flexibility,GRAVY,Aromaticity,AliphaticIndex,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y",
               fmt = "%s")

if __name__ == "__main__":
    main()
