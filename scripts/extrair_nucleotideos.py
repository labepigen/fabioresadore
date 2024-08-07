"""
Este script extrai determinadas sequencias obtidas a partir de 
um conjunto de resultados do BLAST, e salva em um arquivo fasta

Autor: Fábio Resadore
fresadore.bio@hotmail.com   
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


isoenzimas = [
    '6pgdh',
    'g6pdh',
    'mpi',
    'icd',
    'pgm',
    'gpi'
]




todos_genes = {}

arquivo = f'../genes_anotados_nucl/genomas_doc.fas'
with open(arquivo) as handle:
    todos_genes = SeqIO.to_dict(
        SeqIO.parse(handle, "fasta"))


genes_ids = {}

for isoenz in isoenzimas:
    genes_ids[isoenz] = []
    arquivo = f'blasts/extraidos/{isoenz}/todos_genomas.out'

    file = open(arquivo)
    file_lines = file.readlines()
    file.close()

    fasta_iso = []
    for line in file_lines:
        line_ = line.strip().split('\t')
        geneId_ = line_[1]

        # print(todos_genes[geneId_])

        fasta_iso.append(todos_genes[geneId_])
        SeqIO.write(
            fasta_iso, f"blasts/extraidos/{isoenz}_genes.fasta", "fasta")

    # exit()
        # genes_ids[isoenz].append(geneId_)


# print(len(todos_genes))
