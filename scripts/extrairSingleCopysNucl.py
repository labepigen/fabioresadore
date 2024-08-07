"""
Este script utiliza os resultados do orthofinder, extraindo
as sequencias nucletidicas dos single copy genes, salvando
em um arquivo fasta concatenado, para estudos de filogenômica.

Autor: Fábio Resadore
fresadore.bio@hotmail.com   
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os

arquivo_single_copy = '../../Results_Aug25_5/Orthogroups/Orthogroups_SingleCopyOrthologues.txt'
single_orthogroups_dir = '../../Results_Aug25_5/Single_Copy_Orthologue_Sequences/'

cepas = [
    "1365",
    "1545",
    "2335",
    "2371",
    "2490",
    "2689",
    "2690",
    "3283",
    "3315",
    "3356",
    "3358",
    "3371",
    "3481"
]

cepas_nucl_prefix = [
    "resadoreS1365",
    "LshaS1545",
    "Lguy2335",
    "Lguy2371",
    "Lsp2490",
    "resadoreS2689",
    "resadoreS2690",
    "Lsp3283",
    "Llai3315",
    "Lbra3356",
    "Lsp3358",
    "Lguy3371",
    "LshaS3481"
]

# print(len(cepas_nucl_prefix))

file = open(arquivo_single_copy)
file_lines = file.readlines()
file.close()

single_copy_list = {}  # lista dos genes presentes em cada singleCopy orthogroup

for line in file_lines:
    single_copy_list[line.strip()] = []
    # single_copy_list.append(line.strip())

# print(len(single_copy_list))
# print(single_copy_list[1])


"""file = open(single_orthogroups_dir + f'{single_copy_list[1]}.fa')
file_lines = file.readlines()
file.close()

print(file_lines[0])"""


def createNuclSingleCopysFile():
    fastas_iocl_nucl = {}  # fastas_iocl_nucl[cepa][gene_id] = fasta_object

    for cepa in cepas:
        arquivo = f'../../genes_anotados_nucl/IOCL_{cepa}_anotacao.fasta'
        with open(arquivo) as handle:
            fastas_iocl_nucl[cepa] = SeqIO.to_dict(
                SeqIO.parse(handle, "fasta"))

    fasta_infatum = {}

    arquivo_inf = f'TriTrypDB-66_LinfantumJPCM5_AnnotatedCDSs.fasta'
    with open(arquivo_inf) as handle:
        fasta_infatum = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))

    # print(fastas_iocl_nucl['1365']['resadoreS1365000005000'].seq)
    # exit()
    for single in single_copy_list:
        new_fasta = []
        # print(single)
        with open(single_orthogroups_dir + f'{single}.fa') as handle:
            genes_ = {}
            inf_gene_ = ''

            for record in SeqIO.parse(handle, "fasta"):
                # print(record)
                # exit()
                orthogroup_ = single
                id_ = str(record.id)
                cepa_ = id_[-16:-12]
                if cepa_ in cepas or 'LINF_' in str(record.id):
                    if 'LINF_' in id_:
                        inf_gene_ = id_[:-3]
                    else:
                        # print(id_[-11:-2])
                        genes_[cepa_] = id_[-11:-2]
                        # genes_.append()
                # exit()
            for i, cepa in enumerate(cepas):
                try:
                    new_fasta.append(
                        fastas_iocl_nucl[cepa][cepas_nucl_prefix[i]+genes_[cepa]])
                except:
                    print(orthogroup_ + ': ' + cepa)

            new_fasta.append(fasta_infatum[inf_gene_])
            # print(genes_)
        # print(len(new_fasta))
        if len(new_fasta) == 14:
            SeqIO.write(new_fasta, f"orthogroups/{orthogroup_}.fasta", "fasta")
        # exit()


def concatenarSingleCopyNucl():
    dir_path = f'orthogroups/aligned'

    orthogroups_files = []

    for file_path in os.listdir(dir_path):
        if os.path.isfile(os.path.join(dir_path, file_path)):
            orthogroups_files.append(str(file_path))

    orthogroups_sequences = {}

    for orthogroup_file in orthogroups_files:
        arquivo = f'{dir_path}/{orthogroup_file}'
        with open(arquivo) as handle:
            orthogroups_sequences[orthogroup_file[:-6]] = SeqIO.to_dict(
                SeqIO.parse(handle, "fasta"))

    # print(orthogroups_sequences['OG0001456']['resadoreS1365000707400'])

    cepas_ids = []
    sequencias_concatenadas = {}

    for id in orthogroups_sequences['OG0001456']:

        if 'LINF_' in str(id):
            id_ = 'LINF_'
        else:
            id_ = id[:-9].split("_")[-1]
        # print(id_)
        # exit()
        cepas_ids.append(id_)
        sequencias_concatenadas[id_] = ''
        # print(id)
    # exit()

    # print(cepas_ids)
    # exit()
    for i, orthogroup_id in enumerate(orthogroups_sequences):
        print(i)
        for cepa_id in cepas_ids:
            # print(cepa_id)
            for ortho_seq in orthogroups_sequences[orthogroup_id]:
                if cepa_id in ortho_seq:
                    # print(ortho_seq)
                    # exit()
                    sequencia_ = orthogroups_sequences[orthogroup_id][ortho_seq].seq
                    # print(sequencia_)
                    # exit()
                    sequencias_concatenadas[cepa_id] += sequencia_
        # exit()

        # break
        # print(len(sequencias_concatenadas['resadoreS1365']))
        # exit()

    fasta_concatenados = []

    for sequencia_concat in sequencias_concatenadas:
        seq_ = sequencias_concatenadas[sequencia_concat]
        fasta_ = SeqRecord(
            Seq(seq_),
            id=sequencia_concat,
            name=sequencia_concat,
            description='single_copy_concat'
        )

        fasta_concatenados.append(fasta_)

    SeqIO.write(fasta_concatenados, f"single_copy_concat.fasta", "fasta")


createNuclSingleCopysFile()

# ALINHAS AS SEQUENCIAS E SÓ DEPOIS RODAR A PROXIMA FUNCAO
# for i in *.fasta; do mafft --adjustdirection --quiet "$i" > aligned/"$i"; done

# concatenarSingleCopyNucl()
