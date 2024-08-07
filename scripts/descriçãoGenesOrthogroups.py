"""
Este script retorna a descrição dos genes de um determinado 
conjunto de ortólogos, a partir de um conjunto de genomas anotados.    

Autor: Fábio Resadore
fresadore.bio@hotmail.com   
"""


orthogroups = {}


def listandoOrthogroups():

    arquivo = f"../../Results_Aug25_5/Orthogroups/Orthogroups.txt"

    file = open(arquivo)
    file_lines = file.readlines()
    file.close()

    for line in file_lines:
        line_ = line.split(':')
        orthgroup_ = line_[0]
        genes_ = line_[1].strip().split(' ')

        orthogroups[orthgroup_] = genes_


listandoOrthogroups()

gene_description = {}

# isso irá fazer com que os genomas com essas Tags não sejam utilizados
cepas_retirar = ['LPAL13', 'LINF_', 'Lperuviana_SAMEA2767111_',
                 'Llainsoni_LIN2019_', 'Lnaiffi_LnCL223_', 'Lguyanensis_LgCL085_genome', 'LbrM']

cepas = set()

text_file = ''


all_descriptions = set()


# esse é o arquivo que contem todos os genomas anotados, concatenados
arquivo = f"all_proteins.fasta"

file = open(arquivo)
file_lines = file.readlines()
file.close()

for line in file_lines:
    if '>' in line:
        all_descriptions.add(line)

# print(len(all_descriptions))

for ortho in orthogroups:
    text_file_ = f'{ortho}\n'
    for gene in orthogroups[ortho]:
        gene_id = ''
        ja_tem = True

        for retirar in cepas_retirar:
            if retirar in gene:
                ja_tem = False

        if ja_tem:
            if 'resadoreS' in gene or 'CLIOC' in gene or '_S' in gene:

                gene_id = gene[-16:]

            for desc in all_descriptions:
                if gene_id in desc:
                    desc_ = desc.split('.1 ')[1].strip()
                    text_file_ += f'{gene_id}\t{desc_}\n'
                    break

    text_file += f'{text_file_}\n'


# arquivo que serã salvo com as descrições
output = f"descricao_genes.txt"

file = open(output, "w")
file.write(text_file)
file.close()
