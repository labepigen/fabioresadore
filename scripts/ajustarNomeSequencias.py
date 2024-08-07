"""
Este script renomeia e corrige o nome de sequências específicas    

Autor: Fábio Resadore
fresadore.bio@hotmail.com   
"""

isoenzimas = [
    '6pgdh',
    'g6pdh',
    'mpi',
    'icd',
    'pgm',
    'gpi'
]

genes_ids = [
    'Lguy3371',
    'Lguy2371',
    'Lguy2335',
    'Lsp3358',
    'LPAL13',
    'LbrM',
    'LperSAMEA',
    'LshaS3481',
    'LshaS1545',
    'Lbra3356',
    'LguyLgCL',
    'resadoreS2690',
    'resadoreS2689',
    'resadoreS1365',
    'LnaifLnCL223',
    'Lsp3283',
    'LlaiLIN',
    'Llai3315',
    'Lsp2490',
    'LINF'
]

gene_renomear = [
    'L. guyanensis - IOCL 3371',
    'L. guyaneisis - IOCL 2371',
    'L. guyanensis - IOCL 2335',
    'L. sp. (Lbraz/Lguya) - IOCL 3358',
    'L. panamensis - TriTrypDB GCA_000340495.1 ',
    'L. braziliensis - TriTrypDB GCA_000002845.2 ',
    'L. peruviana - GenBank ERR599203',
    'L. shawi - IOCL 3481',
    'L. shawi - IOCL 1545',
    'L. braziliensis - IOCL 3356',
    'L. guyanensis - GenBank ERX180458',
    'L. lindenbergi - IOCL 2690',
    'L. utingensis - IOCL 2689',
    'L. naiffi - IOCL 1365',
    'L. naiffi - GenBank ERX180449',
    'L. sp. (Lbraz/Lnaif) - IOCL 3283',
    'L. lainsoni - GenBank SRX4999894',
    'L. lainsoni - IOCL 3315',
    'L. sp. (Llain/Lnaif) - IOCL 2490',
    'L. infantum'

]

print(len(genes_ids))
print(len(gene_renomear))
# exit()

for isoenz in isoenzimas:

    arquivo = f'blasts/extraidos/align/{isoenz}_genes.fasta'

    new_file = ''

    file = open(arquivo)
    file_lines = file.readlines()
    file.close()

    lines_ = []

    for line in file_lines:
        lines_.append(line)
        """line_ = ''
        for i, cepa_nome in  enumerate(genes_ids):
            if cepa_nome in line:
                #print(line)
                line_ = f'>{gene_renomear[i]}\n'
            #print(line)
            else:
                line_ += line
        new_file += line_"""
    for line in lines_:
        line_ = line
        if '>' in line_:
            for i, cepa_nome in enumerate(genes_ids):
                if cepa_nome in line_:
                    line_ = f'>{gene_renomear[i]}\n'

        # print(line_)
        new_file += line_

    # print(len(new_file))
    # exit()
    file = open(f'{arquivo}.fasta', 'w')
    file.write(new_file)
    file.close()

    # exit()
