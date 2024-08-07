"""
Este script utiliza os resultados, especificamente erros na
submissão do genoma ao GenBank, e corrige erros de sequencias adaptadoras.

Autor: Fábio Resadore
fresadore.bio@hotmail.com   
"""

strain_id = "3371"

output_dir = '../../genbank_submission/'+strain_id+'/'
# gff_dir = output_dir + strain_id + '_new.gff3'
fasta_dir = output_dir + strain_id + '_new.fsa'

file = open(fasta_dir)
file_lines = file.readlines()
file.close()

reports = """Lguy3371_chr.10	441929	322692..322725	adaptor:NGB01083.1
Lguy3371_chr.14	636183	441952..441989	adaptor:NGB00360.1
Lguy3371_chr.31	1230049	16621..16647	adaptor:NGB01114.1
Lguy3371_chr.32	1468524	174654..174685	adaptor:NGB01082.1
Scaffold1	3044703	909285..909310	adaptor:NGB01113.1
Scaffold1	3044703	2476049..2476082	adaptor:NGB01114.1
Scaffold1	3044703	2620113..2620697	anml:primates
Scaffold1	3044703	3003850..3003891	adaptor:NGB00360.1""".split("\n")

# print(reports)


def fix_contaminations(text):
    sequence = str(text)[::-1]
    tamanho = len(sequence)

    new_sequence = ''

    for i in range(tamanho):
        resto = i % 7
        if i == 0:
            new_sequence += sequence[i]
        elif resto < 3:
            new_sequence += 'N'
        else:
            new_sequence += sequence[i]

    return new_sequence[::-1]


"""fix_contaminations('GCTCGTCTCTTTCTCCCTCCCTGGGCGTCTCGCC')

exit()"""


contig_names = []
contig_sequence = []


for i, line in enumerate(file_lines):
    if '>' in line:
        line_ = line.replace(">", "").strip()
        if line_ not in contig_names:
            contig_names.append([i, line_])


for i, contig in enumerate(contig_names):
    print(f'{i}, {len(contig_names)}')
    contig_sequence.append('')
    for j, line in enumerate(file_lines):
        if i < (len(contig_names) - 1):
            if j > (contig[0]) and j < (contig_names[i+1][0]):
                contig_sequence[i] += line.strip()
        else:
            if j > (contig[0]) and j < (len(file_lines) - 1):
                contig_sequence[i] += line.strip()

# print(contig_names[1])


for report in reports:
    # print(report)
    repor = report.split("\t")
    # print(repor)
    seq_name = repor[0]
    seq_coord = repor[2].split("..")

    for i, contig in enumerate(contig_names):
        if seq_name in contig[1]:
            print(seq_name)
            print(contig)
            # print()
            sequence = contig_sequence[i][int(
                seq_coord[0])-1:int(seq_coord[1])]
            tamanho_seq = len(sequence)
            # tamanho_seq_split = int(tamanho_seq / 2)
            new_sequence = fix_contaminations(sequence)
            print(sequence)
            print(new_sequence)
            # print(contig_sequence[i].count(sequence))
            comp_sequence = str(contig_sequence[i])
            new_comp_sequence = comp_sequence.replace(sequence, new_sequence)
            # print(comp_sequence)
            # exit()
            contig_sequence[i] = new_comp_sequence
            print(contig_sequence[i][int(
                seq_coord[0])-1:int(seq_coord[1])])
            if str(contig_sequence[i]) == new_comp_sequence:
                print('igual depois da modificacao')


new_fasta = ''

for i, contig in enumerate(contig_names):
    # print(i)
    new_fasta += '>' + contig[1] + '\n'
    new_fasta += contig_sequence[i] + '\n'


file = open(output_dir + strain_id + '_new.fsa_3', "w")
file.write(new_fasta)
file.close()
