"""
Este script faz uma limpeza no GFF

Autor: Fábio Resadore
fresadore.bio@hotmail.com   
"""

# import json
strain_id = "1545"
file = open('resadoreS'+strain_id+'-refLbraz/pseudo.out.gff3')
file_lines = file.readlines()
file.close()

lines = []
for i in file_lines:
    lines.append(i.split("	"))

#print(lines[50])
# extrair todos os IDs
id_suffix = "Lsha_S1545_"
# id = id_suffix + lines[563][8].split(id_suffix)[1].split(";")[0].split(":")[0].split(".")[0]


def get_gene_ids(lines_in, suffix):
    gene_ids_ = []
    for i in lines_in:
        first_char = i[0][0]
        if first_char != "#":
            att_len = len(i[8].split(suffix))
            if att_len >= 2:
                id = suffix + \
                    i[8].split(suffix)[1].split(";")[0].split(
                        ":")[0].split(".")[0].strip()
                if id not in gene_ids_ and "CTG" not in id:
                    gene_ids_.append(id)

    genes_ids_txt = ""

    for i in gene_ids_:
        genes_ids_txt += i + "\n"

    file = open(strain_id+'_genes_ids.txt', "w")
    file.write(genes_ids_txt)
    file.close()

    return gene_ids_


gene_ids = get_gene_ids(lines, id_suffix)
# print("IDs dos genes encontrados")
print(f"Quantidade de genes IDs encontrados: {len(gene_ids)}")


# -------------------------------


# QUAIS TYPES DE SEQUENCIAS PARA CADA ID

def gene_types_retrieve(gene_ids_in, lines_in):
    genes_types_ids = {}

    for i in gene_ids_in:
        genes_types_ids[i] = []
        for j in lines_in:
            first_char = j[0][0]
            if first_char != "#":
                gene_attr = j[8]
                if i in gene_attr:
                    # if j[2] not in genes_types_ids[i]:
                    genes_types_ids[i].append(j[2])

    uniques_type_config = []

    """lista = ['gene',   'exon']
    lista.sort()
    if lista == ['exon', 'gene']:
        print('IGUAL')

    print(lista)"""

    for i in gene_ids_in:
        types = genes_types_ids[i]
        types.sort()
        if len(uniques_type_config) <= 0:
            uniques_type_config.append(types)
        else:
            if types not in uniques_type_config:
                uniques_type_config.append(types)

    print(uniques_type_config)
    print(len(uniques_type_config))


# gene_types_retrieve(gene_ids, lines)
# print("Types presentes no GFF encontrados")

# -------------------------------


def gene_struc_check(feature_upper, featue_lower, gene_id):
    if featue_lower[0] < feature_upper[0] or featue_lower[1] > feature_upper[1]:
        print(
            f"{gene_id}: erro na hierarquia dos genes {feature_upper} e {featue_lower}")


def gene_pos_check(gene, gene_id):
    if gene[0] >= gene[1]:
        print(f"{gene_id}: gene {gene} inicia depois ou durante o seu fim")


def sobrepos_check(genes, gene_id):
    for i, gene_i in enumerate(genes):
        for j, gene_j in enumerate(genes):
            if i != j:
                if gene_i[0] <= gene_j[1] and gene_j[0] <= gene_i[1]:
                    print(f"{gene_id}: gene {gene_i} e gene {gene_j} se sobrepoem")


# check what gene ID have multiple features

# verificar ID por ID, se algum feature se repete,
# caso sim, salve o ID do gene que isso acontece


def check_multiple_features(lines_in, gene_ids_in):
    multiple_features = []
    #features_not_expected = ['protein_match', 'polypeptide']
    features_not_expected = ['protein_match']
    for gene in gene_ids_in:
        features = []
        multipe = False
        for line in lines_in:
            if line[0][0] != "#":
                if ("ID=" + gene) in line[8]:
                    feature = line[2]
                    if feature not in features_not_expected:
                        if feature in features:
                            multipe = True
                        
                        features.append(feature)
        if multipe:
            multiple_features.append([gene, features])
            # print(f"{gene} multi features ({features})")

    genes_ids_txt = ""

    for i in multiple_features:
        genes_ids_txt += f"{i[0]}   {i[1]}\n"

    file = open(strain_id+'_multiples_features_ids.txt', "w")
    file.write(genes_ids_txt)
    file.close()


check_multiple_features(lines, gene_ids)
print("Multiplas features checadas!")
# -------------------------------


# testar os multiplos CDSs
for i in gene_ids:

    gene_i = i
    characteristics = []
    gene_begin = 0
    gene_end = 0
    mrnas = []
    cdss = []
    relatorio = ""

    for j in lines:
        first_char = j[0][0]
        if first_char != "#":
            # print(j)
            if gene_i in j[8]:
                characteristics.append(j)
                if j[2] == 'gene' or j[2] == 'peseudogene':
                    gene_begin = int(j[3])
                    gene_end = int(j[4])
                elif j[2] == "mRNA" or j[2] == "pseudogenic_transcript":
                    mrnas.append([int(j[3]), int(j[4])])
                elif j[2] == "CDS" or j[2] == "pseudogenic_exon":
                    cdss.append([int(j[3]), int(j[4])])

    """
    gene_begin = 0
    gene_end = 2000
    mrnas = [[0,500],[501,1000],[1000,1500]]
    cdss = [[0,501],[501,1000],[1001,2001]]
    """
    # checar erros de estruturas

    if len(mrnas) >= 1:
        # for l in mrnas:

        for l, k in enumerate(mrnas):
            if int(k[0]) < int(gene_begin):
                relatorio += "gene " + gene_i + ": mRNA inicia antes do inicio do gene\n"
            if int(k[1]) > int(gene_end):
                relatorio += "gene " + gene_i + ": mRNA encerra após o fim do gene\n"
            if int(k[1]) <= int(k[0]):
                relatorio += "gene " + gene_i + ": mRNA encerra antes ou no inicio dele mesmo\n"

            if l > 0:
                # print(k)
                if int(k[0]) <= int(mrnas[l-1][1]):
                    relatorio += "gene " + gene_i + ": mRNAs sobrepondo entre eles\n"
        sobrepos_check(mrnas, gene_i)

    if len(cdss) >= 1:
        for l, k in enumerate(cdss):
            if int(k[0]) < int(gene_begin):
                relatorio += "gene " + gene_i + ": CDS inicia antes do inicio do gene\n"
            if int(k[1]) > int(gene_end):
                relatorio += "gene " + gene_i + ": CDS encerra após o fim do gene\n"
            if int(k[1]) <= int(k[0]):
                relatorio += "gene " + gene_i + ": CDS encerra antes ou no inicio dele mesmo\n"

            sobrepos_check(cdss, gene_i)

print(relatorio)
print("Checagem da estrutura de features finalizado")

# -------------------------------

# -------------------------------
# -------------------------------
# -------------------------------
# -------------------------------
# -------------------------------
