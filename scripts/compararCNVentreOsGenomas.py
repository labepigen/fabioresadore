"""
Este script retorna uma tabela CSV da cobertura do conjunto de genes
de um grupo de genomas sequenciados e anotados.    

Autor: Fábio Resadore
fresadore.bio@hotmail.com   
"""

import sqlite3

con = sqlite3.connect("genes_cov.db")
cur = con.cursor()


try:
    cur.execute(
        "CREATE TABLE genes_bibli(genome, gene_name, orthogroup, cov_median)")
except:
    pass

'''
res = cur.execute("SELECT genome FROM genes_bibli")
res.fetchall() is None

cur.execute("""
    INSERT INTO movie VALUES
        ('Monty Python and the Holy Grail', 1975, 8.2),
        ('And Now for Something Completely Different', 1971, 7.5)
""")
'''


# conjunto de genomas que serão utilizados nessa análise
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

orthogroups = {}


def listandoOrthogroups():

    # arquivo com a lista dos orthogroups
    arquivo = f"Orthogroups2.txt"

    file = open(arquivo)
    file_lines = file.readlines()
    file.close()

    for line in file_lines:
        line_ = line.split(':')
        orthgroup_ = line_[0]
        genes_ = line_[1].strip().split(' ')

        orthogroups[orthgroup_] = genes_


def createDB():

    for genome in cepas:
        # if genome == '3371':

        # arquivo com a mediana de cobertura dos genes anotados de um determinado genoma
        arquivo = f"{genome}/{genome}_depth_medians.txt"

        file = open(arquivo)
        file_lines = file.readlines()
        file.close()

        cov_all_ = file_lines[1].split('\t')
        cur.execute(f"""
                INSERT INTO genes_bibli VALUES
                    ('{genome}', '{cov_all_[0]}', 'all', {cov_all_[1].strip()})
                """)
        con.commit()

        for i, line in enumerate(file_lines):
            if i > 1:
                line_ = line.split('\t')
                gene_ = str(line_[0]).strip()
                cov_median_ = int(line_[1].strip())
                orthogroup_ = ''

                for orthogroup in orthogroups:
                    if (gene_+'.1') in orthogroups[orthogroup]:
                        orthogroup_ = orthogroup
                        break

                cur.execute(f"""
                        INSERT INTO genes_bibli VALUES
                            ('{genome}', '{gene_}', '{orthogroup_}', {cov_median_})
                    """)
                con.commit()


res = cur.execute("DELETE FROM genes_bibli WHERE orthogroup=''")
con.commit()


def gerarTabela():
    table = 'orthogroup\tdesc\t1365\t1545\t2335\t2371\t2490\t2689\t2690\t3283\t3315\t3356\t3358\t3371\t3481\n'

    res = cur.execute("SELECT * FROM genes_bibli")
    database = res.fetchall()

    orthogroup_list = []
    orthogroup_description = {}

    genomes_cov = {}

    for line in database:
        # print(line)
        orthogroup_ = line[2]
        genome_ = line[0]
        gene_name_ = line[1]
        cov_ = line[3]

        if gene_name_ == 'All_genes':
            genomes_cov[genome_] = 0
            genomes_cov[genome_] = cov_
        else:
            if orthogroup_ == 'OG0008160':
                print(line)

            if orthogroup_ not in orthogroups.keys():
                orthogroups[orthogroup_] = {}

                for cepa in cepas:
                    orthogroups[orthogroup_][cepa] = []

            if orthogroup_ not in orthogroup_list:
                orthogroup_list.append(orthogroup_)

            orthogroups[orthogroup_][genome_].append([gene_name_, cov_])

    arquivo_ortho_desc = f"Orthogroup_descriptions_single.txt"

    file = open(arquivo_ortho_desc)
    file_lines = file.readlines()
    file.close()

    for line in file_lines:
        line_ = line.strip().split(':')

        orthogroup_description[line_[0].strip()] = line_[1].strip()

    for orthogroup in orthogroup_list:

        line_ = f'{orthogroup}\t{orthogroup_description[orthogroup]}'
        for cepa in cepas:
            cov_ = 0
            if orthogroups[orthogroup][cepa] != []:
                if len(orthogroups[orthogroup][cepa]) > 1:
                    cov_med = 0
                    for gen_ in orthogroups[orthogroup][cepa]:
                        cov_med += gen_[1]
                    cov_ = cov_med/len(orthogroups[orthogroup][cepa])
                else:
                    gen__ = orthogroups[orthogroup][cepa][0]
                    cov_ = orthogroups[orthogroup][cepa][0][1]

            if cov_ != 0:
                cov_ = int(cov_ / genomes_cov[cepa])

            line_ += f'\t{int(cov_)}'
        line_ += '\n'
        table += line_

    print(table)

    # arquivo CSV do resultado
    file = open('tabela_coverages.txt', "w")
    file.write(table)
    file.close()


# listandoOrthogroups()

# utilize essa função apenas na primeira vez, depois comente novamente:
# createDB()


gerarTabela()


res = cur.execute(
    "SELECT cov_median FROM genes_bibli WHERE genome='3356'")
db = res.fetchall()
print(db)
print(len(db))

con.close()
