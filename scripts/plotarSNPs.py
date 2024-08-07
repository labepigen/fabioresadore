"""
Este script utiliza os resultados do software SnpEff, e 
plota em gráficos os resultados de diversidade genética.

Autor: Fábio Resadore
fresadore.bio@hotmail.com   
"""

import seaborn as sns
import matplotlib.pyplot as plt

cepas = [
    '1365',
    '1545',
    '2335',
    '2371',
    '2490',
    '2689',
    '2690',
    '3283',
    '3315',
    '3356',
    '3358',
    '3371',
    '3481',
]

chr_names = [
        '01',
        '02',
        '03',
        '04',
        '05',
        '06',
        '07',
        '08',
        '09',
        '10',
        '11',
        '12',
        '13',
        '14',
        '15',
        '16',
        '17',
        '18',
        '19',
        '20',
        '21',
        '22',
        '23',
        '24',
        '25',
        '26',
        '27',
        '28',
        '29',
        '30',
        '31',
        '32',
        '33',
        '34',
        '35'
    ]


def extraindoSNPsCromossomos():
    snps_count = {}  # snps_count[cepa][chr] = [total, densidade, tamanho]

    arquivo = f'snps_cromossomos.txt'

    file = open(arquivo)
    file_lines = file.readlines()
    file.close()

    cepa_ = ''
    chr_count_ = 1
    for line in file_lines:
        line_ = line.strip()

        if line_[0] == 'S':
            cepa_ = line_[1:]
            snps_count[cepa_] = {}
            chr_count_ = 1
        else:
            total_ = int(line_.split('\t')[2].replace(',',''))
            densidade_ = int(line_.split('\t')[3].replace(',',''))
            comprimento_ = int(line_.split('\t')[1].replace(',',''))
            snps_count[cepa_][f'chr.{chr_count_}'] = [total_, densidade_, comprimento_]
            chr_count_ += 1

    return snps_count


def plotandoSNPsQuantidade():
    snps_valores_ = extraindoSNPsCromossomos()

    plot_matrix_total_ = []

    for cepa in cepas:
        chr_snps_ = []
        for chr in range(35):
            chr_ = f'chr.{chr+1}'
            chr_snps_.append(snps_valores_[cepa][chr_][0])
        plot_matrix_total_.append(chr_snps_)
    
    plt.figure(figsize=(14.5, 4.5))
    hm_plot = sns.heatmap(plot_matrix_total_, linewidth=0.5,
                          cmap='Wistia', xticklabels=chr_names, yticklabels=cepas)
    # plt.pcolormesh( plot_matrix_ , cmap = 'Wistia' )
    # ticks, labels = plt.xticks()
    # plt.xticks(ticks, labels=chr_names)
    plt.title("Quantidade total de SNVs em cada cromossomo",
              loc="center", fontsize=18)
    plt.xlabel("Cromossomo")
    plt.ylabel("Cepa")
    plt.show()
    
    plot_matrix_densidade_ = []

    for cepa in cepas:
        chr_snps_ = []
        for chr in range(35):
            chr_ = f'chr.{chr+1}'
            tamanho_ = snps_valores_[cepa][chr_][2]
            snps_ = snps_valores_[cepa][chr_][0]
            densidade_ = snps_ / tamanho_
            
            chr_snps_.append(densidade_)
        plot_matrix_densidade_.append(chr_snps_)
    
    plt.figure(figsize=(14.5, 4.5))
    hm_plot = sns.heatmap(plot_matrix_densidade_, linewidth=0.5,
                          cmap='Wistia', xticklabels=chr_names, yticklabels=cepas)
    # plt.pcolormesh( plot_matrix_ , cmap = 'Wistia' )
    # ticks, labels = plt.xticks()
    # plt.xticks(ticks, labels=chr_names)
    plt.title("Densidade de SNVs em cada cromossomo",
              loc="center", fontsize=18)
    plt.xlabel("Cromossomo")
    plt.ylabel("Cepa")
    plt.show()


plotandoSNPsQuantidade()


def mapearCromossomosGenes():
    genes_cromossomos = {}  # genes_cromossomos[cepa][gene] = "chr"
    cepa_cromossomes_names = {}  # cepa_cromossomes_names[cepa] = [chrs]

    for cepa in cepas:
        genes_cromossomos[cepa] = {}
        cepa_cromossomes_names[cepa] = []

        arquivo = f'../../../genbank_submission/{cepa}/{cepa}_new.gff3'

        file = open(arquivo)
        file_lines = file.readlines()
        file.close()

        for line in file_lines:
            if "##" not in line:
                line_ = line.strip().split('\t')
                gene_ = line_[8].split("ID=")[1].split(
                    ';')[0].split('.')[0].split(':')[0][-9:]
                if gene_ not in genes_cromossomos[cepa].keys():
                    chr_ = line_[0].strip()
                    if '20.1' in chr_ or '20.2' in chr_:
                        chr_ = chr_[:-2]
                    genes_cromossomos[cepa][gene_] = chr_
                    if chr_ not in cepa_cromossomes_names[cepa]:
                        cepa_cromossomes_names[cepa].append(chr_)
        # print(cepa_cromossomes_names)
        # exit()
    # print(cepa_cromossomes_names)
    return [genes_cromossomos, cepa_cromossomes_names]


def extrairCromossomoImpacto():
    mapeamento_ = mapearCromossomosGenes()
    snp_gene_impact = {}  # snp_gene_impact[cepa][chr] = 123

    for cepa in cepas:
        snp_gene_impact[cepa] = {}

        for chr in mapeamento_[1][cepa]:

            snp_gene_impact[cepa][chr] = 0
        # print(snp_gene_impact)
        # exit()
        arquivo = f'{cepa}/snpEff_genes.txt'

        file = open(arquivo)
        file_lines = file.readlines()
        file.close()

        # contabilizando  o numero de impactos moderado + alto por cromossomo
        for i, line in enumerate(file_lines):
            if i > 1:
                line_ = line.strip().split('\t')
                gene_ = line_[0].split(' ')[0].split('.')[0].split(':')[0][-9:]
                low_ = int(line_[6])
                moderate_ = int(line_[7])
                high_ = int(line_[5])
                # print(gene_)
                #print(gene_)
                if gene_ in mapeamento_[0][cepa].keys():
                    chr_ = mapeamento_[0][cepa][gene_]
                    """if '20.1' in chr_:
                        print(chr_[:-2])
                        exit()"""
                    snp_gene_impact[cepa][chr_] += moderate_ + high_
        # print(snp_gene_impact)
        # exit()
    return snp_gene_impact


def plotandoGraficos():
    mapeamento_ = mapearCromossomosGenes()
    impactos_ = extrairCromossomoImpacto()

    

    # plotando box plot por cepa
    plot_matrix_ = []  # [cepa[chr[numero_impactos]]]

    for cepa in cepas:
        chr_impacts_ = []
        for i, chr in enumerate(mapeamento_[1][cepa]):
            if i > 0:
                chr_impacts_.append(impactos_[cepa][chr])
        plot_matrix_.append(chr_impacts_)

    plt.figure(figsize=(14.5, 4.5))
    hm_plot = sns.heatmap(plot_matrix_, linewidth=0.5,
                          cmap='Wistia', xticklabels=chr_names, yticklabels=cepas)
    # plt.pcolormesh( plot_matrix_ , cmap = 'Wistia' )
    # ticks, labels = plt.xticks()
    # plt.xticks(ticks, labels=chr_names)
    plt.title("Quantidade de SNVs com possíveis impactos estruturais",
              loc="center", fontsize=18)
    plt.xlabel("Cromossomo")
    plt.ylabel("Cepa")
    plt.show()
    # exit()
    # plotando por cromossomo

    plot_matrix_ = []  # [chr[cepa[numero_impactos]]]

    for chr in range(36):
        if chr > 0:
            cepa_impacts_ = []
            for cepa in cepas:
                chr_name_ = mapeamento_[1][cepa][chr]
                cepa_impacts_.append(impactos_[cepa][chr_name_])
            plot_matrix_.append(cepa_impacts_)

    plt.figure(figsize=(14, 4))
    bplots = plt.boxplot(plot_matrix_,  vert=1, patch_artist=False)
    plt.title("Quantidade de SNVs com impactos estruturais para cada Cromossomo",
              loc="center", fontsize=18)
    plt.xlabel("Cromossomo")
    ticks, labels = plt.xticks()
    plt.xticks(ticks, labels=chr_names)
    plt.ylabel("Quantidade de SNVs")
    plt.show()


plotandoGraficos()
