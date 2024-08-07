"""
Este script corrige erros no arquivo GFF gerado pelo COMPANION

Autor: Fábio Resadore
fresadore.bio@hotmail.com   
"""

cepas = ['2335']
strain_id = "2335"

arquivo = f'../../genbank_submission/{strain_id}/{strain_id}_new.gff3'

file = open(arquivo)
file_lines = file.readlines()
file.close()


def checandoPseudogenes():

    pseudogene_ids = []

    for line in file_lines:
        if '##' not in line:
            line_ = line.strip().split('\t')
            feature_ = line_[2]
            id_ = line_[8].split('ID=')[1].split(';')[0]
            if ':pseudogene' in id_:
                if id_ not in pseudogene_ids:
                    # pseudogene_ids[id_] = []

                    pseudogene_ids.append(id_)

    # print((pseudogene_ids))

    pseudogene_features = {}

    for pseudogene in pseudogene_ids:
        id_ = pseudogene.split(':')[0]

        for line in file_lines:
            if '##' not in line and id_ in line:
                line_ = line.strip().split('\t')
                feature_ = line_[2]

                if id_ not in pseudogene_features.keys():
                    pseudogene_features[id_] = []

                pseudogene_features[id_].append(feature_)

    # print(pseudogene_features)

    gene_check = {}

    for pseudogene in pseudogene_features:
        if 'gene' not in pseudogene_features[pseudogene] or 'mRNA' not in pseudogene_features[pseudogene]:
            gene_check[pseudogene] = 'faltando'

    # print(gene_check)

    if len(gene_check) < 1:
        print('Pseudogenes com features corretas')


def checandoGenes(lines_):
    genes_ids = []

    for line in lines_:
        if '##' not in line:
            line_ = line.strip().split('\t')
            feature_ = line_[2]
            id_ = line_[8].split('ID=')[1].split(';')[0]
            if ':pseudogene' not in id_:
                id_ = id_.split(':')[0].split('.')[0]
                if id_ not in genes_ids:
                    genes_ids.append(id_)

    # print(len(genes_ids))

    gene_features = {}

    for pseudogene in genes_ids:
        id_ = pseudogene.split(':')[0]

        for line in lines_:
            if '##' not in line and id_ in line:
                line_ = line.strip().split('\t')
                feature_ = line_[2]

                if id_ not in gene_features.keys():
                    gene_features[id_] = []

                gene_features[id_].append(feature_)

    # print(gene_features['LshaS1545000005000'])

    gene_check = {}
    gene_remove = []

    for gene in gene_features:
        if 'gene' not in gene_features[gene]:
            gene_check[gene] = 'gene'
            # print(gene)
            gene_remove.append(gene)
        elif 'mRNA' not in gene_features[gene] and 'CDS' not in gene_features[gene]:
            gene_check[gene] = 'mRNA and CDS'
            # print(gene)
            gene_remove.append(gene)

    print(f'Genes incorretos: {len(gene_remove)}')
    # print(gene_check['resadoreS1365000782500'])
    # exit()
    return gene_remove


def arrumandoGff():

    for strain_id in cepas:
        print(f'Cepa: {strain_id}')
        arquivo = f'../../genbank_submission/{strain_id}/{strain_id}_new.gff3'

        file = open(arquivo)
        file_lines = file.readlines()
        file.close()

        genes_remove_ = checandoGenes(file_lines)

        new_gff = ''

        for line in file_lines:
            if '##' in line:
                new_gff += line
            else:
                line_ = line.strip().split('\t')
                # feature_ = line_[2]
                id_ = line_[8].split('ID=')[1].split(
                    ';')[0].split(':')[0].split('.')[0]

                if id_ not in genes_remove_:
                    new_gff += line

        # print(len(new_gff))
        # exit()

        new_gff = (new_gff
                   .replace(' product=', ';product=')
                   .replace(' description=', ';description=')
                   .replace(f'S{strain_id}', f'{strain_id}')
                   )

        output = f'../../genbank_submission/{strain_id}/{strain_id}_new_clean.gff3'

        file = open(output, 'w')
        file.write(new_gff)
        file.close()


arrumandoGff()
