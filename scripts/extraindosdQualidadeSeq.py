"""
Este script utiliza os resutlados do software fastp, e
extrai resultados específicos, sanvando em um arquivo CSV

Autor: Fábio Resadore
fresadore.bio@hotmail.com   
"""

from bs4 import BeautifulSoup
import json
# cepas = ['1365', '1545', '2689', '2690', '3481']
cepas = ['1365', '1545', '2335', '2371', '2490', '2689',
         '2690', '3283', '3315', '3356', '3358', '3371', '3481']


cepa = '3315'
# strain_id = "1545"

dados_seq = {}
for cepa in cepas:
    dados_seq[cepa] = {}

    json_file = f'S{cepa}/fastp.json'
    json_ = open(json_file)
    dados_ = json.load(json_)

    dados_seq[cepa] = dados_['summary']

tabela_ = ''
tabela_ += f'IOCL\tbefore_bases\tbefore_reads\tbefore_q20\tbefore_q30\tafter_bases\tafter_reads\tafter_q20\tafter_q30\n'

for cepa in cepas:
    line_ = f'{cepa}\t'

    before_ = dados_seq[cepa]['before_filtering']
    after_ = dados_seq[cepa]['after_filtering']

    line_ += f"{before_['total_bases']}\t"
    line_ += f"{before_['total_reads']}\t"
    line_ += f"{before_['q20_bases']}\t"
    line_ += f"{before_['q30_bases']}\t"
    line_ += f"{after_['total_bases']}\t"
    line_ += f"{after_['total_reads']}\t"
    line_ += f"{after_['q20_bases']}\t"
    line_ += f"{after_['q30_bases']}\t"

    tabela_ += f'{line_}\n'

out = f'tabela_qualidade.txt'
file = open(out, 'w')
file.write(tabela_)
file.close()
# print(dados_seq)
# arquivo = f'../../genbank_submission/{strain_id}/{strain_id}_new.gff3'
