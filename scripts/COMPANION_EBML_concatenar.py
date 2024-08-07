"""
Este script utiliza os arquivo EMBL gerados pelo COMPANION, 
e concatena todos em um único arquivo EMBL

Autor: Fábio Resadore
fresadore.bio@hotmail.com   
"""


strain_prefix = "3481"
sp_abrev = "Lsha"

caminho = "anotacoes/resadoreS"+strain_prefix+"-refLbraz/"
arquivo = caminho+"Lsha_S"+strain_prefix
arquivos = caminho+"Lsha_S"+strain_prefix+"_LbrM."
embl = []
embl_text = ""
embl_texts = []
cromossomos_list = ["01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20.1","20.2","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35"]

cepa = "IOCL"+strain_prefix
especie = "Leishmania shawi"
nome_amostra = "Lshaw_IOCL"+strain_prefix
locus_tag = "IOCL3481"

novo_embl = caminho+cepa+"_all.embl"


# <chr_leng> = 2343432
# <nome_amostra> = Lnaif_IOCL1365
## <chr_numero> = 1
# <especie> = Leishmania naiffi
## resadoreS1365_ >>>> IOCL1365_
# <nome_amostra+chr_numero> = <nome_amostra> + "." + chrm
# <entry_name> = Lnaif_IOCL1365.chr01


cabecalho = "ID   XXX; XXX; linear; genomic DNA; XXX; XXX; <chr_leng> BP.\nXX\nAC * _<entry_name> \nXX\nPR   Project:PRJEB64591;\nXX\nDT   25-JUL-2023 (Rel. 133, Created)\nXX\nDE   <especie>, <nome_amostra+chr_numero>.\nXX\nKW   .\nXX\nOS   <especie>\nXX\nOC   cellular organisms; Eukaryota; Discoba; Euglenozoa; Kinetoplastea;\nOC   Metakinetoplastina; Trypanosomatida; Trypanosomatidae; Leishmaniinae;\nOC   Leishmania; Viannia.\nXX\nRN   [1]\nRP   1-<chr_leng>\nRG   XXX\nRT   ;\nRL   Submitted (25-JUL-2023) to the INSDC.\nXX\n"

def calcular_introns(texto_intron):
    line_type = ""

    if "complement" in texto_intron:
        line_type = "complement"
    else:
        line_type = "non complement"

    intron = texto_intron.split("join")[1]
    intron = intron.replace("(", "")
    intron = intron.replace(")", "")
    introns = intron.split(",")

    for i in range(len(introns)):
        introns[i] = introns[i].split("..")

    distancias = []

    for i in range(len(introns)-1):
        distancias.append(int(introns[i+1][0]) - int(introns[i][1]))

    new_intron = []
    new_intron.append([introns[0][0], "0"])

    for i in range(len(distancias)):
        if distancias[i] >= 10:
            new_intron[len(new_intron)-1][1] = introns[i][1]
            new_intron.append([introns[i+1][0], "0"])
        else:
            new_intron[len(new_intron)-1][1] = introns[i+1][1]

    new_intron[-1][1] = introns[-1][1]

    new_intron_text = "FT   CDS             "

    if line_type == "complement":
        new_intron_text += "complement("
    
    if len(new_intron) > 1:
        new_intron_text += "join("

    for i in range(len(new_intron)):


        new_intron_text += new_intron[i][0]
        new_intron_text += ".."
        new_intron_text += new_intron[i][1]
        new_intron_text += ","

    new_intron_text = new_intron_text[:-1]

    if line_type == "complement":
        new_intron_text += ")"

    if len(new_intron) > 1:
        new_intron_text += ")"

    return new_intron_text



### primeiro arquivo, n cromossomo
f = open(arquivo+".embl", "r")

for x in f:
    embl.append(x)
f.close()

for line in range(17, len(embl)):
    if embl[line] != "\n":
        if "/colour" not in embl[line]:
            if "/systematic_id" not in embl[line]:
                if "/ortholog" not in embl[line]:
                    if "/translation" not in embl[line]:
                        if "join(" in embl[line]:
                            embl_text = embl_text + calcular_introns(embl[line].strip())  + "\n"
                        else:
                            embl_text = embl_text+embl[line]

embl_text = cabecalho + embl_text
#print(embl_text)

cromossomo_tamanho = embl[19].split("..")[1].strip()

embl_text = embl_text.replace("resadoreS"+strain_prefix+"_", cepa+"_")
embl_text = embl_text.replace("<chr_leng>", cromossomo_tamanho)
embl_text = embl_text.replace("<especie>", especie)
embl_text = embl_text.replace("<nome_amostra+chr_numero>", nome_amostra)
embl_text = embl_text.replace("<entry_name>", sp_abrev+"_IOCL"+strain_prefix+"_scaffold")

embl_texts.append(embl_text)

# arquivos de cromossomos
for chrm in cromossomos_list:
    pass
    #print(chrm)
    embl_text = ""
    embl = []

    f = open(arquivos+chrm+".embl", "r")

    for x in f:
        embl.append(x)
    f.close()

    for line in range(17, len(embl)):
        if embl[line] != "\n":
            if "/colour" not in embl[line]:
                if "/systematic_id" not in embl[line]:
                    if "/ortholog" not in embl[line]:
                        if "/translation" not in embl[line]:
                            if "join(" in embl[line]:
                                embl_text = embl_text + calcular_introns(embl[line].strip())  + "\n"
                            else:
                                embl_text = embl_text+embl[line]

    embl_text = cabecalho + embl_text

    cromossomo_tamanho = embl[19].split("..")[1].strip()

    embl_text = embl_text.replace("resadoreS"+strain_prefix+"_", cepa+"_")
    embl_text = embl_text.replace("<chr_leng>", cromossomo_tamanho)
    embl_text = embl_text.replace("<especie>", especie)
    embl_text = embl_text.replace("<nome_amostra+chr_numero>", nome_amostra+"."+chrm)
    embl_text = embl_text.replace("<entry_name>", sp_abrev+"_IOCL"+strain_prefix+".chr" + chrm)

    embl_texts.append(embl_text)


embl_text = ""
for embl_lista in embl_texts:
    embl_text = embl_text + embl_lista

embl_text = embl_text.replace('/locus_tag="'+sp_abrev+'_S'+ strain_prefix +'_', '/locus_tag="'+locus_tag+'_')

f = open(novo_embl, "w")
f.write(embl_text)
f.close()