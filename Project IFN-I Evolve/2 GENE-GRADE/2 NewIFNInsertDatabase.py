import os
import re
import json
import pymongo
# 定义氨基酸到核酸的翻译字典

translation_table = {
    'A': 'GCT', 'R': 'CGT', 'N': 'AAT', 'D': 'GAT',
    'C': 'TGT', 'Q': 'CAA', 'E': 'GAA', 'G': 'GGT',
    'H': 'CAT', 'I': 'ATT', 'L': 'TTA', 'K': 'AAA',
    'M': 'ATG', 'F': 'TTT', 'P': 'CCT', 'S': 'TCT',
    'T': 'ACT', 'W': 'TGG', 'Y': 'TAT', 'V': 'GTT',
    '*': 'TAA'
}
request_str = 'mongodb://root:123456@localhost:27017/?authMechanism=DEFAULT'
path = "myowndata/newIFN"
# allname=path+"/"+"Alligator_mississippiensis_NW_017712884.1_POS6290216.fasta"
client = pymongo.MongoClient(request_str)
unannoIFNGenes = client.interferon.AllunannoIFN
all59IFNNeighbor=json.load(#23-4-1获取59个新干扰素对应保守位点.py生成文件
    open('myowndata/New59IFNNeighbor.json', 'r', encoding="utf-8"))
def getDataFromFastaAndInsertIntoDB(filename,proteinID):
    allname = path + "/" +filename
    with open(allname, "r",encoding='utf-8') as file:
        # 逐行读取文档内容
        lines = file.readlines()
    # print(len(lines))
    # for line in lines:
    #     print(line.strip())  # 使用strip()方法去除行尾的换行符
    # print(lines[0])
    # print(lines[1])
    #datas
    pattern = r">([^,]+),"
    match = re.search(pattern, lines[0])
    datas = match.group(1)
    # print(datas)
    #chromosome
    chromosome=datas.split("_")[-3]+"_"+datas.split("_")[-2]
    #species
    noused="_"+chromosome+"_"+datas.split("_")[-1]
    species=datas.replace(noused,"").replace("_"," ")
    #pos
    pattern = r'POS(\d+)'
    match = re.search(pattern, datas)
    pos = int(match.group(1))
    # print(pos)
    #start
    pattern = r"start:(\d+),"
    match = re.search(pattern, lines[0])
    start = match.group(1)
    #end
    pattern = r"end:(\d+),"
    match = re.search(pattern, lines[0])
    end = match.group(1)
    #direction
    pattern = r"direction:(.*?),"
    match = re.search(pattern, lines[0])
    direction = match.group(1)
    #proteinSeq
    proteinSeq = lines[1].strip()#去除换行符
    # 反向翻译
    nucleotide_sequence = ""
    for amino_acid in proteinSeq:
        nucleotide_sequence += translation_table.get(amino_acid, '')#找不到对应的键就返回空字符串
    neighbor = []
    for itemNode in all59IFNNeighbor:
        if species==itemNode["species"] and chromosome==itemNode["chromosome"] and pos==itemNode["pos"]:
            neighbor=neighbor+itemNode["neighbors"]
            break
    unannoIFNGenes.insert_one({
        'gene': datas,
        'sequence':nucleotide_sequence,#注意，这个逆转录过程非唯一
        'direction': direction,
        # 'exon':1,
        #'exonPos':[],
        'type':0,
        'species': species,
        'chromosome': chromosome,
        'start': int(start),
        'end': int(end),
        'generation':0,
        'proteinID':proteinID,
        'proteinLength':len(proteinSeq),
        'proteinSequence':proteinSeq,
        'neighbors':neighbor,
        'flag':"unanno",
    })

    # print(datas)
    # print(species)
    # print(chromosome)
    # print(start)
    # print(end)
    # print(direction)
    # print(proteinSeq)
    # print("a")
def insertAllNewIFN():
    unannoIFNGenes.delete_many({})
    items = os.listdir(path)
    num=0
    for item in items:
        # print(item)
        if ".DS" not in item:
            num=num+1
            proteinID="XP_"+str(num)
            getDataFromFastaAndInsertIntoDB(item,proteinID)
    print("成功插入:",num)
if __name__ == "__main__":
    insertAllNewIFN()
    # num=0
    # for x in unannoIFNGenes.find({}):
    #     if x["neighbors"]!=[]:
    #         num=num+1
    # print(num)
    # print(len(all59IFNNeighbor))
    # getDataFromFastaAndInsertIntoDB("Alligator_mississippiensis_NW_017712884.1_POS6290216.fasta")


