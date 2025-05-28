import codecs
import json
import os
'''
主要进行数据合并 将过滤后的比对结果和基因数据合并
将input_path下的所有序列json数据和interferons-03中对应物种的基因数据合并
若是基因数据位于单独的染色体（即该染色体上不存在任何比对片段），那么该基因数据就不显示
'''
'''
'''
species="Anas_platyrhynchos"
ifn_path="myowndata"
input_path="myowndata/third_filter"
output_path="myowndata/last_filter"#存放第三次过滤的数据
def get_dbfile():#获取"myowndata/first_filter下的所有文件名
    datanames = os.listdir(input_path)
    name=[]
    for i in datanames:
        if '.json' in i:
            name.append(i)
    return name
def get_seq(species):
    path=input_path+'/'+species+'.json'
    with open(path) as f:
        data = json.load(f)
    return data
def get_json():
    path=ifn_path+'/'+'interferons-03.json'
    with open(path) as f:
        data = json.load(f)
    return data
def get_species_gene(species):
    datas=get_json()
    for data in datas:
        if species==data["species"]:
            return data["interferon"]  #这个  data["interferon"]包含三类数据 0类 3类
def get_last_data(species):
    seqs=get_seq(species)
    genes=get_species_gene(species)#只包含 0 3类基因
    gene_in=[]#存放已经加入了的gene的索引
    result=[]
    for seq in seqs:#遍历每一个染色体
        new_seq_node = {}
        new_seq_node["length"] = seq["length"]
        new_seq_node["chromosome"] = seq["genome"]
        new_seq_node["chr"] = seq["chr"]
        new_seq_node["children"] = []
        # 处理每一个child
        for child in seq["children"]:
            child_node = {}
            child_node["type"] = "seq"
            child_node["name"] = child["spe_name"] + '_' + child["interferon"]
            child_node["identity"] = child["identity"]
            child_node["positive"] = child["positive"]
            child_node["gap"] = child["gap"]
            child_node["seq_length"] = child["seq_length"]
            child_node["seq_start"] = child["seq_start"]
            child_node["seq_end"] = child["seq_end"]

            child_node["query_seq_content"] = []
            child_node["sbjct_seq_conten"] = []
            child_node["flag_seq_content"] = []
            for char1 in child["query_seq_content"]:
                child_node["query_seq_content"].append(char1)
            for char2 in child["sbjct_seq_content"]:
                child_node["sbjct_seq_conten"].append(char2)
            for char3 in child["flag_seq_content"]:
                child_node["flag_seq_content"].append(char3)
            # 此节点加入
            new_seq_node["children"].append(child_node)
        #遍历每一个基因，找到需要加入的基因
        for index,gene in enumerate(genes):
            #该基因找到了对应染色体
            if seq["genome"] == gene["chromosome"]:
                gene_in.append(index)  # 代表此基因已经加入了
                gene_node={}
                gene_node["type"] = "gene"
                gene_node["name"] = gene["gene"]
                gene_node["inner_type"] = gene["type"]
                gene_node["start"] = gene["start"]
                gene_node["end"] = gene["end"]
                gene_node["exon"] = gene["exon"]
                gene_node["product"] = gene["product"]
                gene_node["direction"] = gene["direction"]
                new_seq_node["children"].append(gene_node)
        result.append(new_seq_node)
    #先不做处理（位于其他染色体的基因不进行显示）
    #处理位于其他染色体的基因（即该基因所处染色体与当前比对序列所处的染色体都不相同）
    # for index, gene in enumerate(genes):
    #     if index not in gene_in: #这个基因所在的染色体还没有加入
    #         new_chr_node = {}
    #        # new_seq_node["length"] = seq["length"]
    #         new_chr_node["chromosome"] = gene["chromosome"]
    #         new_chr_node["chr"] = gene["chromosome"]
    #         new_chr_node["children"] = []
    #         gene_node = {}
    #         gene_node["type"] = "gene"
    #         gene_node["name"] = gene["gene"]
    #         gene_node["inner_type"] = gene["type"]
    #         gene_node["start"] = gene["start"]
    #         gene_node["end"] = gene["end"]
    #         gene_node["exon"] = gene["exon"]
    #         gene_node["product"] = gene["product"]
    #         gene_node["direction"] = gene["direction"]
    #         new_seq_node["children"].append(gene_node)
    json.dump(result, codecs.open(output_path + '/' + species + ".json", 'w+', encoding='utf-8'), ensure_ascii=False)


if __name__ == '__main__':
    names=get_dbfile()
    for name in names:
        spe=name.replace(".json","")
        get_last_data(spe)