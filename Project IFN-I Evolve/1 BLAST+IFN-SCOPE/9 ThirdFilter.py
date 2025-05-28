import codecs
import json
import os
'''
将序列结果进行和基因进行比对，如果该序列位于基因内部，就过滤掉
第三次过滤之后一般就没有重复序列了，但是某些序列可能seq_length很短，如"seq_length":21

'''
'''
修改三个地方ifn_path,input_path,output_path
读取interferons.json中的数据（这里面包含95个物种的干扰素数据 0：I型干扰素 1：其他干扰素 3：常见相邻基因）
interferons.json存储在ifn_path下
将input_path下的json文件(Anas_platyrhynchos.json)文件中的数据过滤（将这些序列和interferons.json中的基因比较）
结果存放在output_path中
'''
species="Anas_platyrhynchos"
ifn_path="myowndata"
input_path="myowndata/second_filter"#存放95个物种文件夹（每个文件夹表示一个物种的第一次数据处理结果）
output_path="myowndata/third_filter"#存放第三次过滤的数据
overlap=20#和I型干扰素的重叠阈值
def get_dbfile():#获取"myowndata/first_filter下的所有文件名
    datanames = os.listdir(input_path)
    name=[]
    for i in datanames:
        if '.json' in i:
            name.append(i)
    return name
def swap(a,b):
    if(a<=b):
        return a,b
    else:
        return b,a
def get_json():
    path=ifn_path+'/'+'interferons0.json'
    with open(path) as f:
        data = json.load(f)
    return data
def get_json1(species):
    path=input_path+'/'+species+'.json'
    with open(path) as f:
        data = json.load(f)
    return data
def get_species_gene(species):
    datas=get_json()
    for data in datas:
        if species==data["species"]:
            return data["interferon"]  #这个  data["interferon"]包含三类数据 0类 1类 3类
def get_seq(species):
    datas=get_json1(species)
    return datas
def third(species):
    genes=get_species_gene(species)
    seqs=get_seq(species)
    if (species == "Anas_platyrhynchos"):
        print(genes)
    # print(genes[0])
    # print(genes[1])
    # print(seqs[0])
    # print(seqs[1])
    # print('seq',len(seqs))
    result=[]
    for seq in seqs:#每个seq代表一个染色体数据
        now_children=seq["children"][:]
        #遍历该物种的每一个基因（type=0 1 3的基因 ）
        for gene in genes:
            new_children = []
            if seq["genome"]==gene["chromosome"]:
                gene_start1=gene["start"]
                gene_end1 = gene["end"]
                gene_start,gene_end=swap(gene_start1,gene_end1)
                # gg_type=gene["type"]
                for child in now_children:
                    seq_start1=child["seq_start"]
                    seq_end1 = child["seq_end"]
                    seq_start,seq_end=swap(seq_start1, seq_end1)
                    #这个if是不允许交错的if
                    # if (gene_start>seq_start and gene_start<seq_end )or(gene_end>seq_start and gene_end<seq_end)or(gene_start==seq_start)or(gene_end==seq_end)or(gene_start<seq_start and gene_end>seq_end):
                    #     pass
                    # #这个if设置了重叠度overlap
                    if (seq_start < gene_start and seq_end> gene_start and seq_end<gene_end and seq_end - gene_start > overlap) or (
                            seq_start > gene_start and seq_end < gene_end) or (
                            seq_start > gene_start and seq_start<gene_end and seq_end>gene_end and seq_end - gene_end > overlap):
                        pass
                    else:
                        #这个序列和基因不存在交错也不位于基因内部，直接加入
                        new_node={}
                        new_node.update(child)
                        new_children.append(new_node)
                now_children=new_children[:]
        new_seq={}
        new_seq.update(seq)
        new_seq["children"]=now_children
        #存在序列的染色体才加入，如果该染色体上所有序列都被过滤掉了，那么就这个染色体就不加入result
        if len(new_seq["children"])!=0:
            result.append(new_seq)
    json.dump(result, codecs.open(output_path + '/' + species + ".json", 'w+', encoding='utf-8'), ensure_ascii=False)


if __name__ == '__main__':
    names=get_dbfile()
    # print(names)
    for name in names:
        spe=name.replace(".json","")
        third(spe)