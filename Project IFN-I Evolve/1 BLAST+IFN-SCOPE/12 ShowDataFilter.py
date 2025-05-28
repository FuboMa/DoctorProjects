import codecs
import json
'''
1 过滤数据
按照seq的长度过滤
更新blast_result.json
2  将每个type="seq"的序列的节点
新增一个属性show_name="LOC104391698"
这个属性的属性值是从name中获得的Chaetura_pelagica_LOC104391698
3  获得detail页面的基因列表选择数据
'''
input_path="myowndata/showdata"
output_path="myowndata/showdata"
t_length=120#长度小于这个值的seq全部过滤，刚好=t_length这个值的是保留了的
def get_json():
    path=input_path+'/'+'blast_result.json'
    data = json.load(open(path, encoding="utf-8"))
    return data
def filter(thresh_length):
    datas=get_json()
    result=[]
    for data in datas:
        r_node={}
        species=data["species"]
        full_name=data["full_name"]
        chr_node=[]
        #每个child代表一个染色体
        for child in data["chr_node"]:
            item={}
            length=child["length"]
            chromosome = child["chromosome"]
            chr=child["chr"]
            children=[]
            for seq_node in child["children"]:
                if seq_node["type"]=='seq':
                    #小于长度阈值直接不加入
                    seq_node["show_name"]=seq_node["name"].split("_")[-1]
                    if seq_node["seq_length"]<thresh_length:
                        continue
                children.append(seq_node)
            if(len(children)==0):
                continue
            item["length"]=length
            item["chromosome"] = chromosome
            item["chr"] = chr
            item["children"] = children
            chr_node.append(item)
        r_node["species"]=species
        r_node["full_name"] = full_name
        r_node["chr_node"] = chr_node
        result.append(r_node)
    json.dump(result, codecs.open(output_path + '/' + "blast_result.json", 'w+', encoding='utf-8'), ensure_ascii=False)

def get_detail_show_genename_list():
    datas = get_json()
    result = []
    for data in datas:
        for chr_node in data["chr_node"]:
            for child in chr_node["children"]:
                if(child["type"]=="seq"):
                    name_node={}
                    ssname=child["name"].split("_")[-1]
                    name_node["value"]=ssname
                    name_node["label"] = ssname
                    if name_node not in result:
                        result.append(name_node)
    json.dump(result, codecs.open(output_path + '/' + "detail_list.json", 'w+', encoding='utf-8'), ensure_ascii=False)

if __name__ == '__main__':
    filter(t_length)
    get_detail_show_genename_list()