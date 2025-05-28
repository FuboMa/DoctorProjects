import codecs
import json
import os
'''
将input_path下的所有数据融合成一个json文件blast_result.json
'''
all_species_path='myowndata'#存放allspecies.txt，将中文名字加入Anas platyrhynchos_绿头鸭
species="Anas_platyrhynchos"
input_path="myowndata/last_filter"
output_path="myowndata/showdata"#存放前端展示的数据
def get_allspecies():
    file_data=[]
    fullname=all_species_path+'/'+'allspecies.txt'
    with open(fullname, "r", encoding="utf-8") as f:
        for line in f.readlines():
            if(line!='\n'):
                line=line.replace("\n","")
                file_data.append(line)
    return file_data
def get_dbfile():#获取"myowndata/first_filter下的所有文件名
    datanames = os.listdir(input_path)
    name=[]
    for i in datanames:
        if '.json' in i:
            name.append(i)
    return name
def get_json(name):
    path=input_path+'/'+name
    with open(path) as f:
        data = json.load(f)
    return data
def join_one():
    result=[]
    names=get_dbfile()
    allspecies=get_allspecies()
    for name in names:
        node={}
        data=get_json(name)
        spe=name.replace(".json","").replace("_",' ')
        node["species"] = spe
        for index,item in enumerate(allspecies):
            if spe in item:
                node["full_name"] = item
        node["chr_node"] = data[:]
        result.append(node)
    json.dump(result, codecs.open(output_path + '/' + "blast_result.json", 'w+', encoding='utf-8'), ensure_ascii=False)


if __name__ == '__main__':
    join_one()
    # names=get_dbfile()
    # for name in names:
    #     spe=name.replace(".json","")
    #     print(spe)

