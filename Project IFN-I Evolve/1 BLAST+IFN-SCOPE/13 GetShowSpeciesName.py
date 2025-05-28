import codecs
import json
import os
'''
将input_path下的所有数据融合成一个json文件
'''

species="Anas_platyrhynchos"
input_path="myowndata/last_filter"
output_path="myowndata/showdata"#存放前端展示的数据
all_species_path='myowndata'##存放allspecies.txt
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
if __name__ == '__main__':
    names=get_dbfile()
    result=[]
    allspecies = get_allspecies()
    for name in names:
        node={}
        spe = name.replace(".json", "").replace("_",' ')
        for index,item in enumerate(allspecies):
            if spe in item:
                node["value"] = item
                node["label"] = item
        result.append(node)
    json.dump(result, codecs.open(output_path + '/' + "species.json", 'w+', encoding='utf-8'), ensure_ascii=False)
