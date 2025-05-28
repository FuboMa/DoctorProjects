import codecs
import json
import os
'''
只需要修改input_path和output_path

第二次过滤,
先将同一个数据库文件（Anas_platyrhynchos）中的所有干扰素的比对结果的一次过滤后的json文件整合到一个json中
整个文件是一个中间文件，名字为Anas_platyrhynchos1.json
再去除Anas_platyrhynchos1.json的重复序列
在output_path文件夹下生成Anas_platyrhynchos.json


将input_path文件夹下的所有物种文件全部处理一遍
过滤结果存放在output_path文件夹下
'''
species="Anas_platyrhynchos"
input_path="myowndata/first_filter"#存放95个物种文件夹（每个文件夹表示一个物种的第一次数据处理结果）
output_path="myowndata/second_filter"#存放第二次过滤的数据
def swap(a,b):
    if(a<=b):
        return a,b
    else:
        return b,a
def create_dir_not_exist(path):
    if not os.path.exists(path):
        os.mkdir(path)
def get_species(species):#获取"myowndata/first_filter/Anas_platyrhynchos"下的所有json文件
    datanames = os.listdir(input_path+'/'+species)
    name=[]
    for i in datanames:
        name.append(i)
    return name
def get_dbfile():#获取"myowndata/first_filter下的所有文件名
    datanames = os.listdir(input_path)
    name=[]
    for i in datanames:
        name.append(i)
    return name

def get_json(species,interferon):
    path=input_path+'/'+species+'/'+interferon#这个interferon是包含.json后缀的
    with open(path) as f:
        data = json.load(f)
    return data
def get_one_json(species):
    path=output_path+'/'+species+'/'+species+"1.json"
    with open(path) as f:
        data = json.load(f)
    return data
#将一个物种下的所有json文件合并成一个json
def in_one(species):
    names = get_species(species)
    result=[]
    for name in names:#name表示Anas_platyrhynchos下面的各个json文件
        datas=get_json(species, name)
        name=name.replace(".json","")
        interferon=name.split("_")[-1]
        spe_name=name.replace("_"+interferon,"")
        #node.update(child)
        #datas是一个列表，其中每个数据是一个字典，代表一个染色体数据
        for data in datas:
            node = {}
            # node["interferon"]=interferon
            # node["spe_name"] = spe_name
            node.update(data)
            for index, child in enumerate(node["children"]):
                node["children"][index]["interferon"]=interferon
                node["children"][index]["spe_name"] = spe_name
            result.append(node)
    ppath=output_path + '/' +species
    create_dir_not_exist(ppath)
    #生成的Anas_platyrhynchos1.json文件是合并的文件还没有进行二次过滤
    json.dump(result, codecs.open(ppath+'/'+species+"1.json", 'w+', encoding='utf-8'), ensure_ascii=False)

def second(species):
    #result中的每个元素中'genome'属性都是唯一的
    result=[]
    datas=get_one_json(species)
    for data in datas:
        if len(result)==0:
            result.append(data)
        else:
            #是否在result中找到相同的染色体
            flag=0
            for old in result:#old表示已经存在result中的染色体
                if old["genome"]==data["genome"]:
                    flag=1
                    #result中已经存在此染色体的相关数据
                    for new_node in data["children"]:
                        fflag = 0  # 该序列片段是否和result中相同染色体中片段存在交错
                        identity = new_node["identity"].replace("%", "")
                        identity = int(identity)
                        a = new_node["seq_start"]
                        b = new_node["seq_end"]
                        seq_start, seq_end = swap(a, b)  # 第一个返回值是小的数字
                        seq_length=new_node["seq_length"]
                        for old_index, old_node in enumerate(old["children"]):
                            sa = old_node["seq_start"]
                            sb = old_node["seq_end"]
                            small,big=swap(sa, sb)
                            i=old_node["identity"].replace("%","")
                            i=int(i)
                            s=old_node["seq_length"]
                            # if (seq_start > small and seq_start < big) or (seq_end > small and seq_end < big) or (seq_start == small) or (seq_end == big):
                            #存在交错，可能会进行更新
                            if(seq_start > small and seq_start < big)or(seq_end > small and seq_end < big)or\
                                    (small>seq_start and small<seq_end)or(big>seq_start and big<seq_end) or(small==seq_start)or(big==seq_end):
                                fflag=1
                                #保留identity大的那个片段
                                if identity > i:
                                    old["children"][old_index].update(new_node)
                                    break  # 更新了之后，直接结束整个for循环
                                if identity==i:
                                    #identity相等保留seq_length长的那一个
                                    if seq_length>s:
                                        old["children"][old_index].update(new_node)
                                    break  # 更新了之后，直接结束整个for循环

                        if fflag == 0:  # 没有交错直接将该段序列加入
                            old["children"].append(new_node)
                    break#找到一次直接break
            #result不存在此染色体，直接加入
            if flag==0:
                result.append(data)
    json.dump(result, codecs.open(output_path + '/' + species + ".json", 'w+', encoding='utf-8'), ensure_ascii=False)
if __name__ == '__main__':
    names=get_dbfile()
    print(names)
    for name in names:
        in_one(name)
        second(name)
