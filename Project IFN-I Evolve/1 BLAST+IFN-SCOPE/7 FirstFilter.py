import codecs
import json
import os
'''
只需要修改input_path和output_path

第一次过滤,去除单个比对结果中的重复序列
将input_path文件夹下的所有物种文件全部处理一遍
过滤结果存放在output_path文件夹下
'''
species="Anas_platyrhynchos"
input_path="myowndata/blast_json"#存放95个物种文件夹（每个文件夹表示一个数据库的比对结果）
output_path="myowndata/first_filter"#存放第一次过滤的数据
def create_dir_not_exist(path):
    if not os.path.exists(path):
        os.mkdir(path)
def get_json(species,interferon):
    path=input_path+'/'+species+'/'+interferon#这个interferon是包含.json后缀的
    with open(path) as f:
        data = json.load(f)
    return data
def swap(a,b):
    if(a<=b):
        return a,b
    else:
        return b,a

def first(species):
    path=input_path+'/'+species
    interferons = os.listdir(path)
    #遍历Anas_platyrhynchos文件夹下的所有json文件
    #每一个interferon表示一个json文件
    for interferon in interferons:
        #此时interferon是带json后缀的
        datas=get_json(species,interferon)
        #interferon=interferon.replace(".json","")
        print(interferon)
        # print(len(datas))
        # print(datas)
        result=[]#result用于存储过滤后的染色体节点，其中每一个元素表示一个染色体节点
        for data in datas:
        # datas表示一个json文件中的内容，它是一个列表，其中每个元素data表示一个染色体
            node={}
            node["chr"]=data["chr"]
            node["genome"] = data["genome"]
            node["length"] = data["length"]
            node["children"] = []#存放不重复的序列
            #result.append(node)
            #每个data表示一个染色体的数据
            print(data)
            #处理data["children"],获得不重复的序列列表node["children"]
            for child in data["children"]:
                #每个child表示一段相似序列
                if child:#child不是空字典
                    if len(node["children"])==0:
                        #node["children"]列表为空，直接加入child
                        node["children"].append(child)
                    else:
                        identity = child["identity"].replace("%","")
                        identity=int(identity)
                        a = child["seq_start"]
                        b = child["seq_end"]
                        seq_start, seq_end = swap(a, b)  # 第一个返回值是小的数字
                        flag=0#是否遍历完了整个node["children"]也没有找到交错序列
                        for index, n in enumerate(node["children"]):
                            a = n["seq_start"]
                            b = n["seq_end"]
                            small,big=swap(a, b)
                            i=n["identity"].replace("%","")
                            i=int(i)
                            if (seq_start>small and seq_start<big) or (seq_end>small and seq_end<big):
                                #存在交错
                                if identity>i:
                                    #并且目前序列的identity大于已经加入node序列的i，那么更新序列
                                    # 需要更新node中的字典
                                    node["children"][index].update(child)
                                    break #更新了之后，直接结束整个for循环
                            if index==len(node["children"])-1:
                                flag=1#遍历完了整个node["children"]也没有交错
                        if flag==1:#没有交错直接加入
                            node["children"].append(child)
            result.append(node)

        #是否存在文件夹myowndata/first_filter/Anas_platyrhynchos
        path = output_path + '/' + species
        create_dir_not_exist(path)
        json.dump(result, codecs.open(path+'/'+interferon, 'w+', encoding='utf-8'), ensure_ascii=False)



def get_species():#获取"data/genes"下的所有json文件
    datanames = os.listdir(input_path)
    name=[]
    for i in datanames:
        name.append(i)
    return name
def all2():
    names=get_species()
    print(names)
    for name in names:
        first(name)
if __name__ == '__main__':
    all2()
    #first(species)




