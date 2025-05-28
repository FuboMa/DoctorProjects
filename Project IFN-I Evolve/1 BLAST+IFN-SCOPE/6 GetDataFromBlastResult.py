import codecs
import json
import os
'''
将每一个.TXT文件变为json
修改两个路径
out_txt_path #比对结果文件存放位置（里面是各个物种的比对结果）
json_out_path#存放提取的json数据
'''
out_txt_path=r'C:\Users\Administrator.DESKTOP-N60LKMN\Desktop\blast_out'#比对结果文件存放位置（里面是各个物种的比对结果）
json_out_path= r'/myowndata/blast_json'
genome_location=[]#存放比对结果中的序列位于哪些染色体上
children_location=[]#存放每段相似序列开始的位置
#是否有多个转录本
def double_seq(file_data):
    #若只有一个转录本，其比对结果中只会出现一个GeneID
    index_list=[]
    for index, data in enumerate(file_data):
        if 'GeneID' in data:
            index_list.append(index)
    if len(index_list)==1:
        #只有一个转录本
        return 0
    else:
        # 至少存在两条转录本，返回第二条转录本在结果文件中的位置
        return index_list[1]


#去除多条转录本的数据,将比对结果中的第二条转录本以及之后的所有数据全部删去
def error_process2(file_data):
    index=double_seq(file_data)
    if index==0:#只有一个转录本数据不需要截断
        return file_data
    else:
        new_data=[]
        for i in range(0,index-1):
                new_data.append(file_data[i])
        return new_data
def mkdir(path):
    folder = os.path.exists(path)
    if not folder:  # 判断是否存在文件夹如果不存在则创建为文件夹
        print("创建了文件夹",path)
        os.makedirs(path)  # makedirs 创建文件时如果路径不存在会创建这个路径
#去除txt文档中的空格，将数据读入file_data
#name指示比对数据库，gene指示比对的基因
def get_useful_data(name,gene):
    file_data=[]
    fullname=out_txt_path+'/'+name+'/'+gene+'.txt'
    with open(fullname, "r", encoding="utf-8") as f:
        for line in f.readlines():
            if(line!='\n'):
                line=line.replace("\n","")
                file_data.append(line)
    return file_data
def get_children_num(file_data):
    children_location.clear()
    for i in range(len(file_data)):
        if 'Expect' in file_data[i] :
            children_location.append(i)
#genome_location=[]#存放比对结果中的序列位于哪些染色体上
#其中最后一个数据用于指示Lambda      K        H        a         alpha位于哪一行
def get_genome_num(file_data):
    genome_location.clear()
    for i in range(len(file_data)):
        if '>' in file_data[i] :
            genome_location.append(i)
        if 'Lambda' in file_data[i]:
            genome_location.append(i)
            return
#有些文档的染色体开始时有四行数据，有些有三行
def error_process(file_data):
    del_list=[]
    for i in range(len(genome_location)-1):#最后一个是元素代表的位置是Lambda      K        H        a         alpha，放入
        index=genome_location[i]
        if 'Length' not in file_data[index + 2]:
            del_list.append(index + 2)
    file_data = [file_data[i] for i in range(0, len(file_data), 1) if i not in del_list]
    if len(del_list)!=0:#file_data改变了，重新计算
        get_genome_num(file_data)
        get_children_num(file_data)

    return file_data

def make_one_node(num1,num2,file_data):
    json_node={}#里面会有三个键：genome（str），length（int），children（list）
    data1=file_data[num1]
    data2=file_data[num1+2]
    # print(data2)
    chr=""
    if  "NC_" in data1:
    #     p=data1
    #     p=p.replace(",","").split(" ")
    #     # print(p)
    #     chr_index=p.index('chromosome')
    #     chr=p[chr_index+1]
        p=file_data[num1]+file_data[num1+1]+file_data[num1+2]
        p = p.replace(",", "").split(" ")
        chr_index = p.index('chromosome')
        chr = p[chr_index + 1]
    json_node["chr"] = chr
    #data1:>NW_008659838.1 Acanthisitta chloris isolate BGI_N310 unplaced genomic scaffold,
    genome=data1.split(' ')[0].replace(">","")
    #data2:Length=84547829
    length=int(data2.split('=')[1])
    json_node["genome"]=genome
    json_node["length"] = length
    #children_location = []  # 存放每段相似序列开始的位置
    children = []
    for index in range(len(children_location)):#最后一个特殊处理
        c_num=children_location[index]
        children_node = {}
        if c_num>num1 and c_num<num2:
            c_dd=file_data[c_num+1]
            c_dd=c_dd.split(',')
            #c_dd:[' Identities = 50/156 (32%)', ' Positives = 91/156 (58%)', ' Gaps = 9/156 (6%)']
            identity=c_dd[0].split('(')[1].replace(')','')#32%
            positive=c_dd[1].split('(')[1].replace(')','')#58%
            gap=c_dd[2].split('(')[1].replace(')','')#6%
            seq_length=int(c_dd[0].split('/')[1].split(' ')[0])#156
            children_node["identity"]=identity
            children_node["positive"] = positive
            children_node["gap"] = gap
            children_node["seq_length"] = seq_length
            query_seq_content=''#用于比对的序列
            sbjct_seq_content=''#比对到的序列
            flag_seq_content = ''  # 比对到的序列
            start_index = children_location[index] + 3

            if index == len(children_location) - 1:
                end_index = num2 - 1
            if index!=len(children_location)-1:
                c_num2=children_location[index+1]-3
                if c_num2 not in genome_location:#代表这不是是新染色体第一个gene
                    end_index=children_location[index+1]-1
                else:#代表这是是新染色体第一个gene
                    end_index = children_location[index + 1] - 4

            # if index == len(children_location) - 1:
            #     end_index = num2 - 1
            # if index!=len(children_location)-1:
            #     # start_index=children_location[index]+3
            #     end_index=children_location[index+1]-1
            seq_start=int(file_data[start_index+2].split('  ')[1])
            # print(seq_start)
            # print(num2)
            # print("c_num")
            # print(c_num)
            # print("end_index")
            # print(end_index)
            # print(file_data[end_index])
            # print(file_data[end_index-1])
            # seq_end=int(file_data[end_index].split('  ')[3])
            seq_end = int(file_data[end_index].split()[3])
            if seq_start>seq_end:
                temp=seq_start
                seq_start=seq_end
                seq_end=temp
            children_node["seq_start"] = seq_start
            children_node["seq_end"] = seq_end
            query_index=start_index
            flag_index=start_index+1
            sbjct_index=start_index+2
            while(query_index<=end_index and sbjct_index<=end_index):
                #query_s:QGMREDF---SLLNQMSTFSLVPCLKDRTNFSFPKEAMQGTQLQRENTTVMVHEMLQQIF  84
                query_s=file_data[query_index][17:]
                query_s=query_s.split('  ')[0]
                sbjct_s=file_data[sbjct_index][17:]
                sbjct_s=sbjct_s.split('  ')[0]
                flag_s=file_data[flag_index][17:]
                # query_s:QGMREDF---SLLNQMSTFSLVPCLKDRTNFSFPKEAMQGTQLQRENTTVMVHEMLQQIF  84
                # query_s=file_data[query_index].replace(" ","").replace("Query","")
                # query_s = ''.join([i for i in query_s if not i.isdigit()])
                # sbjct_s = file_data[sbjct_index].replace(" ","").replace("Sbjct","")
                # sbjct_s = ''.join([i for i in sbjct_s if not i.isdigit()])
                while(len(flag_s)<len(query_s)):
                    flag_s=flag_s+' '
                query_seq_content = query_seq_content+query_s
                sbjct_seq_content = sbjct_seq_content+sbjct_s
                flag_seq_content=flag_seq_content+flag_s
                query_index=query_index+3
                sbjct_index=sbjct_index+3
                flag_index=flag_index+3
            children_node["query_seq_content"] = query_seq_content
            children_node["sbjct_seq_content"] = sbjct_seq_content
            children_node["flag_seq_content"] = flag_seq_content
        children.append(children_node)
    json_node["children"] = children
    return json_node
    # print(json_node)

def make_all(file_data):
    all_json=[]
    # print(genome_location)
    # print(children_location)
    for i in range(len(genome_location)-1):
        genome_node=make_one_node(genome_location[i],genome_location[i+1],file_data)
        all_json.append(genome_node)
    return all_json
def run_one(db,gene):
    o_data=get_useful_data(db, gene)
    file_data = error_process2(o_data)#去除多条转录本的数据
    get_genome_num(file_data)
    get_children_num(file_data)
    file_data=error_process(file_data)
    all_json=make_all(file_data)
    mkdir(json_out_path + '/' + db)#若不存在db文件夹就创建
    json.dump(all_json, codecs.open(json_out_path + '/' + db+'/'+gene+'.json',
                                             'w+', encoding='utf-8'), ensure_ascii=False)

def all():
    db_names = os.listdir(out_txt_path)
    for db_name in db_names:
        path=out_txt_path+'/'+db_name
        genes=os.listdir(path)
        for gene in genes:
            gene=gene.replace(".txt","")
            run_one(db_name,gene)

if __name__ == '__main__':
    # db='Aythya_fuligula'
    # gene='Acanthisitta_chloris_LOC103797348'
    # print('1')
    # run_one(db,gene)
    all()








