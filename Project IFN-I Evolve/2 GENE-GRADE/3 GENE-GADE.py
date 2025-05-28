
import codecs
import os
import re

import pymongo
import networkx as nx
import matplotlib.pyplot as plt
import json
import pickle

from matplotlib.lines import Line2D

list_name = 'myowndata/list_node.json'
datas_path='myowndata/datas_node.json'
joined_list_path='myowndata/joined_list.json'
result_path='myowndata/improved_KDC.json'
request_str = 'mongodb://root:123456@localhost:27017/?authMechanism=DEFAULT'
client = pymongo.MongoClient(request_str)
AlloldAndnewIFN = client.interferon.AlloldAndnewIFN

maxLen = 10
pseudo_flag=0#1表示包含假基因
lowqualit_flag=1#1表示包含lowquality
# allGenes = client.interferon.New2
allGenes = client.interferon.AllGeneDB2

NewIFNsPath=r"myowndata/newIFN"
def getNewIFN():
    IFNs=[]
    fileName = os.listdir(NewIFNsPath)
    speList=[]
    chrList=[]
    speAndchr=[]
    # print("当前文件夹下有:",len(fileName),"个文件")
    for node in fileName:
        if ".DS" not in node:
        # print(node)
            species=node.split("_N")[0].replace("_"," ")
            pos=node.split("POS")[1].replace(".fasta","")
            chromosomeNum=node.split("_N")[1].split("_P")[0].split("_")[1]
            flag=0
            if 'NC' in node:
                flag=1
                chromosome='NC_'+chromosomeNum
            else:
                chromosome = 'NW_' + chromosomeNum
            # print(chromosome)
            # print(species)
            # print(pos)
            if species not in speList:
                speList.append(species)
            if chromosome not in chrList:
                chrList.append(chromosome)
            allName=species+"_"+chromosome
            if allName not in speAndchr:
                speAndchr.append(allName)
            if flag==1:#只加入NC
                IFNs.append({
                    "species":species,
                    "chromosome":chromosome,
                    "pos":int(pos),
                })
    return IFNs,speList,chrList,speAndchr
def findPos(species,chromosome):
    Pos=[]
    IFNs,_,_ ,_= getNewIFN()
    for IFN in IFNs:
        if IFN["species"]==species and IFN["chromosome"]==chromosome:
            Pos.append(IFN["pos"])
    return Pos
def get_json(tpath):
    genes = json.load(open(tpath, "r"))
    return genes
def joinedIFN1():
    # excludeSpecies = ['Xenopus laevis', 'Rana temporaria', 'Bufo gargarizans', 'Bufo bufo', 'Protopterus annectens',
    #                   'Latimeria chalumnae', 'Danio rerio', 'Oryzias latipes', 'Oncorhynchus mykiss', 'Carcharodon carcharias', 'Petromyzon marinus']
    alldatas=[]
    datas = []
    cnt = 0
    allIFN_index = []#干扰素节点在AfterjoinedIFN.json中所在的位置
    IFN_num=0
    for x in allGenes.find({}):
        cnt += 1
        if cnt % 100000 == 0:
            print(cnt)
        if lowqualit_flag==0:
            if 'lowquality' in x['flag']:
                continue
        if pseudo_flag==0:
            if 'pseudo' in x['flag']:
                continue
        if x['type'] ==2 or x['type'] ==3:
            type = 2
        else:
            type=x['type']
        if x['type'] == 0 and 'lowquality' in x['flag']:
            continue
        if x['type'] == 1 or x['chromosome'][:2] == 'NW' :
        # if x['type'] == 1:
            continue
        if x['type'] == 0:
            geneName = '*'.join([x['species'], x['gene'],
                                 str(x['start'])])
        # if x['chromosome'] in chrList:
        pos=int((x['start']+x['end'])/2)
        alldatas.append({'gene': x['gene'], 'species': x['species'], 'chromosome': x['chromosome'],
                      'start': x['start'], 'end': x['end'], 'type': type,'flag': x['flag'],'pos':pos,})
    json.dump(alldatas, codecs.open(
        'myowndata/tempAllData.json', 'w+', encoding="utf-8"), ensure_ascii=False)
    # print(len(alldatas))
    # return alldatas
    #合并相邻干扰素节点 每个干扰素节点多两个属性ifn 和 num
    # i = 0
    # print("开始融合相邻的干扰素节点...")
    # while i<len(alldatas):
    #     #if alldatas[i]['type'] == 2 and alldatas[i]['chromosome'][:2] != 'NW':
    #     if alldatas[i]['type'] == 2:
    #        datas.append(alldatas[i])
    #     #if alldatas[i]['type'] == 0 and alldatas[i]['chromosome'][:2] != 'NW':
    #     if alldatas[i]['type'] == 0:
    #         IFN_num=IFN_num+1
    #         last_index=len(datas)-1
    #         geneName = '*'.join([alldatas[i]['species'], alldatas[i]['gene'],
    #                              str(alldatas[i]['start'])])
    #         if datas[last_index]["type"]==0 and datas[last_index]["chromosome"]==alldatas[i]["chromosome"]\
    #                 and datas[last_index]["species"]==alldatas[i]["species"]:
    #             datas[last_index]["num"]=datas[last_index]["num"]+1
    #             datas[last_index]["ifn"].append(geneName)
    #         else:
    #             alldatas[i]["num"]=1
    #             alldatas[i]["ifn"]=[]
    #             alldatas[i]["ifn"].append(geneName)
    #             datas.append(alldatas[i])
    #             allIFN_index.append(len(datas)-1)
    #     i += 1
    #     # if i % 100000 == 0:
    #     #     print(i)
    # json.dump(datas, codecs.open(
    #     'myowndata/AfterjoinedIFN.json', 'w+', encoding="utf-8"), ensure_ascii=False)
    # json.dump(allIFN_index, codecs.open(
    #     'myowndata/IFNNode_Location.json', 'w+', encoding="utf-8"), ensure_ascii=False)
    # print("融合相邻的干扰素节点完成!!!")
    # print("干扰素数量:",IFN_num)
    # print("融合后的干扰素节点数量:",len(allIFN_index))
def insertDatas():
    IFNs, _, chrList, _ = getNewIFN()
    allDatas = get_json("myowndata/tempAllData.json")
    # print(len(IFNs))
    odllenalldatas=len(allDatas)
    print("插入前：",odllenalldatas)
    flag=0
    for indexIFN,item in enumerate(IFNs):
        print(indexIFN+1,"/",len(IFNs))
        spe=item["species"]
        chr=item["chromosome"]
        # print(chr)
        pos = item["pos"]
        # print(pos)
        # break
        for index,data in enumerate(allDatas):
            if chr==data["chromosome"]:
                # print("ok")
                # break
                if pos>=data["pos"]:
                    if allDatas[index+1]["chromosome"]!=chr:
                        allDatas.insert(index + 1, {'gene': spe + "_" + chr + "_" + str(pos), 'species': spe,
                                                    'chromosome': chr,
                                                    'start': pos, 'end': pos, 'type': 0, 'flag': 'normal',
                                                    'pos': pos, })
                        flag = 1
                        break
                    else:
                        if pos<allDatas[index+1]["pos"]:
                            allDatas.insert(index + 1, {'gene': spe + "_" + chr + "_" + str(pos), 'species': spe,
                                                        'chromosome': chr,
                                                        'start': pos, 'end': pos, 'type': 0, 'flag': 'normal',
                                                        'pos': pos, })
                            flag=1
                            break
        if flag==1:
            print("成功插入")
            # break
    print("插入后：", len(allDatas))
    print("成功插入：",len(allDatas)-odllenalldatas,"个干扰素")
    json.dump(allDatas, codecs.open(
        'myowndata/tempAllData_newIFN_inserted.json', 'w+', encoding="utf-8"), ensure_ascii=False)
def joinedIFN2():
    alldatas=get_json("myowndata/tempAllData_newIFN_inserted.json")
    i = 0
    datas = []
    allIFN_index = []  # 干扰素节点在AfterjoinedIFN.json中所在的位置
    IFN_num = 0
    print("开始融合相邻的干扰素节点...")
    while i<len(alldatas):
        #if alldatas[i]['type'] == 2 and alldatas[i]['chromosome'][:2] != 'NW':
        if alldatas[i]['type'] == 2:
           datas.append(alldatas[i])
        #if alldatas[i]['type'] == 0 and alldatas[i]['chromosome'][:2] != 'NW':
        if alldatas[i]['type'] == 0:
            IFN_num=IFN_num+1
            last_index=len(datas)-1
            geneName = '*'.join([alldatas[i]['species'], alldatas[i]['gene'],
                                 str(alldatas[i]['start'])])
            if datas[last_index]["type"]==0 and datas[last_index]["chromosome"]==alldatas[i]["chromosome"]\
                    and datas[last_index]["species"]==alldatas[i]["species"]:
                datas[last_index]["num"]=datas[last_index]["num"]+1
                datas[last_index]["ifn"].append(geneName)
            else:
                alldatas[i]["num"]=1
                alldatas[i]["ifn"]=[]
                alldatas[i]["ifn"].append(geneName)
                datas.append(alldatas[i])
                allIFN_index.append(len(datas)-1)
        i += 1
        # if i % 100000 == 0:
        #     print(i)
    json.dump(datas, codecs.open(
        'myowndata/AfterjoinedIFN.json', 'w+', encoding="utf-8"), ensure_ascii=False)
    json.dump(allIFN_index, codecs.open(
        'myowndata/IFNNode_Location.json', 'w+', encoding="utf-8"), ensure_ascii=False)
    print("融合相邻的干扰素节点完成!!!")
    print("干扰素数量:",IFN_num)
    print("融合后的干扰素节点数量:",len(allIFN_index))
def getShortList():
    datas = get_json("myowndata/AfterjoinedIFN.json")
    IFN_Location=get_json("myowndata/IFNNode_Location.json")
    #while i < len(datas):
    shortList=[]
    # while f <len(datas):
    for i in IFN_Location:
        fullNode={}
        #if datas[i]['type'] == 2 and datas[i]['chromosome'][:2] != 'NW':
        #if datas[i]['type'] == 2:
        j = i
        jCnt = 0
        l_ifn_flag=[]
        l_index_list=[]
        l_gene_list=[]
        while (j > 0 and jCnt < maxLen and datas[j - 1]['chromosome'] == datas[i]['chromosome']):
            j -= 1
            if datas[j]['type'] == 2:
                jCnt += 1
                lNode={"gene":datas[j]['gene']}
                l_gene_list.insert(0,lNode)
            # if datas[j]['type'] == 0:
            #     # jCnt += 1
            #     l_index_list.insert(0,j)
            #     l_ifn_flag.insert(0,j)
        k = i
        kCnt = 0
        r_ifn_flag = []
        r_index_list = []
        r_gene_list = []
        while (k < len(datas) - 1 and kCnt < maxLen and datas[k + 1]['chromosome'] == datas[i]['chromosome']):
            k += 1
            if datas[k]['type'] == 2:
                kCnt += 1
                rNode = {"gene": datas[k]['gene']}
                r_gene_list.append(rNode)
        fullNode["IFN"]=datas[i]
        fullNode["lGene"] =l_gene_list
        fullNode["rGene"] = r_gene_list
        tempIFN=[{"gene":datas[i]["gene"]}]
        wholeList=l_gene_list+tempIFN+r_gene_list#整个短链11个基因节点
        fullNode["shortList"]=wholeList
        shortList.append(fullNode)
        # i += 1
        # if i % 100000 == 0:
        #     print(i)

    json.dump(shortList, codecs.open(
        'myowndata/shortList.json', 'w+', encoding="utf-8"),indent=2, ensure_ascii=False)
    print("短链组制作完成!!!")
    print("短链组数量:", len(shortList))
def merge_same_name_nodes(graph):
    same_name_nodes = {}  # 存储同名节点
    for node in graph.nodes():
        name = graph.nodes[node]['name']
        if name in same_name_nodes:
            same_name_nodes[name].append(node)
        else:
            same_name_nodes[name] = [node]

    for name, nodes in same_name_nodes.items():
        if len(nodes) > 1:
            keep_node = nodes[0]
            for node in nodes[1:]:
                # 合并属性
                graph.nodes[keep_node].update(graph.nodes[node])
                # 更新边
                for neighbor in graph.neighbors(node):
                    graph.add_edge(keep_node, neighbor)
                graph.remove_node(node)

    return graph

# def joinedShortListold():#->获得vlook2.json
#     G = nx.Graph()
#     nodes = []
#     edges=[]
#     oldShortList = get_json("myowndata/shortList.json")
#     #剔除含有excludes中基因的链表
#     newShortList = []
#     # for ojlist in oldShortList:
#     #     #joined_list.append(ojlist)
#     #     newShortList.append(ojlist)
#     #     for item in ojlist["shortList"]:
#     #         if item['gene'] in excludes:
#     #             newShortList.pop()
#     #             break
#     for ojlist in oldShortList:
#         flag=1
#         for item in ojlist["shortList"]:
#             if item['gene'] in excludes:
#                 flag=0
#         if flag==1:
#             newShortList.append(ojlist)
#     #print(len(newShortList))
#     num=0
#     for item in newShortList:
#         for index, gene in enumerate(item["shortList"]):
#             num=num+1
#             G.add_node(num, name=gene['gene'])
#             if index != 0:  # 当前节点gene["gene"]和前一个节点连接
#                 G.add_edge(num-1,num)
#     merged_graph = merge_same_name_nodes(G)
#     connected_subgraphs = list(nx.connected_components(merged_graph))
#     # 输出每个连通子图中度最大的节点和其 name 属性
#     num=0
#     for subgraph_nodes in connected_subgraphs:
#         subgraph = G.subgraph(subgraph_nodes)
#         max_degree_node = max(subgraph.degree, key=lambda x: x[1])[0]
#         max_degree_node_name = subgraph.nodes[max_degree_node]['name']
#         if subgraph.degree[max_degree_node]>10:
#             num=num+1
#             print("度最大------------")
#             print(f"  Node {max_degree_node}: Name = {max_degree_node_name}, Degree = {subgraph.degree[max_degree_node]}")
#             max_degree_node = max(subgraph.degree, key=lambda x: x[1])[0]
#             node_colors = ['red' if node == max_degree_node else 'green' for node in subgraph.nodes()]
#             edge_colors = ['red' if edge in subgraph.edges(max_degree_node) else '#C0C0C0' for edge in subgraph.edges()]
#             # edge_to_highlight =[if edge in subgraph.edges(max_degree_node) for edge in subgraph.edges()]
#             node_labels = {node: subgraph.nodes[node]['name'] for node in subgraph.nodes()}
#             k=1.3
#             pos = nx.spring_layout(subgraph,seed=42,k=k)  # 设置节点位置
#             # pos = nx.random_layout(subgraph)
#             # pos = nx.kamada_kawai_layout(subgraph)
#             node_size = 600  # 调整节点的大小
#             fig_size = (24, 24)  # 调整图的大小
#             plt.figure(figsize=fig_size)
#             nx.draw( subgraph, pos, with_labels=False, node_color=node_colors, node_size=node_size,edge_color=edge_colors)
#             nx.draw_networkx_labels(subgraph, pos, labels=node_labels)  # 添加节点名称
#             # nx.draw_networkx_edges(subgraph, pos, edgelist=edge_to_highlight = (3, 4), edge_color='red', width=2)
#             plt.savefig("myowndata/article/graph_"+str(num)+".pdf", format="pdf")

excludes = []#,"HACD4","UBAP2","MOB3B"#影响joinedShortList()执行结果
def getLocijsonData(lociGene,graph,interferon0):
    subgraph_nodes_list = list(graph.nodes())
    subgraph_edges_list = list(graph.edges())
    newNodes = []
    newEdges = []
    # loci3=["UBAP2","HACD4","MOB3B"]
    for item in subgraph_nodes_list:
        name = graph.nodes[item]['name']  # 节点名
        # print(name)
        # break
        if name == lociGene:
            newNodes.append({'name': name, 'category': 1})  # 1表示度最大的节点
        elif name in interferon0:
            newNodes.append({'name': name, 'category': 0})  # 0表示干扰素
        else:
            newNodes.append({'name': name, 'category': 2})  # 2表示其他基因
    for edge in subgraph_edges_list:
        name1=graph.nodes[edge[0]]['name']
        name2 = graph.nodes[edge[1]]['name']
        newEdges.append({
            'source': name1,
            'target': name2,
        })
    LASTResult = []
    LASTResult.append({'data': newNodes, 'links': newEdges})
    json.dump(LASTResult, codecs.open(
        'myowndata/' + lociGene + ".json", 'w+', encoding="utf-8"), ensure_ascii=False, indent=2)
    print("------",lociGene,".json数据获取成功")
def joinedShortList():#->获得vlook2.json
    interferonAllData = get_json("myowndata/interferons0-0821-noPsedo.json")
    interferon0 = []
    for Data in interferonAllData:
        for ifn in Data["interferon"]:
            interferon0.append(ifn["gene"])
    print("IFN0数量:", len(interferon0))
    G = nx.Graph()
    nodes = []
    edges=[]
    oldShortList = get_json("myowndata/shortList.json")
    #剔除含有excludes中基因的链表
    newShortList = []
    # for ojlist in oldShortList:
    #     #joined_list.append(ojlist)
    #     newShortList.append(ojlist)
    #     for item in ojlist["shortList"]:
    #         if item['gene'] in excludes:
    #             newShortList.pop()
    #             break
    for ojlist in oldShortList:
        flag=1
        for item in ojlist["shortList"]:
            if item['gene'] in excludes:
                flag=0
        if flag==1:
            newShortList.append(ojlist)
    #print(len(newShortList))
    num=0
    for item in newShortList:
        # newNodes = []
        # newEdges = []
        IFNName=item["IFN"]["gene"]
        for index, gene in enumerate(item["shortList"]):
            num=num+1
            if gene['gene']==IFNName:
                G.add_node(num, name=gene['gene'].upper(),type="interferon gene")
                # newNodes.append({'name': gene['gene'], 'category': 0}) #干扰素
            else:
                G.add_node(num, name=gene['gene'].upper(),type="normal gene")
                # newNodes.append({'name': gene['gene'], 'category': 1})#普通基因
            if index != 0:  # 当前节点gene["gene"]和前一个节点连接
                G.add_edge(num-1,num)
                # newEdges.append({
                #     'source': gene["gene"],
                #     'target': item["shortList"][index - 1]["gene"],
                # })
    merged_graph = merge_same_name_nodes(G)
    connected_subgraphs = list(nx.connected_components(merged_graph))
    num=0
    showLabel=0#是否显示标签，1表示显示，0表示不显示
    color_IFN = '#3498DB'#干扰素节点的颜色
    color_gene = "#27AE60"#普通基因的颜色
    color_maxDegree = "#FA8072"#度最大节点的颜色
    color_label = "black"#标签的颜色
    color_normal_edges='#C0C0C0'#普通边的颜色
    color_related_MaxDegreeNode_edges="#F08080"#与最大度节点相连的边的颜色
    size_label=10#标签大小
    size_IFN=18#干扰素节点的大小
    size_gene=6#普通基因的大小
    size_maxDegree=150#度最大节点的大小
    edge_widths=0.1#边的宽度
    legentSize =15#图例大小
    k = 0.9#引力大小
    for subgraph_nodes in connected_subgraphs:
        subgraph = merged_graph.subgraph(subgraph_nodes)
        # max_degree_node = max(subgraph.degree, key=lambda x: x[1])[0]
        # max_degree_node_name = subgraph.nodes[max_degree_node]['name']
        degree_sequence = sorted(subgraph.degree, key=lambda x: x[1], reverse=True)
        max_degree_node = degree_sequence[0][0]#节点
        max_degree=degree_sequence[0][1]#度
        max_degree_node_name = subgraph.nodes[max_degree_node]['name']#节点名

        sec_degree_node = degree_sequence[1][0]  # 节点
        sec_degree = degree_sequence[1][1]  # 度
        sec_degree_node_name = subgraph.nodes[sec_degree_node]['name']  # 节点名
        if subgraph.degree[max_degree_node]>15:
            #获取网站需要的数据
            getLocijsonData(max_degree_node_name, subgraph, interferon0)
            print("度最大节点:", max_degree_node, "节点名:", max_degree_node_name, "度:", max_degree)
            # print("度第二大节点:", sec_degree_node, "节点名:", sec_degree_node_name, "度:", sec_degree)
            # print(degree_sequence)
            degreeDetailList=[]
            for item in degree_sequence:
                degreeDetailList.append({
                    # "degreeNumber":item[0],
                    "nodeName":subgraph.nodes[item[0]]['name'],
                    "degree":item[1],
                })
            # print(degreeDetailList)
            #----获取每个图的度列表
            json.dump(degreeDetailList, codecs.open(
                'myowndata/article/' + max_degree_node_name + ".json", 'w+', encoding="utf-8"), ensure_ascii=False, indent=2)

            num=num+1
            # print("度最大------------")
            # print(f"  Node {max_degree_node}: Name = {max_degree_node_name}, Degree = {subgraph.degree[max_degree_node]}")
            # max_degree_node = max(subgraph.degree, key=lambda x: x[1])[0]
            node_colors = [color_maxDegree if node == max_degree_node else color_IFN if subgraph.nodes[node]['type'] == "interferon gene" else color_gene for node in subgraph.nodes()]
            edge_colors = [color_related_MaxDegreeNode_edges if edge in subgraph.edges(max_degree_node) else color_normal_edges for edge in subgraph.edges()]
            node_labels = {max_degree_node: subgraph.nodes[max_degree_node]['name']}
            edge_to_highlight = []
            for edge in subgraph.edges():
                if edge in subgraph.edges(max_degree_node):
                    edge_to_highlight.append(edge)
            # node_sizes = [1200 if node == max_degree_node else 10 for node in subgraph.nodes()]
            node_sizes = [size_maxDegree if node == max_degree_node else size_IFN if subgraph.nodes[node]['type'] == "interferon gene" else size_gene for node in subgraph.nodes()]
            pos = nx.spring_layout(subgraph,seed=42,k=k)  #利用force-directed算法生成节点布局
            # pos = nx.random_layout(subgraph)# 生成随机节点布局
            # pos = nx.kamada_kawai_layout(subgraph)#利用kamada_kawai_layout算法布局
            # pos = nx.circular_layout(subgraph)  # 生成圆形节点布局
            # pos = nx.shell_layout(subgraph)  # 生成同心圆节点布局
            # pos = nx.spectral_layout(subgraph)  # 利用图拉普拉斯特征向量生成节点布局
            fig_size = (9, 9)  # 调整图的大小
            plt.figure(figsize=fig_size)
            nx.draw(subgraph, pos, with_labels=False, node_color=node_colors, node_size=node_sizes,
                    edge_color=edge_colors,width=edge_widths)
            # 调整标签大小
            if showLabel==1:
                nx.draw_networkx_labels(subgraph, pos, labels=node_labels, font_size=size_label,font_color=color_label)
            # 与最大度节点关联的边置于最上层
            nx.draw_networkx_edges(subgraph, pos, edgelist=edge_to_highlight, edge_color=color_related_MaxDegreeNode_edges, width=edge_widths)
            # 最大度节置于最上层
            nx.draw_networkx_nodes(subgraph, pos, nodelist=[max_degree_node], node_color=color_maxDegree, node_size=size_maxDegree)
            # plt.savefig("myowndata/article/graph_"+str(num)+".pdf", format="pdf")
            fff=0
            if max_degree_node_name=="HACD4" and fff==1:
                legend_elements = [
                    Line2D([0], [0], marker='o', color='w', markerfacecolor=color_IFN, markersize=legentSize,
                           label='Interferon gene'),
                    Line2D([0], [0], marker='o', color='w', markerfacecolor=color_gene, markersize=legentSize, label='Normal gene'),
                    Line2D([0], [0], marker='o', color='w', markerfacecolor=color_maxDegree, markersize=legentSize,
                           label='Max Degree gene')
                ]
                plt.legend(handles=legend_elements, loc='upper left',fontsize=legentSize)  # 添加图例到右上角
            # plt.savefig("myowndata/article/graph_" + str(num) + ".pdf", format="pdf")max_degree_node_name
            plt.savefig("myowndata/article/" + max_degree_node_name + ".pdf", format="pdf")
    # 在函数结尾的位置添加以下代码来生成图例
    # legend_labels = {
    #     "interferon gene": "Interferon Gene",
    #     "normal gene": "Normal Gene",
    #     "max_degree": "Max Degree Gene"
    # }
    #
    # legend_colors = {
    #     "interferon gene": color_IFN,
    #     "normal gene": color_gene,
    #     "max_degree": color_maxDegree
    # }
    #
    # plt.figure(figsize=(8, 6))  # 调整图例大小
    # handles = [plt.Line2D([0], [0], marker='o', color='w', label=legend_labels[label],
    #                       markersize=10, markerfacecolor=legend_colors[label]) for label in legend_labels]
    #
    # plt.legend(handles=handles, title="Node Types", loc="upper right")
    # plt.axis('off')  # 不显示坐标轴
    # plt.tight_layout()
    # plt.savefig("myowndata/article/legend.pdf", format="pdf")  # 保存图例


def filterIFN():
    loci3 = ["UBAP2", "HACD4", "MOB3B"]
    for lociGene in loci3:
        LASTResult=[]
        data=get_json( 'myowndata/article/' + lociGene + '.json')
        for node in data:
            gene=node["nodeName"].replace(" ","_")
            if "_" in gene:
                numPOS = gene.split("_")[-1]
                gene = gene.replace(numPOS, "POS" + numPOS)
            # print(gene)
            regex = re.compile(gene, re.IGNORECASE)
            sequencesAll = AlloldAndnewIFN.find({"gene": regex})
            temp = []
            for item in sequencesAll:
                temp.append(item)
            if len(temp) == 0:
                LASTResult.append(node)
        json.dump(LASTResult, codecs.open(
            'myowndata/article/new' + lociGene + ".json", 'w+', encoding="utf-8"), ensure_ascii=False, indent=2)
def getStat():
    locigenes = ["HACD4","UBAP2","MOB3B"]#"UBE2R2"
    oldShortList = get_json("myowndata/shortList.json")
    aroundIFNNum = 0
    allIFNnum=0
    allIFNname=[]#所有干扰素名字
    aroundLociIFN=[]#位点周围干扰素
    outLociIFN = []#位点外干扰素
    inLociIFN=[]##位点内干扰素
    for item in oldShortList:#获取所有干扰素名字和数量
        allIFNnum=allIFNnum+item["IFN"]["num"]
        allIFNname=allIFNname+item["IFN"]["ifn"]
    for loci in locigenes:
        for list in oldShortList:
            tempNode={"gene":loci}
            if tempNode in list["shortList"]:
                aroundIFNNum=aroundIFNNum+list["IFN"]["num"]
                aroundLociIFN=aroundLociIFN+list["IFN"]["ifn"]
    for ifn in allIFNname:
        if ifn not in aroundLociIFN:
            outLociIFN.append(ifn)
        else:
            inLociIFN.append(ifn)
    f = open("myowndata/三个位点外的干扰素.txt", "w", encoding='utf-8')
    for i, node in enumerate(outLociIFN):
        strs = node.split("*")[0] + "*" + node.split("*")[1]
        f.write(strs)
        if i != len(outLociIFN) - 1:
            f.write('\n')
    f.close()
    f = open("myowndata/三个位点内的干扰素.txt", "w", encoding='utf-8')
    for i, node in enumerate(inLociIFN):
        strs = node.split("*")[0] + "*" + node.split("*")[1]
        f.write(strs)
        if i != len(inLociIFN) - 1:
            f.write('\n')
    f.close()
    print("总的干扰素数量:",allIFNnum)
    print("位于位点周围的干扰素数量:", aroundIFNNum)
    print("位于外干扰素数量:", allIFNnum-aroundIFNNum)
    print("位于位点周围的干扰素数量所占比例:", aroundIFNNum/allIFNnum)
    print("位点外干扰素数量所占比例:", 1-aroundIFNNum/allIFNnum)
    print("请在myowndata/三个位点外的干扰素.txt中查看位点外干扰素")
    print("请在myowndata/三个位点内的干扰素.txt中查看位点内干扰素")

if __name__ == '__main__':
    # joinedIFN1()#从数据库获取所有基因数据，这个不用执行
    # 修改NewIFNsPath=r"myowndata/NewIFNs"下的内容即可，将新干扰素放入这个文件中，执行下面的函数
    # insertDatas()#生成插入后的新文件
    # joinedIFN2()#融合相邻干扰素节点
    # getShortList()
    joinedShortList()  # 根据getShortList()的结果输出当前度最大的节点
    getStat()  # 根据getShortList()的结果获得统计数据
    filterIFN()
















