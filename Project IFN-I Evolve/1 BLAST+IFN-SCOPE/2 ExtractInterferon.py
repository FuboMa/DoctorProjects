import codecs
import time
import json
import re
#from pymongo import MongoClient


def extractInterferon(species=None):
    def get_json(spe):
        genes_path = 'data/genes'
        with open(genes_path + '/' + spe + '.json') as f:
            data = json.load(f)
        return data

    def updategenes(spe):
        datas = get_json(spe)
        result = []
        for data in datas:
            #针对小家鼠中IFi206系列基因
            #'interferon-activable'
            # if data["type"] == 0:
            #     if "interferon activated".lower() in data["product"].lower():
            #         data["type"]=1
            #     if "interferon-activable" in data["product"]:#针对Orcinus orca虎鲸的LOC101278136
            #         data["type"]=1
            # #针对某些新注释物种中的基因，这些基因名字是LOC，但是product是MOB3B
            # if data["type"] == 2:
            #     if "545645" == data["geneID"]:#针对小家鼠的Gm13283
            #         data["type"]=0
            #     if 'MOB Kinase Activator 3B'.lower() in data["product"].lower():
            #         data["gene"] = 'MOB3B'
            #         data["type"] = 3
            #         # print(data["product"])
            #     if 'Focadhesin'.lower() in data["product"].lower():
            #         data["gene"] = 'FOCAD'
            #         data["type"] = 3
            #     if '3-Hydroxyacyl-CoA Dehydratase 4'.lower() in data["product"].lower():
            #         data["gene"] = 'HACD4'
            #         data["type"] = 3
            #     if 'Ubiquitin Associated Protein 2'.lower() in data["product"].lower():
            #         data["gene"] = 'UBAP2'
            #         data["type"] = 3
            #     if 'Ubiquitin-Associated Protein 2'.lower() in data["product"].lower():
            #         data["gene"] = 'UBAP2'
            #         data["type"] = 3
            #     if 'biquitin Conjugating Enzyme E2 R2'.lower() in data["product"].lower():
            #         data["gene"] = 'UBE2R2'
            #         data["type"] = 3

            result.append(data)
        json.dump(  # 存储该物种的所有基因(有没有干扰素都存了)
            result,
            codecs.open("data/genes/" + spe + ".json", "w", encoding="utf-8"),
            ensure_ascii=False,
        )

    def addGene(gene):
        nonlocal exonPos
        #gene是否在要找的gene组(tital)
        chk3 = gene in tital
        temp = 0           #   chk为1代表存在干扰素，chk2为1时，表示I型干扰素
        flag = 'normal'
        if chk and chk2:     #存在干扰素，不存在regulator、gamma、receptor...
            temp = 0
        elif chk:            #存在干扰素，存在regulator、gamma、receptor...中的一种或几种
            exonPos = []
            temp = 1
            IFNRelated.append({
                "gene":gene,
                "type":temp,
                "cateDcp":cateDcp,
                "info":tempResult,
            })
        elif chk3:           #不存在干扰素，但此gene在要找的gene组中
            exonPos = []
            temp = 3
        else:               #不存在干扰素，此gene也不在要找的gene组
            exonPos = []
            temp = 2
        if chk_pseudo:
            flag='pseudo'
        if chk_lowquality:
            flag='lowquality'
        if chk_pseudo and chk_lowquality:
            flag='pseudo and lowquality'
        exonPos = sorted(exonPos)#除了第一种情况之外，其他情况下的exonPos都为空
        genes.append(
            {
                "gene": gene,
                "geneID": geneID,
                "type": temp,
                "chromosome": chromosome,
                "chr":chromo,
                "start": int(start),
                "end": int(end),
                "exon": int(exon),
                "exonPos": exonPos,
                "product": product,
                "direction": direction,
                "proteinID": proteinID,
                "flag": flag,
            }
        )
    #菜单可以筛选不同“长度”与位点的基因集合。长度指从特定基因（HACD4等菜单中所列基因以及干扰素相关基因）向两边扩展的基因数量。
    def add(gene, length):#lenth用于控制此基因周围的基因数目
        nonlocal result#声明全局变量
        nonlocal cnt    #声明全局变量                      1、 #存在干扰素，不存在regulator、gamma、receptor...，即此基因存在I型干扰素
        if gene["type"] == 0 or gene["type"] == 3:#2、或者不存在干扰素，但此gene在要找的gene组中
            result.append(gene)
            print(gene["gene"])
            cnt = 0                                #满足条件1/2直接将此基因加入，并且cnt设置为0
        elif len(result) < length:
            result.append(gene)
        elif len(result) == length:
            result.pop(0)  #删除第一个元素
            result.append(gene)
        elif cnt < length:     #cnt表示0/2类型的基因周围的基因数目
            result.append(gene)
            cnt += 1
        else:                    #当周围基因数量达到要求后，结果存入results中，并将cnt初始化
            results.append(result)
            result = [gene]
            cnt = 0
        return result

    startTime = time.time()
    #要提取的基因
    tital = [
        "UBE2R2",
        "UBAP2",
        "HACD4",
        "FOCAD",
        "MOB3B",
        "MTAP",
        "TBC1D2",
        "gh1",
        "cd79b",
        "plekhm1",
        "mapta",
        #  'LINGO2', 'NOL6',
    ]
    if species == None:
        species = json.load(open("data/sortedSpecies.json", encoding="utf-8"))
    for x in species[:]:#x=    ["Anas platyrhynchos","绿头鸭"]
        print(x)
        genes = []#存储所有gene的信息
        IFNRelated = []#存储所有type=1类型基因的相关判定信息
        chromosomes = {}#物种的所有染色体信息
        chk = False
        chk2 = True
        chk3 = False
        chk_pseudo = False
        chk_lowquality = False
        previousGene = ""#存储前一个基因的名字
        exon = 0#当前基因的外显子数目
        exonPos = []#当前基因外显子的所处位置[(起点，终点)]
        product = ""#当前基因的产物信息
        chromosome = ""#当前基因所在的染色体
        start = -1
        end = -1
        direction = ""#当前基因位于参考序列的正链(+)或负链(-)上
        result = []
        proteinID = ""#当前基因的产物蛋白质ID
        geneID = ""#当前基因的ID
        gene = ""#存储当前基因的名字
        tempResult = []#存储该基因gff描述信息
        cateDcp = ""
#gff文件内容
#第一列：参考序列，是chromosome or scaffold的编号   即temp[0]
#第二列：注释信息的来源，一般为数据库例或者注释的机构，如果未知，用“."代替
#第三列：注释信息的类型，比如gene、mRNA、exon、CDS、UTR等
#第四列：第三列的注释类型在参考序列上的起始位置
#第五列：第三列的注释类型在参考序列上的终止位置
#第六列：得分，是注释信息可能性的说明，可以是序列相似性比对时的E - values值或者基因预测是的P - values值，“.”表示为空
#第七列：该基因或转录本位于参考序列的正链(+)或负链(-)上
#第八列：这列注释信息仅对第三列为CDS的类型有效，表示起始编码的位置，有效值为0、1、2，0表示该编码框的第一个密码子第一个碱基位于其5末端；1表示该编码框的第一个密码子的第一个碱基位于该编码区外；2表示该编码框的第一个密码子的第一、二个碱基位于该编码区外
#第九列：包含众多注释信息，以多个键值对组成的注释信息描述，不同属性之间以分号相隔，信息比较多我们一一解释：
# ID - -注释信息的编号，在一个GFF文件中必须唯一
# Name - -注释信息的名称，可以重复；
# Alias - -别名
# Parent - -指明feature所从属的上一级ID。用于将exons聚集成transcript，将transripts聚集成gene
# Note - -备注
# Dbxref - -数据库索引
#         with open("D:/data/" + x[0] + ".gff") as f:#x[0]="Anas platyrhynchos"
        with open("gff/"+x[0] + ".gff") as f:  # x[0]="Anas platyrhynchos"
            justify_chrName_flag=0
            if x[0]=="Mauremys reevesii":
                justify_chrName_flag=1
            for line in f.readlines():
                temp = line.split("\t")
    #1、-------------------- 该行是注释信息，滤注释信息--------------------
                if len(temp) < 6:
                    continue
                # elif temp[2] == 'CDS':
                #     continue
    #2、--------------------该行是染色体信息--------------------
                elif temp[2] == "region":
                    #re.search()并不要求必须从字符串的开头进行匹配，也就是说，正则表达式可以是字符串的一部分。
                    #. 代表一个任意的字符,这个字符默认包含所有的不包括换行符在内的所有字符 *代表倍数,只对该符号前一个字符有效
                    #元字符?号也是重复类字符, 但他表示可选  ()字符代表分组,代表一个整体
                    chromo = re.search(r"chromosome=(.*?);", line)
                    mmm_Name=re.search(r"Name=(.*?);", line)
                    if chromo is None:
                        chromo = "UnKnown"
                        if justify_chrName_flag == 1:  # 需要修正乌龟的染色体名字
                            chromo = mmm_Name.group(1)
                    else:#这一行数据中chromosome有对应的值，那么将其存入染色体数组中
                        chromo = chromo.group(1)#group()是按照特定子组数字返回结果，如上面匹配的时候()给正则匹配结果分组了，
                        # group(1)对应正则表达式对象的特定子组1，即第一个()中内容
                    chromosomes[temp[0]] = {#增加一个新对象  #temp[0]是chromosome or scaffold的编号
                        "chromosome": chromo,#chromo代表染色体的名字
                        "start": int(temp[3]),
                        "end": int(temp[4]),
                    }
                    continue
    # 3、--------------------该行是基因信息(处理gene,exon,cds,product)-------------多行表示同一个gene，他们名字相同，中间存在许多exon和cds，
                                                                          # 但是这些cds只对应一个proteinid
                else:
                    # pseudogene：假基因也叫伪基因，他是基因家族在进化过程中形成的无功能的残留物。它与正常基因相似，但丧失正常功能的DNA序列
                    # 假基因可视为基因组中与编码基因序列非常相似的非功能性基因组 DNA 拷贝，一般情况都不被转录，且没有明确生理意义
                    if temp[2] == "gene" or temp[2] == "pseudogene":
                        gene = re.search(r"gene=(.*?);", line)
                        if gene is None:
                            # 元字符[] 表示一个范围,相当于指定匹配一个范围类的字符,
                            gene = re.search(r"GeneID:(.*?)[;,]", line)
                        gene = gene.group(1)#获取基因ID，如此处：gene=104077434
                    if gene == None:#未获取到ID跳过此行
                        continue
                    # 存在多行名字相同的gene，多行表示同一个geneA
                    #进入第一个if，代表一个gene已经完全处理完毕，要处理新的gene了。这个if中的语句都是初始化的语句
                    if gene != previousGene:
                        if previousGene != "":#处理第一个基因
                            addGene(previousGene)#此行代码执行，表示所有名字相同的geneA都处理完毕，准备加入genes数组
                        previousGene = gene
                        geneID = re.search(r"GeneID:(.*?)[;,]", line).group(1)
                        chk = False
                        chk2 = True
                        chk3 = False
                        chk_pseudo = False
                        chk_lowquality = False
                        exon = 0#这个基因的外显子数目，若当前外显子的位置不在exonPos中，那么exon数目+1
                        exonPos = []#记录所有外显子的位置
                        chromosome = temp[0]
                        start = temp[3]
                        end = temp[4]
                        product = ""
                        proteinID = ""
                        tempResult = []  # 存储该基因gff描述信息
                        cateDcp=""#存储具体划分原因
                        if temp[6] == "+" or temp[6] == "-":#该基因或转录本位于参考序列的正链(+)或负链(-)上
                            direction = temp[6]
                    if temp[2] == "CDS":#CDS区是编码区，编码蛋白的。mRNA是包括编码区和非编码区
                        tempProteinID = re.search(r"protein_id=(.+?)(?:;|$)", line)
                        if tempProteinID is not None:
                            proteinID = tempProteinID.group(1)#一个gene对应一个proteinID，即许多行的proteinID是相同的，跟gene一样
                    if temp[2] == "exon":#若当前外显子的位置不在exonPos中，那么exon数目+1
                        tempExon = [int(temp[3]), int(temp[4])]
                        if tempExon not in exonPos:
                            exonPos.append(tempExon)
                            exon += 1
                    chk = chk or "interferon" in line#存在干扰素
                    # chk2为1时当且仅当，此行不存在regulator、gamma、receptor...，含有这些文字代表II型干扰素
                    if "silencer" not in line and "biological_region" not in line and "enhancer" not in line:
                        chk2 = chk2 and "regulator" not in line#不存在regulator
                    chk2 = chk2 and "gamma" not in line#不存在gamma
                    chk2 = chk2 and "receptor" not in line#不存在receptor
                    chk2 = chk2 and "responsive" not in line#不存在responsive
                    chk2 = chk2 and "induced" not in line#不存在induced
                    chk2 = chk2 and "inducible" not in line#不存在inducible
                    chk2 = chk2 and "factor" not in line#不存在factor
                    # if "silencer" not in line and "biological_region" not in line:
                    chk2 = chk2 and "stimulator" not in line#不存在stimulator
                    chk2 = chk2 and "stimulated" not in line#不存在stimulated
                    chk2 = chk2 and "lambda" not in line#不存在lambda

                    if "regulator" in line and "regulator" not in cateDcp:
                        cateDcp=cateDcp+"regulator"+"/"
                        tempResult.append(line.replace("\n", "\t") + "\n")
                    if "receptor" in line and "receptor" not in cateDcp:
                        cateDcp=cateDcp+"receptor"+"/"
                        tempResult.append(line.replace("\n", "\t") + "\n")
                    if "responsive" in line and "responsive" not in cateDcp:
                        cateDcp=cateDcp+"responsive"+"/"
                        tempResult.append(line.replace("\n", "\t") + "\n")
                    if "induced" in line and "induced" not in cateDcp:
                        cateDcp=cateDcp+"induced"+"/"
                        tempResult.append(line.replace("\n", "\t") + "\n")
                    if "inducible" in line and "inducible" not in cateDcp:
                        cateDcp=cateDcp+"inducible"+"/"
                        tempResult.append(line.replace("\n", "\t") + "\n")
                    if "factor" in line and "factor" not in cateDcp:
                        cateDcp=cateDcp+"factor"+"/"
                        tempResult.append(line.replace("\n", "\t") + "\n")
                    if "stimulator" in line and "stimulator" not in cateDcp:
                        cateDcp=cateDcp+"stimulator"+"/"
                        tempResult.append(line.replace("\n", "\t") + "\n")
                    if "stimulated" in line and "stimulated" not in cateDcp:
                        cateDcp=cateDcp+"stimulated"+"/"
                        tempResult.append(line.replace("\n", "\t") + "\n")
                    if "lambda" in line and "lambda" not in cateDcp:
                        cateDcp=cateDcp+"lambda"+"/"
                        tempResult.append(line.replace("\n","\t")+"\n")



                    # # chk2 = chk2 and 'precursor' not in line
                    # chk2 = chk2 and "LOW QUALITY" not in line#不存在LOW QUALITY
                    # chk2 = chk2 and "pseudogene" not in line#不存在pseudogene
                    if "LOW QUALITY" in line:
                        chk_lowquality=True
                    if "pseudogene" in line:
                        chk_pseudo=True
                # chk2为1时，表示I型干扰素
                # ---对基因名字进行修改，将其修改为缩写---
                    if "methylthioadenosine phosphorylase" in line:
                        previousGene = gene = "MTAP"

                    if "MOB Kinase Activator 3B" in line:
                        previousGene = gene = "MOB3B"
                    if "Focadhesin" in line:
                        previousGene = gene = "FOCAD"
                    if "3-Hydroxyacyl-CoA Dehydratase 4" in line:
                        previousGene = gene = "HACD4"
                    if "Ubiquitin Associated Protein 2" in line:
                        previousGene = gene = "UBAP2"
                    if "Ubiquitin Conjugating Enzyme E2 R2" in line:
                        previousGene = gene = "UBE2R2"
                    description = re.search(r"description=(.*?);", line)#对于gene产物的描述
                    if description is not None:
                        product = description.group(1)
                    if product == "":
                        temp = re.search(r"product=(.*?);", line)
                        if temp is not None:
                            product = temp.group(1)
        addGene(previousGene)#处理最后一个gene
# x[0]="Anas platyrhynchos"
        json.dump(#存储该物种的所有染色体
            chromosomes,
            codecs.open("data/chromosomes/" + x[0] + ".json", "w", encoding="utf-8"),
            ensure_ascii=False,
        )
        json.dump(#存储该物种的所有基因(有没有干扰素都存了)
            genes,
            codecs.open("data/genes/" + x[0] + ".json", "w", encoding="utf-8"),
            ensure_ascii=False,
        )
        json.dump(  # 存储该物种的所有基因(有没有干扰素都存了)
            IFNRelated,
            codecs.open("myowndata/detail_info/" + x[0] + ".json", "w", encoding="utf-8"),
            ensure_ascii=False,indent=2,
        )

        updategenes(x[0])#更新genes中被错误的划分到type=2的基因
        #加载刚刚提取的各个物种的gene信息
        genes = json.load(
            codecs.open("data/genes/" + x[0] + ".json", "r", encoding="utf-8")
        )
        print(time.time() - startTime)
        result = []
        results = []
        cnt = 0#控制满足条件的基因A的周围的基因B的数目
        for length in [2, 3, 4, 5]:#length指从特定基因（HACD4等菜单中所列基因以及干扰素相关基因）向两边扩展的基因数量
            result = []
            results = []
            cnt = 0
            for gene in genes:
                add(gene, length)
            json.dump(
                results,
                codecs.open( #这个文件里面存储的是某个物种的基因(即满足条件的基因A和周围的基因B的)
                    "data/interferonLocus/collections/"
                    + str(length) #表示周围基因的数目
                    + "_"
                    + x[0]
                    + ".json",
                    "w",
                    encoding="utf-8",
                ),
                ensure_ascii=False,
            )


if __name__ == "__main__":

    # client = MongoClient('localhost', 27017)
    # collection = client.interferon.chromosomes
    # specises = ['Mus_musculus', 'Oryctolagus_cuniculus', 'Tupaia_chinensis', 'Galeopterus_variegatus',
    #             'Suncus_etruscus', 'Manis_pentadactyla', 'Rhinolophus_ferrumequinum', 'Panthera_uncia', 'Equus_quagga',
    #             'Orcinus_orca', 'Choloepus_didactylus', 'Dasypus_novemcinctus', 'Elephantulus_edwardii',
    #             'Orycteropus_afer_afer', 'Elephas_maximus_indicus', 'Trichechus_manatus_latirostris']
    # for specise in specises:
    #     specise=specise.replace("_"," ")
    #     extractInterferon([[specise]])#如果没传入参数，那么就加载"data/species.json"中的所有物种作为参数
    #所有物种干的gff文件存在"D:/data文件中
    # 第一波：
    # extractInterferon([["Felis catus", "猫"]])
    # extractInterferon([["Canis lupus familiaris", "狗"]])
    # extractInterferon([["Panthera leo", "狮子"]])
    # extractInterferon([["Ailuropoda melanoleuca", "大熊猫"]])
    # extractInterferon([["Mustela putorius furo", "貂"]])
    #第二波：
    # extractInterferon([["Puma yagouaroundi", "南美洲短尾猫"]])
    # extractInterferon([["Puma concolor", "美洲狮"]])
    # extractInterferon([["Panthera pardus", "豹"]])
    # extractInterferon([["Acinonyx jubatus", "猎豹"]])
    # extractInterferon([["Neofelis nebulosa", "云豹"]])
    # extractInterferon([["Lynx rufus", "山猫"]])
    # extractInterferon([["Prionailurus viverrinus", "豹猫"]])
    # extractInterferon([["Panthera onca", "美洲豹"]])



    species = [
        ["Canis lupus", "灰狼"],
        ["Canis familiaris", "狗"],
        ["Alopex lagopus", "北极狐"],
        ["Canis lupus dingo", "澳洲野犬"],
    ]
    extractInterferon([["Canis lupus", "灰狼"]])
    extractInterferon([["Canis familiaris", "狗"]])
    extractInterferon([["Alopex lagopus", "北极狐"]])
    extractInterferon([["Canis lupus dingo", "澳洲野犬"]])

    #猫：Felis catus
    #狗：Canis lupus familiaris
    #狮子：Panthera leo
    ##大熊猫：Ailuropoda melanoleuca
    #貂：Mustela putorius furo
