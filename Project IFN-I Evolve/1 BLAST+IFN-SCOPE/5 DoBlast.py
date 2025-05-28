import os
'''
获取root下的txt文件，里面的物种都是存在数据库的物种
interferon_path中存放所有要比较的干扰素（737个）
db_path存放制作好的数据库
out_path存放比较的结果
log_path存放比对日志
'''
root="/home/yuyangchao/species87.txt"#这里面存的是要比对的物种
# root=r"C:\Users\Administrator.DESKTOP-N60LKMN\Desktop\yyc/species.txt"#这里面存的是要比对的物种


interferon_path=r"/home/yuyangchao/species_3_94"
db_path=r"/home/yuyangchao/db"
out_path=r"/home/yuyangchao/blast_0807_87spe_94IFN"
log_path=r"/home/yuyangchao/log"#生成blast_failed.txt,存放比对失败的干扰素和数据库



def get_interferons():
    datanames = os.listdir(interferon_path)
    result=[]
    for name in datanames:
        if '.faa' in name:
            name=name.replace(".faa","")
            result.append(name)
    return result

def get_dbfile():
    specises = []
    with open(root, 'r', encoding='utf-8') as fp:
        lines = fp.readlines()
        for line in lines:
            line = line.replace("\n", '')
            specises.append(line)
    ename = []
    for name in specises:
        name = name.split("_")[0].replace(" ", "_")
        ename.append(name)
    return ename
def blast(gene,db):
    #输入gene:'Acanthisitta_chloris_LOC103797348'
    #输入db:'Anas_platyrhynchos'
    query_all=interferon_path+'/'+gene+'.faa'
    db_all=db_path+r'/'+db
    opath=out_path+r'/'+db
    if (os.path.exists(opath)==0):
        os.mkdir(opath)
    out_all=opath+r'/'+gene+'.txt'
    #print(''.join(['tblastn -query ', query_all, ' -db ',db_all, ' -out ', out_all, ' -evalue 1e-5 -outfmt 0']))
    os.system(''.join(['tblastn -query ', query_all, ' -db ',db_all, ' -out ', out_all, ' -evalue 1e-5 -outfmt 0']))
#tblastn -query /home/yuyangchao/737/Xenopus_laevis_LOC121399023.faa -db /home/yuyangchao/db/chromosome_1.part2_NC_056727.1 -out /home/yuyangchao/blast_out/Xenopus_laevis_LOC121399023.txt -evalue 1e-5 -outfmt 0
#tblastn -query C:\Users\Administrator.DESKTOP-N60LKMN\Desktop\yyc/737/Xenopus_laevis_LOC108702185.faa -db C:\Users\Administrator.DESKTOP-N60LKMN\Desktop\yyc\db/chromosome_9.part0_NC_056743.1 -out C:\Users\Administrator.DESKTOP-N60LKMN\Desktop\yyc/Xenopus_laevis_LOC108702185.txt -evalue 1e-5 -outfmt 0

if __name__ == '__main__':
    interferons=get_interferons()
    db_names=get_dbfile()
    failed=[]
    db_num=len(db_names)
    interferons_num=len(interferons)
    # print(db_names)
    print("the num of db:",db_num)
    print("the num of ifn:", interferons_num)
    for db_index,db in enumerate(db_names):
        for gene_index,gene in enumerate(interferons):
            try:
                print("db_processing:",db_index+1,'/',db_num,"  ","ifn_processing:",gene_index+1,'/',interferons_num)
                blast(gene,db)
                print("blast",gene,"成功")
            except:
                allfailname=gene+'***'+db
                failed.append(allfailname)
    f = open(log_path + "/blast_failed_0807_87spe_94IFN.txt", "w", encoding='utf-8')
    for node in failed:
        f.write(node + '\n')
    f.close()
