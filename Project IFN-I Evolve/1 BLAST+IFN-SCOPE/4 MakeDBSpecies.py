import os
from pathlib import Path
'''
自动获取fna_path下所有的fna文件将其制作为数据库（生成8个文件）保存至db_path下
若在db_path中已经存在数据库了，那么就不再制作
制作失败的物种会生成日志文件makedb_failed_species.txt保存到log_path下
'''
fna_path=r"/home/yuyangchao/out"
db_path="/home/yuyangchao/db"
log_path="/home/yuyangchao/log"#log文件：makedb_failed_species.txt
# fna_path=r"G:\干扰素数据\西非肺鱼"
# db_path=r"C:\Users\Administrator.DESKTOP-N60LKMN\Desktop\yyc\out"
# log_path=r"C:\Users\Administrator.DESKTOP-N60LKMN\Desktop\yyc\out"#log文件：makedb_failed_species.txt
def makeDB(genome,dbpath=db_path):
    db_path = dbpath
    db = Path(db_path)
    for file in db.rglob('*.nsq'):  # 这样遍历出来的file全部都是绝对路径
        species = file.stem
        if(species==genome):
            print(genome+"的数据库已存在")
            return
    genome = genome.replace(' ', '_')
    os.system(''.join(['makeblastdb -in ',fna_path+'/', genome,
              '.fna -dbtype nucl -out ',db_path+'/', genome]))
# makeblastdb -in C:\Users\Administrator.DESKTOP-N60LKMN\Desktop\yyc\out\chromosome_9.part0_NC_056743.1.fasta -dbtype nucl -out C:\Users\Administrator.DESKTOP-N60LKMN\Desktop\yyc\db\chromosome_9.part0_NC_056743.1
# tblastn -query C:\Users\Administrator.DESKTOP-N60LKMN\Desktop\yyc\737\Acanthisitta_chloris_LOC103797348.faa -db C:\Users\Administrator.DESKTOP-N60LKMN\Desktop\yyc\db\chromosome_9.part0_NC_056743.1 -out C:\Users\Administrator.DESKTOP-N60LKMN\Desktop\yyc\test_out\1.txt -evalue 1e-5 -outfmt 0
if __name__ == '__main__':
    specises = []
    root = "./species.txt"
    with open(root, 'r', encoding='utf-8') as fp:
        lines = fp.readlines()
        for line in lines:
            line = line.replace("\n", '')
            # line=line.replace(" ","_")
            specises.append(line)
    ename = []
    for name in specises:
        name = name.split("_")[0].replace(" ", "_")
        if name!="Protopterus_annectens":#去除西非肺鱼
            ename.append(name)
    print("ename")
    print(ename)
    print(len(ename))
    failed=[]
    for name in ename:
        try:
            makeDB(name)
        except:
            failed.append(name)
    f = open(log_path+"/makedb_failed_species.txt", "w", encoding='utf-8')
    for node in failed:
        f.write(node + '\n')
    f.close()


