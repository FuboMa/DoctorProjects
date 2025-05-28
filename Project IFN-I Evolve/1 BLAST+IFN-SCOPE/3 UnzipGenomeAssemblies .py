import tarfile
import zipfile
from pathlib import Path
import os
import gzip
import shutil
# !/usr/bin/env python
# -*- coding: utf-8 -*-
root="./species.txt"
# inital_path='D:\severtest\pack'
inital_path='/home/yuyangchao/pack'
des_path=r'/home/yuyangchao/middle'
fna_path=r'/home/yuyangchao/out'
# #存放全基因压缩数据的文件
# inital_path='H:\全基因核苷酸序列'
# #存放解压之后数据的文件（这个文件是一个中间文件，代码运行完成后可以直接删除）
# des_path='H:\解压后的全基因核苷酸序列'
# #存放.fna数据的文件
# fna_path='H:\用于制作数据库的数据'
# #第一次解压
def un_tar(name):#这个name是带下划线的
    full_name=inital_path+"/"+name+".tar"
    tar = tarfile.open(full_name)
    names = tar.getnames()
    dst=des_path
    if os.path.isdir(dst):
        pass
    else:
        os.mkdir(dst)
    for name in names:
        tar.extract(name, dst)
    tar.close()
#删除第一次解压产生的txt文档
def deleteTxtFile(des_path):
  file_names = os.listdir(des_path)
  for name in file_names:
      if '.txt' in name:
        os.remove(des_path +'/'+ name)
#对第一次解压之后的文件夹改名（将ncbi-2022-11-7改名为物种名）
def changeName(des_path,species):
    file_names = os.listdir(des_path)
    for name in file_names:
        if 'ncbi' in name:
            old_name=des_path + '/'+name
            new_name=des_path + '/'+species
            os.rename(old_name, new_name)
#第二次解压
def un_gz(name,path=des_path):
    # 获取文件的名称，去掉后缀名
    dir_path=path+'/'+name
    files = os.listdir(dir_path)
    for name in files:
        if 'gz' in name:
            file_name=dir_path+'/'+name
            f_name = file_name.replace(".gz", "")
    # 开始解压
    g_file = gzip.GzipFile(file_name)
    # 读取解压后的文件，并写入去掉后缀名的同名文件（即得到解压后的文件）
    open(f_name, "wb+").write(g_file.read())
    g_file.close()
#第二次改名，将GCF_010993605.1_kPetMar1.pri_genomic.fna改名为物种.fna
def rename(species):
    file_path=des_path+'/'+species
    file_names = os.listdir(file_path)
    for name in file_names:
        if 'fna' in name and 'gz' not in name:
            old_name=file_path + '/'+name
            new_name=file_path + '/'+species+'.fna'
            os.rename(old_name, new_name)
#移动物种.fna文件到另一个文件夹




def remove_file(old_path, new_path):
    print(old_path)
    print(new_path)
    filelist = os.listdir(old_path)      #列出该目录下的所有文件,listdir返回的文件列表是不包含路径的。
    print(filelist)
    for file in filelist:
        src = os.path.join(old_path, file)
        dst = os.path.join(new_path, file)
        print('src:', src)                 # 原文件路径下的文件
        print('dst:', dst)                 # 移动到新的路径下的文件
        shutil.move(src, dst)
def move(species):
    old_path=des_path+'/'+species+'/'+species+'.fna'
    new_path=fna_path+'/'+species+'.fna'
    shutil.move(old_path, new_path)
def all(species):
    try:
        un_tar(species)
    except:
        print(species+"的第一次解压失败")
        return
    try:
        deleteTxtFile(des_path)
    except:
        print("删除第一次解压"+species+'产生的txt文档失败')
        return
    try:
        changeName(des_path,species)
    except:
        print('改变第一次解压'+species+'产生的文件夹名字(ncbi-2020-11)失败')
        return
    try:
        un_gz(species)
    except:
        print(species + "的第二次解压失败")
        return
    try:
        rename(species)
    except:
        print('修改GCF_010993605.1_kPetMar1.pri_genomic.fna的名字为'+ species+'失败' )
        return
    try:
        move(species)
    except:
        print('移动'+species+'.fna文件至最终文件夹失败')
        return
    print(species+".fna数据准备成功")
if __name__=="__main__":
    #species='Acanthisitta_chloris'#这个name是带下划线的
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
    print("ename")
    print(ename)
    print(len(ename))
    # for name in ename:
    #     if name=="Anser_cygnoides_domesticus" or name=="Corvus_cornix_cornix" or name=="Molothrus_ater" or name=="Pan_paniscus" or name=="Protopterus_annectens":
    #         continue
    #     all(name)
    # all(species)
