#物种.txt中不能留有空行，不然会报错
'''
下载NCBI上物种过的gff文件
'''
import os
from time import sleep
from selenium import webdriver
from selenium.webdriver import ActionChains
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import Select#操作select标签
from selenium.webdriver.support.wait import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium import webdriver
specises=[]
ename=[]#存放英文
# cname=[]#存放中文
uncompeleted=[]#存放下载失败的物种
with open("gff下载物种.txt",'r') as fp:
    lines=fp.readlines()
    for line in lines:
        line=line.replace("\n",'')
        e=line.split("_")[0]
        ename.append(e)
print(ename)
#specises=['Mus_musculus', 'Oryctolagus_cuniculus', 'Tupaia_chinensis', 'Galeopterus_variegatus', 'Suncus_etruscus', 'Manis_pentadactyla', 'Rhinolophus_ferrumequinum', 'Panthera_uncia', 'Equus_quagga', 'Orcinus_orca', 'Choloepus_didactylus', 'Dasypus_novemcinctus', 'Elephantulus_edwardii', 'Orycteropus_afer_afer', 'Elephas_maximus_indicus', 'Trichechus_manatus_latirostris']
#ename=specises
options = webdriver.ChromeOptions()
out_path = r'./new_data'  # 下到new_data中
prefs = {'profile.default_content_settings.popups': 0, 'download.default_directory': out_path}
options.add_experimental_option('prefs', prefs)
# options.add_argument("--headless")  # => 为Chrome配置无头模式

#mac上下载到访达中
def download_whole_faa(species):
    #判断是否存在了此压缩包，存在了就不下载了
    file = "/Users/yyc/Downloads"
    file_list= os.listdir(file)
    name=species+'.tar'
    if name in file_list:
        print(species+'  已存在！！！')
        return 2
    else:
        print("开始下载"+species+'...')
    url = 'https://www.ncbi.nlm.nih.gov/'
    driver = webdriver.Chrome(executable_path='./chromedriver', options=options)
    driver.get(url)
    #---------选择下拉框中的基因-------------
    # 等待select元素出现
    WebDriverWait(driver,20).until(EC.visibility_of_element_located((By.ID,'database')))
    # 1.找到select元素
    select_ele = driver.find_element(By.ID,'database')
    # 2.初始化select类
    s = Select(select_ele)
    #3.根据value来选
    s.select_by_value('assembly')

    # ---------向输入框中输入要查的值-------------
    input_ele = driver.find_element(By.ID, 'term')
    keys = species
    input_ele.send_keys(keys)
    # -----------点击按钮-------------
    btn_ele = driver.find_element(By.ID, 'search')
    btn_ele.click()
    # 跳转操作页面
    driver.switch_to.window(driver.window_handles[-1])
    #点击download
    # js = 'document.getElementByXpath("//*[@id="download-asm"]/a").click()'
    # driver.execute_script(js)
    try:
        a = driver.find_element_by_xpath( '//*[@id="download-asm"]/a')
        size=a.size
        ActionChains(driver).move_to_element_with_offset(a,size["width"]/2,size["height"]/2).click().perform()
    except:
        driver.close()
        print("NCBI上不存在" + species + "对应序列")
        return 0
    #-----选择protein_faa
    # 1.找到select元素
    try:
        select_ele = driver.find_element(By.ID,'dl_assembly_file_types')
        # 2.初始化select类
        s = Select(select_ele)
        #3.根据value来选
        s.select_by_value('GENOME_GFF')
    except:
        driver.close()
        print("NCBI上不存在" + species + "对应序列")
        return 0
   #----开始下载----------
    try:
        js = 'document.getElementById("dl_assembly_download").click()'
        driver.execute_script(js)
    except:
        driver.close()
        print("NCBI上不存在" + species+"对应序列")
        return 0

    # 提前关闭
    i=0
    name = 'genome_assemblies_genome_gff.tar'
    while i<10:
        i=i+1
        all_name = os.listdir(file)
        if name in all_name:
            break
        else:
            sleep(60)
    if i==10:
        print("由于网络原因"+species+"下载失败！！！！！！！！！！！！")
        driver.close()
        return 0
    # sleep(300)
    driver.close()
    return 1
#修改访达中的文件名称
def rename(species):
    file = "/Users/yyc/Downloads"
    all_name = os.listdir(file)
    name='genome_assemblies_genome_gff.tar'
    if name in all_name:
        os.rename(file + '/' + name, file+'/'+species+'.tar')
        print(species+'下载完成！！！')
    else:
        return

def judge():
    file = "/Users/yyc/Downloads"
    file_list = os.listdir(file)
    # print(file_list)
    for s in file_list:
        path=file+'/'+s
        if 'genome_assemblies' in s:
            # print(path)
            os.remove(path)
if __name__=='__main__':
    for index,value in enumerate(ename):
        flag=download_whole_faa(value)
        value=value.replace("_"," ")
        if flag == 1:
            rename(value)
        if flag == 0:  # uncompeleted是全局变量 里面存储下载失败的物种，也有可能是物种不存在faa序列
            uncompeleted.append(value)
        judge()
    #未成功下载的物种写入下载情况.txt
    f = open("gff文件下载情况.txt", "w",encoding='utf-8')
    for node in uncompeleted:
        strs=node
        f.write(strs+ '\n')
    f.close()
    # judge()




