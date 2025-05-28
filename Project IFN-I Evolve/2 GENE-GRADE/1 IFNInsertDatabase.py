import json
import pymongo

request_str = 'mongodb://root:123456@localhost:27017/?authMechanism=DEFAULT'
def insertDatabase():
    # client = pymongo.MongoClient('localhost', 27017)
    # request_str = 'mongodb+srv://yyc:yyc123@cluster0.4ptltou.mongodb.net/test'
    client = pymongo.MongoClient(request_str)
    # allGenes = client.interferon.AllGeneRightName
    allGenes = client.interferon.dog4
    allGenes.delete_many({})
    # species = json.load(open('data/sortedSpecies-87.json', 'r', encoding='utf-8'))
    species = json.load(open('data/dog.json', 'r', encoding='utf-8'))
    for index,x in enumerate(species):
        print(index, "/", len(species))
        sp = x[0]
        print(x, "开始")
        genes = json.load(
            open('data/genes/'+sp+'.json', 'r', encoding='utf-8'))
        for y in genes:
            allGenes.insert_one({
                'species': sp,
                'gene': y['gene'],
                'chromosome': y['chromosome'],
                'start': y['start'],
                'end': y['end'],
                'type':y['type'],
                'product':y['product'],
                'flag':y['flag'],
                'direction':y['direction']
            })



def labelInterferon():
    client = pymongo.MongoClient(request_str)
    # allGenes = client.interferon.AllGeneRightName
    allGenes = client.interferon.dog4
    cnt = 0
    # species = json.load(open('data/sortedSpecies-87.json', 'r', encoding='utf-8'))
    species = json.load(open('data/dog.json', 'r', encoding='utf-8'))
    for index,s in enumerate(species):
        print(index,"/",len(species))
        genes = json.load(
            open('data/genes/'+s[0]+'.json', 'r', encoding='utf-8'))
        for x in genes:
            cnt += 1
            if cnt % 10000 == 0:
                print(cnt)
            if (x['type'] == 0 or x['type'] == 1):
                # print(x['chromosome'], x['start'])
                allGenes.update_one({'chromosome': x['chromosome'], 'start': x['start']}, {
                                    '$set': {'type': x['type']}})


def updateDatabase():
    client = pymongo.MongoClient(request_str)
    # allGenes = client.interferon.AllGeneRightName
    # allGenes = client.interferon.dog4
    # change Klhl9 to KLHL9
    # allGenes.update_one({'gene': 'Klhl9'}, {'$set': {'gene': 'KLHL9'}})


if __name__ == "__main__":

    insertDatabase()
    #所有步骤完成了之后后面两个函数再执行。必须再getSequence和getproteinSequence完成之后
    # labelInterferon()
    updateDatabase()
