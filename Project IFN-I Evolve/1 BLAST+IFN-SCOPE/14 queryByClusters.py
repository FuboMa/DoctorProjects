import re
import json
import os
import time
import traceback
import numpy as np
from sklearn.cluster import MeanShift, estimate_bandwidth

"""
1. This file annotate the positions of IFN genes in a certain species based on tblastn alignment results and outputs them to the corresponding FASTA file.
    
    1.1 Firstly, classify the queries in the tblastn file according to their respective chromosomes.
    
    1.2 Secondly, calculate the positions (site_centers) of the queries in the tblastn file using the MeanShift algorithm. 
        The number and density of these  positions are determined by the bd parameter in the find_query_labels() function.
    
    1.3 Classify all queries based on these positions (labels) and save the information accordingly in the structure: species → chromosome → positions → tblastn_result_infos file.
    
    1.4 Each tblastn_result_infos file will serve as an input test sequence for the IFN-SCOPE model in the future.

2. Parameters to be modified:
    
    `bd`: in find_query_labels(), determines the number and density of clusters. 
    
    `out_path`: the root directory for output.
   
    `tblastn_result_path`: the root directory for tblastn results.
"""

out_path = r'E:\Study\Doctor\Projects\IFN\代码\mfb\IFNSeqs\TESTOUTPUT'
tblastn_result_path = r'E:\Study\Doctor\Projects\IFN\代码\mfb\IFNSeqs\TESTFILE'


def get_filenames(path):
    names = os.listdir(path)
    return names


def init_pipeline(tpath, opath):
    names = get_filenames(tpath)  # Retrieve all JSON files under the tpath directory.
    print('\n')
    localtime = time.strftime("%Y-%m-%d", time.localtime())
    error_path = open(f'{opath}/error_species_queryByClusters_{localtime}.txt',
                      'w')  # Create an error_species file to document erroneous outputs.
    for name in names:
        species = name.split('.')[0]
        IFNs_bySpecies = []
        print('----------------------------------------')
        print('{}:'.format(species))
        folder_path = mkdir_forInit(opath, species)
        try:
            json_result_path = f'{tpath}/{name}'
            with open(json_result_path, encoding='utf-8') as f1:
                json_results = json.load(f1)
            for idx, j_result in enumerate(json_results['chr_node']):
                t_result = str(j_result)
                pipeline(species, j_result, t_result)  # Start pipeline
            print('{} complete'.format(species))
            f1.close()
        except:
            print(species, file=error_path)
            error_path.flush()
            exstr = traceback.format_exc()
            print(exstr)
        print('----------------------------------------')
    error_path.close()


def find_chr_bySites(result):
    """
    In the JSON file, search for the corresponding chromosome based on multiple queries of a
    specific conserved locus to provide data for the web crawler to concatenate strings.
    """
    # if result['chr'] == '':
    #     chromosome = [result['chromosome'], 'Unknown']
    # else:
    chromosome = [result['chromosome'], result['chr']]
    print('Chromosome:{}, Chr:{}'.format(chromosome[0], chromosome[1]))
    return chromosome


def find_query_labels(txt):
    seq_start = r'\'seq_start\': (.*?),'
    # starts = re.findall(seq_start, txt, re.S)
    starts = re.findall(seq_start, txt)
    seq_end = r'\'seq_end\': (.*?),'
    ends = re.findall(seq_end, txt)
    int_starts = []
    int_ends = []
    int_poss = []
    for n, start in enumerate(starts):
        int_start = int(start)
        int_end = int(ends[n])
        int_pos = round((int_start + int_end) / 2)
        int_starts.append(int_start)
        int_ends.append(int_end)
        int_poss.append(int_pos)
    array = np.array(int_poss)
    rearray = array.reshape(-1, 1)
    # n_samples = int(len(int_starts) / 1)
    bd = estimate_bandwidth(rearray, quantile=0.3)
    # if bd < 10:
    #     bd = 10
    # elif bd > 100:
    #     bd = 100
    print(f'   bandwidth = {bd}')
    ms = MeanShift(bandwidth=bd, min_bin_freq=5)
    ms.fit(rearray)
    labels = ms.labels_  # labels：所属分类
    centers = []
    for center in ms.cluster_centers_:
        centers.append(int(center))
    print('   Get query_labels complete')
    return centers, labels, int_starts, int_ends


def find_seq_infos(txt):
    idty = r'\'identity\': \'(.*?)\','
    gap = r'\'gap\': \'(.*?)\','
    direction = r'\'direction\': \'(.*?)\','
    len = r'\'seq_length\': (.*?),'
    sbjct_seq = r'\'sbjct_seq_conten\': (.*?),'
    idtys = re.findall(idty, txt)
    gaps = re.findall(gap, txt)
    dirs = re.findall(direction, txt)
    lens = re.findall(len, txt)
    sbjcts = re.findall(sbjct_seq, txt)
    for n, sbjct in enumerate(sbjcts):
        sbjcts[n] = sbjct.strip('\'')
    return idtys, gaps, dirs, lens, sbjcts


def classify_infos_byLabels(centers, labels, starts, ends, idtys, gaps, dirs, lens, sbjcts):
    infos_bL = [[] for x in range(len(centers))]  # bL = byLabels
    labels_to_indices = [[] for x in range(len(centers))]
    for n, label in enumerate(labels):
        labels_to_indices[label].append(n)
    for k, label in enumerate(labels_to_indices):
        for idx in label:
            info = [starts[idx], ends[idx], idtys[idx], gaps[idx], dirs[idx], lens[idx], sbjcts[idx]]
            infos_bL[k].append(info)
    return infos_bL


def output_seq_infos(opath, species, chromosome, centers, infos_bL):
    for n, center in enumerate(centers):
        folder_path = mkdir(opath, species, chromosome, center)
        with open(f'{folder_path}/tblastn_result_infos.fasta', 'w') as f:
            for infos in infos_bL[n]:
                f.write('>' + 'start:' + f'{infos[0]}' +
                        ', ' + 'end:' + f'{infos[1]}' +
                        ', ' + 'idty:' + f'{infos[2]}' +
                        ', ' + 'idty:' + f'{infos[2]}' +
                        ', ' + 'gap:' + f'{infos[3]}' +
                        ', ' + 'direction:' + f'{infos[4]}' +
                        ', ' + 'len:' + f'{infos[5]}' +
                        ', ' + '\n' + f'{infos[6]}' + '\r\n')
        print('   write fasta files complete')


def mkdir_forInit(path, species):
    folder_path = f'{path}/{species}'
    exist = os.path.exists(folder_path)
    if not exist:  # Check if the directory exists; if it does not, create it as a directory.
        os.makedirs(folder_path)
    return folder_path


def mkdir(opath, species, chromosome, center):
    chro = '_'.join(chromosome)
    folder_path = f'{opath}/{species}/{chro}/{center}'
    exist = os.path.exists(folder_path)
    if not exist:  # Check if the directory exists; if it does not, create it as a directory.
        os.makedirs(folder_path)
    return folder_path


def pipeline(species, j_result, t_result):
    list_chromosome = find_chr_bySites(j_result)
    site_centers, query_labels, seq_starts, seq_ends = find_query_labels(t_result)
    seq_idtys, seq_gaps, seq_dirs, seq_lens, seq_sbjcts = find_seq_infos(t_result)
    seq_infos_byLabels = classify_infos_byLabels(site_centers, query_labels, seq_starts, seq_ends, seq_idtys, seq_gaps,
                                                 seq_dirs, seq_lens, seq_sbjcts)
    output_seq_infos(out_path, species, list_chromosome, site_centers, seq_infos_byLabels)


if __name__ == "__main__":
    init_pipeline(tblastn_result_path, out_path)
