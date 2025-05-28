import json
import re
from IFNSCOPE import query_path

'''
1. The purpose of this document is to perform a secondary filtration on the results output by the IFN-SCOPE model (predicted_IFNs_summary.fasta) to obtain the final IFN results:

    1.1 Step 1: Filter again using annotated type I IFN genes (including pseudogenes and low-quality genes).
    
    1.2 Step 2: If multiple results share the same chromosome and have one of their transcription start or transcription end sites in common, filter out duplicates, retaining only one result.
    
    1.3 Output the filtered results.

2. Parameters to be modified:

    `predicted_IFNs_path`: Path to the FASTA file output by IFNSCOPE.py.
    
    `annoIFNs_path`: Directory containing annotated IFN genes (including pseudogenes and low-quality genes).
    
    `overlap`: Length of overlap to consider when filtering out a prediction if it overlaps with a annotated IFN gene (including pseudogenes and low-quality genes).
'''

predicted_IFNs_path = f'{query_path}/predicted_IFNs_summary.fasta'
annoIFNs_path = r'E:\Study\Doctor\Projects\IFN\代码\mfb\IFNSeqs\shape3-withpsedo.json'

overlap = 100


def filter_byAnnoIFNs(ppath, apath):
    with open(ppath) as f:
        pred_f = f.readlines()
    with open(apath) as f:
        annoIFN_f = json.load(f)
    before_cnt = 0
    after_cnt = 0
    result = []
    for n, line in enumerate(pred_f):
        filter_flag = False
        if line.startswith('>'):
            before_cnt += 1
            str_species = r'>(.*?)_N[CW]_'
            IFN_species = re.findall(str_species, line)[0]
            tmp_list = line.split(IFN_species)
            tmp_str = ''.join(tmp_list)
            str_chr = r'>_(.*?)_POS'
            IFN_chr = re.findall(str_chr, tmp_str)[0]
            str_start = r'start:(\d*),'
            IFN_start = int(re.findall(str_start, line)[0])
            str_end = r'end:(\d*),'
            IFN_end = int(re.findall(str_end, line)[0])
            for species_dict in annoIFN_f:
                anno_speciesname = species_dict['species'].replace(' ', '_')
                if IFN_species == anno_speciesname:
                    for species_chrs in species_dict['chrs']:
                        anno_chr = species_chrs['chr']
                        if IFN_chr == anno_chr:
                            for chr_site in species_chrs['IFN_bySites']:
                                for gene_dict in chr_site['IFN_bySite']:
                                    gene_isIFN = gene_dict['isIFN']
                                    if gene_isIFN == 'Y':
                                        gene_start = gene_dict['start']
                                        gene_end = gene_dict['end']
                                        max_start = max(IFN_start, gene_start)
                                        min_end = min(IFN_end, gene_end)
                                        if min_end - max_start >= overlap:
                                            filter_flag = True
                                            break
            if not filter_flag:
                result.append(line)
                result.append(pred_f[n + 1])
                after_cnt += 1
    print(f'Before filter: {before_cnt}')
    print(f'After filter: {after_cnt}')
    return result


def remove_duplications(result):
    contents = []
    infos = []
    seqs = []
    cur_start = 0
    cur_end = 0
    after_cnt = 0
    for n, line in enumerate(result):
        if line.startswith('>'):
            str_start = r'start:(\d*),'
            IFN_start = int(re.findall(str_start, line)[0])
            str_end = r'end:(\d*),'
            IFN_end = int(re.findall(str_end, line)[0])
            if (IFN_start != cur_start) and (IFN_end != cur_end):
                cur_start, cur_end = IFN_start, IFN_end
                contents.append(line)
                contents.append(result[n + 1])
                contents.append('\n')
                infos.append(line)
                seqs.append(result[n + 1])
                after_cnt += 1
    print(f'After removal: {after_cnt}')
    return contents, infos, seqs


def write_final_result(qpath, contents, infos, seqs):
    final_IFN_path = f'{qpath}/final_IFNs_summary.fasta'
    with open(final_IFN_path, 'w') as f:
        for line in contents:
            f.write(line)
    newIFNs_path = r'E:\Study\Doctor\Projects\IFN\代码\mfb\IFNChecks\NewIFNs'
    for n, info in enumerate(infos):
        filename = info.split(',')[0].strip('>')
        with open(f'{newIFNs_path}/{filename}.fasta', 'w') as f:
            f.write(info)
            f.write(seqs[n])


if __name__ == "__main__":
    filter_result = filter_byAnnoIFNs(predicted_IFNs_path, annoIFNs_path)
    IFN_contents, IFN_infos, IFN_seqs = remove_duplications(filter_result)
    write_final_result(query_path, IFN_contents, IFN_infos, IFN_seqs)
