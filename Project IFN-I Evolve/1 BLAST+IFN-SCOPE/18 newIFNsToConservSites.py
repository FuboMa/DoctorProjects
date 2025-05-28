import json
import os
import re





conserv_sites_path = r'E:\Study\Doctor\Projects\IFN\代码\mfb\IFNSeqs\shape3.json'
out_path = r'E:\Study\Doctor\Projects\IFN\代码\mfb\IFNSeqs\queryByClustersOutput'
newIFN_path = r'E:\Study\Doctor\Projects\IFN\代码\mfb\IFNChecks\NewIFNs'


def get_filenames(path):
    names = os.listdir(path)
    return names


def summarize_IFN_infos(IFN_path):
    IFN_infos = []
    IFN_filenames = get_filenames(IFN_path)
    for filename in IFN_filenames:
        file_path = f'{IFN_path}/{filename}'
        with open(file_path) as f:
            f_lines = f.readlines()
        for line in f_lines:
            if line.startswith('>'):
                info_dict = get_info_dict(line)
                if info_dict not in IFN_infos:
                    IFN_infos.append(info_dict)
    return IFN_infos


def get_info_dict(line):
    tmp_list = line.split('_POS')[0]
    NC_flag = False
    if 'NC' in tmp_list:
        tmp_list = tmp_list.split('_NC_')
        NC_flag = True
    else:
        tmp_list = tmp_list.split('_NW_')
    species_name = ''.join(tmp_list[0]).strip('>')
    if NC_flag:
        chr_name = 'NC_' + tmp_list[1]
    else:
        chr_name = 'NW_' + tmp_list[1]
    str_start = r'start:(\d*),'
    IFN_start = int(re.findall(str_start, line)[0])
    str_end = r'end:(\d*),'
    IFN_end = int(re.findall(str_end, line)[0])
    IFN_pos = round((IFN_start + IFN_end) / 2)
    str_dir = r'direction:(.),'
    IFN_dir = re.findall(str_dir, line)[0]
    info_dict = {'species': species_name, 'chr': chr_name, 'start': IFN_start,
                 'end': IFN_end, 'pos': IFN_pos, 'direction': IFN_dir}
    return info_dict


def new_IFN_to_conserv_sites(cpath, opath, IFN_infos):
    with open(cpath) as f:
        conserv_sites_file = json.load(f)
    cnt1, cnt2, cnt3, cnt4, cnt5, cnt6 = 0, 0, 0, 0, 0, 0
    new_IFN_in_exceptions = []
    for IFN_info in IFN_infos:
        IFN_species = IFN_info['species']
        # tmp_flag = False
        for file_items in conserv_sites_file:
            file_species = file_items['species'].replace(' ', '_')
            if IFN_species == file_species:
                # tmp_flag = True
                cnt1 += 1
                IFN_chr = IFN_info['chr']
                IFN_start = IFN_info['start']
                IFN_end = IFN_info['end']
                IFN_dir = IFN_info['direction']
                IFN_pos = IFN_info['pos']
                IFN_dict = {'gene': 'NewIFN', 'pos': IFN_pos, 'start': IFN_start, 'end': IFN_end,
                            'direction': IFN_dir, "typeflag": "normal", "annoflag": "N", "isIFN": "Y"}
                hasChr_flag = False
                for species_chrs in file_items['chrs']:
                    species_chr = species_chrs['chr']
                    if IFN_chr == species_chr:
                        cnt2 += 1
                        hasChr_flag = True
                        in_site_flag = False
                        # IFN_dir = IFN_info['direction']
                        # IFN_gene_pos = IFN_info['pos']
                        # IFN_dict = {'gene': 'new_IFNs', 'pos': IFN_gene_pos, 'direction': IFN_dir,
                        #             "typeflag": "normal", "annoflag": "N", "isIFN": "Y"}
                        for chr_site in species_chrs['IFN_bySites']:
                            min_gene_pos = 99999999
                            max_gene_pos = 0
                            site_genes = chr_site['IFN_bySite']
                            site_name = chr_site['site']
                            for gene in site_genes:
                                if gene['pos'] < min_gene_pos:
                                    min_gene_pos = gene['pos']
                                if gene['pos'] > max_gene_pos:
                                    max_gene_pos = gene['pos']
                            if min_gene_pos < IFN_pos < max_gene_pos:
                                in_site_flag = True
                                if not site_name == 'Exception':
                                    cnt3 += 1
                                else:
                                    cnt4 += 1
                                    # print(f'cnt4: {IFN_species}_{IFN_chr}_POS{IFN_gene_pos}不属于三保守位点')
                                    new_IFN_in_exceptions.append(f'{IFN_species}_{IFN_chr}_POS{IFN_pos}')
                                # value3.append(IFN_info)
                                chr_site['IFN_cnt'] += 1
                                site_genes.append(IFN_dict)
                                break
                        if not in_site_flag:
                            hasExceptionDict_flag = False
                            for chr_site in species_chrs['IFN_bySites']:
                                site_name = chr_site['site']
                                if site_name == 'Exception':
                                    hasExceptionDict_flag = True
                            if not hasExceptionDict_flag:
                                cnt5 += 1
                                new_IFN_in_exceptions.append(f'{IFN_species}_{IFN_chr}_POS{IFN_pos}')
                                exception_dict = {"site": "Exception", "IFN_cnt": 1, "IFN_bySite": [IFN_dict]}
                                species_chrs['IFN_bySites'].append(exception_dict)
                            else:
                                for chr_site in species_chrs['IFN_bySites']:
                                    site_name = chr_site['site']
                                    if site_name == 'Exception':
                                        # IFN_cnt = int(chr_site['IFN_cnt'])
                                        # IFN_cnt += 1
                                        chr_site['IFN_cnt'] += 1
                                        chr_site['IFN_bySite'].append(IFN_dict)
                if not hasChr_flag:
                    exception_dict = {"site": "Exception", "IFN_cnt": 1, "IFN_bySite": [IFN_dict]}
                    chr_dict = {'chr': IFN_info['chr'], 'IFN_bySites': [exception_dict], "chrStart": 0, "chrEnd": 0,
                                "transChrName": "Unknown"}
                    file_items['chrs'].append(chr_dict)
                    cnt6 += 1
                    new_IFN_in_exceptions.append(f'{IFN_species}_{IFN_chr}_POS{IFN_pos}')
        # if not tmp_flag:
        #     print(f'IFN_species: {IFN_species}, file_species: {file_species}')
    with open(f'{opath}/allIFNs_bySites.json', 'w') as f:
        json.dump(conserv_sites_file, f, indent=2)
    # print(f'共有{cnt1}个未注释IFN,其中：\n'
    #       f'{cnt2}个在已知染色体上,'
    #       f'{cnt6}个在已经创建的新染色体上;\n'
    #       f'{cnt3}个在已知三保守位点上,'
    #       f'{cnt4}个在已知Exception位点上,'
    #       f'{cnt5}个在已经创建的新Exception位点上;\n'
    #       f'因此，共有{cnt4 + cnt5 + cnt6}个IFN不在三保守位点，他们是：\n'
    #       f'{new_IFN_in_exceptions}')


if __name__ == "__main__":
    new_IFN_infos = summarize_IFN_infos(newIFN_path)
    new_IFN_to_conserv_sites(conserv_sites_path, out_path, new_IFN_infos)
