import re
import os
import difflib
import time
import traceback
from Digraph import Digraph
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.support.wait import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC

"""
1. This document consolidates each `queryByClusters.py` output, organizing species → chromosome → clustered position → `tblastn_result_infos` file into a single amino acid sequence, 
and compiles all amino acid sequences into one file for input into IFN-SCOPE model.
    
    1.1 Retrieve the "consensus sequence" from each `tblastn_result_infos` file along with its start, end, and other parameters.
    
    1.2 Using a web crawler, download the sequences corresponding to each query region (± elongation nucleotides) from the NCBI website based on the consensus sequence's start, end, and other parameters,
     translate them into amino acids, and store them in `aa_region.txt`.
    
    1.3 Identify all amino acid sequences in `aa_region.txt` that start with "M" and end with "*", 
     with a length greater than `min_IFN_length`, and store them in species → chromosome → clustered position → `possible_B_infos` file.
    
    1.4 Compile all `possible_IFN_infos` files.
    
2. Parameters to be modified:
    
    (The `conserv_ratio` parameter indicates the minimum proportion of sequences in `fasta_seqs` that the `consensus_seq` must include. 
     For example, if `conserv_ratio = 0.9`, then at least 90% of the sequences must still include the `consensus_seq` when extending to the right. 
     Lowering this value slightly increases the `k_value`. This parameter is currently deprecated as the method of assembling amino acid sequences for a locus is no longer intended to be used.)
    
    `elongation_nucleotides`: Indicates the nucleotide distance ± to be fetched from NCBI around the `assembly_seq`. 
     For example, if `elongation_nucleotides = 600`, nucleotides ranging from `assembly_seq_start - 600` to `assembly_seq_end + 600` need to be fetched and converted into amino acids.
    
    `min_IFN_length`: The minimum length of the output IFN amino acid sequence.
    
    `out_path`: The output directory for the new amino acid sequences.
    
3. Other parameter descriptions:
    `nt_aa_ratio`: The ratio of nucleotides to amino acids, where one amino acid in mRNA corresponds to three codons in tRNA.
"""

nt_aa_ratio = 3
conserv_ratio = 0.95
elongation_nucleotides = 300
min_IFN_len = 120

out_path = r'E:\Study\Doctor\Projects\IFN\代码\mfb\IFNSeqs\TESTOUTPUT'


options = webdriver.ChromeOptions()
options.add_argument('--no-sandbox')
options.add_argument('headless')
options.add_argument('--disable-infobars')
options.add_argument('--disable-dev-shm-usage')
s = Service('/usr/bin/chromedriver')
driver = webdriver.Chrome(service=s, options=options)


def get_filenames(path):
    names = os.listdir(path)
    return names


def init_pipeline(opath, elong_nt):
    localtime = time.strftime("%Y-%m-%d", time.localtime())
    error_path = open(f'{opath}/error_species_getIFNSeqs_{localtime}.txt',
                      'w')  # Create an error_species file to document erroneous outputs.
    species_names = get_filenames(opath)
    for species_name in species_names:
        species_path = f'{opath}/{species_name}'
        if os.path.isdir(species_path):
            print('\n' * 5)
            print('{}:'.format(species_name))
            chromosome_names = get_filenames(species_path)
            for chromosome_name in chromosome_names:
                # chromosome_name: NC_051803.1_W
                chromosome_path = f'{species_path}/{chromosome_name}'
                if os.path.isdir(chromosome_path):
                    chr_name = chromosome_name.rsplit('_', 1)[0]
                    # chr_name: NC_051803.1(actually accession_number)
                    while chr_name.count('_') != 1:
                        chr_name = chr_name.rsplit('_', 1)[0]
                        # Some chromosome_name formats are incorrect, such as NC_052627.1_NC_052627.1 or NC_054387.1_9_10L.
                        # This line of code addresses this issue.
                    print(f'chr_name: {chr_name}')
                    centers = get_filenames(chromosome_path)
                    for center in centers:
                        center_path = f'{chromosome_path}/{center}'
                        print(f'center: {center}')
                        fasta_names = get_filenames(center_path)
                        for fasta_name in fasta_names:
                            if 'tblastn_result_infos' in fasta_name:
                                fasta_path = f'{center_path}/{fasta_name}'
                                try:
                                    print(
                                        '----------------------------------------'
                                    )
                                    pipeline(center_path, fasta_path, chr_name,
                                             elong_nt)
                                except:
                                    print('{} {} {} {}'.format(
                                        species_name, chromosome_name,
                                        center_path, fasta_name),
                                        file=error_path)
                                    error_path.flush()
                                    exstr = traceback.format_exc()
                                    print(exstr)
            print('{} complete'.format(species_name))
    error_path.close()


def get_consensus_seq(fpath):
    # The purpose of this function is to retrieve the "consensus sequence" of all sequences in a given tblastn_result_infos file,
    # and to obtain the start, end, and other relevant values of the "consensus sequence".
    seqs = []
    with open(fpath) as f:
        f_lines = f.readlines()
    consen_seq = ''
    min_start = 999999999
    max_end = 0
    seq_dir = ''
    for line in f_lines:
        start = 0
        end = 0
        if line.startswith('>'):
            str_start = r'start:(\d*),'
            seq_start = int(re.findall(str_start,
                                       line)[0])
            if seq_start < min_start:
                min_start = seq_start
            str_end = r'end:(\d*),'
            seq_end = int(re.findall(str_end, line)[0])
            if seq_end > max_end:
                max_end = seq_end
            if not seq_dir:
                str_dir = r'direction:(.),'
                seq_dir = re.findall(str_dir, line)[0]
        elif line is not '\n':
            line = line.replace('-', '')
            seqs.append(line)
            if not consen_seq:
                consen_seq = line
                # cnt += 1
            else:
                dif = difflib.Differ()
                dif_result = list(dif.compare(line, consen_seq))
                flag = True
                for n, aa in enumerate(dif_result):
                    aa = aa.replace(' ', '')
                    if flag:
                        if len(aa) == 1:
                            start = n
                            flag = False
                    elif len(aa) != 1:
                        end = n - 1
                        break
                    elif n == len(dif_result) - 1:
                        end = n
                        break
                tmp_list = []
                for aa in dif_result[start:end + 1]:
                    aa = aa.replace(' ', '')
                    tmp_list.append(aa)
                consen_seq = ''.join(tmp_list)
    return seqs, consen_seq, min_start, max_end, seq_dir


def get_k_value(seqs, consen_seq):
    # The purpose of this function is to determine k_value, which must be greater than or equal to the length of the consensus_seq.
    conserv_cnt = len(seqs)
    conserv_lower_limit = len(seqs) * conserv_ratio
    new_seq = consen_seq
    k = 0
    while conserv_cnt > conserv_lower_limit:
        k = len(consen_seq)
        for seq in seqs:
            start = seq.find(consen_seq)
            end = start + k - 1
            if end == len(seq) - 1:
                continue
            else:
                new_end = end + 1
                tmp_list = [consen_seq, seq[new_end]]
                new_seq = ''.join(tmp_list)
                break
        for seq in seqs:
            ret = seq.find(new_seq)
            if ret == -1:
                conserv_cnt -= 1
        if consen_seq != new_seq:
            consen_seq = new_seq
        else:
            break
    if k % 2 == 0:
        k -= 1
    return k


def get_kmers(seqs, k):
    kmers = []
    for seq in seqs:
        # rseq = seq.replace('-', '')
        kmer_num = len(seq) - k + 1
        start = 0
        for n in range(kmer_num):
            end = start + k
            kmer = seq[start:end]
            start += 1
            kmers.append(kmer)
    kmers_freq = {}
    for kmer in kmers:
        kmers_freq[kmer] = kmers_freq.get(kmer, 0) + 1
    kmers = list(set(kmers))
    return kmers_freq, kmers


def deBrujinGraph_algorithm(kmers):
    graph = Digraph(len(kmers))
    for n, kmer in enumerate(kmers):
        overlap = kmer[1:]
        for m, kmer2 in enumerate(kmers):
            if kmer == kmer2:
                continue
            elif overlap in kmer2:
                graph.point_edge(n, m)
    DBG = graph.adj_list  # DBG = deBrujinGraph
    return DBG


def DFT_alogrithm(kmers, DBG, cur_idx, end_idxs, path):
    # DFT = Depth First Traversal
    path.append(kmers[cur_idx])
    if cur_idx not in end_idxs:
        for next_idxs in DBG[cur_idx]:
            DFT_alogrithm(kmers, DBG, next_idxs, end_idxs, path)
            path.pop()
    else:
        DBG_assemblies.append(path[::])


def get_assembly_by_DBG(kmers_freq, kmers, DBG):
    global DBG_assemblies
    DBG_assemblies = []
    # -----1. Retrieve nodes with out-degree of 0 (end nodes).-----
    end_idxs = []
    neighbor_kmers = []
    for node in DBG:
        for m in node:
            neighbor_kmers.append(m)
    neighbor_kmers = set(neighbor_kmers)
    # -----2. Retrieve nodes with in-degree of 0 (start nodes).-----
    start_idxs = []
    for n, node in enumerate(DBG):
        if n not in neighbor_kmers:
            start_idxs.append(n)  # The neighbor_kmers list stores all nodes with incoming edges, so those not present in it are potential start nodes.
        if not node:
            end_idxs.append(n)
    # -----3. Traverse all edges of the De Bruijn Graph (DBG) using Depth-First Traversal (DFT), and return all paths (usually only one) to `DBG_assemblies`.-----
    path = []
    for start_idx in start_idxs:
        cur_idx = start_idx
        DFT_alogrithm(kmers, DBG, cur_idx, end_idxs, path)
    # -----4. Iterate through `DBG_assemblies`, concatenate each k-mer into an amino acid sequence and store it in
    # `final_seq`, and return the sequence with the highest sum of k-mer frequencies.-----
    final_seqs = []
    final_seqs_freqSums = []  # Each element corresponds to the sum of frequencies of all k-mers traversed in the paths within the final_seqs.
    # max_cnt = 0
    # min_kmer_cnt = 999999999
    for assembly in DBG_assemblies:
        tmp_sum = 0
        flag = False
        seq = 0
        for kmer in assembly:
            # if cnt <= min_kmer_cnt:
            #     min_kmer_cnt = cnt
            if not flag:
                seq = kmer
                flag = True
            else:
                elongation_aa = kmer[-1]
                tmp_list = [seq, elongation_aa]
                seq = ''.join(tmp_list)
            tmp_sum += kmers_freq[kmer]
        final_seqs_freqSums.append(tmp_sum)
        # if min_kmer_cnt >= max_cnt:
        #     max_cnt = min_kmer_cnt
        #     if final_seqs:
        #         final_seqs.pop()
        final_seqs.append(seq)
    if len(final_seqs) == 1:
        ret_seq = final_seqs[0]
    else:
        max_freqSum = max(final_seqs_freqSums)
        max_idx = final_seqs_freqSums.index(max_freqSum)
        ret_seq = final_seqs[max_idx]
    return ret_seq


def find_aa_region(cpath, chr_name, start, end, direction, elong_nt):
    try:
        url = f'https://www.ncbi.nlm.nih.gov/nuccore/{chr_name}'
        aa_region = ''
        region_start = start - elong_nt
        region_end = end + elong_nt
        driver.get(url)
        try:
            print('   Enter NCBI url success')
            # -----1. Filter the specified region based on the region_start and region_end values and present it in FASTA format.-----
            try:
                button1_id = "EntrezSystem2.PEntrez.Nuccore.Sequence_ResultsPanel.Sequence_SingleItemSupl." \
                             "Sequence_ViewerGenbankSidePanel.Sequence_ViewerChangeRegion.Shutter"
                WebDriverWait(driver, 30).until(
                    EC.visibility_of_element_located((By.ID, button1_id)))
                button1 = driver.find_element(By.ID, button1_id)
                button1.click()  # button1 is a dropdown menu button.
                region_text, region_start, region_end = get_region_text(
                    chr_name, region_start, region_end, direction)
                print('   Crawl success')
                nucleotide_region_withBlank = region_text.split("\n", 1)[1]
                nucleotide_region = nucleotide_region_withBlank.replace(
                    '\n', '')
                if direction is '-':
                    nucleotide_region = reverse_complement(nucleotide_region)
                if len(nucleotide_region) % 3 == 0:
                    print('   Nucleotide region length is {}'.format(
                        len(nucleotide_region)))
                else:
                    print(
                        '   Note: Nucleotide region length cannot be divided by 3!'
                    )
                # -----3. Write into aa_region.txt.-----
                aa_region = translate(nucleotide_region)
                with open(f'{cpath}/aa_region.txt', 'w') as f:
                    f.write(aa_region)
            except:
                exstr = traceback.format_exc()
                print(exstr)
        except:
            exstr = traceback.format_exc()
            print(exstr)
    except:
        print('   Note: time out!')
        print(f'region_start: {region_start}')
        print(f'region_end: {region_end}')
    print('   Find aa_region complete')
    print('----------------------------------------')
    return aa_region, region_start, region_end


def get_region_text(chr_name, region_start, region_end, direction):
    try:
        WebDriverWait(driver, 30).until(
            EC.visibility_of_element_located((By.ID, "crfrom")))
        driver.find_element(By.ID, "crfrom").clear()
        driver.find_element(By.ID, "crfrom").send_keys(region_start)
        driver.find_element(By.ID, "crto").clear()
        driver.find_element(By.ID, "crto").send_keys(region_end)
        button2 = driver.find_element(By.ID, "updateselregion")
        button2.click()  # button2 is an update region button
        button3 = driver.find_elements(By.ID, "ReportShortCut6")
        if button3:
            button3[0].click(
            )  # button3 is the button for switching to FASTA format display;
            # if a new search is required due to elong_nt exceeding the chromosome region, button3 will be empty.
        # -----2. Filter the nucleotide_region sequence from the FASTA format.-----
        try:
            WebDriverWait(driver, 20).until(
                EC.visibility_of_element_located((By.TAG_NAME, "pre")))
            div = driver.find_element(By.TAG_NAME, "pre")
            region_text = div.get_attribute(
                'innerText')  # The text here refers to the FASTA format gene sequence displayed in the center of the screen, stored within the `pre` tag.
            # -----* Confirm whether the region exceeds the chromosome values from the FASTA format, and filter the nucleotide_region sequence accordingly.-----

            desc = region_text.split("\n", 1)[0]
            tmp_str = f'{chr_name}:(.*)'
            tmp_text = re.findall(tmp_str, desc)
            region_start_inText = 0
            region_end_inText = 0
            for item in tmp_text:
                region_inText = item.split(' ', 1)[0]
                region_start_inText = int(region_inText.split('-')[0])
                region_end_inText = int(region_inText.split('-')[1])
                break
            # -----* If the region_start or region_end values exceed the chromosome region, re-fetch the data.-----
            if (region_start_inText != region_start) or (region_end_inText
                                                         != region_end):
                print(
                    '   Note: recrawl because of discrepancies! region_start_inText: {}, region_start: {}, '
                    'region_end_inText: {}, region_end: {}'.format(
                        region_start_inText, region_start, region_end_inText,
                        region_end))
                # if direction is '+':
                #     new_start = region_start + 60
                #     new_end = region_end
                # else:
                #     new_start = region_start
                #     new_end = region_end - 60
                new_start = region_start + 15
                new_end = region_end - 15
                region_text, region_start, region_end = get_region_text(
                    chr_name, new_start, new_end, direction)
            return region_text, region_start, region_end
        except:
            exstr = traceback.format_exc()
            print(exstr)
    except:
        exstr = traceback.format_exc()
        print(exstr)


def reverse_complement(DNA_seq):
    complement_dict = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G",
        "N": "N",
        "a": "t",
        "t": "a",
        "g": "c",
        "c": "g",
    }
    sequence_list = list(DNA_seq)
    sequence_list = [complement_dict[base] for base in sequence_list]
    complement_seq = ''.join(sequence_list)
    reverse_complement_seq = complement_seq[::-1]
    return reverse_complement_seq


def translate(DNA_seq):
    gencode = {
        'ATA': 'I',
        'ATC': 'I',
        'ATT': 'I',
        'ATG': 'M',
        'ACA': 'T',
        'ACC': 'T',
        'ACG': 'T',
        'ACT': 'T',
        'AAC': 'N',
        'AAT': 'N',
        'AAA': 'K',
        'AAG': 'K',
        'AGC': 'S',
        'AGT': 'S',
        'AGA': 'R',
        'AGG': 'R',
        'CTA': 'L',
        'CTC': 'L',
        'CTG': 'L',
        'CTT': 'L',
        'CCA': 'P',
        'CCC': 'P',
        'CCG': 'P',
        'CCT': 'P',
        'CAC': 'H',
        'CAT': 'H',
        'CAA': 'Q',
        'CAG': 'Q',
        'CGA': 'R',
        'CGC': 'R',
        'CGG': 'R',
        'CGT': 'R',
        'GTA': 'V',
        'GTC': 'V',
        'GTG': 'V',
        'GTT': 'V',
        'GCA': 'A',
        'GCC': 'A',
        'GCG': 'A',
        'GCT': 'A',
        'GAC': 'D',
        'GAT': 'D',
        'GAA': 'E',
        'GAG': 'E',
        'GGA': 'G',
        'GGC': 'G',
        'GGG': 'G',
        'GGT': 'G',
        'TCA': 'S',
        'TCC': 'S',
        'TCG': 'S',
        'TCT': 'S',
        'TTC': 'F',
        'TTT': 'F',
        'TTA': 'L',
        'TTG': 'L',
        'TAC': 'Y',
        'TAT': 'Y',
        'TAA': '*',
        'TAG': '*',
        'TGC': 'C',
        'TGT': 'C',
        'TGA': '*',
        'TGG': 'W'
    }
    special_symbols = {
        'R': ['A', 'G'],
        'Y': ['T', 'C'],
        'M': ['A', 'C'],
        'K': ['G', 'T'],
        'S': ['G', 'C'],
        'W': ['A', 'T'],
        'H': ['A', 'T', 'C'],
        'B': ['T', 'G', 'C'],
        'V': ['A', 'G', 'C'],
        'D': ['A', 'G', 'T'],
        'N': ['A', 'G', 'T', 'C'],
    }
    aa_seq = ""
    for start in range(0, len(DNA_seq) - 2, nt_aa_ratio):
        stop = start + nt_aa_ratio
        codon = DNA_seq[start:stop]
        aa = ''
        aa_flag = False
        for base in codon:
            # The following code checks whether a codon contains any bases from special_symbols;
            # however, it only considers the case of having one such base, and does not account for multiple occurrences.
            if base in special_symbols:
                aa_flag = True
                possible_aa = []
                for real_base in special_symbols[base]:
                    real_codon = codon.replace(base, real_base)
                    aa = gencode.get(real_codon.upper(), 'X')
                    if aa not in possible_aa:
                        possible_aa.append(aa)
                if len(possible_aa) != 1:
                    aa = 'X'
                else:
                    for item in possible_aa:
                        aa = item
        if not aa_flag:
            aa = gencode.get(codon.upper(), 'X')
        aa_seq = aa_seq + aa
    return aa_seq


def get_possible_IFN_seqs(cpath, aa_region, region_start, region_end,
                          direction):
    IFN_infos = []
    aa_ofstopC = '\*'  # stopC = stop codon
    aa_ofstartC = 'M'  # startC = start codon
    possible_stopC_sites = sorted(
        [x.start() for x in re.finditer(aa_ofstopC, aa_region)])
    possible_startC_sites = sorted(
        [x.start() for x in re.finditer(aa_ofstartC, aa_region)])
    for startC_site in possible_startC_sites:
        for stopC_site in possible_stopC_sites:
            if startC_site < stopC_site:
                seq = aa_region[startC_site:stopC_site]
                len_seq = len(re.sub('[X]', '', seq))
                if ('*' not in seq) and (len_seq > min_IFN_len):
                    if direction is '+':
                        actual_start = region_start + startC_site * nt_aa_ratio
                        actual_end = region_start + stopC_site * nt_aa_ratio
                    else:
                        actual_start = region_end - stopC_site * nt_aa_ratio
                        actual_end = region_end - startC_site * nt_aa_ratio
                    IFN_info = [actual_start, actual_end, direction, seq]
                    IFN_infos.append(IFN_info)
    IFN_infos_path = f'{cpath}/possible_IFN_infos.fasta'
    with open(IFN_infos_path, 'w') as f:
        for IFN_info in IFN_infos:
            f.write('>' + 'start:' + f'{IFN_info[0]}' + ', ' + 'end:' +
                    f'{IFN_info[1]}' + ', ' + 'direction:' + f'{IFN_info[2]}' +
                    ', ' + '\n' + f'{IFN_info[3]}')
            f.write('\r\n')
    print('   Get IFN_seqs complete')


def pipeline(cpath, fpath, chr_name, elong_nt):
    fasta_seqs, consensus_seq, assembly_seq_start, assembly_seq_end, assembly_seq_dir \
        = get_consensus_seq(fpath)
    # k_value = get_k_value(fasta_seqs, consensus_seq)
    # dict_kmers_freq, list_kmers = get_kmers(fasta_seqs, k_value)
    # deBrujinGraph = deBrujinGraph_algorithm(list_kmers)
    # assembly_seq = get_assembly_by_DBG(dict_kmers_freq, list_kmers,
    #                                    deBrujinGraph)
    possible_aa_region, region_start, region_end = find_aa_region(
        cpath, chr_name, assembly_seq_start, assembly_seq_end,
        assembly_seq_dir, elong_nt)
    get_possible_IFN_seqs(cpath, possible_aa_region, region_start, region_end,
                          assembly_seq_dir)


def summarize_IFN_infos(opath):
    IFN_infos = []
    # with open(cpath) as f:
    #     conserv_site_file = json.load(f)
    species_names = get_filenames(opath)
    for species_name in species_names:
        species_path = f'{opath}/{species_name}'
        if os.path.isdir(species_path):
            # print('----------------------------------------')
            # print('{}:'.format(species_name))
            chromosome_names = get_filenames(species_path)
            for chromosome_name in chromosome_names:
                # chromosome_name: NC_051803.1_W
                # print(chromosome_name)
                chromosome_path = f'{species_path}/{chromosome_name}'
                if os.path.isdir(chromosome_path):
                    chr_name = chromosome_name.rsplit('_', 1)[0]
                    if chr_name.endswith(
                            ('C', 'W'
                             )):
                        chr_name = chr_name.rsplit('_', 1)[0]
                    centers = get_filenames(chromosome_path)
                    for center in centers:
                        center_path = f'{chromosome_path}/{center}'
                        fasta_names = get_filenames(center_path)
                        for fasta_name in fasta_names:
                            if 'possible_IFN' in fasta_name:
                                fasta_path = f'{center_path}/{fasta_name}'
                                with open(fasta_path) as f:
                                    f_lines = f.readlines()
                                for n, line in enumerate(f_lines):
                                    if line.startswith('>'):
                                        str_start = r'start:(\d*),'
                                        IFN_start = int(
                                            re.findall(str_start, line)[0])
                                        str_end = r'end:(\d*),'
                                        IFN_end = int(
                                            re.findall(str_end, line)[0])
                                        IFN_pos = round(
                                            (IFN_start + IFN_end) / 2)
                                        str_dir = r'direction:(.),'
                                        IFN_dir = re.findall(str_dir, line)[0]
                                        IFN_seq = f_lines[n + 1]
                                        if len(IFN_seq) > min_IFN_len:
                                            info_dict = {
                                                'species': species_name,
                                                'chr': chr_name,
                                                'start': IFN_start,
                                                'end': IFN_end,
                                                'pos': IFN_pos,
                                                'direction': IFN_dir,
                                                'seq': IFN_seq
                                            }
                                            if info_dict not in IFN_infos:
                                                IFN_infos.append(info_dict)
    return IFN_infos


def write_infos(opath, IFN_infos):
    with open(f'{opath}/possible_IFNs_summary.fasta', 'w') as f:
        for info in IFN_infos:
            f.write('>' + info['species'] + '_' + info['chr'] + '_POS' +
                    str(info['pos']) + ', ' + 'start:' + str(info['start']) +
                    ', ' + 'end:' + str(info['end']) + ', ' + 'direction:' +
                    info['direction'] + ', ' + '\n' + info['seq'] + '\n')


if __name__ == "__main__":
    # init_pipeline(out_path, elongation_nucleotides)
    possible_IFN_infos = summarize_IFN_infos(out_path)
    write_infos(out_path, possible_IFN_infos)
