import os
import re
from multiprocessing import Pool

def get_files_with_suffixes(directory, suffixes):
    matching_files = []
    for root, _, files in os.walk(directory):
        for file in files:
            if any(file.endswith(suffix) for suffix in suffixes):
                matching_files.append(os.path.join(root, file))
    return matching_files


def run_prodigal(fasta, basename, outdir):
    """
    Function to predict ORFs using Prodigal.

    fasta:
    basename:
    outdir:
    """
    faa_file = os.path.join(outdir, basename + '.faa')
    cmd_para = [
                'prodigal', '-q',
                '-i', fasta,
                '-p', 'meta',
                '-a', faa_file,
                '-d', os.path.join(outdir, basename + '.ffn'),
                '-o', os.path.join(outdir, basename + '.gbk')
                ]
    cmd = ' '.join(cmd_para)
    if os.path.exists(faa_file):
        print("ORFs already predicted for bin:", basename)
    else:
        print("ORFs to be predicted for bin:", basename)
        try:
            os.system(cmd)
        except:
            print("Something wrong with prodigal annotation!")


def kegg_annotation(faa, basename, out_dir, db_dir, ko_dic, threads):
    """
    Function to perform KEGG annotation.
    The function invokes hmmsearch.

    faa:
    basename:
    out_dir:
    db_dir:
    ko_dic:
    threads:
    """
    print('KEGG annotation for {}'.format(basename))
    paras = []
    for knum, info in ko_dic.items():

        output = os.path.join(out_dir, knum + '.' + str(basename) + '.hmmout')
        hmm_db = os.path.join(db_dir, 'profiles', knum + '.hmm')

        if os.path.exists(output):
            continue

        if not os.path.exists(hmm_db):
            continue

        if info[1] == 'full':
            threshold_method = '-T'
            outtype = '--tblout'

        elif info[1] == 'domain':
            threshold_method = '--domT'
            outtype = '--domtblout'

        elif info[1] == 'custom':
            threshold_method = '-E'
            outtype = '--tblout'

        paras.append((threshold_method, info[0], outtype, output, hmm_db, faa))

    print("Number of KEGG processes to be performed:", str(len(paras)))

    process = Pool(threads)
    process.map(hmmsearch, paras)


def hmmsearch(paras):
    """
    Function to invoke hmmsearch software.
    """
    (threshold_method, threshold, outtype, output, hmm_db, faa) = paras
    cmd_para = [
        'hmmsearch',
        threshold_method, threshold,
        '--cpu', '1',
        '-o /dev/null',
        outtype, output,
        hmm_db,
        faa
    ]
    cmd = ' '.join(cmd_para)
    try:
        os.system(cmd)
    except:
        print("Something wrong with KEGG hmmsearch!")


def ko_list_parser(ko_list):
    """
    parse ko_list file into a dict object - based on DiTing

    :param ko_list: path of the file ko_list
    :return: a dictionary mapping knum to threshold and score_type
    :rtype: dict
    """
    ko_dic = {}  # { knum : [threshold, score_type] }
    with open(ko_list) as fi:
        next(fi)  # skip the first line (header)
        for line in fi:
            knum, threshold, score_type = line.split('\t')[0:3]
            if threshold == '-':
                continue
            else:
                ko_dic[knum] = [threshold, score_type]
    return ko_dic


# merge kegg annotations into one file
def merge_ko(hmmout_dir, output):
    with open(output, 'w') as fo:
        fo.write('#sample\tgene_id\tk_number\n')
    for hmmout_file in os.listdir(hmmout_dir):
        if hmmout_file.endswith('.hmmout'):
            kobasename = hmmout_file.rsplit('.', 1)[0]
            basename = kobasename.split('.', 1)[1]
            hmmout_file_path = os.path.join(hmmout_dir, hmmout_file)
            with open(hmmout_file_path, 'r') as fi:
                for line in fi:
                    if not line.startswith('#'):
                        gene_id, accession = line.split()[0:2]
                        lines = line.split()
                        if re.match(r'[0-9]+$', lines[2]):
                            k_number = lines[3]
                        else:
                            k_number = lines[2]
                        with open(output, 'a') as fo:
                            fo.write(basename + '\t' + gene_id + '\t' + k_number + '\n')
    #return ko_merged_dict


