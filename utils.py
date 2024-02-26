import os
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

