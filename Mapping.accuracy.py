'''
Simple mapping accuracy evaluation of two SAM/BAM files
This script was modified from https://github.com/milkschen/leviosam2
'''
import argparse
import pysam
import sys

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-g', '--gold', required=True,
        help='Path to the truth SAM/BAM file. [required]'
    )
    parser.add_argument(
        '-q', '--query', required=True,
        help='Path to the query SAM/BAM file. [required]'
    )
    parser.add_argument(
        '-o', '--output', default='',
        help='Path to the output report. Leave empty to write to stdout. [empty]'
    )
    parser.add_argument(
        '-t', '--tolerance', default=10,
        help='Allowed bases of position difference (diff > `-t` is consider incorrect). [10]'
    )
    parser.add_argument(
        '-m', '--high_mapq', default=10,
        help='Lower bound MAPQ for an alignment to be considered as high-quality. [10]'
    )
    args = parser.parse_args()
    return args


def evaluate_mapping_correctness(fn_gold, fn_query, f_out, tolerance, high_mapq):
    gold = {}
    seqs = []
    with pysam.AlignmentFile(fn_gold, 'r') as f:
        for r in f:
            if r.is_read1:
                gold[r.query_name + '_1'] = (r.reference_name, r.reference_start, r.reference_end)
                seqs.append(r.query_name + '_1')
            elif r.is_read2:
                gold[r.query_name + '_2'] = (r.reference_name, r.reference_start, r.reference_end)
                seqs.append(r.query_name + '_2')

    comp = {}
    with pysam.AlignmentFile(fn_query, 'r') as f:
        for r in f:
            if r.is_read1:
                comp[r.query_name + '_1'] = (r.reference_name, r.reference_start, r.reference_end)
            elif r.is_read2:
                comp[r.query_name + '_2'] = (r.reference_name, r.reference_start, r.reference_end)


    cnt_same = 0
    cnt_fp = 0
    cnt_fn = 0
    cnt_incorrect = 0
    cnt_mapped = 0
    for seq in seqs:
        gr = gold[seq]
        cr = comp[seq]
        if gr[0] == None and cr[0] == None:
            cnt_same += 1
        elif gr[0] == None and cr[0] != None:
            cnt_fp += 1
            cnt_mapped += 1
        elif gr[0] != None and cr[0] == None:
            cnt_fn += 1
        else:
            cnt_mapped += 1
            if abs(gr[1] - cr[1]) <= tolerance:
                cnt_same += 1
            else:
                cnt_incorrect += 1

    print(f'total_reads\t{len(gold)}', file=f_out)
    print(f'n_mapped\t{cnt_mapped}', file=f_out)
    print(f'fp_mapped\t{cnt_fp}', file=f_out)
    print(f'fn_mapped\t{cnt_fn}', file=f_out)
    print(f'num_incorrectedMap\t{cnt_incorrect}', file=f_out)
    print(f'num_correctedMap\t{cnt_same}', file=f_out)

    return


if __name__ == '__main__':
    args = parse_args()
    fn_gold = args.gold
    fn_query = args.query
    tolerance = int(args.tolerance)
    high_mapq = args.high_mapq

    if args.output != '':
        f_out = open(args.output, 'w')
    else:
        f_out = sys.stdout

    evaluate_mapping_correctness(
        fn_gold=fn_gold, fn_query=fn_query, f_out=f_out,
        tolerance=tolerance, high_mapq=high_mapq
    )
