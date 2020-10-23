#!/usr/bin/env python

import os
import sys
import time
import subprocess
import argparse as ap
from pathlib import Path
from shutil import which
import pandas as pd

__author__ = 'Anna Elizarova'
__version__ = '1.1'
__date__ = '23 October 2020'

detection_mir_enh_interaction_methods = ['seed_match_needle', 'miranda', 'triplexator']
MAX_PROC=3


def read_params():
    p = ap.ArgumentParser(description="")
    p.add_argument('-d', '--detection_mir_enh_interaction', type=str, default=None,
                   help='Approach to detect miRNA:enhancer interaction: seed_match_needle / miranda / triplexator')
    p.add_argument('-e', '--enhancers', type=str, default=None,
                   help='Path to a fasta file with sequences of enhancers of interest')
    p.add_argument('-m', '--mature_mirnas', type=str, default=None,
                   help='Path to a fasta file with sequences of mature miRNAs of interest (needed only if '
                        'the selected approach is miranda or triplexator')
    p.add_argument('-s', '--seeds_mirnas_forward', type=str, default=None,
                   help='Path to a .tsv file with sequences of seeds of miRNAs of interest (needed only if '
                        'the selected approach is seed_match_needle')
    p.add_argument('-sr', '--seeds_mirnas_reverse_compl', type=str, default=None,
                   help='Path to a .tsv file with reverse complementary sequences of seeds of miRNAs '
                        'of interest (needed only if the selected approach is seed_match_needle')
    p.add_argument('-g', '--genome', type=str, default=None,
                   help='Path to a folder with fasta files with complete genome of the organism of interest: '
                        'one chromosome per file (needed only if the selected approach is seed_match_needle')
    p.add_argument('-eb', '--enhancers_bed', type=str, default=None,
                   help='Path to a bed file with coordinates of enhancers of interest '
                        '(needed only if the selected approach is seed_match_needle')
    p.add_argument('-ge', '--gene_expression', type=str, default=None,
                   help='Path to a .tsv file with gene expression')
    p.add_argument('-me', '--mirnas_expression', type=str, default=None,
                   help='Path to a .tsv file with expression of miRNAs of interest')
    p.add_argument('-ei', '--enh_gene_interaction', type=str, default=None,
                   help='Path to a .tsv file with enhancers and corresponding genes you believe they regulate')
    p.add_argument('-o', '--output', type=str, default=None,
                   help='Full path to output directory')
    p.add_argument('--nproc', type=int, default=MAX_PROC,
                   help='Maximum number of processors to use. Default is 3 or a lower number of available processors.')
    return p.parse_args()


def check_main_args(args):
    if args.detection_mir_enh_interaction:
        if args.detection_mir_enh_interaction not in detection_mir_enh_interaction_methods:
            sys.exit('ERROR: Please select valid method to detect miRNA:enhancer interaction: '
                     'seed_match_needle / miranda / triplexator '
                     '(argument -d or --detection_mir_enh_interaction).\n')
    else:
        sys.exit('ERROR: Please select method to detect miRNA:enhancer interaction '
                 '(argument -d or --detection_mir_enh_interaction).\n')
    if args.enhancers:
        if not os.path.exists(args.enhancers):
            sys.exit('ERROR: Fasta file with enhancers (' + args.enhancers + ') not found\n')
    else:
        sys.exit('ERROR: Please provide a valid fasta file with enhancers (argument -e or --enhancers).\n')

    if args.output:
        output_parent_dir = Path(args.output).parent
        if not os.path.exists(str(output_parent_dir)):
            sys.exit('ERROR:' + output_parent_dir + ' directory not found\n')
    else:
        sys.exit('ERROR: Please provide a valid full path to output directory (argument -o or --output).\n')

    if args.gene_expression:
        if not os.path.exists(args.gene_expression):
            sys.exit('ERROR: .tsv file with gene expression (' + args.gene_expression + ') not found\n')
    else:
        sys.exit('ERROR: Please provide a valid .tsv file with gene expression (argument -ge or --gene_expression).\n')

    if args.mirnas_expression:
        if not os.path.exists(args.mirnas_expression):
            sys.exit('ERROR: .tsv file with expression of miRNAs of interest '
                     '(' + args.mirnas_expression + ') not found\n')
    else:
        sys.exit('ERROR: Please provide a valid .tsv file with expression of miRNAs of interest'
                 ' (argument -me or --mirnas_expression).\n')

    if args.enh_gene_interaction:
        if not os.path.exists(args.enh_gene_interaction):
            sys.exit('ERROR: .tsv file with enhancer:gene interaction (' +
                     args.enh_gene_interaction + ') not found\n')
    else:
        sys.exit('ERROR: Please provide a valid .tsv file with enhancer:gene interaction'
                 ' (argument -ei or --enh_gene_interaction).\n')


def tool_installed(name):
    """Check whether `name` is on PATH and marked as executable."""
    if which(name) is None:
        print('\nERROR: ' + name + ' was not found in PATH and marked as executable. '
                                   'Please, install it as described in README.\n')
        sys.exit()


# ------------------------------------------------------------------------------
#   STEP 1
# ------------------------------------------------------------------------------

def check_mature_mirnas(args):
    if args.mature_mirnas:
        if not os.path.exists(args.mature_mirnas):
            sys.exit('ERROR: Fasta file with mature miRNAs (' + args.mature_mirnas + ') not found\n')
    else:
        sys.exit('ERROR: Please provide a valid fasta file with mature miRNAs (argument -m or --mature_mirnas).\n')


def run_miranda(args):
    index_cmd = ['miranda', args.mature_mirnas, args.enhancers]
    print(' '.join(index_cmd))
    miranda_output_file = os.path.join(args.output, 'miranda.out')
    with open(miranda_output_file, "w") as f:
        subprocess.run(index_cmd, stdout=f)
    print('Miranda finished its work.')


def parse_miranda_output(args):
    print('Parsing miranda output...')
    target_text = '>>'
    file_path = os.path.join(args.output, 'miranda.out')
    res_df = pd.DataFrame()
    with open(file_path) as lines:
        for line in lines:
            if target_text in line:
                # TODO: use full names
                mirna = 'Mir' + ''.join(filter(lambda x: x.isdigit(), line.split()[0].split('-')[2]))
                enh_chr = line.split()[1].split(':')[0]
                enh_start = line.split()[1].split(':')[1].split('-')[0]
                enh_end = line.split()[1].split(':')[1].split('-')[1]
                max_score = line.split()[4]
                max_energy = line.split()[5]
                enh_align_start_pos = line.split()[9]
                res_new = pd.DataFrame({'mirna': [mirna],
                                        'enh.chr': [enh_chr],
                                        'enh.start': [enh_start],
                                        'enh.end': [enh_end],
                                        'max.score': [max_score],
                                        'max.energy': [max_energy],
                                        'enh.align.start.pos': [enh_align_start_pos]
                                        })
                res_df = res_df.append(res_new, ignore_index=True)
    res_df.to_csv(os.path.join(args.output, 'miranda.enh.tsv'),
                  sep='\t',
                  index=False)


def run_triplexator(args):
    error_rate = 19
    lower_length_bound = 11
    index_cmd = ['triplexator', '-ss', args.mature_mirnas, '-ds', args.enhancers, '-p', args.nproc,
                 '--error-rate', error_rate, '--lower-length-bound', lower_length_bound, '-od', args.output]
    print(' '.join(index_cmd))
    subprocess.run(index_cmd)
    index_cmd = ['sed', '-i', '"s/GT (rel)/GT (rel)	/g"', os.path.join(args.output, "triplex_search.summary")]
    subprocess.run(index_cmd)
    print('Triplexator finished its work.')


def check_seed_match_needle_specific_args(args):
    if args.genome:
        if not os.path.exists(args.genome):
            sys.exit('ERROR: ' + args.genome + 'path does not exist.\n')
        elif not os.listdir(args.genome):
            sys.exit('ERROR: ' + args.genome + 'path is empty. Files with fasta per chromosome of genome of interest'
                                               'are expected.\n')
    else:
        sys.exit('ERROR: Please provide a valid path of directory with fasta files with complete genome of the '
                 'organism of interest: one chromosome per file (argument -g or --genome).\n')

    if args.seeds_mirnas_forward:
        if not os.path.exists(args.seeds_mirnas_forward):
            sys.exit('ERROR: .tsv file with mirna seed forward sequences (' + args.seeds_mirnas_forward + ') '
                                                                                                          'not found\n')
    else:
        sys.exit('ERROR: Please provide a valid .tsv file with mirna seed forward sequences'
                 ' (argument -s or --seeds_mirnas_forward).\n')

    if args.seeds_mirnas_reverse_compl:
        if not os.path.exists(args.seeds_mirnas_reverse_compl):
            sys.exit('ERROR: .tsv file with mirna seed reverse complementary sequences (' +
                     args.seeds_mirnas_reverse_compl + ') not found\n')
    else:
        sys.exit('ERROR: Please provide a valid .tsv file with mirna seed reverse complementary sequences'
                 ' (argument -sr or --seeds_mirnas_reverse_compl).\n')

    if args.enhancers_bed:
        if not os.path.exists(args.enhancers_bed):
            sys.exit('ERROR: .bed file with enhancer coordinates (' +
                     args.enhancers_bed + ') not found\n')
    else:
        sys.exit('ERROR: Please provide a valid .bed file with enhancer coordinates'
                 ' (argument -eb or --enhancers_bed).\n')


def run_seed_match_needle(args):
    print('Searching enhancers which contain exact match of provided seeds...')
    index_cmd = ['bash', 'src/get_enh_with_seeds.sh', args.genome, args.seeds_mirnas_forward, args.seeds_mirnas_reverse_compl,
                 args.output, args.enhancers_bed]
    print(' '.join(index_cmd))
    subprocess.run(index_cmd)
    print("Finished searching for enhancers with miRNAs' seeds sequences.")
    pass


# ------------------------------------------------------------------------------
#   STEP 2
# ------------------------------------------------------------------------------
def calc_corr(args):
    out_file_name = 'mir_enh_gene_trios.tsv'
    index_cmd = ['Rscript', 'src/calc_corr.R', args.detection_mir_enh_interaction, args.output, args.gene_expression,
                 args.mirnas_expression, args.enh_gene_interaction]
    print(' '.join(index_cmd))
    subprocess.run(index_cmd)
    print('Calculation is finished. Please find the results in ' + out_file_name)

# ------------------------------------------------------------------------------
#   STEP 3
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
#   MAIN
# ------------------------------------------------------------------------------

def main():
    if not sys.version_info.major == 3:
        sys.stderr.write('ERROR: Python version: ' + sys.version)
        sys.exit('ERROR: This software uses Python3, please update Python')

    args = read_params()
    check_main_args(args)

    print('The method you selected to detect miRNA:enhancer interaction: ' + args.detection_mir_enh_interaction)
    if not os.path.exists(args.output):
        os.makedirs(args.output)
    print('Output directory: ' + args.output)

    print('\nSTEP 1. Detect miRNA:enhancer interaction...')

    if args.detection_mir_enh_interaction == 'seed_match_needle':
        check_seed_match_needle_specific_args(args)
        tool_installed('bedtools')
        run_seed_match_needle(args)
    elif args.detection_mir_enh_interaction == 'miranda':
        check_mature_mirnas(args)
        tool_installed('miranda')
        run_miranda(args)
        parse_miranda_output(args)
    elif args.detection_mir_enh_interaction == 'triplexator':
        check_mature_mirnas(args)
        tool_installed('triplexator')
        run_triplexator(args)
    else:
        sys.exit('ERROR: Please select valid method to detect miRNA:enhancer interaction: '
                 'seed_match_needle / miranda / triplexator '
                 '(argument -d or --detection_mir_enh_interaction).\n')

    print('\nSTEP 2. Calculating correlations of miRNAs with genes associated with enhancers '
          'that they can potentially bind......')
    calc_corr(args)

    if args.detection_mir_enh_interaction == 'seed_match_needle':
        print('Aligning the rest of mature miRNA to the enhancers which contain exact seed...')


if __name__ == '__main__':
    start_time = time.time()
    main()
    mins_elapsed = round((time.time() - start_time) / 60.0, 2)
    print('Finished ' + __file__ + ', ' + str(mins_elapsed) + ' minutes.')
