#!/usr/bin/env python

import os
import sys
import time
import shutil
import subprocess
import argparse as ap
from pathlib import Path
import pandas as pd
from Bio import SeqIO
from Bio.Emboss.Applications import NeedleCommandline
from shutil import which

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
    p.add_argument('-ms', '--mature_mirnas_separate', type=str, default=None,
                   help='Path to a directory containing directories with mature miRNA sequences per one miRNA.'
                        'Each directory with mature mirna fastas should be names exactly as you name your miRNA '
                        'in other input files. Names of fasta files are not important, but will be used further as '
                        '"mature_mirna" column in result tables'
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


def make_reverse_compl_fasta(in_path, out_path):
    for record in SeqIO.parse(in_path, "fasta"):
        output = open(out_path, 'w+')
        output.write(">%s\n%s\n" % (record.id, record.seq.reverse_complement()))
        output.close()

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
    error_rate = '19'
    lower_length_bound = '11'
    index_cmd = ['triplexator', '-ss', args.mature_mirnas, '-ds', args.enhancers, '-p', str(args.nproc),
                 '--error-rate', error_rate, '--lower-length-bound', lower_length_bound, '-od', args.output]
    print(index_cmd)
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

    if args.mature_mirnas_separate:
        if not os.path.exists(args.mature_mirnas_separate):
            sys.exit('ERROR: ' + args.mature_mirnas_separate + 'path does not exist.\n')
        elif not os.listdir(args.mature_mirnas_separate):
            sys.exit('ERROR: ' + args.mature_mirnas_separate +
                     'path is empty.')
    else:
        sys.exit('ERROR: Please provide a valid path of directory with directories per miRNA with fasta files with '
                 'one sequence of mature miRNA per fasta (argument -ms or --mature_mirnas_separate).\n')


def run_seed_match_needle(args):
    print('Searching enhancers which contain exact match of provided seeds...')
    index_cmd = ['bash', 'src/get_enh_with_seeds.sh', args.genome, args.seeds_mirnas_forward, args.seeds_mirnas_reverse_compl,
                 args.output, args.enhancers_bed]
    print(' '.join(index_cmd))
    subprocess.run(index_cmd)
    print("Finished searching for enhancers with miRNAs' seeds sequences.")


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
def get_active_enh_regions(args):
    print('Extracting sequences of "active" regions of enhancers (the ones that align '
          'to miRNAs in seed region +- 14bp)...')
    index_cmd = ['bash', 'src/get_active_regions.sh', args.output, args.genome]
    print(' '.join(index_cmd))
    subprocess.run(index_cmd)
    print("Finished searching for enhancers with miRNAs' seeds sequences.")


def sort_fasta_by_seed(args):
    active_regions_path = os.path.join(args.output, 'enh_active/active_regions/fasta')
    out_dir = os.path.join(active_regions_path, 'fasta_sorted')

    seeds_dict = dict()

    i = 0
    with open(args.seeds_mirnas_forward) as seeds_forward:
        for seed in seeds_forward:
            mirna = seed.split()[0]
            print(mirna)
            mirna_seed_seq = seed.split()[1]
            seeds_dict[mirna + '_' + str(i)] = mirna_seed_seq
            i += 1

    i = 0
    with open(args.seeds_mirnas_reverse_compl) as seeds_rev_compl:
        for seed in seeds_rev_compl:
            mirna = seed.split()[0]
            print(mirna)
            mirna_seed_seq = seed.split()[1]
            seeds_dict[mirna + '_rev_compl_' + str(i)] = mirna_seed_seq
            i += 1

    if os.path.exists(out_dir) and os.path.isdir(out_dir):
        shutil.rmtree(out_dir)
    os.makedirs(out_dir)

    for filename in os.listdir(active_regions_path):
        file_path = os.path.join(active_regions_path, filename)
        if os.path.isfile(file_path):
            for record in SeqIO.parse(file_path, "fasta"):
                mirna_filename = filename.split('_')[0]
                print(mirna_filename)
                for mirna_name, mirna_seed in seeds_dict.items():
                    mirna_name_from_dict = mirna_name.split('_')[0].split('-')[0]
                    if mirna_name_from_dict == mirna_filename and mirna_seed in record.seq:
                        mir_name = "_".join(mirna_name.split('_')[:-1])
                        output = open(os.path.join(out_dir, mir_name), 'a')
                        output.write(">%s\n%s\n" % (record.id, record.seq))
                        output.close()


def align_needle(args):
    fasta_dir = os.path.join(args.output, 'enh_active/active_regions/fasta')
    fastas_sorted_dir = os.path.join(fasta_dir, 'fasta_sorted/')
    fasta_sorted_one_per_file_dir = os.path.join(fasta_dir, 'one_per_file')
    alignments_dir = os.path.join(args.output, 'alignments_needle')

    # if os.path.exists(fasta_sorted_one_per_file_dir) and os.path.isdir(fasta_sorted_one_per_file_dir):
    #     shutil.rmtree(fasta_sorted_one_per_file_dir)
    # os.makedirs(fasta_sorted_one_per_file_dir)

    print('Preparing all outputs as separate fasta files...')
    for fasta_sorted in os.listdir(fastas_sorted_dir):
        fasta_sorted_one_per_file = os.path.join(fasta_sorted_one_per_file_dir, fasta_sorted)
        fasta_sorted_path = os.path.join(fastas_sorted_dir, fasta_sorted)
        if os.path.exists(fasta_sorted_one_per_file) and os.path.isdir(fasta_sorted_one_per_file):
            shutil.rmtree(fasta_sorted_one_per_file)
        os.makedirs(fasta_sorted_one_per_file)
        for record in SeqIO.parse(fasta_sorted_path, "fasta"):
            new_file = os.path.join(fasta_sorted_one_per_file, fasta_sorted + '_' + record.id.replace(':', '_'))
            output = open(new_file, 'w')
            output.write(">%s\n%s\n" % (record.id, record.seq))
            output.close()

    print('getting reverse complementary sequences for the mature mirnas...')
    mature_seqs_rev_compl_path = os.path.join(args.output, 'mature_seqs_rev_compl')
    if os.path.exists(mature_seqs_rev_compl_path) and os.path.isdir(mature_seqs_rev_compl_path):
        shutil.rmtree(mature_seqs_rev_compl_path)
    os.makedirs(mature_seqs_rev_compl_path)

    for mir_dir in os.listdir(args.mature_mirnas_separate):
        rev_compl_dir_per_mir = os.path.join(mature_seqs_rev_compl_path, mir_dir)
        mature_seqs_path = os.path.join(args.mature_mirnas_separate, mir_dir)
        if os.path.exists(rev_compl_dir_per_mir) and \
                os.path.isdir(rev_compl_dir_per_mir):
            shutil.rmtree(rev_compl_dir_per_mir)
        os.makedirs(rev_compl_dir_per_mir)
        for file in os.listdir(mature_seqs_path):
            make_reverse_compl_fasta(in_path=os.path.join(mature_seqs_path, file),
                                     out_path=os.path.join(rev_compl_dir_per_mir, file))

    print('perform alignments each mirna vs each of possibly active regions')
    if os.path.exists(alignments_dir) and os.path.isdir(alignments_dir):
        shutil.rmtree(alignments_dir)
    os.makedirs(alignments_dir)

    mirna_names = [x for x in os.listdir(fasta_sorted_one_per_file_dir) if 'rev_compl' not in x]
    for mirna in mirna_names:
        print(mirna)
        mature_mirna_dir = os.path.join(args.mature_mirnas_separate, mirna)
        one_mir_mature_seqs_rev_compl_path = os.path.join(mature_seqs_rev_compl_path, mirna)
        print('mature_mirna_dir: ', mature_mirna_dir)
        active_regions_path = os.path.join(fasta_sorted_one_per_file_dir, mirna)
        print('active_regions_path: ', active_regions_path)
        out_alignments_per_mirna_dir = os.path.join(alignments_dir, mirna)
        print('out_alignments_per_mirna_dir: ', out_alignments_per_mirna_dir)

        if os.path.exists(out_alignments_per_mirna_dir) and os.path.isdir(out_alignments_per_mirna_dir):
            shutil.rmtree(out_alignments_per_mirna_dir)
        os.makedirs(out_alignments_per_mirna_dir)

        for active_region in os.listdir(active_regions_path):
            active_region_path = os.path.join(active_regions_path, active_region)
            for one_mature_mirna_file in os.listdir(mature_mirna_dir):
                out_alignment = os.path.join(out_alignments_per_mirna_dir, active_region + '_' +
                                             one_mature_mirna_file + "_needle.txt")
                needle_cline = NeedleCommandline(asequence=os.path.join(mature_mirna_dir, one_mature_mirna_file),
                                                 bsequence=active_region_path,
                                                 gapopen=50, gapextend=0.5,
                                                 datafile='EBLOSUM62',
                                                 outfile=out_alignment)
                stdout, stderr = needle_cline()

                out_alignment = os.path.join(out_alignments_per_mirna_dir, active_region + '_' +
                                             one_mature_mirna_file + "rev_compl_needle.txt")
                needle_cline = NeedleCommandline(asequence=os.path.join(one_mir_mature_seqs_rev_compl_path,
                                                                        one_mature_mirna_file),
                                                 bsequence=active_region_path,
                                                 gapopen=50, gapextend=0.5,
                                                 datafile='EBLOSUM62',
                                                 outfile=out_alignment)
                stdout, stderr = needle_cline()

    def get_num_ident_nucl(out_alignment):
        with open(out_alignment) as f:
            for i, line in enumerate(f):
                if i == 25:
                    num_ident_nucl = int(line.split('/')[0].split()[2])
        return num_ident_nucl

    def get_mirna_length(mirna_path):
        for record in SeqIO.parse(mirna_path, "fasta"):
            return len(record.seq)

    def get_percent_identity(alignment_path, mirna_path):
        res = float(get_num_ident_nucl(alignment_path)) / get_mirna_length(mirna_path)
        return round(res, 2)

    print('extract percent identities from outputs')
    pi_df = pd.DataFrame(columns=['mirna', 'mature_mirna', 'active region', 'PI'])
    for mirna in os.listdir(alignments_dir):
        mature_mirna_path = os.path.join(args.mature_mirnas_separate, mirna)
        alignments_per_mirna_path = os.path.join(alignments_dir, mirna)
        for active_region_aligned in os.listdir(alignments_per_mirna_path):
            active_region_chr_corr = ':'.join(active_region_aligned.replace('rev_compl_', '').split('_')[1:3])
            active_region_aligned_path = os.path.join(alignments_per_mirna_path, active_region_aligned)
            for one_mature_mir in os.listdir(mature_mirna_path):
                PI = get_percent_identity(active_region_aligned_path,
                                          os.path.join(mature_mirna_path, one_mature_mir))
                PI_new = pd.DataFrame({'mirna': [mirna],
                                       'mature_mirna': [one_mature_mir],
                                       'active region': [active_region_chr_corr],
                                       'PI': [PI]
                                       })
                pi_df = pi_df.append(PI_new, ignore_index=True)
    pi_df = pi_df.sort_values(by=['PI'], ascending=False)
    pi_df.to_csv(os.path.join(alignments_dir, 'stats.csv'),
                 sep='\t',
                 index=False)
    return pi_df


def get_best_aligned(args, needle_res):
    threshold = 0.5  # 12/22 seen an example in TargetScan but it gives 45 results; let's take less
    needle_res = needle_res[needle_res['PI'] >= threshold].reset_index()
    chr_coord_df = pd.DataFrame(needle_res['active region'].str.split(':', 1).tolist(),
                                columns=['chr', 'coord'])
    start_end_df = pd.DataFrame(chr_coord_df['coord'].str.split('-', 1).tolist(),
                                columns=['start', 'end'])
    needle_res['chr'] = chr_coord_df['chr']
    needle_res = pd.concat([needle_res, start_end_df], axis=1)[['chr', 'start', 'end', 'mirna', 'mature_mirna', 'PI']]
    needle_res = needle_res.drop_duplicates()
    needle_res[['chr', 'start', 'end']].to_csv(os.path.join(args.output, 'alignments_needle',
                                                            'alignment_best_regions.bed'),
                                               sep='\t',
                                               index=False,
                                               header=False)
    needle_res.to_csv(os.path.join(args.output, 'alignment_stats.tsv'), sep='\t', index=False)
    return needle_res


def get_enh_best_aligned(args):
    index_cmd = ['bash', 'src/get_enh_best_aligned.sh', args.output]
    print(' '.join(index_cmd))
    subprocess.run(index_cmd)


def prepare_final_result(alignment_stats, args):
    output_al_needle = os.path.join(args.output, 'alignments_needle/')
    enh_best_aligned_path = os.path.join(output_al_needle, 'alignment_best_enh.bed')
    enh_best_aligned = pd.read_csv(enh_best_aligned_path, '\t', header=None)
    enh_best_aligned.rename(columns={0: 'chr', 1: 'start', 2: 'end'}, inplace=True)
    enh_best_aligned['enhancer'] = enh_best_aligned['chr'].astype(str) + ':' + enh_best_aligned['start'].astype(
        str) + '-' + enh_best_aligned['end'].astype(str)
    enh_best_aligned = enh_best_aligned[[3, 4, 5, 'enhancer']]
    enh_best_aligned.rename(columns={3: 'chr', 4: 'start', 5: 'end'}, inplace=True)

    alignment_stats['chr'] = alignment_stats.chr.astype('str')
    alignment_stats['start'] = alignment_stats.start.astype('int64')
    alignment_stats['end'] = alignment_stats.end.astype('int64')
    enh_best_aligned['chr'] = enh_best_aligned.chr.astype('str')
    enh_best_aligned['start'] = enh_best_aligned.start.astype('int64')
    enh_best_aligned['end'] = enh_best_aligned.end.astype('int64')

    res = pd.merge(alignment_stats,
                   enh_best_aligned,
                   how='left',
                   on=['chr', 'start', 'end'])

    all_corrs = pd.read_csv(os.path.join(args.output, 'mir_enh_gene_trios.tsv'), '\t')
    res['mirna'] = res['mirna'].astype(str)
    res['enhancer'] = res['enhancer'].astype(str)
    all_corrs['mirna'] = all_corrs['mirna'].astype(str)
    all_corrs['enhancer'] = all_corrs['enhancer'].astype(str)

    res = pd.merge(res,
                   all_corrs,
                   how='left',
                   on=['mirna', 'enhancer'])
    res = res[['mirna', 'mature_mirna', 'Gene.Name', 'enhancer', 'corr(miRNA, gene)',
               'p.value adj', 'PI']].drop_duplicates()
    res.to_csv(os.path.join(args.output, 'mir_enh_gene_trios.tsv'), '\t', index=False)


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
        get_active_enh_regions(args)
        sort_fasta_by_seed(args)
        needle_res = align_needle(args)
        alignment_stats = get_best_aligned(args, needle_res)
        get_enh_best_aligned(args)
        prepare_final_result(alignment_stats, args)


if __name__ == '__main__':
    start_time = time.time()
    main()
    mins_elapsed = round((time.time() - start_time) / 60.0, 2)
    print('Finished ' + __file__ + ', ' + str(mins_elapsed) + ' minutes.')
