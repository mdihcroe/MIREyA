# MIREyA (MIRnas functioning through Enhancer Activation)
the description is not complete yet


This is a pipeline for detection of miRNAs and their gene targets up-regulated through triggering their enhancer in the provided expression dataset (as in [Xiao et al., 2017](https://www.tandfonline.com/doi/full/10.1080/15476286.2015.1112487)). It runs in three different modes:
1. Searching for enhancers containing an exact match of the user-provided seed sequence of a miRNA, then expanding each seed by 14 bp of the corresponding mature miRNA and aligning it to the enhancer with Needle tool and keeping only enhancers with the percent identity (PI) > 0.5 (PI defined as a percent of matches between miRNA and DNA region). 
2. Scanning miRNA sequences against enhancer sequences and detecting potential target sites with MiRanda. 
3. Predicting RNA:DNA triplexes between miRNAs and enhancers with Triplexator tool. We relaxed error-rate and lower-length-bound Triplexator default parameters in order to adjust the algorithm to work with extremely short miRNA sequences (error-rate=19, lower-length-bound=11).

The approaches are interchangeable, also the user can merge the results of all approaches to reflect multiple mechanisms of potential miRNA:DNA binding. 


The 2d and 3d modes of the pipeline require the following tools to be pre-installed:
* [MiRanda](http://cbio.mskcc.org/miRNA2003/miranda.html)
* [Triplexator](http://abacus.qfab.org/tools/triplexator/inspector/install.html)

Besides, python>=3.5 and r-base are expected to be pre-installed. The pipeline was tested in Ubuntu and Ubuntu-based linux systems (Ubuntu>=16.04).

Python packages requirements are stored in requirements.txt file; for R -- r-requirements.txt

Here are the examples how to run MIREyA in each mode. Use the format of the example templates from data/ directory to create input files with your own data.
1. `python src/run_mireya.py -d seed_match_needle   -e data/enhancers.macrophages.Mtb.mm9.fasta -o out/seed_match_needle_out/ -ge data/DE_gene_expression.tsv -me data/DE_mirnas_expression.tsv -ei data/enh.gene.assoc.sign.tsv -m data/DE_mirna_mature_seqs.fa -g data/db -s data/seeds_seq_forward_short -sr data/seeds_seq_reverse_compl_short -eb data/enhancers.macrophages.Mtb.bed -ms data/mature_seqs/`
2. `python src/run_mireya.py -d miranda -e data/enhancers.macrophages.Mtb.mm9.fasta -o out/miranda_out/ -ge data/DE_gene_expression.tsv -me data/DE_mirnas_expression.tsv -ei data/enh.gene.assoc.sign.tsv -m data/DE_mirna_mature_seqs.fa`
3. `python src/run_mireya.py -d triplexator -e data/enhancers.macrophages.Mtb.mm9.fasta -o out/triplexator_out -ge data/DE_gene_expression.tsv -me data/DE_mirnas_expression.tsv -ei data/enh.gene.assoc.sign.tsv -m data/DE_mirna_mature_seqs.fa`

Run `python src/run_mireya.py --help` to see description of each parameter

In case you encounter any problem, feel free to contact me: elizarova@phystech.edu
