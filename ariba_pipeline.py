# import pandas as pd
import os
import itertools
import pandas as pd
import argparse


def pairwise(iterable):
    a, b = itertools.tee(iterable)
    next(b, None)
    return zip(a, b)

# This function could be named better.
def get_reads_and_run_ariba(sra_accession: str, strain_id: int, threads: int = 25):
    if not os.path.exists(str(strain_id)):
        os.system('singularity exec singularity_images/sra-toolkit.simg prefetch {}'.format(sra_accession))
        os.system('singularity exec singularity_images/sra-toolkit.simg fastq-dump --split-files {} -O {}'.format(sra_accession, strain_id))
    else:
        print("Reads already retrieved and split....")
    print("Running ARIBA with {} threads".format(str(threads)))
    if not os.path.exists('{}/ariba_card_out'.format(strain_id)):
        os.system('singularity exec singularity_images/ariba.simg ariba run --threads {} out.card.prepareref {}/{} {}/{} {}/ariba_card_out'.format(str(threads), strain_id, sra_accession + '_1.fastq', strain_id, sra_accession + '_2.fastq', strain_id))
    if not os.path.exists('{}/ariba_resfinder_out'.format(strain_id)):
        print("Running ARIBA for Resfinders")
        os.system('singularity exec singularity_images/ariba.simg ariba run --threads {} resfinder.out {}/{} {}/{} {}/ariba_resfinder_out'.format(str(threads), strain_id, sra_accession + '_1.fastq', strain_id, sra_accession + '_2.fastq', strain_id))

parser = argparse.ArgumentParser()
# CSV with "SRA Accession # in one of the columns to allow for download of reads from NCBI
parser.add_argument('--input', help='input_csv', required=False)

# A folder containing paired reads from different organisms. There should be an even number of reads in this folder. 2 For every organism. The names of the reads should be such
# that if sorted they would be paired together in a manner like so when running the pairing after sorting [1, 2, 3, 4] -> (1,2) (3,4) The easiest way to do this would be to have the exact same
# name for paired reads but a _1 or _2 at the end so that sorting produces the correct output.

parser.add_argument('--reads', help='paired_reads_folder', required=False)
parser.add_argument('--genomes', help='genomes_folder', required=False)
args = parser.parse_args()
# Reading in the accessions and downloading them

if __name__ == '__main__':
    if args.input is None:
        df = pd.read_csv(args.input).head()
        report_paths = []
        for index, row in df.iterrows():
            sra_accession = row['SRA Accession #']
            strain_id = row['Strain']
            report_paths.append('{}/ariba_card_out/report.tsv'.format(strain_id))
            if ',' in sra_accession:
                [get_reads_and_run_ariba(sra_accession, strain_id) for sra_accession in sra_accession.replace(',', "").split(" ")]
            else:
                get_reads_and_run_ariba(sra_accession, strain_id)

        os.system('singularity exec singularity_images/ariba.simg ariba summary e_coli_card {}'.format(' '.join(report_paths)))

    if args.reads is None:
        files = os.listdir(args.reads)[0:2]
        files = sorted(files)
        paired_reads = list(pairwise(files))

        for read_one, read_two in paired_reads:
            if not os.path.exists(read_one[0:5]):
                os.system('singularity exec singularity_images/ariba.simg ariba run --threads 15 out.card.prepareref {} {} {}'.format(args.reads + '/' + read_one, args.reads + '/' + read_two, read_one[0:5]))
            else:
                print(read_one[0:5] + "exists already")
            if not os.path.exists(read_one[0:5] + 'resfinder'):
                os.system('singularity exec singularity_images/ariba.simg ariba run --threads 15 resfinder.out {} {} {}'.format(args.reads + '/' + read_one, args.reads + '/' + read_two, read_one[0:5] + 'resfinder'))
                
    genomes = os.listdir(args.genomes)[0:2]
    for genome in genomes:
        if not os.path.exists(genome):
            os.system('singularity exec singularity_images/insilicoseq-latest.simg iss generate --draft {}/{} --model miseq --cpus 25 --output {}'.format(args.genomes, genome,genome))
        else:
            print("Generated reads already exist")
        if not os.path.exists('{}/ariba_resfinder_out'.format(genome)):
            os.system('singularity exec singularity_images/ariba.simg ariba run --threads {} resfinder.out {}/{} {}/{} {}/ariba_resfinder_out'.format('25', genome, genome + '_R1.fastq', genome, genome + '_R2.fastq', genome))
        if not os.path.exists('{}/ariba_card_out'.format(genome)):
            os.system('singularity exec singularity_images/ariba.simg ariba run --threads {} out.card.prepareref {}/{} {}/{} {}/ariba_card_out'.format('25', genome, genome + '_R1.fastq', genome, genome + '_R2.fastq', genome))
