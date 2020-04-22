# import pandas as pd
import os
import itertools
import pandas as pd


# NOTE TO SELF for the future
# I emphasized development speed over reuse here. The code quality is more something meant for throwaway code I figure this won't be done many more times. If needed it can easily be refactored into more modular components as part of a software solution.
# Ideally, this would also include a command line interface and argparsing. It might be better to write functionally using functions like reduce, the functions from functools and itertools to allow massive scalability for a pipeline.
# GNU parallel might be used to make some of the steps go faster but it depends on the step. IO bound tasks run in parallel could benefit from threads whereas CPU bound tasks could use multiple actual cores for a speed up.
# Downloads need to be done with resepect for the bandwith of this machine and the bandwith of NCBI.


# This function could be named better.
def get_reads_from_sra_accession(sra_accession: str, strain_id: int, threads: int = 25):
    if not os.path.exists(str(strain_id)):
        os.system('singularity exec sra-toolkit.simg prefetch {}'.format(sra_accession))
        os.system('singularity exec sra-toolkit.simg fastq-dump --split-files {} -O {}'.format(sra_accession, strain_id))
    else:
        print("Reads already retrieved and split....")
    print("Running ARIBA with {} threads".format(str(threads)))
    if not os.path.exists('{}/ariba_card_out'.format(strain_id)):
        os.system('singularity exec ariba.simg ariba run --threads {} out.card.prepareref {}/{} {}/{} {}/ariba_card_out'.format(str(threads), strain_id, sra_accession + '_1.fastq', strain_id, sra_accession + '_2.fastq', strain_id))
    if not os.path.exists('{}/ariba_resfinder_out'.format(strain_id)):
        print("Running ARIBA for Resfinders")
        os.system('singularity exec ariba.simg ariba run --threads {} resfinder.out {}/{} {}/{} {}/ariba_resfinder_out'.format(str(threads), strain_id, sra_accession + '_1.fastq', strain_id, sra_accession + '_2.fastq', strain_id))

# Reading in the accessions and downloading them
df = pd.read_csv('e_coli/accessions.csv')
report_paths = []

for index, row in df.iterrows():
    sra_accession = row['SRA Accession #']
    strain_id = row['Strain']
    report_paths.append('{}/ariba_card_out/report.tsv'.format(strain_id))
    if ',' in sra_accession:
        [get_reads_from_sra_accession(sra_accession, strain_id) for sra_accession in sra_accession.replace(',', "").split(" ")]
    else:
        get_reads_from_sra_accession(sra_accession, strain_id)

os.system('singularity exec ariba.simg ariba summary e_coli_card {}'.format(' '.join(report_paths)))

def pairwise(iterable):
    a, b = itertools.tee(iterable)
    next(b, None)
    return zip(a, b)

files = os.listdir("google_drive_data/data")
files.remove("Genome")
files = sorted(files)
paired_reads = list(pairwise(files))

for read_one, read_two in paired_reads:
    if not os.path.exists(read_one[0:5]):
        os.system('singularity exec ariba.simg ariba run --threads 15 out.card.prepareref {} {} {}'.format("google_drive_data/data/" + read_one, "google_drive_data/data/" + read_two, read_one[0:5]))
    else:
        print(read_one[0:5] + "exists already")
    if not os.path.exists(read_one[0:5] + 'resfinder'):
        os.system('singularity exec ariba.simg ariba run --threads 15 resfinder.out {} {} {}'.format("google_drive_data/data/" + read_one, "google_drive_data/data/" + read_two, read_one[0:5] + 'resfinder'))

genomes_without_reads = ['AA67', 'AA77', 'AA31', 'AA39', 'AA99', 'AA18', 'AA91', 'AA76', 'AA22', 'AA94', 'AA93', 'AA27', 'AA101', 'AA80', 'AA33', 'AA61', 'AA79', 'AA46', 'AA62', 'AA52', 'AA60', 'AA92', 'AA87', 'AA14', 'AA104', 'AA55', 'AA23', 'AA103', 'AA95']

for genome in genomes_without_reads:
    if not os.path.exists(genome):
        os.system('singularity exec insilicoseq-latest.simg iss generate --draft google_drive_data/data/Genome/{} --model miseq --cpus 25 --output {}'.format(genome,genome))
    else:
        print("Generated reads already exist")
    if not os.path.exists('{}/ariba_resfinder_out'.format(genome)):
        os.system('singularity exec ariba.simg ariba run --threads {} resfinder.out {}/{} {}/{} {}/ariba_resfinder_out'.format('25', genome, genome + '_R1.fastq', genome, genome + '_R2.fastq', genome))
    if not os.path.exists('{}/ariba_card_out'.format(genome)):
        os.system('singularity exec ariba.simg ariba run --threads {} out.card.prepareref {}/{} {}/{} {}/ariba_card_out'.format('25', genome, genome + '_R1.fastq', genome, genome + '_R2.fastq', genome))

