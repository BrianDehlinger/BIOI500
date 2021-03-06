import pandas as pd
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--output', help='output directory', required = True)
parser.add_argument('--genomes', help='directory where genomes are located', required = True)
parser.add_argument('--organism', help='organism name', required = True)
parser.add_argument('-protein', action='store_true')
args = parser.parse_args()

# Again this could be made een more modular. However, this is good enough for a tool comparison

if __name__ == '__main__':
    output_directory = args.output
    genomes_directory = args.genomes
    organism = args.organism
    amr_finder_organism = 'Escherichia' if organism == 'escherichia_coli' else 'Staphylococcus_aureus'
    genome_file_paths = os.listdir(genomes_directory)
    for genome in genome_file_paths:
        if not os.path.exists(f'{output_directory}/{genome}'):
            if not os.path.exists(f'{output_directory}'):
                os.mkdir(f'{output_directory}')
            print("Making directory for organism")
            os.mkdir(f'{output_directory}/{genome}/')
            os.mkdir(f'{output_directory}/{genome}/resfinder')
            os.mkdir(f'{output_directory}/{genome}/amrfinder')
        os.system(f'singularity exec singularity_images/pointfinder.simg pointfinder.py -i {genomes_directory}/{genome} -o {output_directory}/{genome}/resfinder -p pointfinder_db -s {organism} -m blastn -m_p ncbi-blast-2.10.0+/bin/blastn')
        os.system(f'singularity exec singularity_images/resfinder-local.simg resfinder.py -i {genomes_directory}/{genome} -o {output_directory}/{genome}/resfinder -p resfinder_db')
        os.system(f'amrfinder -n {genomes_directory}/{genome} -O {amr_finder_organism} -o {output_directory}/{genome}/amrfinder/amrfinderoutput.txt --threads 17')
        if args.protein:
            os.system(f'amrfinder -p {genomes_directory}/{genome}/*.gpff -O {amr_finder_organism} -o {output_directory}/{genome}/amrfinder/amrfinderoutput.txt --threads 17')



