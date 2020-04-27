import os

import pandas as pd
import os
import argparse
import resource
import timeit


def timer_resfinder():
    os.system(f'singularity exec singularity_images/pointfinder.simg pointfinder.py -i Genome/AA14 -o test/resfinder -p pointfinder_db -s staphylococcus_aureus -m blastn -m_p ncbi-blast-2.10.0+/bin/blastn')
    os.system(f'singularity exec singularity_images/resfinder-local.simg resfinder.py -i Genome/AA14 -o test/resfinder -p resfinder_db')
    return

def time_test_resfinder():
    os.mkdir('test')
    os.mkdir('test/resfinder')
    time = timeit.timeit(timer_resfinder, number=5)
    os.system('rm -rf test')
    return time, '1'

def amr_execute():
    os.system(f'amrfinder -n Genome/AA14 -O Staphylococcus_aureus --threads 17')

def time_test_amrfinder():
    times = []
    avg_time = timeit.timeit(amr_execute, number=5)
    print(avg_time)
    return avg_time, '1'

def time_test_ariba():
    avg_time = timeit.timeit(time_ariba, number=1)
    print(avg_time)
    return avg_time, '1'

def time_ariba():
    os.system(f'singularity exec singularity_images/ariba.simg ariba run resfinder.out --threads 17 AA14/AA14_R1.fastq AA14/AA14_R2.fastq AA14/ariba_resfinder_out')
    os.system(f'singularity exec singularity_images/ariba.simg ariba run staph_snps_db --threads 17 AA14/AA14_R1.fastq AA14/AA14_R2.fastq AA14/ariba_snp')
    os.system('rm -rf AA14/ariba_resfinder_out')
    os.system('rm -rf AA14/ariba_snp')

if __name__ == '__main__':
    organism = 'staphylococcus_aureus'
    amr_finder_organism = 'Escherichia' if organism == 'escherichia_coli' else 'Staphylococcus_aureus'
    avg_time_resfinder = time_test_resfinder()
    avg_time_amr = time_test_amrfinder()
    avg_time_ariba = time_test_ariba()
    print('Resfinder:', avg_time_resfinder)
    print('AMRFinderPlus:', avg_time_amr)
    print('ARIBA:', avg_time_ariba)
