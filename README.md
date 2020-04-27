# BIOI500 - AMRFinder, Resfinder, ARIBA comparison of Antibiotic Resistance Gene finding


BIOI 500 project - Project is indeed a pipeline but meant to compare three tools. As such, development speed was prioritized over reusability, modularity, and code readability. 

# Installation instructions - AMRFinderPlus
Only Linux environments are supported

Install Miniconda(using your user)

(If you run into issues see: https://github.com/ncbi/amr/wiki/Install-with-bioconda)

```bash
cd ~
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash ./Miniconda3-latest-Linux-x86_64.sh 
source ~/.bashrc

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
source ~/miniconda3/bin/activate

conda install -y pandas
conda install -y -c bioconda ncbi-amrfinderplus

amrfinderplus -u
```

# Run with SUDO - Singularity Image Setup and local directory configuration
(If this fails you will need to install singularity and manually pull each of the singularity images in the setup.py to the specified directory)
```bash
python3 setup.py
```

# Running the pipeline - Two steps
The ARIBA pipeline can either be run on a folder containing multiple fasta genomes or on a folder containing paired sequencing reads.
You can also speficy a CSV with a column containing 'SRA Accession #' to download the reads using sra-toolkit and run the pipeline.

```bash
nohup python3 ariba_pipe.py --input csv_file_with_accession_number --genomes path_to_folder_with_fasta_genomes --reads path_to_folder_with_paired_reads &
```

To run Resfinder and AMRFinderPlus

```bash
nohup python3 pipeline.py --genomes genomes_folder_with_one_type_of_organism --output output_folder_name --organism escherichia_coli OR staphylococcus_aureus &
```
