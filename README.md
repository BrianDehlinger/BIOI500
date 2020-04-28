# BIOI500 - AMRFinder, Resfinder, ARIBA comparison of Antibiotic Resistance Gene finding


BIOI 500 project - Project is indeed a pipeline but meant to compare three tools. As such, development speed was prioritized over reusability, modularity, and code readability. 

# Installation instructions - AMRFinderPlus
Only Linux environments are supported and singularity is used to handle dependencies. 

Install Miniconda(using your user not the root user)

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
curl https://rclone.org/install.sh | sudo bash -s beta
```

# PLEASE SETUP RCLONE using your own google drive information.
Data is available here: https://drive.google.com/open?id=1ilRmDnM2Dc5iZP174rWEDUUJvAtomSMs
Please mount this folder in your Google drive using rclone as Data:/Data or modify the Setup for the pipeline below
to use the correct paths for the corresponding folders.

Follow instructions here:
https://www.ostechnix.com/how-to-mount-google-drive-locally-as-virtual-file-system-in-linux/

# Run with SUDO - Singularity Image Setup and local directory configuration
(If this fails you will need to install singularity and manually pull each of the singularity images in the setup.py to the specified directory)
```bash
python3 setup.py
```

# Running the pipeline - Setup
```bash
mkdir google_drive_data
mkdir google_drive_data/data
mkdir Genome/
mkdir e_coli_genomes
rclone copy Data:/Data/data google_drive_data/data
rclone copy Data:/Data/ecoli_genomes e_coli_genomes
rclone copy Data:/Data/Genome Genome/
rclone copy Data:/Data/ncbi ncbi-blast-2.10.0+/
chmod +x ncbi-blast-2.10.0+/bin/blastn
rclone copy Data:Data/AA14 AA14/
rclone copy Data:Data/protein protein_dataset/
```
# Running the pipeline - Analysis

The first command runs ARIBA on both E. coli and S. aureus dataset utilizing the reads.
The second command runs Resfinder and AMRFinderPlus on the E. coli dataset from the goolgle drive utilizing the assemblies.
The third command runs Resfinder and AMRFinderPlus on the S. aureus dataset from the google drive utilizing the assemblies.
The analyze_amr command is meant to be run on the e_coli dataset results produced from the second command and summarizes resistance genes and phenotype predictions.
The analyze_amr_staph.py is meant to be run with the staph_aureus dataset and is run with only two SNPs and two different MIC thresholds.
```
python3 ariba_pipeline.py
python3 pipeline.py --output e_coli_amr_analysis --genomes e_coli_genomes --organism escherichia_coli
python3 pipeline.py --output staph_amr_analysis --genomes Genome/ --organism staphylococcus_aureus
python3 analyze_amr.py --genomes e_coli_amr_analysis
python3 analyze_amr_staph.py --genomes staph_amr_analysis --mic 8 -ariba_only_two
python3 analyze_amr_staph.py --genomes staph_amr_analysis --mic 64 -ariba_only_two
```

# Running the pipeline - Runtime test
 ```bash
python3 time_test.py
```
