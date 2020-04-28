import os

# I know that using os.system is sloppy. It was extremely fast to write. Ideally, this can be made to use subprocess.Popen. However, this is in fact more readable.
os.system("sudo wget -O- http://neuro.debian.net/lists/xenial.us-ca.full | sudo tee /etc/apt/sources.list.d/neurodebian")
os.system("sudo apt-key adv --recv-keys --keyserver hkp://pool.sks-keyservers.net:80 0xA5D32F012649A5A9")
os.system("sudo apt-get update")
os.system("sudo apt-get install -y singularity-container")
os.system("source ~/.bashrc")
os.system("mkdir singularity_images")
os.chdir("singularity_images")
os.system("singularity pull docker://briandehlinger/resfinder:local")
os.system("singularity pull docker://hadrieng/insilicoseq:latest")
os.system("singularity pull docker://sangerpathogens/ariba")
os.system("singularity pull docker://inutano/sra-toolkit")
os.system("singularity pull docker://briandehlinger/pointfinder")
os.system("mkdir google_drive_data")
os.system("mkdir google_drive_data/data")
