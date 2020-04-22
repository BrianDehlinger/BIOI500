import os

os.system("mkdir singularity_images")
os.chdir("singularity_images")
os.system("singularity pull docker://briandehlinger/resfinder:local")
os.system("singularity pull docker://hadrieng/insilicoseq:latest")
os.system("singularity pull docker://sangerpathogens/ariba")
os.system("singularity pull docker://inutano/sra-toolkit")
os.system("singularity pull docker://briandehlinger/pointfinder")
