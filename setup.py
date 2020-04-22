import os

# I know that using os.system is sloppy. It was extremely fast to write. Ideally, this can be made to use subprocess.Popen. However, this is in fact more readable.

os.system("mkdir singularity_images")
os.chdir("singularity_images")
os.system("singularity pull docker://briandehlinger/resfinder:local")
os.system("singularity pull docker://hadrieng/insilicoseq:latest")
os.system("singularity pull docker://sangerpathogens/ariba")
os.system("singularity pull docker://inutano/sra-toolkit")
os.system("singularity pull docker://briandehlinger/pointfinder")
