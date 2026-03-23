Source:
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/035/845/GCF_009035845.1_ASM903584v1/GCF_009035845.1_ASM903584v1_protein.faa.gz
gunzip GCF_009035845.1_ASM903584v1_protein.faa.gz
diamond makedb --in GCF_009035845.1_ASM903584v1_protein.faa --db GCF_009035845.1.protein.dmnd
