# BacMLST
A python script to predict the MLST of bacteria isolates

# Introduction
MultiLocus Sequence Typing is a useful tool to analyse the genetic differences of bacteria isolates within a species, commonly, 7 housekeeping genes were selected to do this job. Traditionally, MLST is a very time/money consuming job, the researchers should amplified the DNA fragments by PCR, and using Sanger sequencing technology to know their sequences, then align them with known alleles. With the development of bioinformatics, some MLST databases were created, the most known is PubMLST, it allow the user upload their sequence data(DNA fragment sequenced by Sanger or whole genome sequence), the appearence of this kind database/tools extremely accelerated the speed of scientific research.


This tool could predict BacMLST using genome data. The database is downloaded from PubMLST, you can give it a species name such as "Streptococcus_suis", or you don't have to give it, It will predict the species by 16S rRNA sequence.


The reason I write this script is similar with the former tools BacSpecies, I'm a beginner of bioinformatics, started to study python three month ago, so, this kind of tools based on a very simple logic, it's very suitable for trainning.


If anyone want to try to use this tools and give me some advices, I will be very grateful! However, If you are seeking for a tool to solve some urgent problem or work, I recommend PubMLST (https://pubmlst.org), mlst in CGE (https://cge.food.dtu.dk/services/MLST/), or mlst developed by tseemann (https://github.com/tseemann/mlst) these tools are very mature, accurancy and faster than mine. However, my tool can also do the job commonly, just a little slow, and maybe have some problems I haven't found yet, if you give the species, almost 10 second for one genome, if not, it will cost 30 - 60 second.


# External Dependencies
BLAST+


# Usage
Put the python script, MLST_database and reference_database file into the folder contains your sequence file

unzip reference_database.zip and MLST_database before use
```
BacMLST [-i] [-o] [-r] [-t] [-v]
Input and Output:
  -i, --input             Input FASTA file
  -o, --output            Output directory
  -r, --refs              fasta file contain all MLST informations
Parameters:
  --outfmt                You can set the output format, we give two options, 'mix' and 'same'
  -t, --threads           Threads to use for BLAST searches
  -v, --version           Show version number and exit
```
# Quick usage
``` Python
python BacMLST.py -i *.fasta 
```
