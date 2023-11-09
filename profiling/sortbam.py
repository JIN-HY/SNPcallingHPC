import os
from pathlib import Path
import glob
import re
import gzip
from optparse import OptionParser
import numpy as np
import pandas as pd

def slurm_head(i): 
    return f"""#!/bin/bash
#SBATCH --job-name=sort{i}
#SBATCH --ntasks-per-node=4
#SBATCH --error=sort{i}.err
#SBATCH --output=sort{i}.out
#SBATCH --mem=30G
#SBATCH --time=24:00:00

ml bwa samtools
"""
cmds = []

parser = OptionParser()
parser.add_option('--indir', default="mergedbam/",help="input directory")
parser.add_option('--outdir', default="sortedbam/", help="output directory")
#parser.add_option('--pattern1', default="pairedR1_001.fastq.gz", help="filename pattern for forward")
#parser.add_option('--pattern2', default="pairedR2_001.fastq.gz", help="filename pattern for forward")
#parser.add_option('--pattern1', default="trim.R1.fastq.gz", help="filename pattern for forward")
parser.add_option('--info')
parser.add_option('--metadf', default="fq2sample.tsv", help="file name to sample name")
parser.add_option('--n', default=10, help="how many command in a bash")
opts, args = parser.parse_args()
indir = opts.indir
outdir = opts.outdir
#pattern1 = opts.pattern1
#pattern2 = opts.pattern2
metadf = opts.metadf
group = opts.n
info = opts.info

#ref = 'GCA_000003195.3_Sorghum_bicolor_NCBIv3_genomic.fna'

df = pd.read_csv(metadf, sep = "\t")
df['seqname'] = [x.split(";")[0] for x in df.submitted_ftp]
df['seqname'] = [os.path.split(x)[1][:-15] for x in df.seqname]
df.seqname = df.seqname + "."
df['PI'] = [x.split("-")[2] for x in df.run_alias]

def sortcommand(PI):
    return f'samtools sort -o {outdir}/{PI}.sorted.bam {indir}/{PI}.bam\n\nsamtools index {outdir}/{PI}.sorted.bam\n'


#fpes = glob.glob(f'{indir}/*bam', recursive = True)

n = 0

i = 1
fslurm = open("sort"+str(i)+".slurm", "w")
fslurm.write(slurm_head(str(i)))

PIs = set(df.PI)
print(len(PIs))
for PI in PIs:
    cmd = sortcommand(PI)
    fslurm.write(f"{cmd}\n")
    n += 1
    if n > group:
        fslurm.close()
        n = 0
        i += 1
        fslurm = open("sort"+str(i)+".slurm", "w")
        fslurm.write(slurm_head(str(i)))
    #fslurm.close()
