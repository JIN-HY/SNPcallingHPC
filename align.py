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
#SBATCH --job-name=SAPalign
#SBATCH --ntasks-per-node=4
#SBATCH --error=align{i}.err
#SBATCH --output=align{i}.out
#SBATCH --mem=30G
#SBATCH --time=24:00:00

ml bwa samtools
"""
cmds = []

parser = OptionParser()
parser.add_option('--indir', default=".",help="input directory")
parser.add_option('--outdir', default="bam", help="output directory")
parser.add_option('--pattern1', default="pairedR1_001.fastq.gz", help="filename pattern for forward")
parser.add_option('--pattern2', default="pairedR2_001.fastq.gz", help="filename pattern for forward")
#parser.add_option('--pattern1', default="trim.R1.fastq.gz", help="filename pattern for forward")
parser.add_option('--info')
parser.add_option('--metadf', default="fq2sample.tsv", help="file name to sample name")
parser.add_option('--n', default=10, help="how many command in a bash")
opts, args = parser.parse_args()
indir = opts.indir
outdir = opts.outdir
pattern1 = opts.pattern1
pattern2 = opts.pattern2
metadf = opts.metadf
group = opts.n
info = opts.info

ref = 'GCA_000003195.3_Sorghum_bicolor_NCBIv3_genomic.fna'

df = pd.read_csv(metadf, sep = "\t")
df['seqname'] = [x.split(";")[0] for x in df.submitted_ftp]
df['seqname'] = [os.path.split(x)[1][:-15] for x in df.seqname]
df.seqname = df.seqname + "."
df['PI'] = [x.split("-")[2] for x in df.run_alias]

def bwacommand(fn):
    fpath, seqname = os.path.split(fn)
    input1 = fn + pattern1
    RG2=re.split(' |-|_|/',seqname)
    ID = RG2[2] + RG2[5] + RG2[11]
    SM = df.loc[df['seqname']==seqname,'PI'].item()
    with gzip.open(input1, 'rt') as f1:
        RG1=re.split(' |:',f1.readline())
        ID = RG1[0]+':'+RG1[1]
        PU = ID+':'+RG1[2]
    LB = fn[:-6]
    input2 = fn + pattern2
    return rf'bwa mem -t 16 {ref} {input1} {input2} -R "@RG\tID:{ID}\tSM:{SM}\tPL:ILLUMINA\tPU:{PU}\tLB:{LB}" | samtools view -bS -o {outdir}/{seqname}bam'


fpes = glob.glob(f'{indir}/*{pattern1}', recursive = True)

n = 0
i = 1


def command(fpes):
    lentail = len(pattern1)
    for fn1 in fpes:
        fn = fn1[:-lentail]
        cmd = bwacommand(fn)
        yield cmd

cmds = command(fpes)

i = 1
fslurm = open("align"+str(i)+".slurm", "w")
fslurm.write(slurm_head(str(i)))


for cmd in cmds:
    fslurm.write(f"{cmd}\n\n")
    n += 1
    if n > group:
        fslurm.close()
        n = 0
        i += 1
        fslurm = open("align"+str(i)+".slurm", "w")
        fslurm.write(slurm_head(str(i)))
