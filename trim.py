import os
from pathlib import Path
import glob
from optparse import OptionParser
import numpy as np

def slurm_head(i): 
    return """#!/bin/bash
#SBATCH --job-name=SAPtrim
#SBATCH --error=trim{i}.err
#SBATCH --output=trim{i}.out
#SBATCH --mem=10G
#SBATCH --time=11:00:00

ml java trimmomatic
"""
cmds = []

parser = OptionParser()
parser.add_option('--indir', default=".",help="input directory")
parser.add_option('--outdir', default="trim", help="output directory")
parser.add_option('--pattern1', default="R1_001.fastq.gz", help="filename pattern for forward")
parser.add_option('--pattern2', default="R2_001.fastq.gz", help="filename pattern for forward")
#parser.add_option('--pattern1', default="trim.R1.fastq.gz", help="filename pattern for forward")
parser.add_option('--info')
parser.add_option('--n', default=10, help="how many command in a bash")
opts, args = parser.parse_args()
indir = opts.indir
outdir = opts.outdir
pattern1 = opts.pattern1
pattern2 = opts.pattern2
group = opts.n
info = opts.info

def trimcommand(input1, input2, pair1, unpair1, pair2, unpair2):
    return f"java -jar $TM_HOME/trimmomatic.jar PE -threads 4 {input1} {input2} {pair1} {unpair1} {pair2} {unpair2} TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:40"

def notpe(indir):    
    folders = glob.iglob(f"{indir}/*")
    npe = {}
    for f in folders:
        seqs = glob.glob(f"{f}/*.fastq.gz")
        if len(seqs)==2: continue
        else: npe[f]=len(seqs)
    return npe

fpes = glob.glob(f'{indir}/**/*{pattern1}', recursive = True)

n = 0
i = 1


def command(fpes, outdir, pattern1, pattern2):
    lentail = len(pattern1)
    for fn1 in fpes:
        fn = fn1[:-lentail]
        fn2 = f"{fn}{pattern2}"
        input1 = fn1
        input2 = fn2
        fpath, fn = os.path.split(fn)
        pair1 = f"{outdir}/{fn}.paired.{pattern1}"
        unpair1 = f"/scratch/{fn}.unpaired{pattern1}"
        pair2 = f"{outdir}/{fn}.paired.{pattern2}"
        unpair2 = f"/scratch/{fn}.unpaired{pattern2}" # didn't work?ls
        cmd = trimcommand(input1, input2, pair1, unpair1, pair2, unpair2)
        yield cmd

cmds = command(fpes, outdir, pattern1, pattern2)

i = 1
fslurm = open("trim"+str(i)+".slurm", "w")
fslurm.write(slurm_head(str(i)))


for cmd in cmds:
    fslurm.write(f"{cmd}\n\n")
    n += 1
    if n > group:
        fslurm.close()
        n = 0
        i += 1
        fslurm = open("trim"+str(i)+".slurm", "w")
        fslurm.write(slurm_job)
