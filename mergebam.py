import os
from pathlib import Path
import glob
from optparse import OptionParser
import numpy as np
import pandas as pd

def slurm_head(i): 
    return f"""#!/bin/bash
#SBATCH --job-name=mergebam
#SBATCH --ntasks-per-node=1
#SBATCH --error=mergebam{i}.err
#SBATCH --output=mergebam{i}.out
#SBATCH --mem=30G
#SBATCH --time=4:00:00

ml samtools
"""


parser = OptionParser()
parser.add_option('--indir', default=".",help="input directory")
parser.add_option('--outdir', default="bam", help="output directory")
parser.add_option('--pattern1', default=".bam", help="filename pattern for forward")
#parser.add_option('--pattern1', default="trim.R1.fastq.gz", help="filename pattern for forward")
parser.add_option('--info')
parser.add_option('--metadf', default="fq2sample.tsv", help="file name to sample name")

opts, args = parser.parse_args()
indir = opts.indir
outdir = opts.outdir
pattern1 = opts.pattern1
metadf = opts.metadf
info = opts.info


df = pd.read_csv(metadf, sep = "\t")
df['seqname'] = [x.split(";")[0] for x in df.submitted_ftp]
df['seqname'] = [os.path.split(x)[1][:-15] for x in df.seqname]
df.seqname = indir + df.seqname + "."
df['PI'] = [x.split("-")[2] for x in df.run_alias]

d = {P:df.loc[df.PI==P, "seqname"].to_list() for P in df.PI}
for P in d:
    fslurm = open("merge"+P+".slurm", "w")
    fslurm.write(slurm_head(P))
    cmd = f"samtools merge -c -o {outdir}/{P}.bam "+ "bam ".join(d[P]) +"bam"
    fslurm.write(cmd)
    fslurm.close()

