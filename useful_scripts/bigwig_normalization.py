import pyBigWig
import numpy as np

input_dmso = "/central/scratchio/mblanco/albert/workup/bigwigs/AB10-A10_E4F1_Bethyl_A300-832A-T.bam.DMSO.bw"
input_fvp = "/central/scratchio/mblanco/albert/workup/bigwigs/AB10-A10_E4F1_Bethyl_A300-832A-T.bam.FVP.bw"
input_bed = "/central/scratchio/mblanco/albert/workup/bed/AB10-A10_E4F1_Bethyl_A300-832A-T_peaks.bed"

dmso = pyBigWig.open(input_dmso)
fvp = pyBigWig.open(input_fvp)

counter = 0

for line in open(input_bed):

    if "#" in line:
        continue
    else:
        counter += 1
        if counter > 5:
            break
        cols = line.strip().split()
        vals_dmso = dmso.values(cols[0], int(cols[1]), int(cols[2]))
        print("dmso", vals_dmso)
        vals_fvp = fvp.values(cols[0], int(cols[1]), int(cols[2]))
        print("fvp", vals_fvp)