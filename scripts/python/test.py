from tag_bam import *
from threshold_tag_and_split import *

input_bam = "/central/scratchio/mblanco/albert/workup/alignments/RD_CD_1-1.merged.RPM.labeled.bam"

bead_labels = set()

with pysam.AlignmentFile(input_bam, "rb") as bam:

    for read in bam.fetch(until_eof=True):

        bead_labels.add(read.get_tag("RG"))

print(bead_labels)
