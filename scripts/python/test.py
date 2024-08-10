from tag_bam import *
from threshold_tag_and_split import *
# input_bam = "/central/scratchio/mblanco/albert/workup/alignments/RD_CD_1-1.merged.RPM.untagged.bam"
output_bam = "/central/scratchio/mblanco/albert/workup/alignments/RD_CD_1-1.merged.BPM.bam"

# label_bam_file_mixed(input_bam=input_bam, output_bam=output_bam, num_tags=9)
a = "dna" in "dna, rna".lower()
print(a)
