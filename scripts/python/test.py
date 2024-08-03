from tag_bam import *

input_bam = "/home/zyang4/chipdip-spidr/workup/alignments/subset_mixed.merged.DNA.untagged.bam"
output_bam = "/home/zyang4/chipdip-spidr/workup/alignments/subset_mixed.merged.DNA.bam"

label_bam_file_mixed(input_bam=input_bam, output_bam=output_bam, num_tags=9)