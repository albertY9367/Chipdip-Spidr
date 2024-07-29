from tag_bam import *

input_bam = "/home/zyang4/chipdip-spidr/workup/alignments/sample1.merged.BPM.untagged.bam"
output_bam = "/home/zyang4/chipdip-spidr/workup/alignments/sample1.merged.BPM.bam"
num_tags = 7

label_bam_file(input_bam=input_bam, output_bam=output_bam, num_tags=num_tags)