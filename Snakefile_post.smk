import json
import os
import sys
import datetime


FINAL = 


rule all:
    input:
        FINAL
##############################################################################
# Profile clusters
##############################################################################

# Generate all statistics
rule generate_all_statistics:
    input:
        TAG_ALL + TAG_SAMP if merge_samples else TAG_SAMP
    output:
        CLUSTER_STATISTICS + ECDFS + CLUSTER_SIZES
    log:
        os.path.join(DIR_LOGS, "generate_all_statistics.log")
    params:
        dir = os.path.join(DIR_WORKUP, "clusters")
    conda:
        conda_env
    shell:
        '''
        python "{generate_all_statistics}" --directory "{params.dir}" --pattern .bam  \
            --xlim 30 &> "{log}"
        '''

##############################################################################
# Logging and MultiQC
##############################################################################

# Aggregate metrics using multiqc
rule multiqc:
    input:
        TAG_SAMP
    output:
        os.path.join(DIR_WORKUP, "qc/multiqc_report.html")
    log:
        os.path.join(DIR_LOGS, "multiqc.log")
    params:
        dir_qc = os.path.join(DIR_WORKUP, "qc")
    conda:
        conda_env
    shell:
        '''
        multiqc -f -o "{params.dir_qc}" "{DIR_WORKUP}" &> "{log}"
        '''

rule pipeline_counts:
    input:
        SPLIT_FQ + TRIM + BARCODEID + SPLIT_DPM_BPM + TRIM_RD + \
        FQ_TO_BAM + MERGE_BEAD + \
        Bt2_DNA_ALIGN + CHR_DNA + MASKED + MERGE_DNA + TAG_SAMP + \
        (TAG_ALL if merge_samples else []) + (SPLITBAMS if generate_splitbams else [])
    output:
        csv = os.path.join(DIR_WORKUP, "qc/pipeline_counts.csv"),
        pretty = os.path.join(DIR_WORKUP, "pipeline_counts.txt")
    log:
        os.path.join(DIR_LOGS, "pipeline_counts.log")
    conda:
        conda_env
    threads:
        10
    shell:
        '''
        {{
            python "{pipeline_counts}" \
              --samples "{samples}" \
              -w "{DIR_WORKUP}" \
              -o "{output.csv}" \
              -n {threads} | \
            column -t -s $'\t' > "{output.pretty}"
        }} &> "{log}"
        '''

##############################################################################
# Splitbams
##############################################################################

# Generate bam files for individual targets based on assignments from clusterfile
rule thresh_and_split:
    input:
        bam = os.path.join(DIR_WORKUP, "alignments/{sample}.merged.DNA.bam"),
        clusters = os.path.join(DIR_WORKUP, "clusters/{sample}.bam")
    output:
        os.path.join(DIR_WORKUP, "alignments/{sample}.DNA.merged.labeled.bam")
    log:
        os.path.join(DIR_LOGS, "{sample}.splitbams.log")
    params:
        dir_splitbams = os.path.join(DIR_WORKUP, "splitbams")
    conda:
        conda_env
    shell:
        '''
        python "{tag_and_split}" \
         -i "{input.bam}" \
         -c "{input.clusters}" \
         -o "{output}" \
         -d "{params.dir_splitbams}" \
         --min_oligos {min_oligos} \
         --proportion {proportion} \
         --max_size {max_size} \
         --num_tags {num_tags} &> "{log}"
        '''

# Generate summary statistics of individiual bam files
rule generate_splitbam_statistics:
    input:
        SPLITBAMS
    output:
        SPLITBAMS_STATISTICS
    log:
        os.path.join(DIR_LOGS, "splitbam_statistics.log")
    params:
        dir = os.path.join(DIR_WORKUP, "splitbams"),
        samples = [f"'{sample}'" for sample in ALL_SAMPLES]
    conda:
        conda_env
    threads:
        4
    shell:
        '''
        {{
            samples=({params.samples})
            for sample in ${{samples[@]}}; do
                for path in "{params.dir}"/"${{sample}}".DNA.merged.labeled*.bam; do
                    count=$(samtools view -@ {threads} -c "$path")
                    echo -e "${{path}}\t${{count}}" >> "{output}"
                done
            done
        }} &> "{log}"
        '''

rule splitbams_all:
    input:
        SPLITBAMS_STATISTICS
    log:
        SPLITBAMS_ALL_LOG
    params:
        dir = os.path.join(DIR_WORKUP, "splitbams"),
        samples = [f"'{sample}'" for sample in ALL_SAMPLES],
    conda:
        conda_env
    threads:
        4
    shell:
        '''
        {{
            targets=$(
                ls "{params.dir}"/*.DNA.merged.labeled_*.bam |
                sed -E -e 's/.*\.DNA\.merged\.labeled_(.*)\.bam/\\1/' |
                sort -u
            )
            echo "$targets"
            for target in ${{targets[@]}}; do
                echo "merging BAM files for target $target"
                path_merged_bam="{params.dir}"/"${{target}}.bam"
                samtools merge -f -@ {threads} "$path_merged_bam" "{params.dir}"/*.DNA.merged.labeled_"$target".bam
                samtools index -@ {threads} "$path_merged_bam"

                # Generate summary statistics of merged bam files
                count=$(samtools view -@ {threads} -c "$path_merged_bam")
                echo -e "${{path_merged_bam}}\t${{count}}" >> "{input}"
            done
        }} &> "{log}"
        '''

rule generate_bigwigs:
    input:
        SPLITBAMS_ALL_LOG,
        sorted_merged_mask = os.path.join(DIR_WORKUP, "mask_merge.bed")
    log:
        BIGWIGS_LOG
    params:
        chrom_map = f"--chrom_map '{path_chrom_map}'" if path_chrom_map is not None else "",
        dir_bam = os.path.join(DIR_WORKUP, "splitbams"),
        dir_bigwig = os.path.join(DIR_WORKUP, "bigwigs"),
        binsize = binsize
    conda:
        conda_env
    threads:
        10
    shell:
        '''
        {{
            mkdir -p "{params.dir_bigwig}"

            targets=$(
                ls "{params.dir_bam}"/*.DNA.merged.labeled_*.bam |
                sed -E -e 's/.*\.DNA\.merged\.labeled_(.*)\.bam/\\1/' |
                sort -u
            )
            for target in ${{targets[@]}}; do
                echo "Generating bigWig for target $target"
                path_bam="{params.dir_bam}"/"${{target}}".bam
                path_bigwig="{params.dir_bigwig}"/"${{target}}".bw

                # deepTools bamCoverage currently does not support generating empty bigWig
                # files from BAM files with no aligned reads. See
                # https://github.com/deeptools/deepTools/issues/598
                #
                # This situation can occur when there are clusters with only oligo (BPM) reads
                # and no chromatin (DPM) reads.
                #
                # In such cases, create an empty (i.e., 0 byte) bigWig file.
                n_reads=$(samtools view -c "$path_bam")
                if [ $n_reads -eq "0" ]; then
                    echo "- No reads in BAM file for target. Creating empty bigWig file."
                    touch "$path_bigwig"
                    continue
                fi

                # calculate genome size from header of BAM file
                # subtract regions (from selected chromosomes) in the mask
                effective_genome_size=$(
                    python {calculate_effective_genome_size} "$path_bam" \
                      -m {input.sorted_merged_mask} \
                      {params.chrom_map} \
                      -t {threads}
                )
                echo "- Effective genome size: $effective_genome_size"

                bamCoverage \
                  --binSize {params.binsize} \
                  --normalizeUsing RPGC \
                  --effectiveGenomeSize $effective_genome_size \
                  -p {threads} \
                  --bam "$path_bam" \
                  --outFileName "$path_bigwig"
            done
        }} &> "{log}"
        '''