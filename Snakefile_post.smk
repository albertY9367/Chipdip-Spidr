import json
import os
import sys
import datetime

##############################################################################
# Initialize settings
##############################################################################

# Copy config file into logs
v = datetime.datetime.now()
run_date = v.strftime("%Y.%m.%d")

configfile: "config.yaml"

##############################################################################
# Load required settings
##############################################################################

bid_config = config.get("bID")
if bid_config is not None:
    print("Using BarcodeID config:", bid_config, file=sys.stderr)
else:
    print("Missing BarcodeID config (bID) in config.yaml", file=sys.stderr)
    sys.exit()

samples = config.get("samples")
if samples is not None:
    print("Using samples file:", samples, file=sys.stderr)
else:
    print("Missing samples file (samples) in config.yaml", file=sys.stderr)
    sys.exit()

DIR_SCRIPTS = config.get("scripts_dir")
if DIR_SCRIPTS is None:
    print("Scripts directory (scripts_dir) not specificed in config.yaml", file=sys.stderr)
    sys.exit()

try:
    dna, rna = False, False
    choose_dna_rna = config["choose_dna_rna"]
    if "dna" in choose_dna_rna.lower():
        dna = True
    if "rna" in choose_dna_rna.lower():
        rna = True
    if not "dna" in choose_dna_rna.lower() and not "rna" in choose_dna_rna.lower():
        print("DNA or RNA data not specified in config.yaml", file=sys.stderr)
        sys.exit()
except:
    print("DNA or RNA data not specified in config.yaml", file=sys.stderr)
    sys.exit()

##############################################################################
# Load optional settings
##############################################################################

email = config.get("email")
if email is not None and email != "":
    print("If any errors are encountered during the pipeline, an email will be sent to:", email, file=sys.stderr)
else:
    print("Email (email) not specified in config.yaml. Will not send email on error.", file=sys.stderr)

out_dir = config.get("output_dir")
if out_dir is not None:
    print("Using output directory:", out_dir, file=sys.stderr)
else:
    out_dir = os.getcwd()
    print("Defaulting to working directory as output directory:", out_dir, file=sys.stderr)

temp_dir = config.get("temp_dir")
if temp_dir is not None:
    print("Using temporary directory:", temp_dir, file=sys.stderr)
else:
    temp_dir = "/central/scratch/"
    print("Defaulting to temporary directory:", temp_dir, file=sys.stderr)

conda_env = config.get("conda_env")
if conda_env is None:
    conda_env = "envs/chipdip.yaml"
    print("No conda environment specified. Defaulting to envs/chipdip.yaml", file=sys.stderr)
if conda_env.strip().lower().endswith(".yaml") or conda_env.strip().lower().endswith(".yml"):
    print("Will create new conda environment from", conda_env, file=sys.stderr)
else:
    print("Using existing conda environment:", conda_env, file=sys.stderr)

merge_samples = config.get("merge_samples", False)

generate_splitbams = config.get("generate_splitbams", False)
if generate_splitbams:
    min_oligos = config.get("min_oligos", 2)
    proportion = config.get("proportion", 0.8)
    max_size = config.get("max_size", 10000)
    print("Will generate BAM files for individual targets using:", file=sys.stderr)
    print("\tmin_oligos:", min_oligos, file=sys.stderr)
    print("\tproportion:", proportion, file=sys.stderr)
    print("\tmax_size:", max_size, file=sys.stderr)
    
##############################################################################
# Location of scripts
##############################################################################

# general

tag_and_split = os.path.join(DIR_SCRIPTS, "python/threshold_tag_and_split.py")
generate_all_statistics = os.path.join(DIR_SCRIPTS, "python/generate_all_statistics.py")

##############################################################################
# Make output directories
##############################################################################

DIR_WORKUP = os.path.join(out_dir, "workup")
DIR_LOGS = os.path.join(DIR_WORKUP, "logs")

DIR_LOGS_CLUSTER = os.path.join(DIR_LOGS, "cluster")
os.makedirs(DIR_LOGS_CLUSTER, exist_ok=True)
out_created = os.path.exists(DIR_LOGS_CLUSTER)
print("Output logs path created:", out_created, file=sys.stderr)

##############################################################################
# Get sample files
##############################################################################

with open(samples) as f:
    FILES = json.load(f)
ALL_SAMPLES = sorted(FILES.keys())

TAG_SAMP = expand(
    os.path.join(DIR_WORKUP, "clusters/{sample}.bam"),
    sample=ALL_SAMPLES)

TAG_ALL = [os.path.join(DIR_WORKUP, "clusters/all.bam")]

CLUSTER_STATISTICS = [os.path.join(DIR_WORKUP, "clusters/cluster_statistics.txt")]

if dna and rna:
    CLUSTER_SIZES = [os.path.join(DIR_WORKUP, "clusters/DPM_read_distribution.pdf"),
            os.path.join(DIR_WORKUP, "clusters/DPM_cluster_distribution.pdf"),
            os.path.join(DIR_WORKUP, "clusters/RPM_read_distribution.pdf"),
            os.path.join(DIR_WORKUP, "clusters/RPM_cluster_distribution.pdf"),
            os.path.join(DIR_WORKUP, "clusters/BPM_cluster_distribution.pdf"),
            os.path.join(DIR_WORKUP, "clusters/BPM_read_distribution.pdf")]
elif dna:
    CLUSTER_SIZES = [os.path.join(DIR_WORKUP, "clusters/DPM_read_distribution.pdf"),
            os.path.join(DIR_WORKUP, "clusters/DPM_cluster_distribution.pdf"),
            os.path.join(DIR_WORKUP, "clusters/BPM_cluster_distribution.pdf"),
            os.path.join(DIR_WORKUP, "clusters/BPM_read_distribution.pdf")]
elif rna:
    CLUSTER_SIZES = [os.path.join(DIR_WORKUP, "clusters/RPM_read_distribution.pdf"),
            os.path.join(DIR_WORKUP, "clusters/RPM_cluster_distribution.pdf"),
            os.path.join(DIR_WORKUP, "clusters/BPM_cluster_distribution.pdf"),
            os.path.join(DIR_WORKUP, "clusters/BPM_read_distribution.pdf")]

ECDFS = [os.path.join(DIR_WORKUP, "clusters/Max_representation_ecdf.pdf"),
         os.path.join(DIR_WORKUP, "clusters/Max_representation_counts.pdf")]

SPLITBAMS = expand(
    [os.path.join(DIR_WORKUP, "alignments/{sample}.merged.DNA.labeled.bam"),
    os.path.join(DIR_WORKUP, "alignments/{sample}.merged.RPM.labeled.bam"),
    os.path.join(DIR_WORKUP, "alignments/{sample}.merged.labeled.bam")],
    sample=ALL_SAMPLES)

SPLITBAMS_STATISTICS = [os.path.join(DIR_WORKUP, "splitbams/splitbam_statistics.txt")]

SPLITBAMS_ALL_LOG = [os.path.join(DIR_LOGS, "splitbams_all.log")]

BIGWIGS_LOG = [os.path.join(DIR_LOGS, "bigwigs.log")]

PIPELINE_COUNTS = [os.path.join(DIR_WORKUP, "pipeline_counts.txt")]

FINAL = CLUSTER_STATISTICS + ECDFS + CLUSTER_SIZES + SPLITBAMS

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
            --xlim 30 --type "{choose_dna_rna}" &> "{log}"
        '''

##############################################################################
# Splitbams
##############################################################################

# Generate bam files for individual targets based on assignments from clusterfile
rule thresh_and_split:
    input:
        bam_dna = os.path.join(DIR_WORKUP, "alignments/{sample}.merged.DNA.bam"),
        bam_rna = os.path.join(DIR_WORKUP, "alignments/{sample}.merged.RPM.bam"),
        clusters = os.path.join(DIR_WORKUP, "clusters/{sample}.bam")
    output:
        dna = temp(os.path.join(DIR_WORKUP, "alignments/{sample}.merged.DNA.labeled.bam")),
        rna = temp(os.path.join(DIR_WORKUP, "alignments/{sample}.merged.RPM.labeled.bam")),
        bam = os.path.join(DIR_WORKUP, "alignments/{sample}.merged.labeled.bam")
    log:
        os.path.join(DIR_LOGS, "{sample}.splitbams.log")
    params:
        dir_splitbams = os.path.join(DIR_WORKUP, "splitbams")
    conda:
        conda_env
    shell:
        '''
        python "{tag_and_split}" \
         --input_dna "{input.bam_dna}" \
         --input_rna "{input.bam_rna}" \
         -c "{input.clusters}" \
         --output_dna "{output.dna}" \
         --output_rna "{output.rna}" \
         -o "{output.bam}" \
         -d "{params.dir_splitbams}" \
         --min_oligos {min_oligos} \
         --proportion {proportion} \
         --max_size {max_size} &> "{log}"
        '''