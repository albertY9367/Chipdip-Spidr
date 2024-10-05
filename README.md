# Merged Chipdip-Spidr Pipeline

This pipeline merges features of Chipdip and Spidr pipeline. It is able to process DNA, RNA, and mixed DNA/RNA data. For now, the pipeline `Snakefile` produces cluster BAM files, cluster statistics, and split bams. The split bams are split by sample and protein. In the next version, split bams with the same protein but from different samples will be merged. This pipeline functions very similar to `chipdip` and `spidr` pipelines. For a user guide, see [here](https://github.com/GuttmanLab/chipdip-pipeline). 

For example settings, check `settings_dna`, `settings_rna`, and `settings_mixed`. Remember to change the config file selection in `run_pipeline.sh` before running the pipeline.

Below are some notable new features of this pipeline. 

## Cluster BAM Format

Previously, `.cluster` file format was used to store information about reads and the clusters they belong to. Now, this information is contained in cluster `.bam` files, that can be found in `workup/clusters` in a normal workflow. Each read is tagged with

- `RT`: Read type. Can be either BPM, DPM, or RPM.
- `RC`: Barcode of cluster. 

To generate cluster statistics, `scripts/python/generate_all_statistics.py` is used, which merges all the statistics generating scripts from before. To generate split bams, `script/python/threshold_tag_and_split.py` is used.

## Settings

Settings for the pipeline, as before, are specified in `config.yaml`. Besides the Chipdip and Spidr fields, two more fields are added.

1. `choose_dna_rna`: (required) takes a string input that specifies the type of data to be analyzed. If data contains DNA reads only, put `"dna"`; if data contains RNA reads only, put `"rna"`; if data contains both DNA and RNA reads, put `"dna, rna"`.
2. `single_end`: (optional) takes a boolean input that specifies the type of RNA data. If `True`, single-end analysis will be performed; if `False`, pair-end analysis will be performed. Default to `False`.

## Sample Spliting

One crucial step in this pipeline is spliting DPM, RPM, and BPM reads. This is done through `scripts/python/split_dpm_rpm_bpm_fq.py`. Chipdip and Spidr pipelines diverge here on how different read types are identified. In Chipdip pipeline, `config.txt` specifies the genomic sequence of DPM reads, which allows Barcode ID to identify and label DPM reads. In Spidr pipeline, the genomic sequence of RPM reads are not specified, preventing Barcode ID to label the RPM reads. This distinction causes their respective spliting script to function differently. As of now, the merged script functions as follows:

1. DNA only data: functions the same as `split_dpm_bpm_fq.py`.
2. RNA only data: functions the same as `split_rpm_bpm_fq.py`.
3. Mixed data: assumes that Barcode ID **is able to** identify RPM reads. This requires `config.txt` to list the genomic sequences of RPM.

## Other Useful Scripts

Folder `useful_scripts` contains a few hopefully useful scripts.

- `clean_up.py`: cleans the intermediate files in a successfully run `workup` folder. To run it, use command `python3 clean_up.py -f PATH_TO_WORKUP`.
- `compute_matrix.py`: batch computeMatrix function from deeptools. For detailed instruction, check the python file. It goes with `compute_matrix_settings.yaml` and `compute_matrix.csv`.
- `fastq2json.py`: generates `samples.json` by the command `fastq2json.py --fastq_dir <path_to_directory_of_FASTQs>`.
- `separate_splitbams.py`: separate splitbams based on experimental conditions (right now, FVP and DMSO). Use command `python3 separate_splitbams.py -d <directory_of_splitbams> --format <suffix_of_splitbams> --type <experiment_type>`. Input for `--type` needs to be the same as `choose_dna_rna` in `config.yaml`.