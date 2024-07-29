# Merged Chipdip-Spidr Pipeline

This pipeline merges features of Chipdip and Spidr pipeline. It is able to process DNA, RNA, and mixed DNA/RNA data. For now, the pipeline `Snakefile` produces cluster BAM files but does not perform statistics collection. If needed, `scripts/python/generate_all_statistics.py` can be used to generate cluster statistics manually. Below are some notes about using this pipeline. 

## Settings

Settings for the pipeline, as before, are specified in `config.yaml`. Besides the Chipdip and Spidr fields, two more fields are added.

1. `choose_dna_rna`: (required) takes a string input that specifies the type of data to be analyzed. If data contains DNA reads only, put `"dna"`; if data contains RNA reads only, put `"rna"`; if data contains both DNA and RNA reads, put `"dna, rna"`.
2. `single_end`: (optional) takes a boolean input that specifies the type of RNA data. If `True`, single-end analysis will be performed; if `False`, pair-end analysis will be performed. Default to `False`.

## Sample Spliting

One crucial step in this pipeline is spliting DPM, RPM, and BPM reads. This is done through `scripts/python/split_dpm_rpm_bpm_fq.py`. Chipdip and Spidr pipelines diverge here on how different read types are identified. In Chipdip pipeline, `config.txt` specifies the genomic sequence of DPM reads, which allows Barcode ID to identify and label DPM reads. In Spidr pipeline, the genomic sequence of RPM reads are not specified, preventing Barcode ID to label the RPM reads. This distinction causes their respective spliting script to function differently. As of now, the merged script functions as follows:

1. DNA only data: functions the same as `split_dpm_bpm_fq.py`.
2. RNA only data: functions the same as `split_rpm_bpm_fq.py`.
3. Mixed data: assumes that Barcode ID **is able to** identify RPM reads. This requires `config.txt` to list the genomic sequences of RPM.