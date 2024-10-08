import pysam
import argparse
import glob
import os
import numpy as np

def parseArgs():

    parser = argparse.ArgumentParser(
        description="Separate splitbams by experimental condition and read type"
    )
    parser.add_argument(
        "-d",
        "--directory",
        dest="dir",
        type=str,
        required=True,
        help="Directory of splitbams",
    )
    parser.add_argument(
        "--format",
        dest="format",
        type=str,
        required=True,
        help="Suffix of splitbams"
    )
    parser.add_argument(
        "--type",
        dest="type",
        type=str,
        required=True,
        help="Type of reads. Either DNA reads, RNA reads only, or mixed DNA RNA reads."
    )
    """
    parser.add_argument(
        "-c",
        "--conditions",
        type=list,
        required=True,
        help="List of all experimental conditions"
    )
    """
    return parser.parse_args()

def main():

    args = parseArgs()
    print("Finding splitbams in directory " + str(args.dir))
    print("With files ending in " + str(args.format))
    print("Type of reads is " + str(args.type))
    # print("The experimental conditions are " + str(args.conditions))
    dna, rna = False, False
    if "dna" in args.type.lower():
        dna = True
    if "rna" in args.type.lower():
        rna = True
    if ("dna" not in args.type.lower()) and \
        ("rna" not in args.type.lower()):
        print("Type of reads is not acceptable")
        print("It must include string dna and/or rna")
        exit()
    pathname = os.path.join(args.dir, "*" + args.format)
    files = glob.glob(pathname=pathname)
    for file in files:
        separateSplitBams(file, dna, rna)

def separateSplitBams(file, dna, rna):

    prefix = file.split(".")[:-1]
    prefix = ".".join(prefix)

    """
    suffix = []
    for condition in conditions:
        if dna and rna:
            suffix.append("_dpm_" + condition + ".bam")
            suffix.append("_rpm_" + condition + ".bam")
            suffix.append("_bpm_" + condition + ".bam")
        elif dna:
            suffix.append("_dpm_" + condition + ".bam")
            suffix.append("_bpm_" + condition + ".bam")
        elif rna:
            suffix.append("_rpm_" + condition + ".bam")
            suffix.append("_bpm_" + condition + ".bam")
    """

    if dna and rna:
        out_dpm_fvp = prefix + "_dpm_fvp.bam"
        out_rpm_fvp = prefix + "_rpm_fvp.bam"
        out_bpm_fvp = prefix + "_bpm_fvp.bam"
        out_dpm_dmso = prefix + "_dpm_dmso.bam"
        out_rpm_dmso = prefix + "_rpm_dmso.bam"
        out_bpm_dmso = prefix + "_bpm_dmso.bam"

        count_out_dpm_fvp = 0
        count_out_rpm_fvp = 0
        count_out_bpm_fvp = 0
        count_out_dpm_dmso = 0
        count_out_rpm_dmso = 0
        count_out_bpm_dmso = 0

        with pysam.AlignmentFile(file, "rb") as bamfile, \
         pysam.AlignmentFile(out_dpm_fvp, "wb", template=bamfile) as dpm_fvp, \
         pysam.AlignmentFile(out_rpm_fvp, "wb", template=bamfile) as rpm_fvp, \
         pysam.AlignmentFile(out_bpm_fvp, "wb", template=bamfile) as bpm_fvp, \
         pysam.AlignmentFile(out_dpm_dmso, "wb", template=bamfile) as dpm_dmso, \
         pysam.AlignmentFile(out_rpm_dmso, "wb", template=bamfile) as rpm_dmso, \
         pysam.AlignmentFile(out_bpm_dmso, "wb", template=bamfile) as bpm_dmso:
        
            for read in bamfile.fetch(until_eof=True):
                
                barcode = read.get_tag("RC")
                read_type = read.get_tag("RT")

                if "fvp" in barcode.lower():
                    if "dpm" in read_type.lower():
                        dpm_fvp.write(read)
                        count_out_dpm_fvp += 1
                    elif "rpm" in read_type.lower():
                        rpm_fvp.write(read)
                        count_out_rpm_fvp += 1
                    elif "bpm" in read_type.lower():
                        bpm_fvp.write(read)
                        count_out_bpm_fvp += 1
                elif "dmso" in barcode.lower():
                    if "dpm" in read_type.lower():
                        dpm_dmso.write(read)
                        count_out_dpm_dmso += 1
                    elif "rpm" in read_type.lower():
                        rpm_dmso.write(read)
                        count_out_rpm_dmso += 1
                    elif "bpm" in read_type.lower():
                        bpm_dmso.write(read)
                        count_out_bpm_dmso += 1

        print("For file " + str(file))
        print("Number of FVP DPM reads:" + str(count_out_dpm_fvp))
        print("Number of FVP RPM reads:" + str(count_out_rpm_fvp))
        print("Number of FVP BPM reads:" + str(count_out_bpm_fvp))
        print("Number of DMSO DPM reads:" + str(count_out_dpm_dmso))
        print("Number of DMSO RPM reads:" + str(count_out_rpm_dmso))
        print("Number of DMSO BPM reads:" + str(count_out_bpm_dmso))

    elif dna:
        out_dpm_fvp = prefix + "_dpm_fvp.bam"
        out_bpm_fvp = prefix + "_bpm_fvp.bam"
        out_dpm_dmso = prefix + "_dpm_dmso.bam"
        out_bpm_dmso = prefix + "_bpm_dmso.bam"

        count_out_dpm_fvp = 0
        count_out_bpm_fvp = 0
        count_out_dpm_dmso = 0
        count_out_bpm_dmso = 0

        with pysam.AlignmentFile(file, "rb") as bamfile, \
         pysam.AlignmentFile(out_dpm_fvp, "wb", template=bamfile) as dpm_fvp, \
         pysam.AlignmentFile(out_bpm_fvp, "wb", template=bamfile) as bpm_fvp, \
         pysam.AlignmentFile(out_dpm_dmso, "wb", template=bamfile) as dpm_dmso, \
         pysam.AlignmentFile(out_bpm_dmso, "wb", template=bamfile) as bpm_dmso:
        
            for read in bamfile.fetch(until_eof=True):
                
                barcode = read.get_tag("RC")
                read_type = read.get_tag("RT")

                if "fvp" in barcode.lower():
                    if "dpm" in read_type.lower():
                        dpm_fvp.write(read)
                        count_out_dpm_fvp += 1
                    elif "bpm" in read_type.lower():
                        bpm_fvp.write(read)
                        count_out_bpm_fvp += 1
                elif "dmso" in barcode.lower():
                    if "dpm" in read_type.lower():
                        dpm_dmso.write(read)
                        count_out_dpm_dmso += 1
                    elif "bpm" in read_type.lower():
                        bpm_dmso.write(read)
                        count_out_bpm_dmso += 1
                        
        print("For file " + str(file))
        print("Number of FVP DPM reads:" + str(count_out_dpm_fvp))
        print("Number of FVP BPM reads:" + str(count_out_bpm_fvp))
        print("Number of DMSO DPM reads:" + str(count_out_dpm_dmso))
        print("Number of DMSO BPM reads:" + str(count_out_bpm_dmso))

    elif rna:

        out_rpm_fvp = prefix + "_rpm_fvp.bam"
        out_bpm_fvp = prefix + "_bpm_fvp.bam"
        out_rpm_dmso = prefix + "_rpm_dmso.bam"
        out_bpm_dmso = prefix + "_bpm_dmso.bam"

        count_out_rpm_fvp = 0
        count_out_bpm_fvp = 0
        count_out_rpm_dmso = 0
        count_out_bpm_dmso = 0

        with pysam.AlignmentFile(file, "rb") as bamfile, \
         pysam.AlignmentFile(out_rpm_fvp, "wb", template=bamfile) as rpm_fvp, \
         pysam.AlignmentFile(out_bpm_fvp, "wb", template=bamfile) as bpm_fvp, \
         pysam.AlignmentFile(out_rpm_dmso, "wb", template=bamfile) as rpm_dmso, \
         pysam.AlignmentFile(out_bpm_dmso, "wb", template=bamfile) as bpm_dmso:
        
            for read in bamfile.fetch(until_eof=True):
                
                barcode = read.get_tag("RC")
                read_type = read.get_tag("RT")

                if "fvp" in barcode.lower():
                    if "rpm" in read_type.lower():
                        rpm_fvp.write(read)
                        count_out_rpm_fvp += 1
                    elif "bpm" in read_type.lower():
                        bpm_fvp.write(read)
                        count_out_bpm_fvp += 1
                elif "dmso" in barcode.lower():
                    if "rpm" in read_type.lower():
                        rpm_dmso.write(read)
                        count_out_rpm_dmso += 1
                    elif "bpm" in read_type.lower():
                        bpm_dmso.write(read)
                        count_out_bpm_dmso += 1

        print("For file " + str(file))
        print("Number of FVP RPM reads:" + str(count_out_rpm_fvp))
        print("Number of FVP BPM reads:" + str(count_out_bpm_fvp))
        print("Number of DMSO RPM reads:" + str(count_out_rpm_dmso))
        print("Number of DMSO BPM reads:" + str(count_out_bpm_dmso))

if __name__ == "__main__":
    main()