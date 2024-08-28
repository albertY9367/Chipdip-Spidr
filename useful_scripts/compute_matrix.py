import csv
import os
import yaml
import argparse

""" Python script to batch execute computeMatrix

Specify general settings in compute_matrix_settings.yaml
Specify bed files, bigwigs, and output location in compute_matrix.csv
    - Each row produces 1 computeMatrix file
    - In each row:  
        - First entry is bed file
        - Second to second last entries are bigwigs
        - Last entry is output file name

When running script on command line, specify
    '-f': folder of settings file and csv file
"""

parser = argparse.ArgumentParser()
parser.add_argument("-f",
                    action="store",
                    dest="folder",
                    required=True)
args = parser.parse_args()

print("Finding settings and csv file in " + str(args.folder))

csvfile = os.path.join(args.folder, "compute_matrix.csv")
settingfile = os.path.join(args.folder, "compute_matrix_settings.yaml")

with open(settingfile, "r") as file:
    print("Loading settings ...")
    settings = yaml.safe_load(file)
    # is absolute paths used?
    absolute_path = settings["absolute_path"]
    if absolute_path:
        print("Use absolute path")
        pass
    else:
        print("Use relative path")
        bed_folder = settings["bed_folder"]
        print("Use bed folder " + str(bed_folder))
        bigwig_folder = settings["bigwig_folder"]
        print("Use bigwig folder " + str(bigwig_folder))
        output_folder = settings["output_folder"]
        print("Use output folder " + str(output_folder))
    mode = settings["mode"]
    print("Use computeMatrix mode " + str(mode))

# open the bed and bigwigs csv file
with open(csvfile, newline="") as f:
    print("Computing matrices ...")
    reader = csv.reader(f)
    if absolute_path:
        for row in reader:
            # read the bed and bigwigs paths
            bed = row[0]
            bigwigs = row[1:-1]
            bigwigs = " ".join(bigwigs)
            # read output file path
            outfile = row[-1]
            # write out the command
            cmd = "computeMatrix " + mode + " -S " + bigwigs +\
                " -R " + bed + \
                " -o " + outfile
            # run deeptools
            os.system(cmd)
            print("All matrices computed")
    else:
        for row in reader:
            # read the bed and bigwigs paths
            bed = os.path.join(bed_folder, row[0])
            bigwigs = row[1:-1]
            for i in range(len(bigwigs)):
                bigwigs[i] = os.path.join(bigwig_folder, bigwigs[i])
            bigwigs = " ".join(bigwigs)
            # read output file path
            outfile = os.path.join(output_folder, row[-1])
            # write out the command
            cmd = "computeMatrix " + mode + " -S " + bigwigs +\
                " -R " + bed + \
                " -o " + outfile
            # run deeptools
            os.system(cmd)
            print("All matrices computed")