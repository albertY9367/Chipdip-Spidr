import argparse
import os
import shutil

""" Python script to clean up workup folders

Currently deletes the following subfolders
    - alignments_parts
    - fastqs
    - splitfq
    - trimmed
"""

parser = argparse.ArgumentParser()
parser.add_argument("-f",
                    action="store",
                    dest="folder",
                    required=True)
args = parser.parse_args()

print("Cleaning up " + str(args.folder))

trash_list = ["alignments_parts", "fastqs", "splitfq", "trimmed"]

for folder in trash_list:

    path = os.path.join(args.folder, folder)
    if os.path.isdir(path):
        shutil.rmtree(path)
        print("Deleted " + path)