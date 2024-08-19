from tag_bam import *
from threshold_tag_and_split import *
import glob

pathname = "/central/scratchio/mblanco/albert/workup/splitbams/*.bam"
print(glob.glob(pathname=pathname))
print(len(glob.glob(pathname=pathname)))