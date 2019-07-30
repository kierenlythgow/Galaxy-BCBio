import os
import sys
import subprocess
import shutil
import glob
import re
import csv
import yaml
import binascii
import copy

manifest = sys.argv[1]
metadata = sys.argv[2]
genome = sys.argv[3]
aln_type = sys.argv[4]
multiqc = sys.argv[5]
final_tar = sys.argv[6]
rdata = sys.argv[7]
rnaseq_script = "/data/galaxy/tools/bcbio/rnaseq.r" #NEED TO AMEND
cwd = 'rnaseq/work'
data = {}
elements = []
tsv = {}
tsvdata = []

subprocess.check_call(["python", "--version"])
print(os.getcwd())

print("GENOME: " + genome)
print("ALIGNER: " + aln_type)

#Make bcbio directory tree
os.makedirs("rnaseq/config")
os.makedirs("rnaseq/final")
os.makedirs("rnaseq/work")

# Input manifest file
with open(manifest, "r") as csvfile:
    #Remove trailing and interspersed whitespace
    for line in csvfile:
        cleaned = re.sub("[\\s+]", "", line).rstrip(",\n")

        if cleaned:
            elements.append(cleaned)

    #Create dictionary for samples
    csvReader = csv.DictReader(elements)
    for row in csvReader:
        id = row['sample_id']
        data[id] = row

# Input metadata file
with open(metadata, "r") as tsvfile:
    # Check sample id's match
    for line in tsvfile:
        l = line.strip()
        tsvdata.append(l)
    print(tsvdata)
    # Create dictionary for samples
    csvReader = csv.DictReader(tsvdata, delimiter='\t')
    for row in csvReader:
        print(row)
        id = row['Sample']
        tsv[id] = row

# Check sample id's match
for key in data:
    if key in tsv:
        print(key + ' matched')
    else:
        print(key + ' ERROR: Sample does not match')

# Iterate through all samples and fastq pairs
samplelist = []
fq1list = []
fq2list = []
for (k, v) in data.items():
    print("Sample: " + k)

    sample_name = k
    samplelist.append(sample_name)
    dat_1 = v['file1']
    dat_2 = v['file2']

    #Check if fastqs are compressed
    with open(dat_1, 'rb') as f:
        if binascii.hexlify(f.read(2)) == b'1f8b':
            fastq_1 = sample_name + "_1.fastq.gz"
        else:
            fastq_1 = sample_name + "_1.fastq"

        fq1list.append(fastq_1)

    with open(dat_2, 'rb') as f2:
        if binascii.hexlify(f2.read(2)) == b'1f8b':
            fastq_2 = sample_name + "_2.fastq.gz"
        else:
            fastq_2 = sample_name + "_2.fastq"

        fq2list.append(fastq_2)

    # Create symlinks of fastqs to working directory
    os.symlink(dat_1, cwd + '/' + fastq_1)
    os.symlink(dat_2, cwd + '/' + fastq_2)

print(fq1list)
print(fq2list)

bcbio_yaml = """\
details:
- algorithm:
    aligner: hisat2
  analysis: RNA-seq
  genome_build: hg38
fc_name: rnaseq
upload:
  dir: ../final
"""

#TO BE AMENDED
analysis = 'RNA-seq'
#genome = 'GRCh37'
#aligner = 'bwa'

config = yaml.load(bcbio_yaml)

#Create template YAML file
for sample, fq1, fq2 in zip(samplelist, fq1list, fq2list):

     #Create deep copy of sample data
     rep = copy.deepcopy(config['details'][0])

     rep['analysis'] = analysis
     rep['genome_build'] = genome
     rep['algorithm']['aligner'] = aln_type

     rep['description'] = sample
     rep['files'] = [fq1,fq2]

     config['details'].append(rep)

del config['details'][0]

yaml.dump(config, sys.stdout, default_flow_style=False)

with open("rnaseq/config/template_rnaseq.yaml", "w") as f:
        yaml.dump(config, f, default_flow_style=False)

subprocess.check_call(["bcbio_nextgen.py", "../config/template_rnaseq.yaml", "-n 8"], cwd=cwd)

# Glob for timestamped bcbio directory
multiqc_file = glob.glob("rnaseq/final/*rnaseq/multiqc/multiqc_report.html")

# Assert only one result directory exists
try:
    assert len(multiqc_file) == 1
    print('One result directory')
except AssertionError:
    print('ERROR: Directory number not == 1')

print(multiqc_file)
print(multiqc_file[0])

# Copy multiqc html file and tar the multiqc directory
shutil.copy(multiqc_file[0], multiqc)

subprocess.check_call(["tar",
                 "-cvf",
                 final_tar,
                 "rnaseq/final"])

# Run Rscript on rnaseq data
subprocess.check_call(["Rscript",
                       rnaseq_script,
                       "-f", metadata,
                       "-o", rdata,
                       "-b", genome,
                       "-d", "rnaseq/final"])

'''
#Example JSON file
with open("galaxy.json", "r") as read_file:
    data = json.load(read_file)

print(data["name"])

'''


