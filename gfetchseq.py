# gfetchseq [interval_file] [genome] [api_key] [yes/no]
import bioblend
from bioblend import galaxy
from bioblend.galaxy import objects
from bioblend.galaxy.tools.inputs import inputs
import json
import csv
import time
import datetime
import logging
import pandas as pd
from sys import argv, exit, stderr
from os import remove, rename, path
from termcolor import colored
from Bio import SeqIO

# setting up logger
logger = logging.getLogger()
logger.setLevel(logging.INFO)
formatter = logging.Formatter('%(message)s')
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
ch.setFormatter(formatter)
logger.addHandler(ch)

def convert_interval_to_bed(interval_file):
    records_list = []
    with open(interval_file, 'rt') as f:
        content = csv.reader(f, delimiter='\t')
        for line in content:
            records_list.append(line)
    bed_list = []
    # CASE: genName	 genStart	genEnd
    if len(records_list[1]) == 3:
        for record in records_list[1:]:
            line = [record[0], record[1], record[2]]
            bed_list.append(line)
    # CASE: genName	 genStart	genEnd  Strand
    if (len(records_list[1]) == 4) and (records_list[1][3] == '-' or records_list[1][3] == '+'):
        for record in records_list[1:]:
            line = [record[0], record[1], record[2], 'n/a', 0, record[3]]
            bed_list.append(line)
    # CASE: genName	 genStart	genEnd  repName
    if (len(records_list[1]) == 4) and (records_list[1][3] != '-' and records_list[1][3] != '+'):
        for record in records_list[1:]:
            line = [record[0], record[1], record[2], record[3]]
            bed_list.append(line)
    # CASE: genName	 genStart	genEnd  Strand  repName           
    if (len(records_list[1]) == 5) and (records_list[1][3] == '-' or records_list[1][3] == '+'):
        for record in records_list[1:]:    
            line = [record[0], record[1], record[2], record[4], 0, record[3]]
            bed_list.append(line)
    # CASE: genName	 genStart	genEnd  repName   Strand      
    if (len(records_list[1]) == 5) and (records_list[1][4] == '-' or records_list[1][4] == '+'):
        for record in records_list[1:]:    
            line = [record[0], record[1], record[2], record[3], 0, record[4]]
            bed_list.append(line)
    bed_filename = f'{interval_file[:-9]}.bed'
    with open(bed_filename, 'wt') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(bed_list)
    logger.info(f"Converted {colored(interval_file, 'yellow')} to {colored(bed_filename, 'green')}")
    return bed_filename

def upload_bedfile_and_get_dataset(_hist, _file, _type, _genome):
    logger.info(f"Uploading {colored(_file, 'green')} to Galaxy...")
    start = time.time()
    job = gi.tools.upload_file(_file, _hist, file_type=_type, dbkey=_genome, space_to_tab=True)
    job_id = job['jobs'][0]['id']
    while True:
        job_state = gi.jobs.get_state(job_id)
        if job_state == 'ok':
            break
        time.sleep(30)
    end = time.time()
    job_time = f'{str(datetime.timedelta(seconds=int(end - start)))}'
    logger.info(f"File uploaded in {colored(job_time, 'yellow')}")
    logger.info(colored('Retrieving Dataset...' , 'blue'))
    try:
        gio = objects.GalaxyInstance(gurl, gkey)
        dataset = gio.histories.list()[0].get_datasets()[-1]
    except Exception as e:
        logger.info(f'{colored("Problem processing the dataset on Galaxy server", "red")}')
        logger.info(f'{colored(f"Please retry...", "yellow")}')
        exit(0)        

    return dataset

def fetch_sequences_and_get_dataset_id(dataset):
    start = time.time()
    logger.info(f"Fetching DNA sequences from {colored(dataset.genome_build, 'green')} in {colored('fasta', 'green')} format...")
    dataset_name = dataset.name
    dataset_id = dataset.id
    hist_id = dataset.container_id
    tool_id = 'Extract genomic DNA 1'
    tool_input = inputs().set_dataset_param(dataset_name, dataset_id, src="hda")
    job = gi.tools.run_tool(history_id = hist_id, tool_id = tool_id, tool_inputs = tool_input)
    job_id = job['jobs'][0]['id']
    dataset_id = job['outputs'][0]['id']
    while True:
        job_state = gi.jobs.get_state(job_id)
        if job_state == 'ok':
            break
        time.sleep(30)
    end = time.time()
    job_time = f'{str(datetime.timedelta(seconds=int(end - start)))}'
    logger.info(f"Extraction completed in {colored(job_time, 'yellow')}")
    return dataset_id

def download_fasta(dataset_id, fasta_file):
    start = time.time()
    fasta_res = gi.datasets.download_dataset(dataset_id, file_path=fasta_file, use_default_filename=False)
    logger.info(f"Downloading {colored(fasta_res, 'green')} to current folder...")
    end = time.time()
    dl_time = f'{str(datetime.timedelta(seconds=int(end - start)))}'
    logger.info(f"Downloading completed in {colored(dl_time, 'yellow')}")

def fasta_to_lower(fasta_file):
    records = (rec.lower() for rec in SeqIO.parse(fasta_file, "fasta"))
    out_file = f'{fasta_file[:-6]}_lower.fasta'
    count = SeqIO.write(records, out_file, "fasta")
    remove(fasta_file)
    logger.info(f"Converted {colored(count, 'yellow')} records to lower case")
    return out_file
    
def fasta_format_id(bed_file, fasta_file):
    df = pd.read_csv(bed_file, sep='\t', skiprows=(0), header=None)
    index = 0
    records = []

    # CASE: geName, genStart, genEnd
    if len(df.columns) == 3:
        for record in SeqIO.parse(fasta_file, "fasta"):
            g_name = df.loc[index, 0]
            g_start = str(df.loc[index, 1])
            g_end = str(df.loc[index, 2])
            if '+' in record.id:
                g_strand = '+'
            else:
                g_strand = '-'
            record.id = f"{g_name}_{g_start}-{g_end}{g_strand}"
            record.name = f"{[index]}"
            record.description = f"{[index]}"
            records.append(record)
            index += 1

    # CASE: genName	 genStart	genEnd  repName
    if len(df.columns) == 4:
        for record in SeqIO.parse(fasta_file, "fasta"):
            g_name = df.loc[index, 0]
            g_start = str(df.loc[index, 1])
            g_end = str(df.loc[index, 2])
            r_name = df.loc[index, 3]
            if '+' in record.id:
                g_strand = '+'
            else:
                g_strand = '-'
            record.id = f"{g_name}_{g_start}-{g_end}{g_strand}"
            record.name = f"{r_name}"        
            record.description = f"[{index+1}] {r_name}"
            records.append(record)
            index += 1

    # CASE: geName, genStart, genEnd, repName, Score, Strand
    if len(df.columns) == 6:
        for record in SeqIO.parse(fasta_file, "fasta"):
            g_name = df.loc[index, 0]
            g_start = str(df.loc[index, 1])
            g_end = str(df.loc[index, 2])
            r_name = df.loc[index, 3]
            g_strand = df.loc[index, 5]
            record.id = f"{g_name}_{g_start}-{g_end}{g_strand}"
            record.name = f"{r_name}"
            record.description = f"[{index+1}] {r_name}"
            records.append(record)
            index += 1

    out_file = f'{fasta_file[:-6]}_fid.fasta'
    with open(out_file, "w") as output_handle:
        SeqIO.write(records, output_handle, "fasta")
    remove(fasta_file)
    final_file = f'{out_file[:-16]}.fasta'
    rename(out_file, final_file)
    logger.info("Formatted fasta records IDs")

def clean_history(hist_id, dataset_id_list):
    start = time.time()
    logger.info(f'{colored("Removing files from Galaxy...", "blue")}')

    for ds in dataset_id_list:
        gi.histories.delete_dataset(hist_id, ds, purge=True)
    end = time.time()
    purge_time = f'{str(datetime.timedelta(seconds=int(end - start)))}'
    logger.info(f"Purge completed in {colored(purge_time, 'yellow')}")

# PROGRAM START

# check if arguments are healthy
in_file = argv[1]

if len(argv) < 4:
    logger.info(f'{colored("Argument missing", "red")}')
    logger.info(f'{colored("gfetchseq [interval_file] [genome] [api_key] [yes/no]", "green")}')
    exit(0)

if len(argv) > 5:
    logger.info(f'{colored("Too many arguments", "red")}')
    logger.info(f'{colored("gfetchseq [interval_file] [genome] [api_key] [yes/no]", "green")}')
    exit(0)

if in_file[-9:] != '.interval':
    logger.info(f'{colored("Not a valid .interval file extension", "red")}')
    logger.info(f'{colored("example.interval", "green")}')
    exit(0)

if (len(argv) == 5 and argv[4] != 'yes') and (len(argv) == 5 and argv[4] != 'no'):
    logger.info(f'{colored("Invalid last argument", "red")}')
    logger.info(f'Do you want to purge Galaxy history after processing? {colored("yes/no", "green")}')
    exit(0)

if not path.isfile(in_file):
    logger.info(f'{colored("Cannot find interval file in current folder", "red")}')
    exit(0)

# loading list of genomes available in Galaxy
with open('galaxy_dbkeys.csv', 'rt') as f:
    reader = csv.reader(f)
    dbkeys = list(reader)[0]

# check if genome arg is in the list
if argv[2] not in dbkeys:
    logger.info(f'{colored("Genome invalid", "red")}')
    logger.info(f'{colored(f"Check dkeys file {dbkeys[:5]}...", "green")}')
    exit(0)

# check if interval file content is healthy
records_list = []
with open(in_file, 'rt') as f:
    content = csv.reader(f, delimiter='\t')
    index = 0
    for line in content:
        records_list.append(line)
        if index == 3:
            break
        index += 1
if len(records_list[1]) < 3:
    logger.info(f"{colored('Invalid content in interval file', 'red')}")
    logger.info(f'Must have min 3 columns: {colored("genName  genStart  genEnd", "green")}  {colored("Strand  repName", "yellow")}')
    exit(0)

if (len(records_list[1][1]) == 0 or len(records_list[1][2]) == 0):
    logger.info(f"{colored('Coordinates in interval file are invalid', 'red')}")
    logger.info(f'Example: {colored("chr1 100662981 100669420", "green")} {colored(" - L1PA4", "yellow")}')
    exit(0)

full_start = time.time()

logger.info(f'{colored("Connecting to Galaxy server...", "blue")}')

# initialize both galaxy clients (tools and objects)
try:
    gkey = argv[3]
    gurl = 'https://usegalaxy.org'
    gi = galaxy.GalaxyInstance(gurl, gkey)
    guser = gi.users.get_current_user()['username']
except bioblend.ConnectionError as bce:
    logger.info(f'{colored("Cannot connect to Galaxy", "red")}')
    logger.info(f'{colored(f"{bce.body}...", "yellow")}')
    exit(0)

logger.info(f"Connected to user account {colored(guser, 'red')}")

bed_file = convert_interval_to_bed(in_file)

_hist = gi.histories.get_histories()[0]['id']
_type = 'bed'
_genome = argv[2]

dataset = upload_bedfile_and_get_dataset(_hist, bed_file, _type, _genome)

dataset_id = dataset.id

output_id = fetch_sequences_and_get_dataset_id(dataset)

fasta_file = f'{bed_file[:-4]}.fasta'
download_fasta(output_id, fasta_file)

logger.info(colored('Formatting fasta file...', 'blue'))

lower_fasta = fasta_to_lower(fasta_file)
fasta_format_id(bed_file, lower_fasta)

id_list = [dataset_id, output_id]

if argv[4] == 'yes':
    clean_history(_hist, id_list)

full_end = time.time()
full_time = f'{str(datetime.timedelta(seconds=int(full_end - full_start)))}'
logger.info(f"Full process completed in {colored(full_time, 'yellow')}")
