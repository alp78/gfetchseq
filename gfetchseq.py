# gfetchseq [interval_file] [genome] [api_key] [yes/no]
from bioblend import galaxy
from bioblend.galaxy import objects
from bioblend.galaxy.tools.inputs import inputs
import json
import csv
import time
import datetime
import logging
import pandas as pd
from sys import argv
from os import remove, rename
from termcolor import colored
from Bio import SeqIO

def convert_interval_to_bed(interval_file):
    records_list = []
    with open(interval_file, 'rt') as f:
        content = csv.reader(f, delimiter='\t')
        for line in content:
            records_list.append(line)
    bed_list = []
    for record in records_list[1:]:
        line = [record[0], record[1], record[2], record[4], 0, record[3]]
        bed_list.append(line)
    bed_filename = f'{interval_file.rstrip(".interval")}.bed'
    with open(bed_filename, 'wt') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(bed_list)
    logger.info(f"Converted {colored(interval_file, 'green')} to {colored('bed', 'green')} format")
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
    logger.info(f'Retrieving Dataset...')
    dataset = gio.histories.list()[0].get_datasets()[-1]
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

def download_fasta(dataset_id, path):
    start = time.time()
    fasta_res = gi.datasets.download_dataset(dataset_id, file_path=path, use_default_filename=False)
    logger.info(f"Downloading {colored(fasta_res, 'green')} to current folder...")
    end = time.time()
    dl_time = f'{str(datetime.timedelta(seconds=int(end - start)))}'
    logger.info(f"Downloading completed in {colored(dl_time, 'yellow')}")

def fasta_to_lower(source):
    records = (rec.lower() for rec in SeqIO.parse(source, "fasta"))
    out_file = f'{source.rstrip(".fasta")}_lower.fasta'
    count = SeqIO.write(records, out_file, "fasta")
    remove(source)
    logger.info(f'Converted {colored(count, "yellow")} records to lower case')
    return out_file

def fasta_format_id(interval_file, fasta_file):
    df = pd.read_csv(interval_file, sep='\t', skiprows=(0), header=(0))
    index = 0
    records = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        g_name = df.loc[index,'genoName']
        g_start = str(df.loc[index,'genoStart'])
        g_end = str(df.loc[index,'genoEnd'])
        g_strand = df.loc[index,'strand']
        record.id = f"{g_name}_{g_start}-{g_end}{g_strand}"
        record.name = f"[{index+1}] {df.loc[index,'repName']}"
        record.description = f"[{index+1}] {df.loc[index,'repName']}"
        records.append(record)
        index += 1
    out_file = f'{fasta_file.rstrip(".fasta")}_fid.fasta'
    with open(out_file, "w") as output_handle:
        SeqIO.write(records, output_handle, "fasta")
    remove(fasta_file)
    final_file = f'{fasta_file.rstrip("_lower_fid.fasta")}.fasta'
    rename(out_file, final_file)
    logger.info('Formatted fasta records IDs')

def clean_history(hist_id, dataset_id_list):
    start = time.time()
    logger.info('Removing files from Galaxy...')

    for ds in dataset_id_list:
        gi.histories.delete_dataset(hist_id, ds, purge=True)
    end = time.time()
    purge_time = f'{str(datetime.timedelta(seconds=int(end - start)))}'
    logger.info(f"Purge completed in {colored(purge_time, 'yellow')}")



full_start = time.time()

# setting up logger
logger = logging.getLogger()
logger.setLevel(logging.INFO)
formatter = logging.Formatter('%(message)s')
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
ch.setFormatter(formatter)
logger.addHandler(ch)
logger.info('Connecting to Galaxy server...')

# initialize both galaxy clients (tools and objects)
gkey = argv[3]
gurl = 'https://usegalaxy.org'
gi = galaxy.GalaxyInstance(gurl, gkey)
gio = objects.GalaxyInstance(gurl, gkey)
guser = gi.users.get_current_user()['username']
hist_id = ''

logger.info(f"Connected to user account {colored(guser, 'red')}")

in_file = argv[1]
bed_file = convert_interval_to_bed(in_file)

_hist = gi.histories.get_histories()[0]['id']
_type = 'bed'
_genome = argv[2]

dataset = upload_bedfile_and_get_dataset(_hist, bed_file, _type, _genome)

dataset_id = dataset.id

output_id = fetch_sequences_and_get_dataset_id(dataset)

fasta_file = f'{bed_file.rstrip(".bed")}.fasta'
download_fasta(output_id, fasta_file)

logger.info('Formatting fasta file...')

lower_fasta = fasta_to_lower(fasta_file)
fasta_format_id(in_file, lower_fasta)

id_list = [dataset_id, output_id]

if argv[4] == 'yes':
    clean_history(_hist, id_list)

full_end = time.time()
full_time = f'{str(datetime.timedelta(seconds=int(full_end - full_start)))}'
logger.info(f"Full process completed in {colored(full_time, 'yellow')}")
