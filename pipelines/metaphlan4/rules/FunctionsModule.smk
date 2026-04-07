import datetime
import os
import pandas as pd

def getCurrentTime():
    return datetime.datetime.now().strftime("%Y_%b_%d_%H:%M")

def get_fastq_R1(wildcards):
    return metadata[metadata.sample_seq_id == wildcards.sample].R1.tolist()
def get_fastq_R2(wildcards):
    return metadata[metadata.sample_seq_id == wildcards.sample].R2.tolist()

def GatherReadCountFiles(metadata):
    inputString  = []

    barcodes = metadata["sample_barcode"].unique().tolist()
    for barcode in barcodes:
        samples = metadata[metadata["sample_barcode"] == barcode]["sample_seq_id"]
        for sample in samples:
            inputString.append(OUT_DIR+barcode+"/"+sample+"_ReadCount.txt")

    return inputString


def GetReadCountFiles(wildcards):
    ReadCountFiles = OUT_DIR+metadata.sample_barcode+"/"+metadata.sample_seq_id+"_NonHuman_ReadCount.txt"
    return ReadCountFiles



def GetSampleFastqFilesR1(wildcards):
    input_samples = metadata[metadata.sample_barcode == wildcards.barcode].sample_seq_id.tolist()
    input = [OUT_DIR+"{barcode}/" + x + "_R1_NonHuman.fastq.gz" for x in input_samples]
    return input

def GetSampleFastqFilesR2(wildcards):
    input_samples = metadata[metadata.sample_barcode == wildcards.barcode].sample_seq_id.tolist()
    input = [OUT_DIR+"{barcode}/" + x + "_R2_NonHuman.fastq.gz" for x in input_samples]
    return input



def getSubSamples(CountsFile,cutoff):
    ids = []
    
    file = pd.read_table(CountsFile,sep="\t",header=None)
    file.columns = ["sample_barcode","sample_seq_id", "read_count"]
    
    for barcode in file["sample_barcode"].unique().tolist():
        samples = file[file["sample_barcode"] == barcode]
        if sum( samples.read_count ) >= int(cutoff):
            ids.append(barcode)
    return ids






