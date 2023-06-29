import os
from Bio import SeqIO
import subprocess
from nested.core.te_harvest import TE_harvest 


def extract_sequence(transposon, record):
    """
    Extracts the transposon from the BioPython record

    Parameters
    ----------
    transposon : nested.core.te.TE
        the transposon to extract
    record : Bio.SeqIO.SeqRecord
        the sequence to extract from

    Returns
    -------
        string that represents the transposon
    """
    return record.seq[transposon.location[0]:transposon.location[1]]


def create_ltr_database(input_file, database_path):
    """
    Uses the static method of the TE_harvest class from TE-nester to extract transposons from the input sequence
    and create a database of LTRs

    Args:
        input_file str: input fasta file name
        database_path str: path where database is to be created 

    Returns:
        str : path to database
    """
    with open(database_path, "w+") as out_fasta:
        records = SeqIO.parse(input_file, "fasta")
        for record in records:
            transposons = TE_harvest.run(record.id, record.seq)
            for transposon in transposons:
                for t in transposon:
                    seq_count += 1
                    out_fasta.write(">{} {}\n".format(record.id, seq_count))
                    out_fasta.write("{}\n".format(extract_sequence(t, record)))
    return database_path
    


def run_harvest(input_file):
    """ Runs LTR_harvest on input sequence and extracts transposon sequences

    Args:
        input_file str: input fasta file path

    Returns:
        str : path to database
    """
    path = "/tmp/TEster/LTRHarvest"
    os.makedirs(path, exist_ok=True)
    database_path = "{}artificial_database.fa".format(path)

    return create_ltr_database(input_file, database_path)