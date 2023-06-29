from Bio import SeqIO
import os
from nested.core.te import run_ltr_finder


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
    Uses parse_ltr_table from TE-nester to extract transposons from ltr_finder
    output and create a database of LTRs

    Parameters
    ----------
    input_file : str
        query sequence file
    database_path : str
        path and name of the database to be created

    Returns
    -------
    the same database_path for convenience
    """

    seq_count = 0
    with open(database_path, "w+") as out_fasta:
        records = SeqIO.parse(input_file, "fasta")
        for record in records:
            transposons = run_ltr_finder(record.id, record.seq)
            for transposon in transposons:
                seq_count += 1
                out_fasta.write(">{} {}\n".format(record.id, seq_count))
                out_fasta.write("{}\n".format(extract_sequence(transposon, record)))
    return database_path


def run_finder(input_file):
    """
    Creates a file for the LTR finder output and calls the parsing function
    on it

    Parameters
    ----------
    input_file : str
        path to the query file

    Returns
    -------
    path to created TE database file
    """

    path = "/tmp/TEster/LTR_finder/"
    os.makedirs(path, exist_ok=True)
    database_path = "{}artificial_database.fa".format(path)

    return create_ltr_database(input_file, database_path)
