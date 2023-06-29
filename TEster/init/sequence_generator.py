import subprocess
from Bio import SeqIO
import pyfastx
from nested.core.te import TE


class EmptyInputFileException(Exception):
    """ Raised if the given database or file
        to be analysed is empty """

    def __init__(self, message):
        self.message = message


class Element:
    """
    Contains info about the average element in the input database
    parameters in the analysis will be set according to the information
    obtained from Element
    """

    def __init__(self, input_db, length):
        """
        Parameters
        ----------
        input_db : str
            the database to get info from
        length : int
            length of the average element
        """
        self.length = length
        self.score, self.site_presence = self.get_average_info(input_db)

    def get_average_info(self, input_db):
        """
        Obtains info using TE-nester's run_ltr_finder

        Parameters
        ----------
        input_db : str
            input database of elements
        """
        te_count = 0
        site_presence = 0
        score = 0

        records = SeqIO.parse(input_db, "fasta")
        for record in records:
            TEs = TE.run(record.id, record.seq)
            at_least_two = 0
            for TE in TEs:
                te_count += 1
                if TE.score is not None:
                    score += TE.score
                if "nan" not in TE.ppt:
                    at_least_two += 1
                if "nan" not in TE.pbs:
                    at_least_two += 1
                if "nan" not in TE.tsr_left:
                    at_least_two += 1
                if at_least_two >= 2:
                    site_presence += 1
        if te_count > 0:
            return score/te_count, site_presence/te_count
        return 0, 0


def sequence_generator(input_file, input_db, percentage):
    """
    Calculates the average lengths of the transposons to be given to generator
    and the average sequence length in the query file, then generates
    a corresponding sequence

    Parameters
    ----------
    input_file : str
        query sequence file path
    input_db : str
        input database file path

    Returns
    -------
    element : Element
        inof about the average element
    path to generated file : str

    """
    avg_element_length = round(pyfastx.Fasta(input_db).mean)

    # set this to 1Mbp for shorter run time
    avg_sequence_length = 1000000

    # avg_sequence_length = round(pyfastx.Fasta(input_file).mean)
    #print("Analysing properties of elements in database")
    #element = Element(input_db, avg_element_length)

    iterations = round((avg_sequence_length / avg_element_length) * (percentage/100))

    subprocess.run(['nested-generator', '-l', str(avg_sequence_length), '-i',
                    str(iterations), '-d', "/tmp/TEster",
                    input_db, "TEster_generated.fa"])

    return "/tmp/TEster/generated_data/"
