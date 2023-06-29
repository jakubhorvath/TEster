import csv
import os

class ElementCounter:
    """
    Class that holds the amount of FP, FN, TP
    """

    def __init__(self):
        self.true_positives = 0
        self.false_positives = 0
        self.false_negatives = 0


def get_gff_path(path):
    """
    Solves the issues emerging from the unknown
    input sequence name and path.
    returns the path to the required gff file created by nester

    Parameters
    ----------
    path : string
        path to the directory containing the gff file

    Returns
    -------
    string
        full path to the gff file present in the

    """
    sequence_name = os.listdir(path + "/")[0]
    gff_file_name = "{}_genome_browser.gff".format(sequence_name)
    return "{}/{}/{}".format(path, sequence_name, gff_file_name)


def calculate_accuracy(true_positives, false_positives, false_negatives):
    """
    Calculates the accuracy with emphasis on false positives

    Parameters
    ----------
    true_positives : int
        number of elements correctly identified
    false_positives : int
        number of elements incorrectly identified
    false_negatives : int
        number of unidentified elements

    Returns
    -------
    F-score calculated from the input amounts
    """
    precision = 0
    if true_positives + false_positives != 0:
        precision = true_positives/(true_positives + false_positives)

    recall = true_positives / (true_positives + false_negatives)

    if 0.09 * precision + recall != 0:
        f1_score = ((1.09) * precision * recall) / (0.09 * precision + recall)
        return f1_score
    return 0


def find_false_positives(generated_gff, nester_gff, element_counter):
    """
    Detects number of true and false positives in the resulting .gff file
    from nester compared to the one provided by generator

    Parameters
    ----------
    generated_gff : TextIOWrapper
        GFF3 file produced by TE-generator
    nester_gff : TextIOWrapper
        GFF3 file produced by TE-nester
    element_counter : ElementCounter
        Object to load the detected counts into
    """

    for line in nester_gff:
        TP_detected = False
        # if at beginning, empty line or not a TE base
        if len(line) > 0 and line[0] != "#" and line.split()[2] == "nested_repeat":
            generated_gff.seek(0)
            split_line = line.split()

            # beginning and end sites ofdetected TEs
            ltr_base = int(split_line[3])
            ltr_end = int(split_line[4])

            relative_deviation = (ltr_end - ltr_base) * 0.07

            # for each line of TE-nester detected transposons, parse all generated TEs
            for generated_line in generated_gff:

                generated_split = generated_line.split()

                # if at beginning, empty line or not a TE base
                if len(generated_line) > 0 and generated_line[0] != "#" and generated_split[2] == "te_base":

                    # beginning and end sites of generated TEs
                    base = int(generated_split[3])
                    end = int(generated_split[4])

                    # detect if found TE is true or false positive
                    if abs(base - ltr_base) <= relative_deviation and abs(end - ltr_end) <= relative_deviation:
                        element_counter.true_positives += 1
                        TP_detected = True

            # if no true positives were detected
            if not TP_detected:
                element_counter.false_positives += 1


def find_false_negatives(generated_gff, nester_gff, element_counter):
    """
    Detects number of false negatives in the resulting .gff file
    from generator compared to the one created by nester

    Parameters
    ----------
    generated_gff : TextIOWrapper
        GFF3 file produced by TE-generator
    nester_gff : TextIOWrapper
        GFF3 file produced by TE-nester
    element_counter : ElementCounter
        Object to load the detected counts into
    """
    generated_gff.seek(0)

    for gen_line in generated_gff:
        is_false_negative = True

        if len(gen_line) > 0 and gen_line[0] != "#" and gen_line.split()[2] == "te_base":
            nester_gff.seek(0)

            split_gen_row = gen_line.split()

            # beginning and end sites ofdetected TEs
            gen_base = int(split_gen_row[3])
            gen_end = int(split_gen_row[4])

            # the permitted deviation from the actual position of the LTR
            relative_deviation = (gen_end - gen_base) * 0.07

            for nester_line in nester_gff:
                split_nester_row = nester_line.split()

                if len(nester_line) > 0 and nester_line[0] != "#" and split_nester_row[2] == "nested_repeat":
                    nest_base = int(split_nester_row[3])
                    nest_end = int(split_nester_row[4])

                    # detect if found TE is true or false positive
                    if abs(nest_base - gen_base) <= relative_deviation and abs(nest_end - gen_end) <= relative_deviation:
                        is_false_negative = False

            if is_false_negative:
                element_counter.false_negatives += 1


def get_accuracy(generated_gff, nester_gff, parameter="", out_csv=None):
    """
    Calculates the accuracy from the given gff files and outputs them
    to a csv file if one is provided

    Parameters
    ----------
    generated_gff : TextIOWrapper
        GFF3 file produced by TE-generator
    nester_gff : TextIOWrapper
        GFF3 file produced by TE-nester
    parameter : string
        parameter name
    out_csv : TextIOWrapper
        csv output file

    Returns
    -------
    int
        number of true positives
    int
        number of false positives
    int
        number of false negatives
    float
        accuracy calculated from the three above

    """

    nester_gff.seek(0)
    element_counter = ElementCounter()

    find_false_positives(generated_gff, nester_gff, element_counter)

    find_false_negatives(generated_gff, nester_gff, element_counter)

    accuracy = calculate_accuracy(element_counter.true_positives,
                                  element_counter.false_positives,
                                  element_counter.false_negatives)

    if out_csv is not None:
        writer = csv.writer(out_csv, delimiter=',',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
        writer.writerow([parameter,
                         element_counter.true_positives,
                         element_counter.false_positives,
                         element_counter.false_negatives,
                         accuracy])

    return element_counter.true_positives, element_counter.false_positives, element_counter.false_negatives, accuracy
