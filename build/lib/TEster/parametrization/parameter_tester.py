from scipy import stats
import subprocess
import csv
import os
import TEster.utils.tester_utils as tester_utils
from TEster.analysis.gff_parser import get_accuracy, get_gff_path
from TEster.parametrization.parameter import Parameter


def split_gb_results(csv_path, mean) -> (dict, dict):
    """
    Splits the results obtained into good/bad based on median

    Parameters
    ----------
    csv_path : str
        The path to the csv file that contains the results
    mean : float
        The value that splits results into good and bad

    Returns
    -------
    good
        dictionary mapping all parameters to list of good values
    bad
        dictionary mapping all parameters to list of bad values
    """
    line_count = 0
    good = {}
    bad = {}
    for p in tester_utils.param_defaults:
        good[p] = []
        bad[p] = []
    good["Accuracy"] = []
    bad["Accuracy"] = []

    with open(csv_path, "r") as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        for row in csv_reader:
            line_count += 1
            if line_count == 1:
                continue
            if round(float(row[-3]), 3) > mean:
                tester_utils.add_to_dict(row, good)
            else:
                tester_utils.add_to_dict(row, bad)

    return good, bad


def ks_test(parameters, good, bad) -> bool:
    """
    Wraps the Kolmogorov-Smirnov test and sets the values distribution
    to good values

    Parameters
    ----------
    parameters : Parameter
        list of Parameter objects
    good : dict
        dictionary mapping all parameters to list of good values
    bad : dict
        dictionary mapping all parameters to list of bad values

    Returns
    -------
    bool
        True if there is a significant difference between good and bad
        False if the difference is less than 0.2
    """

    differing_distributions = False

    for param in parameters:
        if len(good[param.name]) == 0:
            return False
        param.calculate_kde(good[param.name], good["Accuracy"])
        bad_values = tester_utils.create_values_list(param.name, min(bad[param.name]), max(bad[param.name]))
        if len(param.values) > 0 and len(bad_values) > 0:
            KS_statistic, p_value = stats.ks_2samp(param.values, bad_values)

            if KS_statistic > 0.1:
                differing_distributions = True

    return differing_distributions


def run_nester_iterations(outcsv, generated_path, parameters, iterations) -> int:
    """
    Runs iterations of the parametrisation based on values chosen by
    tester_utils.choose_value.
    Adds them to config.yml and runs nested-nester

    Parameters
    ----------
    outcsv : TextIOWrapper
        file wrapped for writing results into
    generated_path : str
        path to the generated sequence file
    parameters : list
        list of Parameter objects

    Returns
    -------
    int
        sum of all the accuracies obtained to be used in calculating the mean
    """
    accuracy_sum = 0
    generated_file = "{}TEster_generated.fa".format(generated_path)

    for i in range(iterations):
        print("Iteration", i)
        param_values = []
        for param in parameters:
            chosen_value = param.choose_value()

            tester_utils.edit_config(param.name, chosen_value)
            param_values.append(chosen_value)

        print("With parameters:", [i for i in tester_utils.param_defaults])
        print("on values:", param_values)
        subprocess.run(["nested-nester", "-d", "/tmp/TEster/nester_results", generated_file])
        nester_gff_path = get_gff_path("/tmp/TEster/nester_results/data/")
        generated_gff_path = "{}data/GENERATED_1/GENERATED_1.gff".format(generated_path)

        with open(nester_gff_path, "r") as nester_gff, open(generated_gff_path, "r") as generated:
            TP, FP, FN, accuracy = get_accuracy(generated, nester_gff)

        accuracy_sum += accuracy

        param_values.extend((accuracy, FP, FN))
        outcsv.writerow(param_values)

    return accuracy_sum


def run_analysis(generated_path, iterations, element, out_dir=".", parameters=[], run_number=1):
    """
    Runs nester multiple times on distributed parameter values
    Recursively narrowing down the distributions until good and bad results
    differ very little. This final distribution is then optimal

    Parameters
    ----------
    generated_path : str
        path to the generated sequence file
    element : Element
        class containing information about the elements in the input database
    out_dir : str
        path to output the results provided by the user
    parameters : list
        empty in the first level of recursion, used later
        to narrow the distribution
    run_number : int
        indicates which run is taking place

    Returns
    -------
    dict
        good configurations
    dict
        bad configurations
    """
    os.makedirs(out_dir, exist_ok=True)

    if len(parameters) == 0:
        parameters = [Parameter(p_name, element) for p_name in tester_utils.param_defaults]

    print("Initiating run number: ", run_number)
    with open("{}/counts_run{}.csv".format(out_dir, run_number), "w+") as csv_file:
        outcsv = tester_utils.prepare_csv(csv_file)
        accuracy_sum = run_nester_iterations(outcsv, generated_path, parameters, iterations)

    good, bad = split_gb_results("{}/counts_run{}.csv".format(out_dir, run_number), round(accuracy_sum/iterations, 3))

    differing_distributions = ks_test(parameters, good, bad)

    if differing_distributions:
        return run_analysis(generated_path, iterations, element, out_dir, parameters, run_number+1)
    else:
        return good, bad
