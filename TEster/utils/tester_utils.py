import csv
import ruamel.yaml
import yaml
import subprocess
from numpy import arange
from nested.config.config import config_path
from TEster.toolchain.TEster import tool_used

if tool_used == "ltr_finder":
    param_defaults = {"o": 3, "t": 1, "e": 0, "m": 2, "u": -2, "D": 20000, 
                  "d": 1000, "L": 3500, "l": 100, "p": 20, "g": 50, "G": 2,
                  "T": 4, "S": 6, "M": 0}
elif tool_used == "ltr_harvest":
    param_defaults = {"minlenltr": 100, "maxlenltr": 1000, "mindistltr":1000,
                  "maxdistltr": 15000, "similar": 85, "mintsd": 4, "maxtsd": 20,
                  "vic": 60, "xdrop": 5, "mat": 2, "mis": -2, "ins": -3, "del": -3} 
else:
    pass #TODO raise exception


def set_config_to_final(good_values, sequence_path, out_dir):
    """
    Detects the best accuracy and writes the parameters that correspond to this
    accuracy into the config.yml file

    Parameters
    ----------
    good_values : dict
        dictionary containing parameter values recognized as good
    sequence_path : str
        path to the query file
    out_dir : str
        path to the output directory
    """

    # finds index of configuration with best result
    best_position = good_values["Accuracy"].index(max(good_values["Accuracy"]))

    # writes that configuration to config.yml
    for p_name in param_defaults:
        edit_config(p_name, good_values[p_name][best_position])

    # runs TE-nester
    subprocess.run(["nested-nester", "-d", out_dir, sequence_path])


def edit_config(parameter, value):
    """
    Edits the config.yml file for the particular parameter to the given value

    Parameters
    ----------
    parameter : str
        the parameter to change
    value : int/float
        the value to assign the parameter
    """
    yaml_dumper = ruamel.yaml.YAML()

    with open(config_path, "r") as config_file:
        config = yaml.safe_load(config_file)

    config['ltr']['args'][parameter] = value

    with open(config_path, "w") as output:
        yaml_dumper.dump(config, output)


def reset_config():
    """
    Resets the parameters of ltr_finder in config.yml to their default values
    """
    for p_name in param_defaults:
        edit_config(p_name, param_defaults[p_name])


def prepare_csv(file):
    """
    Opens the csv buffer and writes the first row

    Parameters
    ----------
    file : TextIOWrapper
        the csv file wrapper

    Returns
    -------
    outcsv
        the csv writer object created
    """
    outcsv = csv.writer(file, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

    first_row = []

    for p_name in param_defaults:
        first_row.append(p_name)

    first_row.extend(("Accuracy", "False Positives", "False Negatives"))
    outcsv.writerow(first_row)

    return outcsv


def add_to_dict(row, dict):
    """
    Converts row of parameter values into a dictionary

    Parameters
    ----------
    row : list
        row read from csv file
    dict : dict
        dictionary of the parameter:values_list type
    """
    for p, value in zip(param_defaults, row[:-1]):
        if p == "M":
            dict[p].append(round(float(value), 1))
        else:
            dict[p].append(int(value))
    dict["Accuracy"].append(round(float(row[-3]), 3))


def create_values_list(param_name, min, max):
    """
    Creates list based on max and min value and the corresponding parameter

    Parameters
    ----------
    param_name : str
        name of the parameter to create the list for
    min : int/float
        minimum value to be contained in the list
    max : int/float
        maximum value to be contained in the list

    Returns
    -------
    the created list
    """
    if min == max:
        return [min]
    if param_name == "M":
        return ([i for i in arange(min, max+0.1, 0.1)])
    elif max > 100000:
        return ([i for i in range(min, max+1, 10)])
    else:
        return ([i for i in range(min, max+1)])
