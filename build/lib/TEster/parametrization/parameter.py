from statistics import mean, stdev
from scipy import stats
from TEster.utils.tester_utils import create_values_list
from numpy.random import choice


class Parameter:
    """
    Class representing the attributes of the called parameter.
    self.values = all the values available for the current distribution
    """

    def __init__(self, param, element):
        """
        Parameters
        ----------
        param : str
            the string name of the parameter to be created
        element : sequence_generator.Element
            info of the average element in input database
        """
        switcher = {
            "o": (0, 30, "o"),  # gap open penalty
            "t": (0, 25, "t"),  # gap extension penaltu
            "e": (0, 20, "e"),  # gap end penalty
            "m": (0, 25, "m"),  # match score
            "u": (-20, 0, "u"),  # unmatch score
            "D": (5000, 35000, "D"),  # Max 5' 3' LTR distance 140000
            "d": (200, 1800, "d"),  # Min 5' 3' LTR distance
            "L": (2000, 14000, "L"),  # Max 5' 3' LTR length
            "l": (10, 400, "l"),  # Min 5' 3' LTR length
            "p": (10, 100, "p"),  # Min length of exact match pair
            "g": (10, 200, "g"),  # Max joined pair gap
            "G": (0, 20, "G"),  # Max gap between RT sub-domains
            "T": (2, 10, "T"),  # Min subdomains found in RT domain
            "S": (0, 10, "S"),  # Output score limit
            "M": (0, 1, "M")  # min LTR similarity threshold
        }
        self.min, self.max, self.name = switcher.get(param)
        self.set_values_initial(create_values_list(self.name, self.min, self.max))

    def set_values_initial(self, values):
        """
        Sets the new initial values and updates distributions, mean
        and the probability distribution function of each value

        Parameters
        ----------
        values : list
            list of values to assign to the Parameter object

        """
        self.values = values
        self.mean = mean(self.values)
        if len(self.values) == 1:
            self.distribution = [1]
        else:
            self.deviation = stdev(self.values)
            self.distribution = [stats.norm.pdf(i, self.mean, self.deviation) for i in self.values]

    def set_values(self, values):
        """
        Sets the new values based on the limited range obtained during the analysis
        using functionality from the initial setter function

        Parameters
        ----------
        values : list
            list of new values to assign to the Parameter object
        """
        self.set_values_initial(create_values_list(self.name, min(values), max(values)))

    def choose_value(self):
        """
        Takes a parameter, based on its values and distributions
        chooses one number based on probabilities listed in its distribution
        list

        Returns
        -------
            returns the chosen value
        """
        chosen = choice(self.values, 1, self.distribution)[0]
        if self.name == "M":
            return round(float(chosen), 1)
        else:
            return int(chosen)

    def calculate_kde(self, values, accuracies):
        """
        Calculates the Kernel Density Estimation from the values
        using accuracies as weight assignments and assigns it as
        distribution to the parameter

        Parameters
        ----------
        values : list
            input values to calculate the kde from
        accuracties : list
            accuracies achieved upon the given values
        """
        min_accuracy = min(accuracies)

        self.values = create_values_list(self.name, min(values), max(values))
        if len(self.values) == 1:
            self.distribution = [1]
        else:
            kernel = stats.gaussian_kde(values, bw_method='scott', weights=[(i-min_accuracy) for i in accuracies])
            self.distribution = list(kernel(self.values))
