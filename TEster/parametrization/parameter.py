from statistics import mean, stdev
from scipy import stats
from TEster.utils.tester_utils import create_values_list
from numpy.random import choice


class Parameter:
    """
    Class representing the attributes of the called parameter.
    self.values = all the values available for the current distribution
    """

    def __init__(self):
        """
        Parameters
        ----------a distribution of values is highly likely to resemble the simulated
        param : str
            the string name of the parameter to be created
        element : sequence_generator.Element
            info of the average element in input database
        """

        
    def set_values_initial(self, values, expected_value):
        """
        Sets the new initial values and updates distributions, mean
        and the probability distribution function of each value

        Parameters
        ----------
        values : list
            list of values to assign to the Parameter object

        """
        self.values = values
        if len(self.values) == 1:
            self.distribution = [1]
        else:
            self.deviation = stdev(self.values)
            self.distribution = [stats.norm.pdf(self.values, expected_value, self.deviation)]

    def set_values(self, values):
        """
        Sets the new values based on the limited range obtained during the analysis
        using functionality from the initial setter function

        Parameters
        ----------a distribution of values is highly likely to resemble the simulated
        values : list
            list of new values to assign to the Parameter object
        """
        values = create_values_list(self.name, min(values), max(values))
        self.set_values_initial(values, mean(values))

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
        distribution to the parametera distribution of values is highly likely to resemble the simulated

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
            accuracy = [(i-min_accuracy) for i in accuracies]
            if max(accuracy) == 0:
                self.set_values(values)
            else:
                kernel = stats.gaussian_kde(values, bw_method='scott', weights=accuracy)
                self.distribution = list(kernel(self.values))


class LtrFinderParameter(Parameter):
    def __init__(self, param):
        super().__init__()
        switcher = {
            "o": (0, 17, "o", 3),           # gap open penalty
            "t": (0, 15, "t", 1),           # gap extension penaltu
            "e": (0, 15, "e", 1),           # gap end penalty
            "m": (0, 15, "m", 2),           # match score
            "u": (-10, 0, "u", -2),         # unmatch score
            "D": (5000, 35000, "D", 20000), # Max 5' 3' LTR distance 140000
            "d": (200, 1800, "d", 1000),    # Min 5' 3' LTR distance
            "L": (2000, 14000, "L", 3500),  # Max 5' 3' LTR length
            "l": (10, 400, "l", 100),       # Min 5' 3' LTR length
            "p": (10, 100, "p", 20),        # Min length of exact match pair
            "g": (10, 200, "g", 50),        # Max joined pair gap
            "G": (0, 20, "G", 2),           # Max gap between RT sub-domains
            "T": (2, 10, "T", 4),           # Min subdomains found in RT domain
            "S": (0, 10, "S", 6),           # Output score limit
            "M": (0, 1, "M", 0.2),          # min LTR similarity threshold
            "r": (1, 18, "r", 14),          # PBS detecting threshold, min tRNA match
            "E": (0, 1, "E", 0.7)           # LTR must have edge signal (2 of PBS, PPT, TSR)
        }
        self.min, self.max, self.name, self.default = switcher.get(param)
        super().set_values_initial(create_values_list(self.name, self.min, self.max), self.default)



class LtrHarvestParameter(Parameter):
    def __init__(self, param):
        super().__init__()
        switcher = {
            "minlenltr": (0, 160, "minlenltr", 100),
            "maxlenltr": (400, 2000, "maxlenltr", 1000),
            "mindistltr": (100, 1200, "mindistltr", 1000),
            "maxdistltr": (5000, 18000, "maxdistltr", 15000),
            "similar": (0, 100, "similar", 85),
            "mintsd": (0, 10, "mintsd", 4),
            "maxtsd": (0, 40, "maxtsd", 20),
            "vic": (0, 100, "vic", 60),
            "xdrop": (0, 10, "xdrop", 5),
            "mat": (0, 10, "mat", 2),
            "mis": (-10, 0, "mis", -2),
            "ins": (-15, 0, "ins", -3),
            "del": (-15, 0, "del", -3)
        }
        self.min, self.max, self.name, self.default = switcher.get(param)
        super().set_values_initial(create_values_list(self.name, self.min, self.max), self.default)
