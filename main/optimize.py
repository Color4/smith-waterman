
from main.false_positives import calculate_pos_rate, tpr
from scipy.optimize import minimize
from main.utils import read_scoring_matrix
import pandas as pd
import numpy as np


def summed_tpr(matrix, positives, negatives, gap, extension):
    """
    calculates the summed true positive rate for 0, 0.1, 0.2 false positive rates
    :param matrix: np matrix
    :param positives: list of positive tupules
    :param negatives: list of negative tupules
    :param gap: gap penalty
    :param extension: gap extension penalty
    :return: summed true positive rate
    """

    # my smith-waterman algorithm requires a labeled 24X24 pandas dataframe. reshape from numpy array
    # coerced by scikit learn
    matrix = np.reshape(matrix, (24, 24))
    names = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X', '*']
    matrix = pd.DataFrame(data=matrix, columns=names, index=names)

    summed_tpr = 0

    # calculate the true positives for 3 false positive thresholds
    for rate in [0, 10, 20]:
        summed_tpr += calculate_pos_rate(positives, negatives, gap, extension, matrix, rate, False, tpr)

    print(-summed_tpr)

    # return negative tpr because we are maximizing
    return -summed_tpr


def optimize(matrix, positive_pairs, negative_pairs, gap, extension):
    """
    using sklearn's batch gradient descent and the nelder mead minimization algorithm, optimize the starting
    BLOSUM62 matrix so that we have the best possible TPR

    :param matrix: starting matrix
    :param positive_pairs: sequence tupules of true positives
    :param negative_pairs: sequence tupules of false positives
    :param gap: gap penalty
    :param extension: gap extension
    :return: the new optimized matrix
    """

    matrix = read_scoring_matrix(matrix)
    # print(list(matrix))
    m = minimize(fun=summed_tpr, x0=matrix, args=(positive_pairs, negative_pairs, gap, extension), method="nelder-mead",
                 options={"maxiter":5, "display":True})

    return m