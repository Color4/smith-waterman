
from main.false_positives import calculate_pos_rate
from scipy.optimize import minimize
from main.utils import read_scoring_matrix
import pandas as pd
import numpy as np


def summed_tpr(matrix, positives, negatives, gap, extension):

    matrix = np.reshape(matrix, (24, 24))
    names = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X', '*']
    matrix = pd.DataFrame(data=matrix, columns=names)
    # matrix.set_index(matrix)

    # print(matrix)

    tpr = 0
    for rate in [0, 10, 20]:
        # print(matrix)
        tpr += calculate_pos_rate(positives, negatives, gap, extension, matrix, rate, False, tpr)

    return -tpr


def optimize(matrix, positive_pairs, negative_pairs, gap, extension):

    matrix = read_scoring_matrix(matrix)
    # print(list(matrix))
    m = minimize(fun=summed_tpr, x0=matrix, args=(positive_pairs, negative_pairs, gap, extension), method="nelder-mead")
    print(m)