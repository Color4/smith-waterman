
from main.false_positives import calculate_pos_rate, tpr
from scipy.optimize import minimize
from main.utils import read_scoring_matrix
import pandas as pd
import numpy as np


def summed_tpr(matrix, positives, negatives, gap, extension):

    matrix = np.reshape(matrix, (24, 24))
    names = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X', '*']
    matrix = pd.DataFrame(data=matrix, columns=names, index=names)

    print("tpr")

    summed_tpr = 0
    for rate in [0, 10, 20]:
        # print(matrix)
        summed_tpr += calculate_pos_rate(positives, negatives, gap, extension, matrix, rate, False, tpr)

    return -summed_tpr


def optimize(matrix, positive_pairs, negative_pairs, gap, extension):

    matrix = read_scoring_matrix(matrix)
    # print(list(matrix))
    m = minimize(fun=summed_tpr, x0=matrix, args=(positive_pairs, negative_pairs, gap, extension), method="nelder-mead")
    print(m)