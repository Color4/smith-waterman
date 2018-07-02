# main method to generate smith-waterman alignments
# of amino acid sequence and assess different parameters
# 21 Feb 17
# author <christa.caggiano@ucsf.edu>


from main.utils import *
import itertools
from main.false_positives import calculate_pos_rate, fpr
from glob import glob
import matplotlib.pyplot as plt
from main.optimize import optimize

def pick_gap_penalties(positive_pairs, negative_pairs, gap, extension, scoring_matrix):
    """
    controls picking the best combination of gap penalties
    :param positive_pairs: list of tupules containing true positives
    :param negative_pairs: list of tupules containing false positives
    :param gap: int to consider gap penalties
    :param extension: int to consider range of gap extension penalties
    :param scoring_matrix: path to a scoring matrix
    :return: best gap penalty/extension
    """

    blosum50 = read_scoring_matrix(scoring_matrix)  # read in scoring matrix

    penalties = list(itertools.product(list(range(1, gap)), list(range(1, extension))))  # find all combinations of gap/extension penalties
    best_fpr = (1, (0, 0))

    for gap, extension in penalties:

        # calculate the false positive rate given a true positive rate of 0.7
        calc_fpr = calculate_pos_rate(positive_pairs, negative_pairs, gap, extension, blosum50, 70, False, fpr)
        print(calc_fpr, gap, extension)

        # if the current combination gives a better fpr, set it as the best fpr
        if calc_fpr < best_fpr[0]:
            best_fpr = calc_fpr, (gap, extension)

    return best_fpr[1]


def test_scoring_matrix(scoring_matrix_list, positive_pairs, negative_pairs, gap, extension, true_positive_rates, normalize):
    """
    given a list of scoring matrices, and fixed penalties, pick matrix with the best average fpr
    :param scoring_matrix_list: list of scoring matrices path
    :param positive_pairs: list of true positive tupules
    :param negative_pairs: list of false positive tupules
    :param gap: set gap penalty
    :param extension: set extension penalty
    :param true_positive_rates: range of true positive rates to consider
    :param normalize: whether the scores should be normalized. should always be false.
    :return: a dictionary with the false positive scores of each scoring matrix as keys
    """
    fpr_matrix = {}

    for m in scoring_matrix_list:

        matrix = read_scoring_matrix(m)  # read each scoring matrix

        fpr_matrix[m] = []
        for tpr in true_positive_rates:

            # calculate the fpr for every tpr given
            fpr_matrix[m].append(calculate_pos_rate(positive_pairs, negative_pairs, gap, extension, matrix, tpr, normalize, fpr))

    print(fpr_matrix)
    return fpr_matrix


def plot_roc(true_positive_rates, false_positive_rates):
    """
    show plot of rocs with a range of tpr and fprs
    :param true_positive_rates: list of true positive rates
    :param false_positive_rates: dictionary of false positive rates (as generated in test_scoring_matrix)
    :return: None
    """

    for matrix in false_positive_rates:

        # generate line plot of fpr vs tpr for each matrix
        plt.plot(false_positives[matrix], true_positive_rates, label=str(matrix))

    # plot a square [0, 1] x [0, 1] ROC curve
    plt.legend()
    plt.xlabel("False Positive Rate")
    plt.xlim(0, 1.0)
    plt.ylabel("True Positive Rate")
    plt.ylim(0, 1.0)
    plt.show()

def pick_best_scoring(false_positive_rates):
    """
    finds the scoring matrix with the best average fpr
    :param false_positive_rates:
    :return: file path of matrix with best fpr
    """

    # calculate the mean fpr for each matrix
    mean_fpr = [(sum(false_positive_rates[m])/len(false_positive_rates[m]), m) for m in false_positive_rates]

    print(mean_fpr)
    return min(mean_fpr)[1]

if __name__ == "__main__":

    # read in positive and negative sequences as string tupules
    positive_pairs = read_sequence_pairs("../Negpairs.txt")
    negative_pairs = read_sequence_pairs("../Pospairs.txt")

    # pick the best gap penalties in the range (1,20) and best extension penalty combination in range (1, 5)
    # gap, extension = pick_gap_penalties(positive_pairs, negative_pairs, 20, 5, "../scoring/BLOSUM50")

    scoring_matrices = glob("../scoring/*")  # get all scoring matrices
    true_positive_rates = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

    # range of fps for each scoring matrix
    # false_positives = test_scoring_matrix(scoring_matrices, positive_pairs, negative_pairs, 1, 1, true_positive_rates, False)

    # plot rocs for every scoring matrix and a range of tpr
    # plot_roc([x/100 for x in true_positive_rates], false_positives)

    # pick the matrix with the best average fpr
    false_positives = {'../scoring/PAM250': [0.14000000000000001, 0.28000000000000003, 0.41999999999999998, 0.62, 0.78000000000000003, 0.88, 0.93999999999999995, 0.93999999999999995, 0.97999999999999998, 0.97999999999999998, 1.0], '../scoring/MATIO': [0.14000000000000001, 0.22, 0.38, 0.54000000000000004, 0.73999999999999999, 0.81999999999999995, 0.90000000000000002, 0.90000000000000002, 0.93999999999999995, 0.95999999999999996, 1.0], '../scoring/BLOSUM62': [0.12, 0.32000000000000001, 0.41999999999999998, 0.66000000000000003, 0.80000000000000004, 0.90000000000000002, 0.92000000000000004, 0.93999999999999995, 0.95999999999999996, 1.0, 1.0], '../scoring/PAM100': [0.14000000000000001, 0.29999999999999999, 0.41999999999999998, 0.69999999999999996, 0.73999999999999999, 0.88, 0.93999999999999995, 0.93999999999999995, 0.95999999999999996, 1.0, 1.0], '../scoring/BLOSUM50': [0.12, 0.32000000000000001, 0.35999999999999999, 0.64000000000000001, 0.76000000000000001, 0.88, 0.92000000000000004, 0.93999999999999995, 0.95999999999999996, 1.0, 1.0]}
    # best_matrix = pick_best_scoring(false_positives)
    best_matrix = "../scoring/BLOSUM62"

    # given the best average matrix, generate roc for matrix, compare to same matrix normalized by shortest sequence distance
    # normalized = {}
    # normalized[best_matrix] = false_positives[best_matrix]
    # normalized["normalized"] = list(test_scoring_matrix([best_matrix], positive_pairs,
    #                                                     negative_pairs, 1, 3, true_positive_rates, True).values())[0]
    # print(normalized)
    #
    # plot_roc([x/100 for x in true_positive_rates], normalized)

    # given the best matrix, optimize matrix to produce an ideal matrix for our combination of positive/negative sequences
    optimize("../scoring/BLOSUM62", positive_pairs, negative_pairs, 1, 3)