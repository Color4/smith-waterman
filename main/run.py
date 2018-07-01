# main method to generate smith-waterman alignments
# of amino acid sequence and assess different parameters
# 21 Feb 17
# author <christa.caggiano@ucsf.edu>


from main.utils import *
from main.smith_waterman import align
import itertools
from main.false_positives import calculate_pos_rate
from main.optimize import optimize

from glob import glob
import matplotlib.pyplot as plt


def pick_gap_penalties(positive_pairs, negative_pairs, gap, extension, scoring_matrix):

    blosum50 = read_scoring_matrix(scoring_matrix)
    penalties = list(itertools.product(list(range(1, gap)), list(range(1, extension))))
    best_fpr = (1, (0, 0))

    for gap, extension in penalties:

        fpr = calculate_pos_rate(positive_pairs, negative_pairs, gap, extension, blosum50, 70, False, fpr)
        if fpr < best_fpr[0]:
            best_fpr = fpr, (gap, extension)

    return best_fpr[1]


def test_scoring_matrix(scoring_matrix_list, positive_pairs, negative_pairs, gap, extension, true_positive_rates, normalize):

    fpr = {}

    for m in scoring_matrix_list:

        matrix = read_scoring_matrix(m)
        fpr[m] = []
        for tpr in true_positive_rates:
            fpr[m].append(calculate_pos_rate(positive_pairs, negative_pairs, gap, extension, matrix, tpr, normalize, fpr))

    return fpr


def plot_roc(true_positive_rates, false_positive_rates):

    for matrix in false_positive_rates:
        plt.plot(false_positive_rates[matrix], true_positive_rates, label=str(matrix).replace("../scoring/", ""))

    plt.legend()
    plt.xlabel("False Positive Rate")
    plt.xlim(0, 1.0)
    plt.ylabel("True Positive Rate")
    plt.ylim(0, 1.0)
    plt.show()

def pick_best_scoring(false_positive_rates):

    mean_fpr = [(sum(false_positive_rates[m])/len(false_positive_rates[m]), m) for m in false_positive_rates]

    return min(mean_fpr)[1]

if __name__ == "__main__":

    positive_pairs = read_sequence_pairs("../test_pos.txt")
    negative_pairs = read_sequence_pairs("../test_neg.txt")

    # gap, extension = pick_gap_penalties(positive_pairs, negative_pairs, 20, 5, "../scoring/BLOSUM50")

    scoring_matrices = glob("../scoring/*")
    true_positive_rates = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

    # false_positives = test_scoring_matrix(scoring_matrices, positive_pairs, negative_pairs, 1, 3, true_positive_rates, False)
    # plot_roc([x/100 for x in true_positive_rates], false_positives)
    #
    # best_matrix = pick_best_scoring(false_positives)
    # # best_matrix = "../scoring/BLOSUM62"
    #
    # normalized = {}
    # normalized[best_matrix] = false_positives[best_matrix]
    # normalized["normalized"] = list(test_scoring_matrix([best_matrix], positive_pairs,
    #                                                     negative_pairs, 1, 3, true_positive_rates, True).values())[0]
    # plot_roc([x/100 for x in true_positive_rates], normalized)

    optimize("../scoring/BLOSUM62", positive_pairs, negative_pairs, 1, 3)