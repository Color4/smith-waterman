# smith waterman algorithm
# 20 Feb 2018
# author <christa.caggiano@ucsf.edu>

import pandas as pd
import numpy as np


def align(seq1, seq2, gap_opening, gap_extension, scoring_matrix):
    """
    controls the smith-waterman alignment algorithm
    :param seq1: sequence 1 to be aligned
    :param seq2: sequence 2 to be aligned
    :param gap_opening: penalty for opening a gap
    :param gap_extension: penalty for extending a gap
    :param scoring_matrix: scoring matrix to generate alignment
    :return: names of aligned sequences and alignment score
    """

    alignment_matrix = np.zeros((len(seq1), len(seq2)))  # intialize to zero

    for i in range(1, len(seq1)):  # for each sequence to be aligned
        for j in range(1, len(seq2)):
            left = find_max_gap(alignment_matrix[:i, j-1], gap_opening, gap_extension)  # find leftmost value
            above = find_max_gap(alignment_matrix[i-1, :j], gap_opening, gap_extension)  # find above value
            diag = int(scoring_matrix[seq1[i]][seq2[j]]) + int(alignment_matrix[i-1, j-1]) # find diagonal

            max_score = max(left, above, diag, 0)  # pick the max of the 3 previous options, or 0

            alignment_matrix[i, j] = max_score  # set the newest part of the matrix as the max

    return align_sequences(alignment_matrix, seq1, seq2)  # return alignments


def read_scoring_matrix(matrix_file):
    """
    process the scoring matrix file
    :param matrix_file: location of a scoring matrix
    :return: return dataframe of matrix
    """

    scoring_matrix = []

    with open(matrix_file) as matrix:
        for line in matrix:
            if "#" not in line:
                scoring_matrix.append(line.split())

    return make_scoring_df(scoring_matrix)


def make_scoring_df(scoring_matrix_list):
    """
    scoring matrix to be turned into a dataframe
    :param scoring_matrix_list:
    :return: df of scoring matrix
    """

    headers = scoring_matrix_list[0]
    scoring_matrix_list = scoring_matrix_list[1:]

    return pd.DataFrame(scoring_matrix_list, columns=headers, index=headers)


def find_max_gap(scores, gap_opening, gap_extension):
    """
    finds a length dependant gap
    :param scores:
    :param gap_opening:
    :param gap_extension:
    :return: max of gap penalty options (ie best alignment)
    """
    gap = gap_opening + gap_extension*(len(scores))
    gap_array = np.arange(gap,gap_opening, -gap_extension)
    scores_with_penalties = scores-gap_array
    return np.amax(scores_with_penalties)  # max of the gap scores with penalties


def align_sequences(alignment_matrix, seq1, seq2):
    """
    iterate over all the possible tracebacks, and pick the one with the highest score
    :param alignment_matrix: possible alignment scores
    :param seq1:
    :param seq2:
    :return: highest score of all paths
    """

    alignment1 = ""
    alignment2 = ""

    max_value = alignment_matrix.max()  # find the max value of the alignment matrix to begin traceback
    score = 0
    i, j = np.unravel_index(alignment_matrix.argmax(), alignment_matrix.shape)
    alignment1 += seq1[i]
    alignment2 += seq2[j]


    while max_value != 0:  # continue until coming across a 0

        score += max_value
        left = alignment_matrix[i][j - 1]
        above = alignment_matrix[i - 1][j]
        diag = alignment_matrix[i - 1][j - 1]

        # out of all the options (left, diagonal, or above) pick the one that produces the highest alignment score
        if left > above and left > diag:
            max_value = left
            j -= 1
            alignment2 += seq2[j]
            alignment1 += "-"
        elif above > diag and above > left:
            max_value = above
            i -= 1
            alignment1 += seq1[i]
            alignment2 += "-"
        else:
            max_value = diag
            i -= 1
            j -= 1
            alignment1 += seq1[i]
            alignment2 += seq2[j]

    return "".join(x for x in reversed(alignment1)), "".join(x for x in reversed(alignment2)), score

