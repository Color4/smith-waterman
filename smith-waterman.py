# smith waterman algorithm
# 20 Feb 2018
# author <christa.caggiano@ucsf.edu>

import pandas as pd
import numpy as np


def read_scoring_matrix(matrix_file):

    scoring_matrix = []

    with open(matrix_file) as matrix:
        for line in matrix:
            if "#" not in line:
                scoring_matrix.append(line.split())

    return make_scoring_df(scoring_matrix)


def make_scoring_df(scoring_matrix_list):

    headers = scoring_matrix_list[0]
    scoring_matrix_list = scoring_matrix_list[1:]

    return pd.DataFrame(scoring_matrix_list, columns=headers, index=headers)

def generate_alignment_matrix(seq1, seq2, gap_opening, gap_extension, scoring_matrix):

    alignment_matrix = np.zeros((len(seq1), len(seq2)))

    for i in range(1, len(seq1)):
        for j in range(1, len(seq2)):

            left = find_max_gap(alignment_matrix[:, j-1], gap_opening, gap_extension)
            above = find_max_gap(alignment_matrix[i-1, :], gap_opening, gap_opening)
            diag = int(scoring_matrix[seq1[i]][seq2[j]]) + int(alignment_matrix[i-1, j-1])

            max_score = max(left, above, diag, 0)

            alignment_matrix[i, j] = max_score

    return(alignment_matrix)

def find_max_gap(scores, gap_opening, gap_extension):

    max = 0
    gap = gap_opening

    for score in reversed(scores):

        if (score - gap) >= max:
            max = score - gap

        gap += gap_extension

    return max

def generate_traceback(alignment_matrix, seq1, seq2):

    alignment1 = ""
    alignment2 = ""

    max = alignment_matrix.max()
    i, j = np.unravel_index(alignment_matrix.argmax(), alignment_matrix.shape)

    alignment1 += seq1[i]
    alignment2 += seq2[i]

    while



if __name__ == "__main__":

    scoring_matrix = read_scoring_matrix("BLOSUM50")

    seq1 = "HEAGAWGHEE"
    seq2 = "PAWHEAE"

    gap_opening = 1
    gap_extension = 1

    alignment_matrix = generate_alignment_matrix(seq1, seq2, gap_opening, gap_extension, scoring_matrix)
    print(alignment_matrix)
    generate_traceback(alignment_matrix)