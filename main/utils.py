from Bio import SeqIO
from main.smith_waterman import align
import pandas as pd


def read_scoring_matrix(matrix_file):
    """
    process the scoring matrix file
    :param matrix_file: location of a scoring matrix
    :return: return dataframe of matrix
    """
    with open(matrix_file, "r") as matrix:
        scoring_matrix = [line.split() for line in matrix if "#" not in line]
    return make_scoring_df(scoring_matrix)


def read_sequence_pairs(input):

    with open(input, "r") as input:
        pair_list = [read_sequence(tuple(line.split(" "))) for line in input]

    return pair_list


def read_sequence(pair):
    """
    uses biopython package to read fasta sequences
    :param line: line containing fasta file location
    :return: seq object
    """

    seq1, seq2 = pair
    seq2 = seq2.strip("\n")

    return SeqIO.read("../" + seq1, "fasta").seq, SeqIO.read("../" + seq2, "fasta").seq


def make_scoring_df(scoring_matrix_list):
    """
    scoring matrix to be turned into a dataframe
    :param scoring_matrix_list:
    :return: df of scoring matrix
    """

    headers = scoring_matrix_list[0]
    scoring_matrix_list = scoring_matrix_list[1:]

    return pd.DataFrame(scoring_matrix_list, columns=headers, index=headers)


def generate_alignment(pair, gap, extension, matrix, normalize):
    """
    takes two sequences, a gap/extension penalty, and a scoring matrix, and generates alignment
    helper function to abstract away from smith_waterman.py in other parts of the code
    :param pair: sequence pair
    :param gap: gap penalty
    :param extension: gap extension penalty
    :param matrix: scoring matrix
    :param normalize: whether the alignment score should be normalized by the shortest sequence
    :return: alignments, and the score, normalized if applicable
    """

    align1, align2, score = align(pair[0], pair[1], gap, extension, matrix)

    normalization_constant = min(len(pair[0]), len(pair[1])) if normalize else 1

    return align1, align2, score/normalization_constant
