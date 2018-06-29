# main method to generate smith-waterman alignments
# of amino acid sequence and assess different parameters
# 21 Feb 17
# author <christa.caggiano@ucsf.edu>


from .smith_waterman import read_scoring_matrix,  align
from .false_positives import pick_best_gap, pick_best_scoring, calculate_tpr
from Bio import SeqIO
import collections

def read_sequence(line):
    """
    uses biopython package to read fasta sequences
    :param line: line containing fasta file location
    :return: seq object
    """
    return SeqIO.read(line.split(" ")[0], "fasta"), SeqIO.read(line.split(" ")[1].strip("\n"), "fasta")


def generate_alignment_scores(sequence_pairs, scoring_matrix, gap_extension, gap_opening, ouput_name, normalize):
    """
    calls smith waterman alignment. Reads in amino acid sequences from fasta files and generates
    a file of alignment scores

    :param sequence_pairs: file containing locations of fasta files to be aligned
    :param scoring_matrix:
    :param gap_extension:
    :param gap_opening:
    :param ouput_name: file to output alignments
    :return: none
    """

    with open(sequence_pairs, "r") as f, open(ouput_name, "w") as o:
        for line in f:
            seq1, seq2 = read_sequence(line)
            alignment1, alignment2, score = align(seq1.seq, seq2.seq, gap_opening, gap_extension, scoring_matrix)
            if not normalize:
                print(seq1.id, seq2.id, score, file=o)
            else:
                print(seq1.id, seq2.id, score/min(len(seq1), len(seq2)), file=o)


def assess_gap_penalties(scoring_matrix_name):
    """
    call alignments for a range of gap opening/extension penalties using the BLOSUM50 scoring matrix
    generates a file of the alignment scores for each condition

    :return: none
    """

    scoring_matrix = read_scoring_matrix(scoring_matrix_name)
    for i in range(1, 20):
        for j in range(1, 5):
            for pair in ["Pospairs.txt", "Negpairs.txt"]:
                generate_alignment_scores(pair, scoring_matrix, j, i,
                                          pair.replace(".txt", str(j) + str(i) + "_aligned_scores.txt"), False)


def assess_scoring_matrices(gap_opening, gap_extension):
    """
    run smith-waterman alignment on a several different scoring matrices

    :param gap_opening: best performing gap opening penalty determined from assess_gap_penalties()
    :param gap_extension: best performing gap extension penalty
    :return: none
    """
    scoring_matrices = ["scoring/BLOSUM50", "scoring/BLOSUM62", "scoring/MATIO", "scoring/PAM100", "scoring/PAM250"]
    for matrix in scoring_matrices:
        matrix_df = read_scoring_matrix(matrix)
        for pair in ["Pospairs.txt", "Negpairs.txt"]:
            generate_alignment_scores(pair, matrix_df, gap_extension, gap_opening,
                                      str(matrix) + pair.replace(".txt", "output.txt"), False)


def generate_normalized_scores(scoring_matrix_name, gap_opening, gap_extension):
    """
    normalizes sequence alignment scores based on the largest sequence
    :param scoring_matrix_name: matrix to generate alignments
    :param gap_opening:
    :param gap_extension:
    :return: none
    """
    scoring_matrix = read_scoring_matrix(scoring_matrix_name)
    for pair in ["Pospairs.txt", "Negpairs.txt"]:
        generate_alignment_scores(pair, scoring_matrix, gap_extension, gap_opening,
                                  str(scoring_matrix_name) + pair.replace(".txt", "output.txt"), True)


def most_common_aa(sequence_pairs):
    """
    find the most common amino acid in the entire set of sequences
    :param sequence_pairs: file of sequence pairs
    :return: none
    """
    total_sequence = ""
    with open(sequence_pairs) as f:
        for line in f:
            seq1, seq2 = read_sequence(line)
            total_sequence += seq1.seq
            total_sequence += seq2.seq
    return collections.Counter(total_sequence).most_common()


def optimize_matrix(scoring_matrix, gap_extension, gap_opening):
    """
    optimize matrix by rewarding amino acids that are more common to the positive group
    :param scoring_matrix: scoring matrix to use
    :param gap_extension:
    :param gap_opening:
    :return: none
    """
    # neg_most_common = most_common_aa("Negpairs.txt")
    optimized = False
    matrix = read_scoring_matrix(scoring_matrix)
    pos_most_common = most_common_aa("Pospairs.txt")

    while not optimized:
        for aa in pos_most_common:
            matrix[aa[0]][aa[0]] *= 10

            generate_alignment_scores("Pospairs.txt", matrix, gap_extension, gap_opening, "pos_optimized", False)
            generate_alignment_scores("Negpairs.txt", matrix, gap_extension, gap_opening, "neg_optimized", False)

            print(calculate_tpr("pos_optimized", "neg_optimized"))
            if calculate_tpr("pos_optimized", "neg_optimized") == 4:
                optimized = True


# controls function calls

assess_gap_penalties("scoring/BLOSUM50")  # calculate alignments for a variety of  gap penalties
fpr, gap_opening, gap_extension = pick_best_gap()  # out of all alignments, pick best gap penalties

assess_scoring_matrices(int(gap_opening), int(gap_extension))  # pick the best scoring matrix out of those provided
pick_best_scoring()  # pick scoring matrix that produces best fpr

generate_normalized_scores("scoring/BLOSUM50", int(gap_opening), int(gap_extension))  # test normalization

optimize_matrix("scoring/MATIO", int(gap_opening), int(gap_extension))  # try to optimize matrix

