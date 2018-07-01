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




def make_scoring_df(scoring_matrix_list):
    """
    scoring matrix to be turned into a dataframe
    :param scoring_matrix_list:
    :return: df of scoring matrix
    """

    headers = scoring_matrix_list[0]
    scoring_matrix_list = scoring_matrix_list[1:]

    return pd.DataFrame(scoring_matrix_list, columns=headers, index=headers)



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

