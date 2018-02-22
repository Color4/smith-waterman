# smith waterman algorithm
# 20 Feb 2018
# author <Christa.caggiano@ucsf.edu>

import pandas as pd

def read_scoring_matrix(matrix_file):

    scoring_matrix = []

    with open(matrix_file) as matrix:
        for line in matrix:
            if "#" not in line:
                scoring_matrix.append(line.split())

    make_scoring_df(scoring_matrix)

def make_scoring_df(scoring_matrix_list):

    headers = scoring_matrix_list[0]
    scoring_matrix_list = scoring_matrix_list[1:]

    df = pd.DataFrame(scoring_matrix_list, columns=headers)

    df.index(headers)

    print(df)


if __name__ == "__main__":

    read_scoring_matrix("BLOSUM50")