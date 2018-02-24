import glob
import numpy as np


def calculate_percentile(array, percentile):
    return np.percentile(array, percentile)


def read_scores(file):
    with open(file) as f:
        pos_scores = []
        for line in f:
            pos_scores.append(line.split(" ")[2])
        pos_scores = np.asarray(pos_scores, dtype=float)

    return (pos_scores)


def get_gap(file_name):

    whole_numbers = list(range(20))
    whole_numbers = ', '.join(str(x) for x in whole_numbers)

    gap = ""

    for char in file_name:
        if char in whole_numbers:
            gap += char

    gap_extension = gap[0]
    gap_opening = gap[1:]

    return gap_extension, gap_opening


def pick_best_gap():

    pos_list = glob.glob("Pospairs*_aligned_scores.txt")

    with open("false_positive_rates.txt", "w") as o:
        print("gap extension", "gap_opening", "true_pos", "false_pos_num", "false_post_rate", file=o)
        all_false_pos_rates = []
        for file in pos_list:
            gap_extension, gap_opening = get_gap(file)

            pos_scores = read_scores(file)
            neg_scores = read_scores(file.replace("Pospairs", "Negpairs"))

            true_pos_val = calculate_percentile(pos_scores, 30)
            false_pos_rate = np.sum(neg_scores > true_pos_val)
            all_false_pos_rates.append((false_pos_rate/len(neg_scores), gap_opening, gap_extension))


            print(gap_extension, gap_opening, true_pos_val, false_pos_rate, false_pos_rate/len(neg_scores), file=o)

    return min(all_false_pos_rates)


def pick_best_scoring():

    pos_list = glob.glob("scoring/*pos_output.txt")
    with open("scoring_matrix_output.txt", "w") as o:
        for file in pos_list:

            pos_scores = read_scores(file)
            neg_scores = read_scores(file.replace("pos", "neg"))

            true_pos_score = list(range(0, 100, 10))
            false_pos_score = []

            for score in true_pos_score:
                true_pos_val = calculate_percentile(pos_scores, score)
                false_pos_score.append((np.sum(neg_scores > true_pos_val))/len(neg_scores))
            print(str(file), false_pos_score, file=o)


def calculate_tpr(pos_alignments, neg_alignments):

    false_pos_rates = [0, 10, 20, 30]
    neg_scores = read_scores(neg_alignments)
    pos_scores = read_scores(pos_alignments)

    num_passing = 0

    for rate in false_pos_rates:
        fpr = calculate_percentile(neg_scores, rate)
        if np.sum(pos_scores < fpr) == 0:
            num_passing += 1

    return num_passing