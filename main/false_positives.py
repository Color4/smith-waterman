import numpy as np
from main.utils import generate_alignment


def fpr(positive, negative, tpr):
    """
    calculate the false positive rate given a true positive
    :param positive: positive score values
    :param negative: false positive score values
    :param tpr: tpr rate
    :return: fpr rate
    """
    threshold = np.percentile(np.asarray(positive), 100-tpr)  # uses numpy percentile (ex 70th tpr is the 30th percentile)
    total_false_positives = sum(negative > threshold)

    return total_false_positives/len(negative)


def tpr(positive, negative, fpr):
    """
    same as fpr, but calculates true positive given false positive rates
    :param positive: alignment scores for true positives
    :param negative: alignment scores for false positives
    :param fpr: fpr rate
    :return: trp
    """
    threshold = np.percentile(np.asarray(negative), 100 - fpr)
    total_true_positives = sum(positive > threshold)

    return total_true_positives / len(positive)


def calculate_pos_rate(positive_pairs, negative_pairs, gap, extension, matrix, rate, normalize, type):

    """
    controls calculation of true and false positive rates. Generates alignment scores, and noramlizes if applicable,
    and then calculates either tpr or fpr
    :param positive_pairs: tuple of sequences
    :param negative_pairs: tuple of sequences
    :param gap: gap penalty
    :param extension: extension penalty
    :param matrix: scoring matrix
    :param rate: rate of either tpr or fpr
    :param normalize: whether scores should be normalized
    :param type: fpr or tpr method to be calculated
    :return: either the fpr or tpr
    """

    positive_scores = []
    negative_scores = []

    for positive, negative in zip(positive_pairs, negative_pairs):

        # calculates the smith-waterman alignment scores
        positive_scores.append(generate_alignment(positive, gap, extension, matrix, normalize)[2])
        negative_scores.append(generate_alignment(negative, gap, extension, matrix, normalize)[2])

    return type(positive_scores, negative_scores, rate)

