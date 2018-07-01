import numpy as np
from main.utils import generate_alignment


def fpr(positive, negative, tpr):
    threshold = np.percentile(np.asarray(positive), 100-tpr)
    total_false_positives = sum(negative > threshold)

    return total_false_positives/len(negative)


def tpr(positive, negative, fpr):
    threshold = np.percentile(np.asarray(negative), 100 - fpr)
    total_true_positives = sum(positive > threshold)

    return total_true_positives / len(positive)


def calculate_pos_rate(positive_pairs, negative_pairs, gap, extension, matrix, rate, normalize, type):

    positive_scores = []
    negative_scores = []

    for positive, negative in zip(positive_pairs, negative_pairs):

        positive_scores.append(generate_alignment(positive, gap, extension, matrix, normalize)[2])
        negative_scores.append(generate_alignment(negative, gap, extension, matrix, normalize)[2])

    return type(positive_scores, negative_scores, rate)

