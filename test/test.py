import unittest
from main.smith_waterman import align
from main.utils import read_scoring_matrix, generate_alignment
from main.false_positives import fpr

class TestSmithWaterman(unittest.TestCase):
    def test_same(self):
        matrix = read_scoring_matrix("../scoring/BLOSUM50")
        score = align("G", "G", 1, 1, matrix)[-1]
        self.assertEqual(score, 0)

    def test_overlap(self):
        matrix = read_scoring_matrix("../scoring/BLOSUM50")
        score = align("G", "GC", 1, 1, matrix)[-1]
        self.assertEqual(score, 0)

class TestFalsePostive(unittest.TestCase):

    def test_fpr(self):
        positive = [0, 1, 2, 3, 4, 5]
        negative = [0, 1, 2, 3, 4, 5]

        self.assertAlmostEqual(fpr(positive, negative, 70), 0.666666666667)

    def test_fpr_zeros(self):
        positive = [0, 1, 2, 3, 4, 5]
        negative = [0, 0, 0, 0, 0, 0]

        self.assertEqual(fpr(positive, negative, 70), 0.0)

    def test_fpr_low_percentile(self):
        positive = [0, 1, 2, 3, 4, 5]
        negative = [0, 1, 2, 3, 4, 5]

