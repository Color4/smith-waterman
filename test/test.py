import unittest
from main.smith_waterman import align
from main.utils import read_scoring_matrix

class TestSmithWaterman(unittest.TestCase):
    def test_same(self):
        matrix = read_scoring_matrix("scoring/BLOSUM50")
        score = align("G", "G", 1, 1, matrix)[-1]
        self.assertEqual(score, 0)

    def test_overlap(self):
        matrix = read_scoring_matrix("scoring/BLOSUM50")
        score = align("G", "GC", 1, 1, matrix)[-1]
        self.assertEqual(score, 0)
