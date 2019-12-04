import helper_functions
import re
import numpy as np
import math
from random import choice

"""
1) Basic dynamic programming that runs in quadratic time and space [up to 50 marks].
2) Dynamic programming that runs in linear space [up to 65 marks].
3)A Heuristic procedure that runs in sub-quadratic time (similar to FASTA and BLAST) [up to 85 marks].
They should output three items: the score of the best local alignment found by the algorithm plus two lists of indices,
one for each input sequences, that realise the matches/mismatches in the alignment.
"""

def get_score(alpha: str, scoring: list, char_s: str, char_t: str) -> int:
    alpha += '_'
    return scoring[alpha.index(char_s)][alpha.index(char_t)]

def check_score(alphabet, scoring_matrix, seq_s, seq_t, alignment_s, alignment_t):
    score = 0
    for i in range(alignment_s[0], alignment_s[-1]):
        if i not in alignment_s:
            score += get_score(alphabet, scoring_matrix, seq_s[i], '_')

    for i in range(alignment_t[0], alignment_t[-1]):
        if i not in alignment_t:
            score += get_score(alphabet, scoring_matrix, '_', seq_t[i])

    while alignment_s and alignment_t:
        score += get_score(alphabet, scoring_matrix, seq_s[alignment_s[0]], seq_t[alignment_t[0]])
        alignment_s = alignment_s[1:]
        alignment_t = alignment_t[1:]
    return score


class SmithWaterman:
    """
    1) Basic dynamic programming that runs in quadratic time and space [up to 50 marks].
     - for local alignment
    Implementation of the Smith-Waterman algorithm: https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm
    :param seq1: sequence of chars, str
    :param seq2: sequence of chars, str
    :return: 2 arr's of each chars alignment
    """

    def __init__(self, seq1, seq2, scoring_matrix, alphabet):
        # Scores of pairings
        self.scoring_matrix = scoring_matrix

        # Set of unique characters (same order as in scoring matrix)
        self.alphabet = alphabet

        # Sequences
        self.seq1 = seq1
        self.seq2 = seq2

        # Setup cost matrix
        self.cost_matrix = helper_functions.create_cost_matrix(seq1, seq2)

        # Setup both backtrack and cost matrix initial values
        self.cost_matrix, self.backtrack_matrix = helper_functions.matrix_setup(self.cost_matrix, local=True)

    # Scoring function
    def score(self, a, b):
        # Get index in scoring matrix for chars a & b
        if a == '-':
            a_index = len(self.scoring_matrix[0])-1
        else:
            a_index = self.alphabet.index(a)
        if b == '-':
            b_index = len(self.scoring_matrix)-1
        else:
            b_index = self.alphabet.index(b)

        # Return score from matrix
        return self.scoring_matrix[b_index][a_index]

    # Align 2 sequences
    def align(self):
        # Max score tracker
        max_score = -float('inf')
        max_index = []  # can have >1 greatest local alignment

        # Iterate over scoring matrix and generate scoring (start at 1,1 and work from there) (O(n^2) method)
        for y in range(1, len(self.seq2)+1):  # y -> seq2
            for x in range(1, len(self.seq1)+1):  # x -> seq1
                vals = [
                    # seq2[y-1], seq1[x-1] as matrix has empty row & col at start
                    self.cost_matrix[y-1][x-1] + self.score(self.seq2[y-1], self.seq1[x-1]),  # diagonal
                    self.cost_matrix[y-1][x] + self.score(self.seq2[y-1], '-'),  # up
                    self.cost_matrix[y][x - 1] + self.score('-', self.seq1[x-1]),  # left
                    0]  # 0 for local alignment
                # Update scoring matrix
                self.cost_matrix[y][x] = max(vals)
                # Get index of max
                index = vals.index(max(vals))
                # Update backtrack matrix
                if index == 0:
                    self.backtrack_matrix[y][x] = 'D'
                elif index == 1:
                    self.backtrack_matrix[y][x] = 'U'
                elif index == 2:
                    self.backtrack_matrix[y][x] = 'L'
                # Check if new greatest score seen
                if max(vals) > max_score:
                    max_score = max(vals)
                    max_index = [y, x]

        # Return max score and best local alignment chars + indices
        return [max_score, helper_functions.backtrack(self.backtrack_matrix, max_index, self.seq1, self.seq2)]


class NeedlemanWunsch():
    """
    NeedlmanWunsch algorithm for global alignment (used in Hirschberg))
    :param seq1: sequence of chars, str
    :param seq2: sequence of chars, str
    :return: 2 arr's of each chars alignment
    """

    def __init__(self, seq1, seq2, scoring_matrix, alphabet):
        # Calculating score of pairings
        self.scoring_matrix = scoring_matrix

        # Set of unique characters (same order as in scoring matrix)
        self.alphabet = alphabet

        # Sequences
        self.seq1 = seq1
        self.seq2 = seq2

        # Setup cost matrix
        self.cost_matrix = helper_functions.create_cost_matrix(seq1, seq2)

        # Setup both backtrack and scoring matrix for global alignment
        self.cost_matrix, self.backtrack_matrix = helper_functions.matrix_setup(self.cost_matrix, local=False,
                                                                                scoring_matrix=self.scoring_matrix,
                                                                                alphabet=self.alphabet,
                                                                                seq1=self.seq1,
                                                                                seq2=self.seq2)

    # Scoring function
    def score(self, a, b):
        # Get index in scoring matrix for chars a & b
        if a == '-':
            a_index = len(self.scoring_matrix[0]) - 1
        else:
            a_index = self.alphabet.index(a)
        if b == '-':
            b_index = len(self.scoring_matrix) - 1
        else:
            b_index = self.alphabet.index(b)

        # Return score from matrix
        return self.scoring_matrix[b_index][a_index]

    # Align 2 sequences
    def align(self):
        # Global align -> always at bottom right index
        max_index = [len(self.backtrack_matrix)-1, len(self.backtrack_matrix[0])-1]

        # Iterate over cost matrix and generate scoring (start at 1,1 and work from there) (O(n^2) method)
        for y in range(1, len(self.seq2)+1):  # y -> seq2
            for x in range(1, len(self.seq1)+1):  # x -> seq1
                vals = [
                    # seq[y-1], seq[x-1] as matrix has empty row & col at start
                    self.cost_matrix[y-1][x-1] + self.score(self.seq2[y-1], self.seq1[x-1]),  # diagonal
                    self.cost_matrix[y-1][x] + self.score(self.seq2[y-1], '-'),  # up
                    self.cost_matrix[y][x - 1] + self.score('-', self.seq1[x-1]),  # left
                ]
                # Update cost matrix
                self.cost_matrix[y][x] = max(vals)

                # Get index of max
                index = vals.index(max(vals))
                # Update backtrack matrix
                if index == 0:
                    self.backtrack_matrix[y][x] = 'D'
                elif index == 1:
                    self.backtrack_matrix[y][x] = 'U'
                elif index == 2:
                    self.backtrack_matrix[y][x] = 'L'

        # Score = value in bottom right cell (NW is global alignment)
        max_score = self.cost_matrix[len(self.cost_matrix)-1][len(self.cost_matrix[0])-1]

        # Find all max alignments
        return [max_score, helper_functions.backtrack(self.backtrack_matrix, max_index, self.seq1, self.seq2)]


class Hirschberg():
    """
    2) Dynamic programming that runs in linear space [up to 65 marks].
     - for local alignment
     Implementation of Hirschberg's algorithm: https://en.wikipedia.org/wiki/Hirschberg%27s_algorithm
    :param seq1: sequence of chars, str
    :param seq2: sequence of chars, str
    :return: 2d array of each chars alignment
    """

    def __init__(self, seq1, seq2, scoring_matrix, alphabet):
        # Calculating score of pairings
        self.scoring_matrix = scoring_matrix

        # Set of unique characters (same order as in scoring matrix)
        self.alphabet = alphabet

        # Sequences
        self.seq1 = seq1
        self.seq2 = seq2

        # Setup cost matrix
        self.cost_matrix = helper_functions.create_cost_matrix(seq1, seq2)

        # Setup both backtrack and cost matrix initial values
        self.cost_matrix, self.backtrack_matrix = helper_functions.matrix_setup(self.cost_matrix, local=True)

    # Scoring function
    def score(self, a, b):
        # Get index in scoring matrix for chars a & b
        if a == '-':
            a_index = len(self.scoring_matrix[0]) - 1
        else:
            a_index = self.alphabet.index(a)
        if b == '-':
            b_index = len(self.scoring_matrix) - 1
        else:
            b_index = self.alphabet.index(b)

        # Return score from matrix
        return self.scoring_matrix[b_index][a_index]

    # Last Row method (returns last row of scoring matrix - linear space complexity)
    def last_row(self, seq1, seq2):
        max_val, max_index = -float('inf'), [1, 1]

        # Init rows to 0s (as local alignment)
        prev_row = [0 for _ in range(len(seq1) + 1)]
        current_row = [0 for _ in range(len(seq1) + 1)]

        # Init first row
        for j in range(1, len(seq1) + 1):
            prev_row[j] = max(0, prev_row[j-1] + self.score('-', seq1[j-1]))  # insert/left

        # Loop over remaining rows and calc scores
        for i in range(1, len(seq2) + 1):
            # Get first value in new row
            current_row[0] = max(0, prev_row[0] + self.score(seq2[i-1], '-'))  # del/up

            # Evaluate each value in row
            for j in range(1, len(seq1) + 1):
                score_sub = prev_row[j - 1] + self.score(seq2[i - 1], seq1[j - 1])  # diagonal
                score_ins = current_row[j - 1] + self.score('-', seq1[j-1])  # left
                score_del = prev_row[j] + self.score(seq2[i-1], '-')  # up

                # Local alignment -> max(vals, 0)
                current_row[j] = max(0, score_sub, score_del, score_ins)

                # Update max_val / index if score > max
                if current_row[j] > max_val:
                    max_val = current_row[j]
                    max_index = [i, j]  # y, x

            # Update prev row & clear current row
            prev_row = current_row
            current_row = [0 for _ in range(len(seq1) + 1)]

        return prev_row, max_val, max_index

    # Calls recursive function with vals and returs
    def run(self):
        # Get max index from forward pass
        _, max_val, max_index = self.last_row(self.seq1, self.seq2)

        # Get min index from backward pass
        _, _, min_index = self.last_row(self.seq1[::-1], self.seq2[::-1])

        # Subtract lengths from min index (s.t. actual start position)
        min_index[1] = len(self.seq1) - min_index[1]
        min_index[0] = len(self.seq2) - min_index[0]

        # Get local alignment
        return self.align(self.seq1[min_index[1]:max_index[1]], self.seq2[min_index[0]:max_index[0]], min_index[1], min_index[0])

    # Hirschberg algorithm (ref. https://en.wikipedia.org/wiki/Hirschberg%27s_algorithm)
    def align(self, seq1, seq2, seq1_offset, seq2_offset):
        out1_chars, out2_chars = [], []
        out1_indices, out2_indices = [], []
        max_score = 0

        # Empty seq1
        if len(seq1) == 0:
            # Add -'s to remaining alignment s.t. valid
            for i in range(len(seq2)):
                out1_chars.append('-')
                out2_chars.append(seq2[i])
            # Score produced alignment
            for i in range(len(out1_chars)):
                max_score += self.score(out1_chars[i], out2_chars[i])

        # Empty seq2
        elif len(seq2) == 0:
            # Add -'s to remaining alignment s.t. valid
            for i in range(len(seq1)):
                out1_chars.append(seq1[i])
                out2_chars.append('-')
            # Score produced alignment
            for i in range(len(out1_chars)):
                max_score += self.score(out1_chars[i], out2_chars[i])

        # Apply SW for optimal local alignment
        elif len(seq1) == 1 or len(seq2) == 1:
            NW = NeedlemanWunsch(seq1, seq2, self.scoring_matrix, self.alphabet)
            needleman_output = NW.align()
            max_score, out1_indices, out2_indices, out1_chars, out2_chars = needleman_output[0], needleman_output[1][0], \
                                                                        needleman_output[1][1], needleman_output[1][2],\
                                                                        needleman_output[1][3]
            out1_indices = [x + seq1_offset for x in out1_indices]
            out2_indices = [x + seq2_offset for x in out2_indices]

        else:
            # Get midpoint of Seq2
            seq2_mid = len(seq2) // 2

            # Get scoring of lhs (in linear space)
            r_left, _, _ = self.last_row(seq1, seq2[:seq2_mid])

            # Get scoring of rhs (in linear space) [reversed]
            r_right, _, _ = self.last_row(seq1[::-1], seq2[seq2_mid:][::-1])
            r_right.reverse()  # flip back again for calculating total score

            # Sum values and find argmax
            row = [l + r for l, r in zip(r_left, r_right)]
            maxidx, maxval = max(enumerate(row), key=lambda a: a[1])

            # Partition seq1 at argmax
            seq1_mid = maxidx

            # Recursively call align on each half
            max_score_left, aligned_1_left_indices, aligned_2_left_indices, aligned_1_left_chars, aligned_2_left_chars = self.align(seq1[:seq1_mid], seq2[:seq2_mid], seq1_offset, seq2_offset)
            max_score_right, aligned_1_right_indices, aligned_2_right_indices, aligned_1_right_chars, aligned_2_right_chars = self.align(seq1[seq1_mid:], seq2[seq2_mid:], seq1_offset+seq1_mid, seq2_offset+seq2_mid)

            # Add results of recursive calls to out vars
            out1_chars = aligned_1_left_chars + aligned_1_right_chars
            out2_chars = aligned_2_left_chars + aligned_2_right_chars
            out1_indices = aligned_1_left_indices + aligned_1_right_indices
            out2_indices = aligned_2_left_indices + aligned_2_right_indices
            max_score = max_score_left + max_score_right

        return max_score, out1_indices, out2_indices, out1_chars, out2_chars


class BandedSmithWaterman:
    def __init__(self, seq1, seq2, scoring_matrix, alphabet, width, seed):
        # Scores of pairings
        self.scoring_matrix = scoring_matrix

        # Set of unique characters (same order as in scoring matrix)
        self.alphabet = alphabet

        # Width to explore
        self.width = width

        # Seed of interest (ungapped alignment (s1, e1). (s2, e2))
        self.seed = seed

        # Create region for checking whether inside/outside region
        self.region = self.create_region_set(seed, seq1, seq2)

        # Find start intercept with axis
        x, y = seed[0][0], seed[1][0]  # seq1 -> x, seq2 -> y
        x_intercept, y_intercept = x - min(x, y), y - min(x, y)

        # Banded region is width space away from cells on diagonal in dirs: left, right, up, down (that exist)
        # Get starts of seq1 & seq2 (leftmost cell and topmost cell)
        seq1_start = max(0, x_intercept - width)
        seq2_start = max(0, y_intercept - width)
        # Get ends of seq1 & seq2 (rightmost and bottommost cell)
        seq1_end = len(seq1)-seq2_start
        seq2_end = len(seq2)-seq1_start

        # Can crop sequences s.t. ignore search area outside of region produced by diagonal
        self.seq1 = seq1[seq1_start:seq1_end]
        self.seq2 = seq2[seq2_start:seq2_end]

        # Offset for checking coors against & providing correct indices on output
        self.seq1_offset = seq1_start
        self.seq2_offset = seq2_start

        # Setup cost matrix
        self.cost_matrix = helper_functions.create_cost_matrix(self.seq1, self.seq2)

        # Setup both backtrack and cost matrix initial values
        self.cost_matrix, self.backtrack_matrix = helper_functions.matrix_setup(self.cost_matrix, local=True)

    def create_region_set(self, seed, seq1, seq2):
        """
        Given the starting seed, iterate along its diagonal and store all coors that belong to the banded search
        in a set
        :param seed: arr of 2 tuples containing s1, e1 & s2, e2 indices
        :param seq1: sequence 1, str
        :param seq2: sequence 2, str
        :return: set of coordinates that are valid for the given seed and width
        """
        # Set to store all (x, y) coors of valid regions in the grid
        region = set()

        # Iterate up + left add all along diagonal + any width away to set
        x, y = seed[0][0], seed[1][0]
        while x >= 0 and y >= 0:
            for i in range(-self.width, self.width+1):
                for j in range(-self.width, self.width+1):
                    if 0 <= x+i < len(seq1) and 0 <= y+j < len(seq2):
                        region.add((x+i, y+j))
            x -= 1
            y -= 1
        # Iterate down + right add all along diagonal + any width away to set
        x, y = seed[0][0]+1, seed[1][0]+1

        while x < len(seq1) and y < len(seq2):
            for i in range(-self.width, self.width + 1):
                for j in range(-self.width, self.width + 1):
                    if 0 <= x + i < len(seq1) and 0 <= y + j < len(seq2):
                        region.add((x + i, y + j))
            x += 1
            y += 1

        # # Debug
        # for y in range(len(seq2)):
        #     for x in range(len(seq1)):
        #         if (x, y) in region:
        #             print("X", end=" ")
        #         else:
        #             print("-", end=" ")
        #     print("")

        return region

    def in_region(self, x, y):
        """
        Given x & y coordinates return T/F dependent on whether the location lies within the banded region.
        :param x:
        :param y:
        :return:
        """
        return (x, y) in self.region

    # Scoring function
    def score(self, a, b):
        # Get index in scoring matrix for chars a & b
        if a == '-':
            a_index = len(self.scoring_matrix[0])-1
        else:
            a_index = self.alphabet.index(a)
        if b == '-':
            b_index = len(self.scoring_matrix)-1
        else:
            b_index = self.alphabet.index(b)

        # Return score from matrix
        return self.scoring_matrix[b_index][a_index]

    # Align 2 sequences
    def align(self):
        # Max score tracker
        max_score = -float('inf')
        max_index = [0, 0]  # init to 0,0 (0 score)

        # Iterate over scoring matrix and generate scoring (start at 1,1 and work from there) (O(n^2) method)
        for y in range(1, len(self.seq2)+1):  # y -> seq2
            for x in range(1, len(self.seq1)+1):  # x -> seq1
                # Check if current scoring cell lies in diagonal
                if self.in_region(x + self.seq1_offset, y + self.seq2_offset):
                    # If in diagonal, score as normal (cell not in diagonal all have score = 0)
                    vals = [
                        # seq2[y-1], seq1[x-1] as matrix has empty row & col at start
                        self.cost_matrix[y - 1][x - 1] + self.score(self.seq2[y - 1], self.seq1[x - 1]),  # diagonal
                        self.cost_matrix[y - 1][x] + self.score(self.seq2[y - 1], '-'),  # up
                        self.cost_matrix[y][x - 1] + self.score('-', self.seq1[x - 1]),  # left
                        0]  # 0 for local alignment
                    # Update scoring matrix
                    self.cost_matrix[y][x] = max(vals)
                    # Get index of max
                    index = vals.index(max(vals))
                    # Update backtrack matrix if score it come from is a valid cell
                    if index == 0 and self.in_region(x-1 + self.seq1_offset, y-1 + self.seq2_offset):
                        self.backtrack_matrix[y][x] = 'D'
                    elif index == 1 and self.in_region(x + self.seq1_offset, y-1 + self.seq2_offset):
                        self.backtrack_matrix[y][x] = 'U'
                    elif index == 2 and self.in_region(x-1 + self.seq1_offset, y + self.seq2_offset):
                        self.backtrack_matrix[y][x] = 'L'
                    # Check if new greatest score seen (score for vals outside diagonals score = 0)
                    if max(vals) > max_score:
                        max_score = max(vals)
                        max_index = [y, x]
                else:
                    # If cell doesn't lie in diagonal -> score = 0 (local alignment still)
                    self.cost_matrix[y][x] = 0

        # Return max score and best local alignment chars + indices
        return [max_score, helper_functions.backtrack(self.backtrack_matrix, max_index, self.seq1, self.seq2)]



def dynprog(alphabet, scoring_matrix, sequence1, sequence2):
    SW = SmithWaterman(sequence1, sequence2, scoring_matrix, alphabet)
    results = SW.align()
    return results[0], results[1][0], results[1][1], results[1][2], results[1][3]


def dynproglin(alphabet, scoring_matrix, sequence1, sequence2):
    HB = Hirschberg(sequence1, sequence2, scoring_matrix, alphabet)
    results = HB.run()
    return results[0], results[1], results[2], results[3], results[4]


def heuralign(alphabet, scoring_matrix, sequence1, sequence2):
    FA = FASTA(sequence1, sequence2, scoring_matrix, alphabet)
    results = FA.align()
    return results[0], results[1][0], results[1][1], results[1][2], results[1][3]



if __name__ == "__main__":
    # # Debug input 1
    # alphabet = "ABC"
    # scoring_matrix = [[1, -1, -2, -1], [-1, 2, -4, -1], [-2, -4, 3, -2], [-1, -1, -2, 0]]
    # sequence1 = "AABBAACA"
    # sequence2 = "CBACCCBA"
    # # Debug input 2
    # alphabet = "ABCD"
    # scoring_matrix = [
    #         [ 1,-5,-5,-5,-1],
    #         [-5, 1,-5,-5,-1],
    #         [-5,-5, 5,-5,-4],
    #         [-5,-5,-5, 6,-4],
    #         [-1,-1,-4,-4,-9]]
    # sequence1 = "AAAAACCDDCCDDAAAAACC"
    # sequence2 = "CCAAADDAAAACCAAADDCCAAAA"
    # # Debug input 3
    # alphabet = "ABCD"
    # scoring_matrix = [
    #         [ 1,-5,-5,-5,-1],
    #         [-5, 1,-5,-5,-1],
    #         [-5,-5, 5,-5,-4],
    #         [-5,-5,-5, 6,-4],
    #         [-1,-1,-4,-4,-9]]
    sequence1 = "AACAAADAAAACAADAADAAA"
    sequence2 = "CDCDDD"    sequence1 = "AACAAADAAAACAADAADAAA"
    sequence2 = "CDCDDD"    sequence1 = "AACAAADAAAACAADAADAAA"
    sequence2 = "CDCDDD"
    # Debug input 4
    alphabet = "ABCD"
    scoring_matrix = [
            [ 1,-5,-5,-5,-1],
            [-5, 1,-5,-5,-1],
            [-5,-5, 5,-5,-4],
            [-5,-5,-5, 6,-4],
            [-1,-1,-4,-4,-9]]
    x = 10
    sequence1 = "".join(choice(list(alphabet)) for i in range(x))
    sequence2 = "".join(choice(list(alphabet)) for i in range(x))


    #  Part 3 - < O(n^2) heuristic procedure, similar to FASTA and BLAST (time)
    score, out1_indices, out2_indices, out1_chars, out2_chars = heuralign(alphabet, scoring_matrix, sequence1, sequence2)


    # Output - print results
    print("------------")
    print("Score: {0}".format(score))
    print("Indices: {0} | {1}".format(out1_indices, out2_indices))
    score = check_score(alphabet + '_', scoring_matrix, sequence1, sequence2, out1_indices, out2_indices)
    print('CHECKING SCORE: {} \n'.format(score))
