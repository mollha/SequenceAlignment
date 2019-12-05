import helper_functions
import re
import numpy as np
import math

def backtrack(paths: list, max_indices: tuple) -> tuple:
    i, j = max_indices
    alignment_s, alignment_t = [], []

    while True:
        path = paths[i][j]
        if path == 'D':
            alignment_s.append(i - 1)
            alignment_t.append(j - 1)
            i, j, = i - 1, j - 1
        elif path == 'U':
            i = i - 1
        elif path == 'L':
            j = j - 1
        else:
            break
    alignment_s.reverse()
    alignment_t.reverse()
    return alignment_s, alignment_t

# def backtrack(backtrack_matrix, index, seq1, seq2):
#     """
#     Iterate over max_indexes and find all best local alignments
#     :param backtrack_matrix: backtrack matrix
#     :param index: arr of y, x coor of starting point for max alignment
#     :return: 2d arr of all alignments
#     """
#     # Start at max index & backtrack
#     out1_chars, out2_chars = [], []
#     out1_indices, out2_indices = [], []
#     # Stored as [y, x]
#     x = index[1]
#     y = index[0]
#
#     # NB: seq1[x-1] & seq2[y-1] as backtrack & scoring matrix len = len(seq) + 1 (empty @ beginning)
#     while backtrack_matrix[y][x] is not None and x >= 0 and y >= 0:
#         if backtrack_matrix[y][x] == 'D':
#             out1_chars.append(seq1[x-1])
#             out2_chars.append(seq2[y-1])
#             # Add indices of matches to output
#             out1_indices.append(x-1)
#             out2_indices.append(y-1)
#             x -= 1
#             y -= 1
#         elif backtrack_matrix[y][x] == 'L':
#             out1_chars.append(seq1[x - 1])
#             out2_chars.append('_')
#             x -= 1
#         else:
#             out1_chars.append('_')
#             out2_chars.append(seq2[y - 1])
#             y -= 1
#
#     # Return alignment -> indices1, indices2, chars1, chars2
#     return list(reversed(out1_indices)), list(reversed(out2_indices)),  list(reversed(out1_chars)), list(reversed(out2_chars))


def get_score(alpha: str, scoring: list, char_s: str, char_t: str) -> int:
    alpha += '_'
    return scoring[alpha.index(char_s)][alpha.index(char_t)]

def check_score(alphabet, scoring_matrix, seq_s, seq_t, alignment_s, alignment_t):
    score = 0
    sequence_s = ''
    sequence_t = ''
    for i in range(alignment_s[0], alignment_s[-1]):
        if i not in alignment_s:
            sequence_s += seq_s[i]
            sequence_t += '_'
            score += get_score(alphabet, scoring_matrix, seq_s[i], '_')

    for i in range(alignment_t[0], alignment_t[-1]):
        if i not in alignment_t:
            sequence_s += '_'
            sequence_t += seq_t[i]
            score += get_score(alphabet, scoring_matrix, '_', seq_t[i])

    while alignment_s and alignment_t:
        score += get_score(alphabet, scoring_matrix, seq_s[alignment_s[0]], seq_t[alignment_t[0]])
        sequence_s += seq_s[alignment_s[0]]
        sequence_t += seq_t[alignment_t[0]]
        alignment_s = alignment_s[1:]
        alignment_t = alignment_t[1:]
    print(sequence_s)
    print(sequence_t)
    return score

"""
1) Basic dynamic programming that runs in quadratic time and space [up to 50 marks].
2) Dynamic programming that runs in linear space [up to 65 marks].
3)A Heuristic procedure that runs in sub-quadratic time (similar to FASTA and BLAST) [up to 85 marks].
They should output three items: the score of the best local alignment found by the algorithm plus two lists of indices,
one for each input sequences, that realise the matches/mismatches in the alignment.
"""

def banded_SW(alphabet, scoring_matrix, seq1, seq2, width, seed):
    # Set to store all (x, y) coors of valid regions in the grid
    region = set()
    # Iterate up + left add all along diagonal + any width away to set
    x, y = seed[0][0], seed[1][0]
    while x >= 0 and y >= 0:
        for i in range(-width, width + 1):
            for j in range(-width, width + 1):
                if 0 <= x + i < len(seq1) and 0 <= y + j < len(seq2):
                    region.add((x + i, y + j))
        x -= 1
        y -= 1
    # Iterate down + right add all along diagonal + any width away to set
    x, y = seed[0][0] + 1, seed[1][0] + 1

    while x < len(seq1) and y < len(seq2):
        for i in range(-width, width + 1):
            for j in range(-width, width + 1):
                if 0 <= x + i < len(seq1) and 0 <= y + j < len(seq2):
                    region.add((x + i, y + j))
        x += 1
        y += 1

    # region has been created
    x, y = seed[0][0], seed[1][0]  # seq1 -> x, seq2 -> y
    x_intercept, y_intercept = x - min(x, y), y - min(x, y)
    # Banded region is width space away from cells on diagonal in dirs: left, right, up, down (that exist)
    # Get starts of seq1 & seq2 (leftmost cell and topmost cell)
    seq1_start = max(0, x_intercept - width)
    seq2_start = max(0, y_intercept - width)
    # Get ends of seq1 & seq2 (rightmost and bottommost cell)
    seq1_end = len(seq1) - seq2_start
    seq2_end = len(seq2) - seq1_start

    seq1 = seq1[seq1_start:seq1_end]
    seq2 = seq2[seq2_start:seq2_end]
    seq1_offset = seq1_start
    seq2_offset = seq2_start

    # initialises cost and backtrack (paths) matrix here

    values = [[0 for _ in range(len(seq1) + 1)] for _ in range(len(seq2) + 1)]
    backtrack_matrix = [['R' for _ in range(len(values[0]))] for _ in range(len(values))]
    backtrack_matrix[0] = ['L' for _ in range(len(backtrack_matrix[0]))]
    for i in range(len(backtrack_matrix)):
        backtrack_matrix[i][0] = 'U'
    # Set 0,0 to None (always terminate here)
    backtrack_matrix[0][0] = None
    values[0] = [0 for _ in range(len(values[0]))]
    for i in range(len(values)):
        values[i][0] = 0

    # Max score tracker
    max_score = -float('inf')
    max_index = [0, 0]  # init to 0,0 (0 score)

    # Iterate over scoring matrix and generate scoring (start at 1,1 and work from there)
    for y in range(1, len(seq2) + 1):  # y -> seq2
        # Break flag -> set true once in banded region, then if out -> can break (as dont need to consider further)
        flag = False
        for x in range(1, len(seq1) + 1):  # x -> seq1
            # Check if current scoring cell lies in diagonal
            if (x + seq1_offset, y + seq2_offset) in region:
                # Set flag (as reached banded region)
                flag = True
                # If in diagonal, score as normal (cell not in diagonal all have score = 0)
                vals = [
                    # seq2[y-1], seq1[x-1] as matrix has empty row & col at start
                    values[y - 1][x - 1] + get_score(alphabet, scoring_matrix, seq2[y - 1], seq1[x - 1]),  # diagonal
                    values[y - 1][x] + get_score(alphabet, scoring_matrix, seq2[y - 1], '_'),  # up
                    values[y][x - 1] + get_score(alphabet, scoring_matrix, '_', seq1[x - 1]),  # left
                    0]  # 0 for local alignment
                # Update scoring matrix
                values[y][x] = max(vals)
                # Get index of max
                index = vals.index(max(vals))
                # Update backtrack matrix if score it come from is a valid cell
                if index == 0 and (x - 1 + seq1_offset, y - 1 + seq2_offset) in region:
                    backtrack_matrix[y][x] = 'D'
                elif index == 1 and (x + seq1_offset, y - 1 + seq2_offset) in region:
                    backtrack_matrix[y][x] = 'U'
                elif index == 2 and (x - 1 + seq1_offset, y + seq2_offset) in region:
                    backtrack_matrix[y][x] = 'L'
                # Check if new greatest score seen (score for vals outside diagonals score = 0)
                if max(vals) > max_score:
                    max_score = max(vals)
                    max_index = [y, x]
            else:
                # If cell doesn't lie in diagonal -> score = 0 (local alignment still)
                values[y][x] = 0
                # If flag = True, have passed over banded region and back into region outside of band -> can break
                if flag:
                    break
    alignment_s, alignment_t = backtrack(backtrack_matrix, max_index)
    return max_score, alignment_s, alignment_t

class FASTA:
    """
    3) Heuristic procedure that runs in sub-quadratic time (similar to FASTA and BLAST) [up to 85 marks].
    - for local alignment
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

        # Word length
        self.word_length = max(3, math.ceil(min(len(self.seq1), len(self.seq2)) * 0.015))
        # self.word_length = 3

        # Setup cost matrix
        self.values = helper_functions.create_cost_matrix(seq1, seq2)

        # Setup both backtrack and cost matrix initial values
        self.values, self.backtrack_matrix = helper_functions.matrix_setup(self.values, local=True,
                                                                                scoring_matrix=self.scoring_matrix,
                                                                                alphabet=self.alphabet,
                                                                                seq1=self.seq1,
                                                                                seq2=self.seq2)

    def seed(self, seq1, seq2):
        """
        Given 2 sequences and the word_length param, find the start and ends for all matching subwords in seq2 that
        are also in seq1
        :param seq1: str, seq1
        :param seq2: str, seq2
        :return: list of indexes of matching subwords
        """
        # Found words
        words = []
        # Get all words of length word_length and check if present in seq2
        for i in range(0, len(seq1)-self.word_length+1):
            # Get substring of word length
            word = seq1[i:i+self.word_length]
            # Get start & end indexes of matches
            matches = [(m.start(0), m.end(0)) for m in re.finditer(word, seq2)]
            if matches:
                for match in matches:
                    # Store in format (seq1 - start, end), (seq2 - start, end)
                    words.append([(i, i+self.word_length), match])
        return words

    def get_diagonals(self, seeds):
        """
        Given the word indices, find ones that lie on the same diagonal.
        They lie on the same diagonal if the difference in start index is the same
        :param seeds: list of tuples of  matches in words
        :return:
        """
        # Store the difference as the key, and the indices as the values
        diagonals = {}
        # Iterate over all seeds
        for item in seeds:
            # Get difference in starting index
            diff = item[0][0] - item[1][0]
            # Add item to dict s.t. all items w/ same diff in starting index have same key in dict
            try:
                diagonals[diff].append(item)
            except KeyError:
                diagonals[diff] = [item]
        return diagonals

    def get_threshold(self):
        """
        Iterate over scoring matrix and determine the threshold value
        :return: int, threshold hold value
        """
        # Calculate weighted average according to symbol count within the sequences
        char_freq = {i: (self.seq1+self.seq2).count(i)/len(self.seq1+self.seq2) for i in self.alphabet}

        # For all pairs, get scores and weight by char_freq
        pairs = [[i, j] for i in self.alphabet for j in self.alphabet]
        pair_weights = [char_freq[x]+char_freq[y] for x, y in pairs]

        # Expected score
        threshold = np.average([get_score(alphabet, scoring_matrix, x[0], x[1]) for x in pairs], weights=pair_weights)

        return threshold

    def extend(self, diagonals):
        """
        Given a dict of all seeds along the same diagonal, extend them if they lie upon the
        same diagonal and return the new start-end indices
        :param diagonals: dict, output of get_diagonals
        :return: 2d arr of tuples, [[(start1, end1), (start2, end2)], ...]
        """

        # Calculate the expected score of any given string
        threshold = self.get_threshold()

        # Repeat until no more merges occur
        while True:
            # Exit condition -> no more merges have occured
            flag = True

            # Loop over every diagonal
            for diagonal_id in diagonals:

                # Get all subsequences that lie on that diagonal
                seeds = diagonals[diagonal_id]
                new_seeds = []
                # If len(seeds) == 1 [no option to merge]
                if len(seeds) == 1:
                    new_seeds = seeds
                else:
                    # Loop over all seeds on the same diagonal
                    for i in range(len(seeds)-1):
                        # Get subsequences between seeds & score it
                        subseq1, subseq2 = self.seq1[seeds[i][0][1]:seeds[i+1][0][0]], \
                                           self.seq2[seeds[i][1][1]:seeds[i+1][1][0]]
                        score = 0
                        for j in range(len(subseq1)):
                            score += get_score(alphabet, scoring_matrix, subseq1[j], subseq2[j])

                        # If score >= expected value -> merge
                        if score >= threshold:
                            s = [(seeds[i][0][0], seeds[i+1][0][1]), (seeds[i][1][0], seeds[i+1][1][1])]
                            if s not in seeds:
                                # print("Original {0} | {1}".format(seeds[i], seeds[i+1]))
                                # print("Merged: {0}".format(s))
                                new_seeds.append(s)
                                flag = False
                        else:
                            new_seeds.append(seeds[i])
                            if i == len(seeds)-2:  # Add final seed into pool even if dont match
                                new_seeds.append(seeds[i+1])

                # Update seeds to contain merges
                diagonals[diagonal_id] = new_seeds

            # If no merges have occured -> break
            if flag:
                break

        return diagonals

    def get_best(self, diagonals):
        """
        Given the (now merged diagonals), get the best n seeds.
        :param diagonals: dict of seeds along same diagonals.
        :return:
        """
        # Fixed val as otherwise time grows ++ w/ increase in sequence length
        to_return = 1
        top_scores = []

        # Gets top 'to_return' # scores
        for diagonal_id in diagonals:
            for seed in diagonals[diagonal_id]:
                seq1, seq2 = self.seq1[seed[0][0]:seed[0][1]], self.seq2[seed[1][0]:seed[1][1]]
                score = 0
                for i in range(len(seq1)):
                    score += get_score(alphabet, scoring_matrix,seq1[i], seq2[i])
                if len(top_scores) < to_return:
                    top_scores.append([score, seed])
                elif score > top_scores[-1][0]:
                    del top_scores[-1]
                    top_scores.append([score, seed])
                # Sort in desc order by score
                top_scores = sorted(top_scores, key=lambda x: x[0], reverse=True)
        return [x[1] for x in top_scores]

    def banded_smith_waterman(self, best_seeds):
        """
        Run banded smith waterman with the width = avg distance between diagonals
        :param best_seeds: 2d arr of (start, end) indices for best seeds
        :return: results for FASTA
        """
        width = 32

        # --- 2) Run BandedSmithWaterman on each pair using the found average distance between diagonals ---
        max_score = -float('inf')
        best_results = None
        for seed in best_seeds:
            # print("Running BSW...")
            results = banded_SW(alphabet, scoring_matrix, self.seq1, self.seq2, width, seed)
            # print("Input Seed {0} | Output - {1}".format(seed, results))
            if results[0] > max_score:
                max_score = results[0]
                best_results = results
        return best_results

    def align(self):
        # --- 1) Seed sequences - find word sequences the sequences have in common ---
        seeds = []
        while not seeds:
            seeds = self.seed(self.seq1, self.seq2)
            if not seeds:
                self.word_length -= 1
                if self.word_length == 0:
                    print("Sequences contain NO matching characters! Cannot seed!")
                    return [0, [[], [], [], []]]
        # print("Got seeds.")

        # --- 2) Identify matching diagonals in seeds ---
        diagonals = self.get_diagonals(seeds)
        # print("Got seeds along same diagonal.")

        # --- 3) Greedily extend words along diagonals if score of subsequence > threshold ---
        diagonals = self.extend(diagonals)
        # print("Seeds extended.")

        # --- 4) Get top n% of seeds ---
        best_seeds = self.get_best(diagonals)
        # print("Got best seeds.")

        # --- 5) Run banded SmithWaterman (with width A) on best seed ---
        return self.banded_smith_waterman(best_seeds)


def heuralign(alphabet, scoring_matrix, sequence1, sequence2):
    FA = FASTA(sequence1, sequence2, scoring_matrix, alphabet)
    return FA.align()


if __name__ == "__main__":
    alphabet = "ABCD"
    scoring_matrix = [
            [ 1,-5,-5,-5,-1],
            [-5, 1,-5,-5,-1],
            [-5,-5, 5,-5,-4],
            [-5,-5,-5, 6,-4],
            [-1,-1,-4,-4,-9]]
    sequence1 = "DDCDDCCCDCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCDDDCDADCDCDCDCD"
    sequence2 = "DDCDDCCCDCBCCCCDDDCDBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBDCDCDCDCD"

    # Strip to ensure no whitespace
    sequence1, sequence2 = sequence1.strip(), sequence2.strip()
    print("Seq 1 - {0} ".format(sequence1))
    print("Seq 2 - {0}".format(sequence2))
    print("------------")

    # Part 2 - O(n) dynamic prog. (space)
    score, out2_indices , out1_indices = heuralign(alphabet, scoring_matrix, sequence1, sequence2)

    # Output - print results
    print("Score: {0}".format(score))
    print("Indices: {0} | {1}".format(out1_indices, out2_indices))
    score = check_score(alphabet + '_', scoring_matrix, sequence1, sequence2, out1_indices, out2_indices)
    print('CHECKING SCORE: {} \n'.format(score))

