import helper_functions
import re
import numpy as np
import math

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

def bandedSW(alphabet, scoring_matrix, seq_s, seq_t, diagonal, search_width):
    # Crop sequences around the diagonal
    seq_s = seq_s[]
    # If diff = 2, then (0, 2), (1, 3), (2, 4) etc until either i or j is at its max
    # therefore (i = len(t) - diff, j = len(t))
    # if diff = 0, then (0, 0), (1, 1) etc until both i and j are at their max
    # If diff = -2, then (2, 0), (3, 1), (4, 2) etc until either i or j is at its max
    # therefore (i = len(s), j = len(s) - diff)


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
        seq1_end = len(seq1) - seq2_start
        seq2_end = len(seq2) - seq1_start

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
            for i in range(-self.width, self.width + 1):
                for j in range(-self.width, self.width + 1):
                    if 0 <= x + i < len(seq1) and 0 <= y + j < len(seq2):
                        region.add((x + i, y + j))
            x -= 1
            y -= 1
        # Iterate down + right add all along diagonal + any width away to set
        x, y = seed[0][0] + 1, seed[1][0] + 1

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
        if a == '_':
            a_index = len(self.scoring_matrix[0]) - 1
        else:
            a_index = self.alphabet.index(a)
        if b == '_':
            b_index = len(self.scoring_matrix) - 1
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
        for y in range(1, len(self.seq2) + 1):  # y -> seq2
            for x in range(1, len(self.seq1) + 1):  # x -> seq1
                # Check if current scoring cell lies in diagonal
                if self.in_region(x + self.seq1_offset, y + self.seq2_offset):
                    # If in diagonal, score as normal (cell not in diagonal all have score = 0)
                    vals = [
                        # seq2[y-1], seq1[x-1] as matrix has empty row & col at start
                        self.cost_matrix[y - 1][x - 1] + self.score(self.seq2[y - 1], self.seq1[x - 1]),  # diagonal
                        self.cost_matrix[y - 1][x] + self.score(self.seq2[y - 1], '_'),  # up
                        self.cost_matrix[y][x - 1] + self.score('_', self.seq1[x - 1]),  # left
                        0]  # 0 for local alignment
                    # Update scoring matrix
                    self.cost_matrix[y][x] = max(vals)
                    # Get index of max
                    index = vals.index(max(vals))
                    # Update backtrack matrix if score it come from is a valid cell
                    if index == 0 and self.in_region(x - 1 + self.seq1_offset, y - 1 + self.seq2_offset):
                        self.backtrack_matrix[y][x] = 'D'
                    elif index == 1 and self.in_region(x + self.seq1_offset, y - 1 + self.seq2_offset):
                        self.backtrack_matrix[y][x] = 'U'
                    elif index == 2 and self.in_region(x - 1 + self.seq1_offset, y + self.seq2_offset):
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


def heuralign(alphabet: str, scoring_matrix: list, seq_s: str, seq_t: str) -> list:
    # Word length
    word_length = 3

    # Merge distance
    merge_distance = 2

    """
        3) Heuristic procedure that runs in sub-quadratic time (similar to FASTA and BLAST) [up to 85 marks].
        - for local alignment
        :param seq1: sequence of chars, str
        :param seq2: sequence of chars, str
        :return: 2 arr's of each chars alignment
        """



    def get_threshold():
        """
        Iterate over scoring matrix and determine the threshold value
        :return: int, threshold hold value
        """
        # Calculate weighted average according to symbol count within the sequences
        char_freq = {i: (seq_s + seq_t).count(i) / len(seq_s + seq_t) for i in alphabet}

        # For all pairs, get scores and weight by char_freq
        pairs = [[i, j] for i in alphabet for j in alphabet]
        pair_weights = [char_freq[x] + char_freq[y] for x, y in pairs]

        # Expected score
        threshold = np.average([get_score(alphabet, scoring_matrix, x[0], x[1]) for x in pairs], weights=pair_weights)

        return threshold

    def extend(diagonals):
        """
        Given a dict of all seeds along the same diagonal, extend them if they lie upon the
        same diagonal and return the new start-end indices
        :param diagonals: dict, output of get_diagonals
        :return: 2d arr of tuples, [[(start1, end1), (start2, end2)], ...]
        """

        # Calculate the expected score of any given string
        threshold = get_threshold()

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
                    for i in range(len(seeds) - 1):
                        # Get subsequences between seeds & score it
                        subseq1, subseq2 = seq_s[seeds[i][0][1]:seeds[i + 1][0][0]], \
                                           seq_t[seeds[i][1][1]:seeds[i + 1][1][0]]
                        score = 0
                        for j in range(len(subseq1)):
                            score += get_score(alphabet, scoring_matrix, subseq1[j], subseq2[j])

                        # If score >= expected value -> merge
                        if score >= threshold:
                            s = [(seeds[i][0][0], seeds[i + 1][0][1]), (seeds[i][1][0], seeds[i + 1][1][1])]
                            if s not in seeds:
                                # print("Original {0} | {1}".format(seeds[i], seeds[i+1]))
                                # print("Merged: {0}".format(s))
                                new_seeds.append(s)
                                flag = False
                        else:
                            new_seeds.append(seeds[i])
                            if i == len(seeds) - 2:  # Add final seed into pool even if dont match
                                new_seeds.append(seeds[i + 1])

                # Update seeds to contain merges
                diagonals[diagonal_id] = new_seeds

            # If no merges have occured -> break
            if flag:
                break

        return diagonals

    def get_best(diagonals):
        """
        Given the (now merged diagonals), get the best n% of seeds.
        :param diagonals: dict of seeds along same diagonals.
        :return:
        """
        to_return = math.ceil(len(diagonals) * 0.1)
        top_scores = []

        for diagonal_id in diagonals:
            for seed in diagonals[diagonal_id]:
                seq1, seq2 = seq_s[seed[0][0]:seed[0][1]], seq_t[seed[1][0]:seed[1][1]]
                score = 0
                for i in range(len(seq1)):
                    score += get_score(alphabet, scoring_matrix, seq_s[i], seq_t[i])
                if len(top_scores) < to_return:
                    top_scores.append([score, seed])
                elif score > top_scores[-1][0]:
                    del top_scores[-1]
                    top_scores.append([score, seed])
                # Sort in desc order by score
                top_scores = sorted(top_scores, key=lambda x: x[0], reverse=True)
        return [x[1] for x in top_scores]

    def banded_smith_waterman(best_seeds):
        """
        Run banded smith waterman with the width = avg distance between diagonals
        :param best_seeds: 2d arr of (start, end) indices for best seeds
        :return: results for FASTA
        """
        # --- 1) Calculate average distance between diagonals as banded region to search in --
        if len(best_seeds) == 1:
            print("Only 1 Seed found. Setting banded region = word length.")
            avg_width = word_length
        else:
            avg_width = []
            for i in range(len(best_seeds) - 1):
                i_diag = best_seeds[i][0][0] - best_seeds[i][1][0]
                for j in range(i + 1, len(best_seeds)):
                    j_diag = best_seeds[j][0][0] - best_seeds[j][1][0]
                    avg_width.append(abs(i_diag - j_diag))
            avg_width = round(sum(avg_width) / len(avg_width))
        # print("Average Distance between best scoring seeds (diagonally) {0}".format(avg_width))

        # --- 2) Run BandedSmithWaterman on each pair using the found average distance between diagonals ---
        max_score = -float('inf')
        best_results = None
        for seed in best_seeds:
            BSW = BandedSmithWaterman(seq_s, seq_t, scoring_matrix, alphabet, avg_width, seed)
            results = BSW.align()
            # print("Input Seed {0} | Output - {1}".format(seed, results))
            if results[0] > max_score:
                max_score = results[0]
                best_results = results

        return best_results

    # TODO: what happens if cant find a single seed between 2 strings?
    # (ref. above) best local alignment is then best match of non-matching chars

    def tune_ktup() -> int:
        min_length = min(len(seq_s), len(seq_t))
        # TODO tune ktup somehow
        return 3

    def get_seeds(ktup_val: int) -> tuple:
        while True:
            seed_dictionary = {}
            seed_list = []
            for index in range(len(seq_s) - ktup_val + 1):
                subword = seq_s[index: index + ktup_val]
                if subword in seed_dictionary:
                    seed_dictionary[subword].append(index)
                else:
                    seed_dictionary[subword] = [index]
            for index in range(len(seq_s) - ktup_val + 1):
                subword = seq_t[index: index + ktup_val]
                score_list = [get_score(alphabet, scoring_matrix, char, char) for char in subword]
                seed_score = sum(score_list)

                if subword in seed_dictionary:
                    seed_list += [(s_index, index, seed_score) for s_index in seed_dictionary[subword]]
            if not seed_list and ktup_val > 1:
                ktup_val -= 1
            else:
                break
        return ktup_val, seed_list

    def get_diagonals(seed_list: list) -> dict:
        seed_diagonals = {}
        for seed in seed_list:
            # seeds are of length ktup, and are in the form (s_index, t_index, score)
            difference = seed[1] - seed[0]    # t_index - s_index
            if difference in seed_diagonals:
                seed_diagonals[difference].append(seed)
            else:
                seed_diagonals[difference] = [seed]
        return seed_diagonals

    def extend_diagonal(diagonal_seeds: list) -> tuple:
        # extend diagonals and merge seeds if necessary
        # seeds start at length ktup, and begin in the form ((s_start, t_start), score)
        # they end as length >= ktup, and are translated to the form ((s_start, t_start), score, length)
        extended_seeds = set()
        total_score = 0
        for seed_tuple in diagonal_seeds:
            s_start, t_start, seed_score = seed_tuple
            i, j = s_start - 1, t_start - 1
            count = 0

            # extend from top
            while i and j:
                addition = get_score(alphabet, scoring_matrix, seq_s[i], seq_t[j])
                if addition < 0:
                    break
                else:
                    seed_score += addition
                    count += 1
                    i, j = i - 1, j - 1
            i, j = s_start + ktup, t_start + ktup
            s_start, t_start = s_start - count, t_start - count
            length = ktup + count
            count = 0
            # extend from bottom
            while i < len(seq_s) and j < len(seq_t):
                addition = get_score(alphabet, scoring_matrix, seq_s[i], seq_t[j])
                if addition < 0:
                    break
                else:
                    seed_score += addition
                    count += 1
                    i, j = i + 1, j + 1
            length += count

            clone_set = extended_seeds.copy()

            # # --- RESOLVE OVERLAP ---
            for existing_seed in extended_seeds:
                lower, upper = existing_seed[0], existing_seed[0] + length
                if (lower <= s_start <= upper) or (lower <= s_start + length <= upper):
                    # Overlap exists
                    if existing_seed[2] < seed_score:
                        total_score -= existing_seed[2]
                        clone_set.remove(existing_seed)
            extended_seeds = clone_set
            extended_seeds.add((s_start, t_start, seed_score, length))
            total_score += seed_score
        return total_score, list(extended_seeds)


    # promesses


    #--- 1) Seed sequences - find word sequences the sequences have in common ---
    ktup, seeds = get_seeds(tune_ktup())

    # if there are NO seeds, then just choose the middle diagonal
    if seeds:
        diagonals = get_diagonals(seeds)
        diagonal_scores = []
        for diagonal_key in diagonals:
            diagonal_scores.append(extend_diagonal(diagonals[diagonal_key]))
        diagonal_scores.sort(key=lambda x: x[0], reverse=True)
        diagonal_scores = diagonal_scores[0: min(3, len(diagonal_scores))]  # get top 3 diagonals
        print(diagonal_scores)
    else:
        print("Sequences contain NO matching characters!!")
        # middle diagonal
    print('HI')




    # seeds = []
    # while not seeds:
    #     """
    #     Given 2 sequences and the word_length param, find the start and ends for all matching subwords in seq2 that
    #     are also in seq1 to make list of indexes of matching subwords
    #     """
    #     # Found words
    #     words = []
    #     # Get all words of length word_length and check if present in seq2
    #     for i in range(0, len(seq_s) - word_length + 1):
    #         # Get substring of word length
    #         word = seq_t[i:i + word_length]
    #         # Get start & end indexes of matches
    #         matches = [(m.start(0), m.end(0)) for m in re.finditer(word, seq_t)]
    #         if matches:
    #             for match in matches:
    #                 # Store in format (seq1 - start, end), (seq2 - start, end)
    #                 words.append([(i, i + word_length), match])
    #     seeds = words
    #     if not seeds:
    #         word_length -= 1
    #         if word_length == 0:
    #             print("Sequences contain NO matching characters!!")
    #             return [0, [[], [], [], []]]
    # print("Seeds: {0}".format(seeds))

    # --- 2) Identify matching diagonals in seeds ---
    diagonals = get_diagonals(seeds)
    # print("Diagonals: {0}".format(diagonals))

    # --- 3) Greedily extend words along diagonals if score of subsequence > threshold ---
    diagonals = extend(diagonals)

    # --- 4) Get top n% of seeds ---
    best_seeds = get_best(diagonals)

    # --- 5) Run banded SmithWaterman (with width A) on best seed ---
    results = banded_smith_waterman(best_seeds)

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
    # sequence1 = "AACAAADAAAACAADAADAAA"
    # sequence2 = "CDCDDD"
    # Debug input 4
    alphabet = "ABCD"
    scoring_matrix = [
        [1, -5, -5, -5, -1],
        [-5, -1, -5, -5, -1],
        [-5, -5, 5, -5, -4],
        [-5, -5, -5, 6, -4],
        [-1, -1, -4, -4, -9]]
    sequence1 = "ABCCBABBCDBCDABAB"
    sequence2 = "BDDCDBBCDAABDBCAA"

    print("Starting:")
    # Strip to ensure no whitespace
    sequence1, sequence2 = sequence1.strip(), sequence2.strip()
    print("Seq 1 - {0} ".format(sequence1))
    print("Seq 2 - {0}".format(sequence2))
    print("------------")

    # TODO: improve runtime of scoring by using dicts!

    #  Part 3 - < O(n^2) heuristic procedure, similar to FASTA and BLAST (time)
    score, out1_indices, out2_indices, out1_chars, out2_chars = heuralign(alphabet, scoring_matrix, sequence1, sequence2)

    # Output - print results
    print("------------")
    print("Score: {0}".format(score))
    print("Indices: {0} | {1}".format(out1_indices, out2_indices))
    score = check_score(alphabet + '_', scoring_matrix, sequence1, sequence2, out1_indices, out2_indices)
    print('CHECKING SCORE: {} \n'.format(score))