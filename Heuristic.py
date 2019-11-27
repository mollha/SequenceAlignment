import helper_functions
import re
import numpy as np
import math

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

"""
1) Basic dynamic programming that runs in quadratic time and space [up to 50 marks].
2) Dynamic programming that runs in linear space [up to 65 marks].
3)A Heuristic procedure that runs in sub-quadratic time (similar to FASTA and BLAST) [up to 85 marks].
They should output three items: the score of the best local alignment found by the algorithm plus two lists of indices,
one for each input sequences, that realise the matches/mismatches in the alignment.
"""

def BandedSmithWaterman(seq1, seq2, scoring_matrix, alphabet, width, seed):

    def create_region_set(seed, seq1, seq2):
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

        # # Debug
        # for y in range(len(seq2)):
        #     for x in range(len(seq1)):
        #         if (x, y) in region:
        #             print("X", end=" ")
        #         else:
        #             print("-", end=" ")
        #     print("")

        return region

    def in_region(x, y):
        """
        Given x & y coordinates return T/F dependent on whether the location lies within the banded region.
        :param x:
        :param y:
        :return:
        """
        return (x, y) in region

    # Align 2 sequences
    def align():
        # Max score tracker
        max_score = -float('inf')
        max_index = [0, 0]  # init to 0,0 (0 score)

        # Iterate over scoring matrix and generate scoring (start at 1,1 and work from there) (O(n^2) method)
        for y in range(1, len(seq2) + 1):  # y -> seq2
            for x in range(1, len(seq1) + 1):  # x -> seq1
                # Check if current scoring cell lies in diagonal
                if in_region(x + seq1_offset, y + seq2_offset):
                    # If in diagonal, score as normal (cell not in diagonal all have score = 0)
                    vals = [
                        # seq2[y-1], seq1[x-1] as matrix has empty row & col at start
                        cost_matrix[y - 1][x - 1] + get_score(alphabet, scoring_matrix, seq2[y - 1], seq1[x - 1]),  # diagonal
                        cost_matrix[y - 1][x] + get_score(alphabet, scoring_matrix, seq2[y - 1], '_'),  # up
                        cost_matrix[y][x - 1] + get_score(alphabet, scoring_matrix, '_', seq1[x - 1]),  # left
                        0]  # 0 for local alignment
                    # Update scoring matrix
                    cost_matrix[y][x] = max(vals)
                    # Get index of max
                    index = vals.index(max(vals))
                    # Update backtrack matrix if score it come from is a valid cell
                    if index == 0 and in_region(x - 1 + seq1_offset, y - 1 + seq2_offset):
                        backtrack_matrix[y][x] = 'D'
                    elif index == 1 and in_region(x + seq1_offset, y - 1 + seq2_offset):
                        backtrack_matrix[y][x] = 'U'
                    elif index == 2 and in_region(x - 1 + seq1_offset, y + seq2_offset):
                        backtrack_matrix[y][x] = 'L'
                    # Check if new greatest score seen (score for vals outside diagonals score = 0)
                    if max(vals) > max_score:
                        max_score = max(vals)
                        max_index = [y, x]
                else:
                    # If cell doesn't lie in diagonal -> score = 0 (local alignment still)
                    cost_matrix[y][x] = 0

        # Return max score and best local alignment chars + indices
        return [max_score, helper_functions.backtrack(backtrack_matrix, max_index, seq1, seq2)]

    # Create region for checking whether inside/outside region
    region = create_region_set(seed, seq1, seq2)

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
    seq1 = seq1[seq1_start:seq1_end]
    seq2 = seq2[seq2_start:seq2_end]

    # Offset for checking coors against & providing correct indices on output
    seq1_offset = seq1_start
    seq2_offset = seq2_start

    # Setup cost matrix
    cost_matrix = helper_functions.create_cost_matrix(seq1, seq2)

    # Setup both backtrack and cost matrix initial values
    cost_matrix, backtrack_matrix = helper_functions.matrix_setup(cost_matrix, local=True)
    return align()


def heuralign(alphabet, scoring_matrix, sequence1, sequence2):

    def seed(word_length, seq1, seq2) -> list:
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
        for i in range(0, len(seq1) - word_length + 1):
            # Get substring of word length
            word = seq1[i:i + word_length]
            # Get start & end indexes of matches
            matches = [(m.start(0), m.end(0)) for m in re.finditer(word, seq2)]
            if matches:
                for match in matches:
                    # Store in format (seq1 - start, end), (seq2 - start, end)
                    words.append([(i, i + word_length), match])
        return words

    def get_diagonals(seeds):
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

    def get_threshold(seq1, seq2):
        """
        Iterate over scoring matrix and determine the threshold value
        :return: int, threshold hold value
        """
        # Calculate weighted average according to symbol count within the sequences
        char_freq = {i: (seq1 + seq2).count(i) / len(seq1 + seq2) for i in alphabet}

        # For all pairs, get scores and weight by char_freq
        pairs = [[i, j] for i in alphabet for j in alphabet]
        pair_weights = [char_freq[x] + char_freq[y] for x, y in pairs]

        # Expected score
        threshold = np.average([get_score(alphabet, scoring_matrix, x[0], x[1]) for x in pairs], weights=pair_weights)

        return threshold

    def extend(seq1, seq2, diagonals):
        """
        Given a dict of all seeds along the same diagonal, extend them if they lie upon the
        same diagonal and return the new start-end indices
        :param diagonals: dict, output of get_diagonals
        :return: 2d arr of tuples, [[(start1, end1), (start2, end2)], ...]
        """

        # Calculate the expected score of any given string
        threshold = get_threshold(seq1, seq2)

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
                        subseq1, subseq2 = seq1[seeds[i][0][1]:seeds[i + 1][0][0]], \
                                           seq2[seeds[i][1][1]:seeds[i + 1][1][0]]
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

    def get_best(diagonals, seq1, seq2):
        """
        Given the (now merged diagonals), get the best n% of seeds.
        :param diagonals: dict of seeds along same diagonals.
        :return:
        """
        to_return = math.ceil(len(diagonals) * 0.1)
        top_scores = []

        for diagonal_id in diagonals:
            for seed in diagonals[diagonal_id]:
                seq1, seq2 = seq1[seed[0][0]:seed[0][1]], seq2[seed[1][0]:seed[1][1]]
                score = 0
                for i in range(len(seq1)):
                    score += get_score(alphabet, scoring_matrix, seq1[i], seq2[i])
                if len(top_scores) < to_return:
                    top_scores.append([score, seed])
                elif score > top_scores[-1][0]:
                    del top_scores[-1]
                    top_scores.append([score, seed])
                # Sort in desc order by score
                top_scores = sorted(top_scores, key=lambda x: x[0], reverse=True)
        return [x[1] for x in top_scores]

    def banded_smith_waterman(best_seeds, word_length):
        """
        Run banded smith waterman with the width = avg distance between diagonals
        :param word_length:
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
            results = BandedSmithWaterman(sequence1, sequence2, scoring_matrix, alphabet, avg_width, seed)
            # print("Input Seed {0} | Output - {1}".format(seed, results))
            if results[0] > max_score:
                max_score = results[0]
                best_results = results

        return best_results

    def align(word_length, seq1, seq2):
        # TODO: what happens if cant find a single seed between 2 strings?
        # (ref. above) best local alignment is then best match of non-matching chars

        # --- 1) Seed sequences - find word sequences the sequences have in common ---
        seeds = []
        while not seeds:
            seeds = seed(word_length, seq1, seq2)
            if not seeds:
                word_length -= 1
                if word_length == 0:
                    print("Sequences contain NO matching characters!!")
                    return [0, [[], [], [], []]]
        # print("Seeds: {0}".format(seeds))

        # --- 2) Identify matching diagonals in seeds ---
        diagonals = get_diagonals(seeds)
        # print("Diagonals: {0}".format(diagonals))

        # --- 3) Greedily extend words along diagonals if score of subsequence > threshold ---
        diagonals = extend(seq1, seq2, diagonals)
        # --- 4) Get top n% of seeds ---
        best_seeds = get_best(diagonals, seq1, seq2)

        # --- 5) Run banded SmithWaterman (with width A) on best seed ---
        return banded_smith_waterman(best_seeds, word_length)

    results = align(word_length=2, seq1=sequence1, seq2=sequence2)
    return results[0], results[1][0], results[1][1], results[1][2], results[1][3]


string_1 ="BDABCBBDCAABDDCDABACDADDCBCABA"
string_2 ="CCAABBBDABACADBCCADCADAAACDCCA"
scoring_matrix = [[ 1,-5,-5,-5,-1],[-5, 1,-5,-5,-1],[-5,-5, 5,-5,-4],[-5,-5,-5, 6,-4],[-1,-1,-4,-4,-9]]
alphabet = "ABCD"

a = heuralign(alphabet, scoring_matrix, string_1, string_2)
print("Score:   ", a[0])
print("Indices: ", a[1],a[2])
score = check_score(alphabet + '_', scoring_matrix, string_1, string_2, a[1],a[2])
print('CHECKING SCORE: {} \n'.format(score))