from random import choice
import helper_functions

"""
1) Basic dynamic programming that runs in quadratic time and space [up to 50 marks].
2) Dynamic programming that runs in linear space [up to 65 marks].
3)A Heuristic procedure that runs in sub-quadratic time (similar to FASTA and BLAST) [up to 85 marks].
They should output three items: the score of the best local alignment found by the algorithm plus two lists of indices,
one for each input sequences, that realise the matches/mismatches in the alignment.
"""

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
        print('seed: ', seed)
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
        # Max score tracker
        max_score = -float('inf')
        max_index = [0, 0]  # init to 0,0 (0 score)

        # Iterate over scoring matrix and generate scoring (start at 1,1 and work from there)
        for y in range(1, len(self.seq2) + 1):  # y -> seq2
            # Break flag -> set true once in banded region, then if out -> can break (as dont need to consider further)
            flag = False
            for x in range(1, len(self.seq1) + 1):  # x -> seq1
                # Check if current scoring cell lies in diagonal
                if self.in_region(x + self.seq1_offset, y + self.seq2_offset):
                    # Set flag (as reached banded region)
                    flag = True
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
                    # If flag = True, have passed over banded region and back into region outside of band -> can break
                    if flag:
                        break

        # Return max score and best local alignment chars + indices
        return [max_score, helper_functions.backtrack(self.backtrack_matrix, max_index, self.seq1, self.seq2)]


def heuralign(alphabet: str, scoring_matrix: list, seq_s: str, seq_t: str) -> tuple:
    band_width = 3

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
            while i >= 0 and j >= 0:
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

    #--- 1) Seed sequences - find word sequences the sequences have in common ---
    ktup, seeds = get_seeds(tune_ktup())

    # if there are NO seeds, then just choose the middle diagonal
    if not seeds:
        BSW = BandedSmithWaterman(seq_s, seq_t, scoring_matrix, alphabet, band_width, [(0, 0)])
        return BSW.align()

    diagonals = get_diagonals(seeds)
    diagonal_scores = []
    print(diagonals)
    for diagonal_key in diagonals:
        diagonal_scores.append(extend_diagonal(diagonals[diagonal_key]))
    diagonal_scores.sort(key=lambda x: x[0], reverse=True)
    diagonal_scores = diagonal_scores[0: min(3, len(diagonal_scores))]  # get top 3 diagonals
    diagonal_scores = [x[1] for x in diagonal_scores]

    highest_score = float('-inf')
    best_alignment = None

    for diagonal_line in diagonal_scores:
        first_seed = diagonal_line[0]
        BSW = BandedSmithWaterman(seq_s, seq_t, scoring_matrix, alphabet, band_width, [(first_seed[0], first_seed[0] + first_seed[3]), (first_seed[1], first_seed[1] + first_seed[3])])
        diff_score, alignments = BSW.align()
        print(BSW.align())
        print(alignments)
        alignment_s, alignment_t, _, _ = alignments
        if diff_score > highest_score:
            best_alignment = diff_score, alignment_s, alignment_t
    return best_alignment


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
        [1, 1, -5, -5, -1],
        [1, 1, -5, -5, -1],
        [-5, -5, 5, -5, -4],
        [-5, -5, -5, 6, -4],
        [-1, -1, -4, -4, -9]]
    x = 5
    sequence1 = "".join(choice(list(alphabet)) for i in range(x))
    sequence2 = "".join(choice(list(alphabet)) for i in range(x))
    sequence1 = "AACAAADAAAACAADAADAAA"
    sequence2 = "CDCDAAACCACACAAAADD"

    print("Starting:")
    # Strip to ensure no whitespace
    sequence1, sequence2 = sequence1.strip(), sequence2.strip()
    print("Seq 1 - {0} ".format(sequence1))
    print("Seq 2 - {0}".format(sequence2))
    print("------------")


    #  Part 3 - < O(n^2) heuristic procedure, similar to FASTA and BLAST (time)
    score, out1_indices, out2_indices = heuralign(alphabet, scoring_matrix, sequence1, sequence2)

    # Output - print results
    print("------------")
    print("Score: {0}".format(score))
    print("Indices: {0} {1}".format(out1_indices, out2_indices))
    score = check_score(alphabet + '_', scoring_matrix, sequence1, sequence2, out1_indices, out2_indices)
    print('CHECKING SCORE: {} \n'.format(score))