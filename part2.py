def backtrack(backtrack_matrix, index, seq1, seq2):
    """
    Iterate over max_indexes and find all best local alignments
    :param backtrack_matrix: backtrack matrix
    :param index: arr of y, x coor of starting point for max alignment
    :return: 2d arr of all alignments
    """
    # Start at max index & backtrack
    out1_chars, out2_chars = [], []
    out1_indices, out2_indices = [], []
    # Stored as [y, x]
    x = index[1]
    y = index[0]

    # NB: seq1[x-1] & seq2[y-1] as backtrack & scoring matrix len = len(seq) + 1 (empty @ beginning)
    while backtrack_matrix[y][x] is not None and x >= 0 and y >= 0:
        if backtrack_matrix[y][x] == 'D':
            out1_chars.append(seq1[x-1])
            out2_chars.append(seq2[y-1])
            # Add indices of matches to output
            out1_indices.append(x-1)
            out2_indices.append(y-1)
            x -= 1
            y -= 1
        elif backtrack_matrix[y][x] == 'L':
            out1_chars.append(seq1[x - 1])
            out2_chars.append('-')
            x -= 1
        else:
            out1_chars.append('-')
            out2_chars.append(seq2[y - 1])
            y -= 1

    # Return alignment -> indices1, indices2, chars1, chars2
    return list(reversed(out1_indices)), list(reversed(out2_indices)),  list(reversed(out1_chars)), list(reversed(out2_chars))



def matrix_setup(cost_matrix, local, scoring_matrix=None, alphabet=None, seq1=None, seq2=None):
    """
    Given the cost_matrix, create the backtrack matrix & init the vals in both
    :param cost_matrix: n x m array
    :param local: bool, whether local alignment or not
    :return: n x m scoring matrix & n x m backtrack matrix (both w/ first row & col initialized)
    """
    # --- 1) Create backtrack matrix & initialize -> len + 1 as top left is blank ---
    backtrack_matrix = [[None for _ in range(len(cost_matrix[0]))] for _ in range(len(cost_matrix))]
    backtrack_matrix[0] = ['L' for _ in range(len(backtrack_matrix[0]))]
    for i in range(len(backtrack_matrix)):
        backtrack_matrix[i][0] = 'U'
    # Set 0,0 to None (always terminate here)
    backtrack_matrix[0][0] = None

    # --- 2) Initialize values in cost matrix ---
    # TODO: surely want to initialize values to the max(0, score(a,b)) rather than just 0?
    # If local alignment, init cost_matrix vals = 0
    if local:
        # Init first row & cols of matrices to 0
        cost_matrix[0] = [0 for _ in range(len(cost_matrix[0]))]
        for i in range(len(cost_matrix)):
            cost_matrix[i][0] = 0

    # If global, init cost matrix values using scoring matrix
    else:
        # Scoring function
        def score(a, b):
            # Get index in scoring matrix for chars a & b
            if a == '-':
                a_index = len(scoring_matrix[0]) - 1
            else:
                a_index = alphabet.index(a)
            if b == '-':
                b_index = len(scoring_matrix) - 1
            else:
                b_index = alphabet.index(b)
            return scoring_matrix[b_index][a_index]

        # 0,0 has score set to 0
        cost_matrix[0][0] = 0
        # Init first row and col of cost matrix using score function
        for i in range(len(seq1)):  # init 1st row
            cost_matrix[0][i+1] = cost_matrix[0][i] + score(seq1[i], '-')

        for i in range(len(seq2)):  # init 1st col
            cost_matrix[i+1][0] = cost_matrix[i][0] + score(seq2[i], '-')

    return cost_matrix, backtrack_matrix


def backtrack(backtrack_matrix, index, seq1, seq2):
    """
    Iterate over max_indexes and find all best local alignments
    :param backtrack_matrix: backtrack matrix
    :param index: arr of y, x coor of starting point for max alignment
    :return: 2d arr of all alignments
    """
    # Start at max index & backtrack
    out1_chars, out2_chars = [], []
    out1_indices, out2_indices = [], []
    # Stored as [y, x]
    x = index[1]
    y = index[0]

    # NB: seq1[x-1] & seq2[y-1] as backtrack & scoring matrix len = len(seq) + 1 (empty @ beginning)
    while backtrack_matrix[y][x] is not None and x >= 0 and y >= 0:
        if backtrack_matrix[y][x] == 'D':
            out1_chars.append(seq1[x-1])
            out2_chars.append(seq2[y-1])
            # Add indices of matches to output
            out1_indices.append(x-1)
            out2_indices.append(y-1)
            x -= 1
            y -= 1
        elif backtrack_matrix[y][x] == 'L':
            out1_chars.append(seq1[x - 1])
            out2_chars.append('-')
            x -= 1
        else:
            out1_chars.append('-')
            out2_chars.append(seq2[y - 1])
            y -= 1

    # Return alignment -> indices1, indices2, chars1, chars2
    return list(reversed(out1_indices)), list(reversed(out2_indices)),  list(reversed(out1_chars)), list(reversed(out2_chars))





"""
1) Basic dynamic programming that runs in quadratic time and space [up to 50 marks].
2) Dynamic programming that runs in linear space [up to 65 marks].
3)A Heuristic procedure that runs in sub-quadratic time (similar to FASTA and BLAST) [up to 85 marks].
They should output three items: the score of the best local alignment found by the algorithm plus two lists of indices,
one for each input sequences, that realise the matches/mismatches in the alignment.
TODO: For second part & introducing new scoring rules -> 'Gap Penalty (special penalty for consecutive “-”)'
      https://www.site.uottawa.ca/~lucia/courses/5126-10/lecturenotes/03-05SequenceSimilarity.pdf (slide 31+)
"""

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
        self.cost_matrix = [[None for _ in range(len(seq1) + 1)] for _ in range(len(seq2) + 1)]

        # Setup both backtrack and scoring matrix for global alignment
        self.cost_matrix, self.backtrack_matrix = matrix_setup(self.cost_matrix, local=False,
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
        max_index = [len(self.backtrack_matrix) - 1, len(self.backtrack_matrix[0]) - 1]

        # Iterate over cost matrix and generate scoring (start at 1,1 and work from there) (O(n^2) method)
        for y in range(1, len(self.seq2) + 1):  # y -> seq2
            for x in range(1, len(self.seq1) + 1):  # x -> seq1
                vals = [
                    # seq[y-1], seq[x-1] as matrix has empty row & col at start
                    self.cost_matrix[y - 1][x - 1] + self.score(self.seq2[y - 1], self.seq1[x - 1]),  # diagonal
                    self.cost_matrix[y - 1][x] + self.score(self.seq2[y - 1], '-'),  # up
                    self.cost_matrix[y][x - 1] + self.score('-', self.seq1[x - 1]),  # left
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
        max_score = self.cost_matrix[len(self.cost_matrix) - 1][len(self.cost_matrix[0]) - 1]

        # Find all max alignments
        return [max_score, backtrack(self.backtrack_matrix, max_index, self.seq1, self.seq2)]


# ----------------
def hirschberg(seq1, seq2, scoring_matrix, alphabet):
# -------------


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
        self.cost_matrix = [[None for _ in range(len(seq1) + 1)] for _ in range(len(seq2) + 1)]

        # Setup both backtrack and cost matrix initial values
        self.cost_matrix, self.backtrack_matrix = matrix_setup(self.cost_matrix, local=True)

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
            prev_row[j] = max(0, prev_row[j - 1] + self.score('-', seq1[j - 1]))  # insert/left

        # Loop over remaining rows and calc scores
        for i in range(1, len(seq2) + 1):
            # Get first value in new row
            current_row[0] = max(0, prev_row[0] + self.score(seq2[i - 1], '-'))  # del/up

            # Evaluate each value in row
            for j in range(1, len(seq1) + 1):
                score_sub = prev_row[j - 1] + self.score(seq2[i - 1], seq1[j - 1])  # diagonal
                score_ins = current_row[j - 1] + self.score('-', seq1[j - 1])  # left
                score_del = prev_row[j] + self.score(seq2[i - 1], '-')  # up

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
        return self.align(self.seq1[min_index[1]:max_index[1]], self.seq2[min_index[0]:max_index[0]], min_index[1],
                          min_index[0])

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
                                                                            needleman_output[1][1], needleman_output[1][
                                                                                2], \
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
            max_score_left, aligned_1_left_indices, aligned_2_left_indices, aligned_1_left_chars, aligned_2_left_chars = self.align(
                seq1[:seq1_mid], seq2[:seq2_mid], seq1_offset, seq2_offset)
            max_score_right, aligned_1_right_indices, aligned_2_right_indices, aligned_1_right_chars, aligned_2_right_chars = self.align(
                seq1[seq1_mid:], seq2[seq2_mid:], seq1_offset + seq1_mid, seq2_offset + seq2_mid)

            # Add results of recursive calls to out vars
            out1_chars = aligned_1_left_chars + aligned_1_right_chars
            out2_chars = aligned_2_left_chars + aligned_2_right_chars
            out1_indices = aligned_1_left_indices + aligned_1_right_indices
            out2_indices = aligned_2_left_indices + aligned_2_right_indices
            max_score = max_score_left + max_score_right

        return max_score, out1_indices, out2_indices, out1_chars, out2_chars

def dynproglin(alphabet, scoring_matrix, sequence1, sequence2):
    HB = Hirschberg(sequence1, sequence2, scoring_matrix, alphabet)
    results = HB.run()
    return results[0], results[1], results[2], results[3], results[4]

a = dynproglin("ABCD", [[1, -5, -5, -5, -1], [-5, 1, -5, -5, -1], [-5, -5, 5, -5, -4], [-5, -5, -5, 6, -4], [-1, -1, -4, -4, -9]],
            "AABCCDCA", "CDCDDD")
print("Score:   ", a[0])
print("Indices: ", a[1],a[2])