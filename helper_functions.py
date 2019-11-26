"""
Helper functions for rest of work
"""


def create_cost_matrix(seq1, seq2):
    # Generate cost matrix: x -> len(seq1)+1, y -> len(seq2)+1
    cost_matrix = [[None for _ in range(len(seq1) + 1)] for _ in range(len(seq2) + 1)]
    return cost_matrix


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


def alignment_pretty_print(out1, out2, seq1, seq2):
    """
    Given alignment, print in format:
        G C A T
        |   | |
        G T A T
    :param alignment: 2d arr of matches in alignment
    :return: Nothing
    """
    # Setup vars
    s1 = []
    hyphens = []
    s2 = []
    for i in range(len(out1)):
        s1.append(out1[i])
        s2.append(out2[i])
        if out1[i] == out2[i]:
            hyphens.append("|")
        else:
            hyphens.append(" ")
    # Print
    print(" ".join(s1))
    print(" ".join(hyphens))
    print(" ".join(s2))


def matrix_pretty_print(matrix, seq1, seq2):
    seq1 = "-" + seq1
    seq2 = "-" + seq2
    print("   ", end="")
    for item in seq1:
        print(item, end="    ")
    print()

    for i in range(len(matrix)):
        print(seq2[i], matrix[i])