import numpy as np
from itertools import product

# http://biorecipes.com/DynProgBasic/code.html
def dynprog(alphabet, scoring_matrix, s, t):

    global ALPHABET, SCORING_MATRIX

    ALPHABET = alphabet + "_"
    SCORING_MATRIX = np.array(scoring_matrix)

    D = build_score_matrix(s, t)




    score, s_align, t_align, s_matches, t_matches = traceback(D, s, t)

    return score, s_matches, t_matches, s_align, t_align

# https://en.wikipedia.org/wiki/Hirschberg's_algorithm
def dynproglin(alphabet, scoring_matrix, s, t):

    global ALPHABET, SCORING_MATRIX

    ALPHABET = alphabet + "_"
    SCORING_MATRIX = np.array(scoring_matrix)

    def recurse(s, t):

        s_align, t_align = "", ""

        if len(s) == 0 or len(t) == 0:

            s_align = ""
            t_align = ""

        elif len(s) == 1 or len(t) == 1:

            result = dynprog(alphabet, scoring_matrix, s, t)
            s_align = result[3]
            t_align = result[4]

        else:

            s_mid = int(len(s) / 2)
            score_l = build_score_matrix(s[:s_mid], t, sublinear=True)
            score_r = build_score_matrix(rev(s[s_mid:]), rev(t), sublinear=True)
            t_mid = np.argmax(score_l + rev(score_r))

            z_l, w_l = recurse(s[:s_mid], t[:t_mid])
            z_r, w_r = recurse(s[s_mid:], t[t_mid:])

            s_align = z_l + z_r
            t_align = w_l + w_r

        return s_align, t_align

    s_align, t_align = recurse(s, t)

    # So it's BACA, BACA. This doesn't exist.
    # If the scoring matrix is the same... let's just double-check again that it is. For my peace of mind, at least.
    # maybe that's part of the complexity?
    # Ah. It's a different matrix. On the second-last row, mind, which indicates to me it's an off-by-one.
    # The alignment, score, and matches are different. Of course it's all of them. But more interestingly, we're appended characters that shouldn't eist.

    print(s_align, t_align)

    score = align_score(s_align, t_align)
    s_matches, t_matches = get_alignment_indices(s_align, t_align)

    # No... it's more than that.
    # We're actually a row short.

    return score, s_matches, t_matches, s_align, t_align

# file:///C:/Users/user/Documents/MEGA/University/Year%203/CCS/Bioinformatics/HeuristicAlign.pdf
def heuralign(alphabet, scoring_matrix, s, t):

    # s (database), |s| >> |t| is the y-axis
    # t (query), is the x-axis

    global SCORING_MATRIX, ALPHABET

    SCORING_MATRIX = np.array(scoring_matrix)
    ALPHABET = alphabet + "_"

    ktup = 3
    lookup = {"".join(i): [] for i in product(alphabet, repeat=ktup)}

    for y in range(len(s) - ktup + 1):

        word = s[y:y + ktup]
        lookup[word].append(y)

    print(lookup)

    matches = []

    for x in range(len(t) - ktup + 1):

        word = t[x:x + ktup]

        print(word)

        for y in lookup[word]:

            extended = extend(s, t, (x, y, ktup))

            matches.append(extended)

    print(matches)

    diagonals = dict()
    k = 9

    for match in matches:

        x, y, length = match

        i = k * ((y - x) // k)

        if i not in diagonals:

            diagonals[i] = [match]

        else:

            diagonals[i].append(match)

    for key, value in diagonals.items():

        start_x = 0
        start_y = 0
        end_x_match = None
        end_y_match = None

        for match in value:

            if match[0][0] < start_x:

                start_x = match[0][0]

            if match[1][1] < start_y:

                start_y = match[1][1]

            if end_y_match is None or match[1] + match[2] > end_y_match:

                end_y_match = end_y_match

            if end_x_match is None or match[0] + match[2] > end_x_match:

                end_x_match = end_x_match

        end_x = end_x_match[0] + end_x_match[2]
        end_y = end_y_match[1] + end_y_match[2]

        s = s[start_y:end_y]
        t = t[start_x:end_x]

        D = build_score_matrix_banded(s, t, k)

        score, s_align, t_align, s_matches, t_matches = traceback(D, s, t)

        return score, s_matches, t_matches, s_align, t_align



def extend(s, t, match):

    x, y, length = match

    while x > 0 and y > 0:

        if cost(t[x - 1], s[y - 1]) > 0:

            x -= 1
            y -= 1
            length += 1

        else:

            break

    while x + length + 1 < len(t) and y + length + 1 < len(s):

        if cost(t[x + 1], s[y + 1]) > 0:

            length += 1

        else:

            break

    return x, y, length


def build_score_matrix_banded(s, t, k):

    size_y = len(s) + 1  # y-axis
    size_x = len(t) + 1  # x-axis
    D = np.zeros((size_y, size_x), dtype="int8")
    width = int(k / 2)

    for y in range(size_y):

        for x in range(y - width, y + width + 1):

            if 0 <= x < size_x and y > 0:

                D[y, x] = max(
                    0,  # For local alignment
                    match(D, y, x, s, y, t),  # The cost of matching two characters
                    insert_gap_into_s(D, y, x, t),  # The cost of matching a gap in s with a character in t
                    insert_gap_into_t(D, y, x, s, y)  # The cost of matching a gap in t with a character in s
                )

    return D

def get_alignment_indices(s_align, t_align):

    s_matches, t_matches = [], []
    s_point, t_point = 0, 0

    for i in range(len(s_align)):

        if s_align[i] != "_" and t_align[i] != "_":

            s_matches.append(s_point)
            t_matches.append(t_point)

            s_point += 1
            t_point += 1

        if s_align[i] != "_" and t_align[i] == "_":

            s_point += 1

        if s_align[i] == "_" and t_align[i] != "_":

            t_point += 1

    return s_matches, t_matches

def build_score_matrix(s, t, sublinear=False):

    size_y = len(s) + 1  # y-axis
    size_x = len(t) + 1  # x-axis
    shape = (2, size_x) if sublinear else (size_y, size_x)
    D = np.zeros(shape, dtype="int8")

    for y in range(1, size_y):  # y

        D_y = 1 if sublinear else y

        for x in range(1, size_x):

            D[D_y, x] = max(
                0,  # For local alignment
                match(D, D_y, x, s, y, t),  # The cost of matching two characters
                insert_gap_into_s(D, D_y, x, t),  # The cost of matching a gap in s with a character in t
                insert_gap_into_t(D, D_y, x, s, y)  # The cost of matching a gap in t with a character in s
            )

        if sublinear and y < size_y - 1: # Copy the 1st row onto the second row unless it's the final iteration

            D[0] = D[1].copy()
            D[1] = 0

    return D[1] if sublinear else D

def match(D, D_y, x, s, s_y, t):   # Conceptually D

    return D[D_y - 1][x - 1] + cost(s[s_y - 1], t[x - 1])

def insert_gap_into_s(D, D_y, x, t):  # s is the y-axis string: conceptually L

    return D[D_y][x - 1] + cost(t[x - 1], "_")

def insert_gap_into_t(D, D_y, x, s, s_y):  # t is the x-axis string: conceptually U

    return D[D_y - 1][x] + cost(s[s_y - 1], "_")

def traceback(D, s, t):

    score = np.amax(D)
    y, x = np.unravel_index(D.argmax(), D.shape)

    # We don't need the traceback, right?

    s_align, t_align = "", ""
    s_matches = []
    t_matches = []

    while y != 0 or x != 0:

        current = D[y][x]

        if current == 0:  # The end of the best local alignment

            return score, s_align, t_align, s_matches[::-1], t_matches[::-1]

        elif current == match(D, y, x, s, y, t): # D

            x -= 1
            y -= 1
            s_align = s[y] + s_align
            t_align = t[x] + t_align

            s_matches.append(y)
            t_matches.append(x)

        elif current == insert_gap_into_s(D, y, x, t):  # L

            x -= 1
            s_align = "_" + s_align
            t_align = t[x] + t_align


        elif current == insert_gap_into_t(D, y, x, s, y):  # U

            y -= 1
            s_align = s[y] + s_align
            t_align = "_" + t_align

        else:

            raise ValueError("Something's fucked!")

    return score, s_align, t_align, s_matches[::-1], t_matches[::-1]

def cost(c1, c2):

    global ALPHABET, SCORING_MATRIX

    print(ALPHABET)
    print(SCORING_MATRIX)

    return SCORING_MATRIX[ALPHABET.index(c1), ALPHABET.index(c2)]

def align_score(s_align, t_align):

    score = 0

    for i in range(len(s_align)):

        score += cost(s_align[i], t_align[i])

    return score

def rev(l):

    return l[::-1]


ALPHABET = "ABC_"
SCORING_MATRIX = np.array([
    [1, -1, -2, -1],
    [-1, 2, -4, -1],
    [-2, -4, 3, -2],
    [-1, -1, -2, 0]
])
s = "AABBAACA"
t = "CBACCCBA"

print(build_score_matrix(s, t))

print(build_score_matrix(s, t, sublinear=True))

print(heuralign(ALPHABET, SCORING_MATRIX,s, t))