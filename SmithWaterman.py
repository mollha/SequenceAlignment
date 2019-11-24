# TODO functions need cleaning

# -------------------------- HELPER FUNCTIONS -----------------------------
def check_score(alphabet, scoring_matrix, seq_s, seq_t, alignment_s, alignment_t):
    score = 0
    align_seq_s, align_seq_t = '', ''
    longest_sequence = max(len(seq_s), len(seq_t))
    for i in range(longest_sequence):
        if i in alignment_s:
            align_seq_s += seq_s[i]
        else:
            align_seq_s += '_'
    for j in range(longest_sequence):
        if j in alignment_t:
            align_seq_t += seq_t[j]
        else:
            align_seq_t += '_'
    for index in range(longest_sequence):
        char_s, char_t = align_seq_s[index], align_seq_t[index]
        score += get_score(alphabet, scoring_matrix, char_s, char_t)

    return score

def get_score(alphabet: str, scoring_matrix: list, char_s: str, char_t: str):
    return scoring_matrix[alphabet.index(char_s)][alphabet.index(char_t)]


def needleman_wunsch(alphabet: str, scoring_matrix: list, seq_s: str, seq_t: str, shift: tuple):
    nw_matrix = []
    paths_matrix = []
    for i in range(0, len(seq_s) + 1):
        nw_matrix.append([])
        paths_matrix.append([])
        for j in range(0, len(seq_t) + 1):
            if not i and not j:
                val = 0
                path_val = 'R'
            elif not i:
                val = nw_matrix[i][j - 1] + get_score(alphabet, scoring_matrix, '_', seq_t[j - 1])
                path_val = 'L'
            elif not j:
                val = nw_matrix[i-1][j] + get_score(alphabet, scoring_matrix, seq_s[i - 1], '_')
                path_val = 'U'
            else:
                diag = nw_matrix[i-1][j - 1] + get_score(alphabet, scoring_matrix, seq_s[i - 1], seq_t[j - 1])
                up = nw_matrix[i-1][j] + get_score(alphabet, scoring_matrix, seq_s[i - 1], '_')
                left = nw_matrix[i][j - 1] + get_score(alphabet, scoring_matrix, '_', seq_t[j - 1])
                val = max(diag, up, left)
                path_val = ['D', 'U', 'L'][[diag, up, left].index(val)]
            nw_matrix[i].append(val)
            paths_matrix[i].append(path_val)
    #
    # print('NW MATRIX')
    # for row in nw_matrix:
    #     print(row)
    # print('\nPATH MATRIX')
    # for row in paths_matrix:
    #     print(row)

    # ---- BACKTRACK -----
    i, j = len(seq_s), len(seq_t)
    alignment_s, alignment_t = [], []

    while True:  # while we haven't yet had to restart
        path = paths_matrix[i][j]
        if path == 'D':
            # print('add to alignment s ', shift[0] + i - 1)
            # print('add to alignment t ', shift[1] + j - 1)
            alignment_s.append(shift[0] + i - 1)
            alignment_t.append(shift[1] + j - 1)
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


def hirschberg(alphabet: str, scoring_matrix: list, seq_s: str, seq_t: str, shift: tuple):
    # shift describes the number of shifts that need to be applied to the sequence indices
    if not seq_s or not seq_t:
        return [], []
    if len(seq_s) == 1 or len(seq_t) == 1:
        return needleman_wunsch(alphabet, scoring_matrix, seq_s, seq_t, shift)

    def get_last_row(seq_1, seq_2):
        row1, row2 = [], []
        for i in range(0, len(seq_1) + 1):
            if row2:
                row1 = row2
                row2 = []

            for j in range(0, len(seq_2) + 1):
                if not i and not j:
                    row1.append(0)
                elif not i:
                    row1.append(row1[j - 1] + get_score(alphabet, scoring_matrix, '_', seq_t[j - 1]))
                elif not j:
                    row2.append(row1[j] + get_score(alphabet, scoring_matrix, seq_s[i - 1], '_'))
                else:
                    # At this stage, we can assume we have completed row1 - we only append to row2
                    diag = row1[j - 1] + get_score(alphabet, scoring_matrix, seq_s[i - 1], seq_t[j - 1])
                    up = row1[j] + get_score(alphabet, scoring_matrix, seq_s[i - 1], '_')
                    left = row2[j - 1] + get_score(alphabet, scoring_matrix, '_', seq_t[j - 1])
                    val = max(diag, up, left)
                    row2.append(val)
        return row2

    if len(seq_s) == 1:
        combined_rows = get_last_row(seq_s, seq_t)
        max_index = combined_rows.index(max(combined_rows))
        return [shift[0]], [shift[1] + max_index]

    # split the first sequence in half
    last_row_1 = get_last_row(seq_s[0:len(seq_s) // 2], seq_t)
    last_row_2 = get_last_row(seq_s[len(seq_s) // 2:][::-1], seq_t[::-1])

    combined_rows = [x + y for x, y in zip(last_row_1, last_row_2[::-1])]  # sum the last rows
    # print('sequence s', seq_s)
    # print('sequence t', seq_t)
    #
    # print(combined_rows)
    max_index = combined_rows.index(max(combined_rows))

    h1_s_ind, h1_t_ind = hirschberg(alphabet, scoring_matrix, seq_s[0:len(seq_s) // 2], seq_t[0:max_index], shift)
    h2_s_ind, h2_t_ind = hirschberg(alphabet, scoring_matrix, seq_s[len(seq_s) // 2:], seq_t[max_index:],
               (shift[0] + (len(seq_s) // 2), shift[1] + max_index))
    return h1_s_ind + h2_s_ind, h1_t_ind + h2_t_ind


# b = hirschberg("AGTC_", [[2,-1,-1,-1,-2],[-1,2,-1,-1,-2],[-1,-1,2,-1,-2],[-1,-1,-1,2,-2],[-2,-2,-2,-2,0]], "AGTACGCA", "TATGC", (0,0))
# print(b)
# -----------------------------------------------------------------------


def dynprog(alphabet: str, scoring_matrix: list, seq_s: str, seq_t: str):
    alphabet += '_'
    high_score = -float('inf')
    max_indices = None

    # Compute scores
    values, paths = [], []
    for i in range(0, len(seq_s) + 1):
        values.append([])
        paths.append([])
        for j in range(0, len(seq_t)+1):
            path_val = ''

            if not i and not j:
                val = 0
                path_val = 'R'
            elif not i:
                val = max(0, values[i][j-1] + get_score(alphabet, scoring_matrix, '_', seq_t[j - 1]))
                if not val:
                    path_val = 'R'
                else:
                    path_val = 'L'
            elif not j:
                val = max(0, values[i-1][j] + get_score(alphabet, scoring_matrix, seq_s[i - 1], '_'))
                if not val:
                    path_val = 'R'
                else:
                    path_val = 'U'
            else:
                diag = values[i-1][j-1] + get_score(alphabet, scoring_matrix, seq_s[i - 1], seq_t[j - 1])
                up = values[i-1][j] + get_score(alphabet, scoring_matrix, seq_s[i - 1], '_')
                left = values[i][j-1] + get_score(alphabet, scoring_matrix, '_', seq_t[j - 1])
                val = max(0, diag, up, left)

                if val == diag:
                    path_val = 'D'
                elif val == up:
                    path_val = 'U'
                elif val == left:
                    path_val = 'L'
                elif val == 0:
                    path_val = 'R'

            values[i].append(val)
            paths[i].append(path_val)
            if val > high_score:
                high_score = val
                max_indices = (i, j)

    # ---- BACKTRACK -----
    i, j = max_indices
    alignment_s, alignment_t = [], []

    while True:       # while we haven't yet had to restart
        path = paths[i][j]
        if path == 'D':
            alignment_s.append(i-1)
            alignment_t.append(j-1)
            i, j, = i - 1, j - 1
        elif path == 'U':
            i = i - 1
        elif path == 'L':
            j = j - 1
        else:
            break
    alignment_s.reverse()
    alignment_t.reverse()
    return [high_score, alignment_s, alignment_t]


def dynproglin(alphabet: str, scoring_matrix: list, seq_s: str, seq_t: str):
    alphabet += '_'

    # Compute Smith-Waterman matrix scores using two columns at a time
    def scan_matrix(seq_1, seq_2):
        col1, col2 = [], []
        high_score = -float('inf')
        max_indices = None

        for j in range(0, len(seq_2) + 1):
            for i in range(0, len(seq_1) + 1):
                if not i and not j:
                    # Append 0 to col1
                    val = 0
                    col1.append(val)
                elif not i:
                    # Append value to col2
                    val = max(0, col1[i] + get_score(alphabet, scoring_matrix, '_', seq_t[j - 1]))
                    col2.append(val)
                elif not j:
                    # Append value to col1
                    val = max(0, col1[i - 1] + get_score(alphabet, scoring_matrix, seq_s[i - 1], '_'))
                    col1.append(val)
                else:
                    # At this stage, we can assume we have completed col1 - we only append to col2
                    diag = col1[i - 1] + get_score(alphabet, scoring_matrix, seq_s[i - 1], seq_t[j - 1])
                    up = col2[i - 1] + get_score(alphabet, scoring_matrix, seq_s[i - 1], '_')
                    left = col1[i] + get_score(alphabet, scoring_matrix, '_', seq_t[j - 1])
                    val = max(0, diag, up, left)
                    col2.append(val)

                if val > high_score:
                    high_score = val
                    max_indices = (i, j)
            if col2:
                # check col2 isn't empty first i.e. it isn't the first pass
                col1 = col2
                col2 = []
        return max_indices, high_score

    # -------------------------- Find the sub-sequences upon which to apply global alignment --------------------------
    (i_index, j_index), high_score = scan_matrix(seq_s, seq_t)
    seq_s, seq_t = seq_s[0:i_index][::-1], seq_t[0:j_index][::-1]
    s_add, t_add = i_index, j_index
    (i_index, j_index), _ = scan_matrix(seq_s, seq_t)
    seq_s, seq_t = seq_s[0:i_index][::-1], seq_t[0:j_index][::-1]
    s_add -= i_index        # s_add describes the value to add to indices of sequence s based on what has been removed
    t_add -= j_index        # s_add describes the value to add to indices of sequence t based on what has been removed
    # -----------------------------------------------------------------------------------------------------------------
    alignment_s, alignment_t = hirschberg(alphabet, scoring_matrix, seq_s, seq_t, (s_add, t_add))
    return [high_score, alignment_s, alignment_t]



# HIRSCHBERG TESTS
a = dynprog ("AGTC", [[2,-1,-1,-1,-2],[-1,2,-1,-1,-2],[-1,-1,2,-1,-2],[-1,-1,-1,2,-2],[-2,-2,-2,-2,0]], "AGTACGCA", "TATGC")
print("Score:   ", a[0])
print("Indices: ", a[1],a[2])

a = dynproglin ("AGTC", [[2,-1,-1,-1,-2],[-1,2,-1,-1,-2],[-1,-1,2,-1,-2],[-1,-1,-1,2,-2],[-2,-2,-2,-2,0]], "AGTACGCA", "TATGC")
print("Score:   ", a[0])
print("Indices: ", a[1],a[2])

#
# string1 = "AABBAACA"
# string2 = "AABBAACA"
# # a = dynprog ("ABC", [[1,-1,-2,-1],[-1,2,-4,-1],[-2,-4,3,-2],[-1,-1,-2,0]], "AABBAACA", "CBACCCBA")
#
# # a = dynprog ("ABC", [[1,-1,-2,-1],[-1,2,-4,-1],[-2,-4,3,-2],[-1,-1,-2,0]], string1, string2)
# # print("Score:   ", a[0])
# # print("Indices: ", a[1],a[2])
#
#
# # a = dynprog("ABCD", [[1, -5, -5, -5, -1], [-5, 1, -5, -5, -1], [-5, -5, 5, -5, -4], [-5, -5, -5, 6, -4], [-1, -1, -4, -4, -9]],
# #             "DDCDDCCCDCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCDDDCDADCDCDCDCD",
# #             "DDCDDCCCDCBCCCCDDDCDBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBDCDCDCDCD")
# # print("Score:   ", a[0])
# # print("Indices: ", a[1],a[2])
#
# a = dynproglin("ABCD", [[1, -5, -5, -5, -1], [-5, 1, -5, -5, -1], [-5, -5, 5, -5, -4], [-5, -5, -5, 6, -4], [-1, -1, -4, -4, -9]],
#             "AABCCDCA", "CADBCDBBDD")
#
# print('\n')
#
# a = dynprog("ABCD", [[1, -5, -5, -5, -1], [-5, 1, -5, -5, -1], [-5, -5, 5, -5, -4], [-5, -5, -5, 6, -4], [-1, -1, -4, -4, -9]],
#             "AABCCDCA", "CADBCDBBDD")
# print("Score:   ", a[0])
# print("Indices: ", a[1], a[2])
# print('\n')
#
# a = dynprog("ABCD", [[1, -5, -5, -5, -1], [-5, 1, -5, -5, -1], [-5, -5, 5, -5, -4], [-5, -5, -5, 6, -4], [-1, -1, -4, -4, -9]],
#             "ACDCCBAA", "DDBBDCBDAC")
# print("Score:   ", a[0])
# print("Indices: ", a[1],a[2])
#
# # a = dynprog("ABCD", [[1, -5, -5, -5, -1], [-5, 1, -5, -5, -1], [-5, -5, 5, -5, -4], [-5, -5, -5, 6, -4], [-1, -1, -4, -4, -9]],
# #             "ACDCCBAA", "CDCDDD")
# # print("Score:   ", a[0])
# # print("Indices: ", a[1],a[2])
#
#
# # a = dynproglin("ABCD", [[1, -5, -5, -5, -1], [-5, 1, -5, -5, -1], [-5, -5, 5, -5, -4], [-5, -5, -5, 6, -4], [-1, -1, -4, -4, -9]],
# #             "AACAAADAAAACAADAADAAA", "CDCDDD")
# # print("Score:   ", a[0])
# # print("Indices: ", a[1],a[2])
#
# # a = dynprog ("ABC", [[1,-1,-2,-1],[-1,2,-4,-1],[-2,-4,3,-2],[-1,-1,-2,0]], string1[::-1], string2[::-1])
# # print("Score:   ", a[0])
# # print("Indices: ", a[1],a[2])
#
# # a = dynproglin ("ABC", [[1,-1,-2,-1],[-1,2,-4,-1],[-2,-4,3,-2],[-1,-1,-2,0]], "AABBAACA", "CBACCCBA")
# # b = dynproglin ("ABC", [[1,-1,-2,-1],[-1,2,-4,-1],[-2,-4,3,-2],[-1,-1,-2,0]], "AABBAAC", "CBAB")
#
# # def heuralign(alphabet: str, scoring_matrix: list, seq_s: str, seq_t: str):
# #     # BLAST
# #
# #     indel = len(alphabet)
# #     alphabet += '_'
# #     threshold = 1  # threshold