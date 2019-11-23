def dynprog(alphabet: str, scoring_matrix: list, seq_s: str, seq_t: str):
    indel = len(alphabet)
    alphabet += '-'
    high_score = -float('inf')
    max_indices = None

    # Compute scores
    values, paths = [], []
    for i in range(0, len(seq_s) + 1):
        values.append([])
        paths.append([])
        for j in range(0, len(seq_t)+1):
            index_s = alphabet.index(seq_s[i - 1])
            index_t = alphabet.index(seq_t[j - 1])
            path_val = ''

            if not i and not j:
                val = 0
                path_val = 'R'
            elif not i:
                val = max(0, values[i][j-1] + scoring_matrix[indel][index_t])
                if not val:
                    path_val = 'R'
                else:
                    path_val = 'L'
            elif not j:
                val = max(0, values[i-1][j] + scoring_matrix[index_s][indel])
                if not val:
                    path_val = 'R'
                else:
                    path_val = 'U'
            else:
                diag = values[i-1][j-1] + scoring_matrix[index_s][index_t]
                up = values[i-1][j] + scoring_matrix[index_s][indel]
                left = values[i][j-1] + scoring_matrix[indel][index_t]
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
    indel = len(alphabet)
    alphabet += '-'

    # Compute SW scores
    def scan_matrix(seq_1, seq_2):
        col1, col2 = [], []
        high_score = -float('inf')
        max_indices = None

        for j in range(0, len(seq_2) + 1):
            for i in range(0, len(seq_1) + 1):
                index_s = alphabet.index(seq_1[i - 1])
                index_t = alphabet.index(seq_2[j - 1])

                if not i and not j:
                    # Append 0 to col1
                    val = 0
                    col1.append(val)
                elif not i:
                    # Append value to col2
                    val = max(0, col1[i] + scoring_matrix[indel][index_t])
                    col2.append(val)
                elif not j:
                    # Append value to col1
                    val = max(0, col1[i - 1] + scoring_matrix[index_s][indel])
                    col1.append(val)
                else:
                    # At this stage, we can assume we have completed col1 - we only append to col2
                    diag = col1[i - 1] + scoring_matrix[index_s][index_t]
                    up = col2[i - 1] + scoring_matrix[index_s][indel]
                    left = col1[i] + scoring_matrix[indel][index_t]
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

    # TODO - need to add the start and end points to the indices alignment at an APPROPRIATE time
    # Alternatively, could add them, then sort at the end

    # Get the start and end-points of a local alignment
    # -------------------------- Find the sub-sequences upon which to apply global alignment --------------------------
    (i_index, j_index), high_score = scan_matrix(seq_s, seq_t)
    seq_s, seq_t = seq_s[0:i_index][::-1], seq_t[0:j_index][::-1]
    s_add, t_add = i_index, j_index
    (i_index, j_index), _ = scan_matrix(seq_s, seq_t)
    seq_s, seq_t = seq_s[0:i_index][::-1], seq_t[0:j_index][::-1]
    s_add -= i_index        # s_add describes the value to add to indices of sequence s based on what has been removed
    t_add -= j_index        # s_add describes the value to add to indices of sequence t based on what has been removed
    # -----------------------------------------------------------------------------------------------------------------

    def hirschberg(seq_s, seq_t):

        def get_last_row(seq_1, seq_2):
            row1, row2 = [], []

            for i in range(0, len(seq_1) + 1):
                for j in range(0, len(seq_2) + 1):
                    index_s = alphabet.index(seq_1[i - 1])
                    index_t = alphabet.index(seq_2[j - 1])

                    if not i and not j:
                        # Append 0 to col1
                        val = 0
                        row1.append(val)
                    elif not i:
                        # Append value to col2
                        val = max(row1[i] + scoring_matrix[indel][index_t])
                        row2.append(val)
                    elif not j:
                        # Append value to col1
                        val = max(row1[i - 1] + scoring_matrix[index_s][indel])
                        row1.append(val)
                    else:
                        # At this stage, we can assume we have completed col1 - we only append to col2
                        diag = row1[i - 1] + scoring_matrix[index_s][index_t]
                        up = row2[i - 1] + scoring_matrix[index_s][indel]
                        left = row1[i] + scoring_matrix[indel][index_t]
                        val = max(diag, up, left)
                        row2.append(val)
                if row2:
                    # check row2 isn't empty first i.e. it isn't the first pass
                    row1 = row2
                    row2 = []
            return row2

        # split the first sequence in half
        last_row_1 = get_last_row(seq_s[0:len(seq_s) // 2], seq_t)
        last_row_2 = get_last_row(seq_s[len(seq_s) // 2:][::-1], seq_t[::-1])

        # sum the last rows
        combined_rows = [x + y for x, y in zip(last_row_1, last_row_2[::-1])]
        max_index = combined_rows.index(max(combined_rows))

        hirschberg(seq_s[0:len(seq_s) // 2], seq_t[0:max_index])
        hirschberg(seq_s[len(seq_s) // 2:], seq_t[max_index:])

string1 = "AABBAACA"
string2 = "AABBAACA"
# a = dynprog ("ABC", [[1,-1,-2,-1],[-1,2,-4,-1],[-2,-4,3,-2],[-1,-1,-2,0]], "AABBAACA", "CBACCCBA")

# a = dynprog ("ABC", [[1,-1,-2,-1],[-1,2,-4,-1],[-2,-4,3,-2],[-1,-1,-2,0]], string1, string2)
# print("Score:   ", a[0])
# print("Indices: ", a[1],a[2])


# a = dynprog("ABCD", [[1, -5, -5, -5, -1], [-5, 1, -5, -5, -1], [-5, -5, 5, -5, -4], [-5, -5, -5, 6, -4], [-1, -1, -4, -4, -9]],
#             "DDCDDCCCDCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCDDDCDADCDCDCDCD",
#             "DDCDDCCCDCBCCCCDDDCDBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBDCDCDCDCD")
# print("Score:   ", a[0])
# print("Indices: ", a[1],a[2])

a = dynproglin("ABCD", [[1, -5, -5, -5, -1], [-5, 1, -5, -5, -1], [-5, -5, 5, -5, -4], [-5, -5, -5, 6, -4], [-1, -1, -4, -4, -9]],
            "AABCCDCA", "CADBCDBBDD")

print('\n')

a = dynprog("ABCD", [[1, -5, -5, -5, -1], [-5, 1, -5, -5, -1], [-5, -5, 5, -5, -4], [-5, -5, -5, 6, -4], [-1, -1, -4, -4, -9]],
            "AABCCDCA", "CADBCDBBDD")
print("Score:   ", a[0])
print("Indices: ", a[1], a[2])
print('\n')

a = dynprog("ABCD", [[1, -5, -5, -5, -1], [-5, 1, -5, -5, -1], [-5, -5, 5, -5, -4], [-5, -5, -5, 6, -4], [-1, -1, -4, -4, -9]],
            "ACDCCBAA", "DDBBDCBDAC")
print("Score:   ", a[0])
print("Indices: ", a[1],a[2])

# a = dynprog("ABCD", [[1, -5, -5, -5, -1], [-5, 1, -5, -5, -1], [-5, -5, 5, -5, -4], [-5, -5, -5, 6, -4], [-1, -1, -4, -4, -9]],
#             "ACDCCBAA", "CDCDDD")
# print("Score:   ", a[0])
# print("Indices: ", a[1],a[2])


# a = dynproglin("ABCD", [[1, -5, -5, -5, -1], [-5, 1, -5, -5, -1], [-5, -5, 5, -5, -4], [-5, -5, -5, 6, -4], [-1, -1, -4, -4, -9]],
#             "AACAAADAAAACAADAADAAA", "CDCDDD")
# print("Score:   ", a[0])
# print("Indices: ", a[1],a[2])

# a = dynprog ("ABC", [[1,-1,-2,-1],[-1,2,-4,-1],[-2,-4,3,-2],[-1,-1,-2,0]], string1[::-1], string2[::-1])
# print("Score:   ", a[0])
# print("Indices: ", a[1],a[2])

# a = dynproglin ("ABC", [[1,-1,-2,-1],[-1,2,-4,-1],[-2,-4,3,-2],[-1,-1,-2,0]], "AABBAACA", "CBACCCBA")
# b = dynproglin ("ABC", [[1,-1,-2,-1],[-1,2,-4,-1],[-2,-4,3,-2],[-1,-1,-2,0]], "AABBAAC", "CBAB")

# def heuralign(alphabet: str, scoring_matrix: list, seq_s: str, seq_t: str):
#     # BLAST
#
#     indel = len(alphabet)
#     alphabet += '-'
#     threshold = 1  # threshold