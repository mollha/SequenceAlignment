# ------------------------ HELPER FUNCTIONS --------------------------
def check_score(alphabet: str, scoring_matrix: list, seq_s: str, seq_t: str, alignment_s: list,
                alignment_t: list) -> int:
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


def last_row(alphabet: str, scoring_matrix: list, seq_1: str, seq_2: str) -> tuple:
    max_val, max_index = -float('inf'), [1, 1]

    # Init rows to 0s (as local alignment)
    prev_row = [0 for _ in range(len(seq_1) + 1)]
    current_row = [0 for _ in range(len(seq_1) + 1)]

    # Init first row
    for j in range(1, len(seq_1) + 1):
        prev_row[j] = max(0, prev_row[j - 1] + get_score(alphabet, scoring_matrix, seq_1[j - 1], '_'))

    # Loop over remaining rows and calc scores
    for i in range(1, len(seq_2) + 1):
        # Get first value in new row
        current_row[0] = max(0, prev_row[0] + get_score(alphabet, scoring_matrix, '_', seq_2[i - 1]))  # del/up

        # Evaluate each value in row
        for j in range(1, len(seq_1) + 1):
            score_sub = prev_row[j - 1] + get_score(alphabet, scoring_matrix, seq_1[j - 1], seq_2[i - 1])
            score_ins = current_row[j - 1] + get_score(alphabet, scoring_matrix, seq_1[j - 1], '_')
            score_del = prev_row[j] + get_score(alphabet, scoring_matrix, '_', seq_2[i - 1])

            # Local alignment -> max(vals, 0)
            current_row[j] = max(0, score_sub, score_del, score_ins)

            # Update max_val / index if score > max
            if current_row[j] > max_val:
                max_val = current_row[j]
                max_index = [i, j]  # y, x

        # Update prev row & clear current row
        prev_row = current_row
        current_row = [0 for _ in range(len(seq_1) + 1)]

    return prev_row, max_val, max_index


def get_score(alphabet: str, scoring_matrix: list, char_s: str, char_t: str) -> int:
    alphabet += '_'
    return scoring_matrix[alphabet.index(char_s)][alphabet.index(char_t)]


def backtrack(paths: list, max_indices: tuple) -> tuple:
    i, j = max_indices
    alignment_s, alignment_t = [], []

    while True:  # while we haven't yet had to restart
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


def needleman_wunsch(alphabet: str, scoring_matrix: list, seq_s: str, seq_t: str) -> tuple:
    values = []
    paths = []
    for i in range(0, len(seq_s) + 1):
        values.append([])
        paths.append([])
        for j in range(0, len(seq_t) + 1):
            if not i and not j:
                val = 0
                path_val = 'R'
            elif not i:
                val = values[i][j - 1] + get_score(alphabet, scoring_matrix, '_', seq_t[j - 1])
                path_val = 'L'
            elif not j:
                val = values[i - 1][j] + get_score(alphabet, scoring_matrix, seq_s[i - 1], '_')
                path_val = 'U'
            else:
                diag = values[i - 1][j - 1] + get_score(alphabet, scoring_matrix, seq_s[i - 1], seq_t[j - 1])
                up = values[i - 1][j] + get_score(alphabet, scoring_matrix, seq_s[i - 1], '_')
                left = values[i][j - 1] + get_score(alphabet, scoring_matrix, '_', seq_t[j - 1])
                val = max(diag, up, left)
                path_val = ['D', 'U', 'L'][[diag, up, left].index(val)]
            values[i].append(val)
            paths[i].append(path_val)
    return backtrack(paths, (len(seq_s), len(seq_t)))


def hirschberg(alphabet: str, scoring_matrix: list, seq_s: str, seq_t: str, offset: tuple) -> tuple:
    if len(seq_s) == 1 or len(seq_t) == 1:
        alignment_s, alignment_t = needleman_wunsch(alphabet, scoring_matrix, seq_s, seq_t)
        return [s + offset[0] for s in alignment_s], [t + offset[1] for t in alignment_t]
    else:
        last_row_1, _, _ = last_row(alphabet, scoring_matrix, seq_s, seq_t[:len(seq_t) // 2])
        last_row_2, _, _ = last_row(alphabet, scoring_matrix, seq_s[::-1], seq_t[len(seq_t) // 2:][::-1])
        last_row_2.reverse()
        index, _ = max(enumerate([x + y for x, y in zip(last_row_1, last_row_2)]), key=lambda x: x[1])

        alignment_s_half1, alignment_t_half1 = hirschberg(alphabet, scoring_matrix, seq_s[:index],
                                                          seq_t[:len(seq_t) // 2], offset)
        alignment_s_half2, alignment_t_half2 = hirschberg(alphabet, scoring_matrix, seq_s[index:],
                                                          seq_t[len(seq_t) // 2:],
                                                          (offset[0] + index, offset[1] + len(seq_t) // 2))

    return alignment_s_half1 + alignment_s_half2, alignment_t_half1 + alignment_t_half2


# ---------------------------------------------------------------------


def dynprog(alphabet: str, scoring_matrix: list, seq_s: str, seq_t: str):
    high_score = -float('inf')
    max_indices = None

    # Compute scores
    values, paths = [], []
    for i in range(0, len(seq_s) + 1):
        values.append([])
        paths.append([])
        for j in range(0, len(seq_t) + 1):
            path_val = ''

            if not i and not j:
                val = 0
                path_val = 'R'
            elif not i:
                val = max(0, values[i][j - 1] + get_score(alphabet, scoring_matrix, '_', seq_t[j - 1]))
                if not val:
                    path_val = 'R'
                else:
                    path_val = 'L'
            elif not j:
                val = max(0, values[i - 1][j] + get_score(alphabet, scoring_matrix, seq_s[i - 1], '_'))
                if not val:
                    path_val = 'R'
                else:
                    path_val = 'U'
            else:
                diag = values[i - 1][j - 1] + get_score(alphabet, scoring_matrix, seq_s[i - 1], seq_t[j - 1])
                up = values[i - 1][j] + get_score(alphabet, scoring_matrix, seq_s[i - 1], '_')
                left = values[i][j - 1] + get_score(alphabet, scoring_matrix, '_', seq_t[j - 1])
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

    alignment_s, alignment_t = backtrack(paths, max_indices)
    return [high_score, alignment_s, alignment_t]


def dynproglin(alphabet: str, scoring_matrix: list, seq_s: str, seq_t: str) -> list:
    _, high_score, max_index = last_row(alphabet, scoring_matrix, seq_s, seq_t)
    # Get min index from backward pass
    _, _, min_index = last_row(alphabet, scoring_matrix, seq_s[::-1], seq_t[::-1])

    # Subtract lengths from min index (s.t. actual start position)
    i, j = len(seq_t) - min_index[0], len(seq_s) - min_index[1]
    alignment_s, alignment_t = hirschberg(alphabet, scoring_matrix, seq_s[j:max_index[1]],
                                          seq_t[i:max_index[0]], (j, i))
    return [high_score, alignment_s, alignment_t]


# ---------------------------------------------------------------------------------------------------------------------

# TESTS
string_1, string_2 = "AACCDDAACC", "CADDACDDAA"
scoring_matrix = [[1, -5, -5, -5, -1], [-5, 1, -5, -5, -1], [-5, -5, 5, -5, -4], [-5, -5, -5, 6, -4],
                  [-1, -1, -4, -4, -9]]
alphabet = "ABCD"

a = dynprog(alphabet, scoring_matrix, string_1, string_2)
print("Score:   ", a[0])
print("Indices: ", a[1], a[2])
score = check_score(alphabet + '_', scoring_matrix, string_1, string_2, a[1], a[2])
print('CHECKING SCORE: {} \n'.format(score))
recent_score = score
#

a = dynproglin(alphabet, scoring_matrix, string_1, string_2)
print("Score:   ", a[0])
print("Indices: ", a[1], a[2])
score = check_score(alphabet + '_', scoring_matrix, string_1, string_2, a[1], a[2])
print('CHECKING SCORE: {} \n'.format(score))
if score != recent_score:
    print(string_1 + ' and ' + string_2 + ' do not have matching alignments...')
