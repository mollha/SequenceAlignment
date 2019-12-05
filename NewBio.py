from random import choice
import time

# ------------------------ HELPER FUNCTIONS --------------------------


# TODO remove this function once checks are complete
def check_score(alpha: str, scoring: list, seq_s: str, seq_t: str, alignment_s: list, alignment_t: list):
    score_val = 0
    for i in range(alignment_s[0], alignment_s[-1]):
        if i not in alignment_s:
            score_val += get_score(alpha, scoring, seq_s[i], '_')

    for i in range(alignment_t[0], alignment_t[-1]):
        if i not in alignment_t:
            score_val += get_score(alpha, scoring, '_', seq_t[i])

    while alignment_s and alignment_t:
        score_val += get_score(alpha, scoring, seq_s[alignment_s[0]], seq_t[alignment_t[0]])
        alignment_s = alignment_s[1:]
        alignment_t = alignment_t[1:]
    return score_val


def get_score(alpha: str, scoring: list, char_s: str, char_t: str) -> int:
    alpha += '_'
    return scoring[alpha.index(char_s)][alpha.index(char_t)]


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


def last_row(alpha: str, scoring: list, seq_s: str, seq_t: str) -> list:
    row1, row2 = [], []
    for i in range(0, len(seq_s) + 1):
        if row2:
            row1 = row2
            row2 = []

        for j in range(0, len(seq_t) + 1):
            if not i and not j:
                row1.append(0)
            elif not i:
                row1.append(row1[j - 1] + get_score(alpha, scoring, '_', seq_t[j - 1]))
            elif not j:
                row2.append(row1[j] + get_score(alpha, scoring, seq_s[i - 1], '_'))
            else:
                diag = row1[j - 1] + get_score(alpha, scoring, seq_s[i - 1], seq_t[j - 1])
                up = row1[j] + get_score(alpha, scoring, seq_s[i - 1], '_')
                left = row2[j - 1] + get_score(alpha, scoring, '_', seq_t[j - 1])
                val = max(diag, up, left)
                row2.append(val)
    return row2


def needleman_wunsch(alpha: str, scoring: list, seq_s: str, seq_t: str) -> tuple:
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
                val = values[i][j - 1] + get_score(alpha, scoring, '_', seq_t[j - 1])
                path_val = 'L'
            elif not j:
                val = values[i - 1][j] + get_score(alpha, scoring, seq_s[i - 1], '_')
                path_val = 'U'
            else:
                diag = values[i - 1][j - 1] + get_score(alpha, scoring, seq_s[i - 1], seq_t[j - 1])
                up = values[i - 1][j] + get_score(alpha, scoring, seq_s[i - 1], '_')
                left = values[i][j - 1] + get_score(alpha, scoring, '_', seq_t[j - 1])
                val = max(diag, up, left)
                path_val = ['D', 'L', 'U'][[diag, left, up].index(val)]
            values[i].append(val)
            paths[i].append(path_val)
    return backtrack(paths, (len(seq_s), len(seq_t)))


def hirschberg(alpha: str, scoring: list, seq_s: str, seq_t: str, offset: tuple, left=True) -> tuple:
    if len(seq_s) == 1 or len(seq_t) == 1:
        alignment_s, alignment_t = needleman_wunsch(alpha, scoring, seq_s, seq_t)
        alignment_s, alignment_t = [s + offset[0] for s in alignment_s], [t + offset[1] for t in alignment_t]
        return alignment_s, alignment_t
    else:

        last_row_1 = last_row(alpha, scoring, seq_s[:len(seq_s) // 2], seq_t)
        last_row_2 = last_row(alpha, scoring, seq_s[len(seq_s) // 2:][::-1], seq_t[::-1])
        last_row_2.reverse()

        zipped = [x + y for x, y in zip(last_row_1, last_row_2)]
        max_val = max(zipped)
        max_positions = [i for i, j in enumerate(zipped) if j == max_val]

        if left:
            index = min(max_positions)
        else:
            index = max(max_positions)

        alignment_s_half1, alignment_t_half1 = hirschberg(alpha, scoring, seq_s[:len(seq_s) // 2],
                                                          seq_t[:index], offset, left=True)
        alignment_s_half2, alignment_t_half2 = hirschberg(alpha, scoring, seq_s[len(seq_s) // 2:],
                                                          seq_t[index:],
                                                          (offset[0] + (len(seq_s) // 2), offset[1] + index), left=False)
    return alignment_s_half1 + alignment_s_half2, alignment_t_half1 + alignment_t_half2



# ---------------------------------------------------------------------


def dynprog(alphabet: str, scoring_matrix: list, seq_s: str, seq_t: str):
    high_score = -float('inf')
    max_indices = None

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

    alignment_s, alignment_t = hirschberg(alphabet, scoring_matrix, seq_s, seq_t, (s_add, t_add))
    return [high_score, alignment_s, alignment_t]





# ---------------------------------------------------------------------------------------------------------------------

# TESTS

string_1 ="BDABCBBDCAABDDCDABACDADDCBCABA"
string_2 ="CCAABBBDABACADBCCADCADAAACDCCA"



scoring_matrix = [[ 1,-5,-5,-5,-1],[-5, 1,-5,-5,-1],[-5,-5, 5,-5,-4],[-5,-5,-5, 6,-4],[-1,-1,-4,-4,-9]]
alphabet = "ABCD"

x = 500
string_1 = "".join(choice(list(alphabet)) for i in range(x))
string_2 = "".join(choice(list(alphabet)) for i in range(x))

stamp1 = time.time()
a = dynprog (alphabet, scoring_matrix, string_1, string_2)
print("Score:   ", a[0])
print("Indices: ", a[1],a[2])
score = check_score(alphabet + '_', scoring_matrix, string_1, string_2, a[1],a[2])
print('CHECKING SCORE: {} \n'.format(score))
recent_score = score
stamp2 = time.time()
print('Seconds for 1st', stamp2-stamp1)

a = dynproglin (alphabet, scoring_matrix, string_1, string_2)
print("Score:   ", a[0])
print("Indices: ", a[1],a[2])
score = check_score(alphabet + '_', scoring_matrix, string_1, string_2, a[1],a[2])
print('CHECKING SCORE: {} \n'.format(score))
if score != recent_score:
    print(string_1 + ' and ' + string_2 + ' do not have matching alignments...')
stamp3 = time.time()
print('Seconds for 1st', stamp3-stamp2)

