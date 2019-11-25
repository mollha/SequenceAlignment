# ------------------------ HELPER FUNCTIONS --------------------------
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


def hirschberg(alphabet, scoring_matrix, seq_s, seq_t):

    def last_row(seq_1, seq_2):
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

    def align(seq1, seq2, seq1_offset, seq2_offset):
        out1_indices, out2_indices = [], []
        max_score = 0


        # Apply SW for optimal local alignment
        if len(seq1) == 1 or len(seq2) == 1:
            needleman_output = needleman_wunsch(alphabet, scoring_matrix, seq1, seq2)
            out1_indices, out2_indices = needleman_output[0], needleman_output[1]
            out1_indices = [x + seq1_offset for x in out1_indices]
            out2_indices = [x + seq2_offset for x in out2_indices]

        else:
            # Get midpoint of Seq2
            seq2_mid = len(seq2) // 2

            # Get scoring of lhs (in linear space)
            r_left, _, _ = last_row(seq1, seq2[:seq2_mid])

            # Get scoring of rhs (in linear space) [reversed]
            r_right, _, _ = last_row(seq1[::-1], seq2[seq2_mid:][::-1])
            r_right.reverse()  # flip back again for calculating total score

            # Sum values and find argmax
            row = [l + r for l, r in zip(r_left, r_right)]
            maxidx, maxval = max(enumerate(row), key=lambda a: a[1])

            # Partition seq1 at argmax
            seq1_mid = maxidx

            # Recursively call align on each half
            max_score_left, aligned_1_left_indices, aligned_2_left_indices = align(
                seq1[:seq1_mid], seq2[:seq2_mid], seq1_offset, seq2_offset)
            max_score_right, aligned_1_right_indices, aligned_2_right_indices = align(
                seq1[seq1_mid:], seq2[seq2_mid:], seq1_offset + seq1_mid, seq2_offset + seq2_mid)

            # Add results of recursive calls to out vars
            out1_indices = aligned_1_left_indices + aligned_1_right_indices
            out2_indices = aligned_2_left_indices + aligned_2_right_indices
            max_score = max_score_left + max_score_right

        return max_score, out1_indices, out2_indices

    def run():
        # Get max index from forward pass
        _, max_val, max_index = last_row(seq_s, seq_t)

        # Get min index from backward pass
        _, _, min_index = last_row(seq_s[::-1], seq_t[::-1])

        # Subtract lengths from min index (s.t. actual start position)
        min_index[1] = len(seq_s) - min_index[1]
        min_index[0] = len(seq_t) - min_index[0]

        # Get local alignment
        return align(seq_s[min_index[1]:max_index[1]], seq_t[min_index[0]:max_index[0]], min_index[1],
                          min_index[0])
    return run()


# ---------------------------------------------------------------------


def dynprog(alphabet: str, scoring_matrix: list, seq_s: str, seq_t: str):
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

    alignment_s, alignment_t = backtrack(paths, max_indices)
    return [high_score, alignment_s, alignment_t]

def dynproglin (alphabet, scoring_matrix, string_1, string_2):
    return hirschberg(alphabet, scoring_matrix, string_1, string_2)

string_1, string_2 = "AACCDDAACC", "CADDACDDAA"
scoring_matrix = [[ 1,-5,-5,-5,-1],[-5, 1,-5,-5,-1],[-5,-5, 5,-5,-4],[-5,-5,-5, 6,-4],[-1,-1,-4,-4,-9]]
alphabet = "ABCD"
#
a = dynprog (alphabet, scoring_matrix, string_1, string_2)
print("Score:   ", a[0])
print("Indices: ", a[1],a[2])
score = check_score(alphabet + '_', scoring_matrix, string_1, string_2, a[1],a[2])
print('CHECKING SCORE: {} \n'.format(score))
recent_score = score
#

a = dynproglin (alphabet, scoring_matrix, string_1, string_2)
print("Score:   ", a[0])
print("Indices: ", a[1],a[2])
score = check_score(alphabet + '_', scoring_matrix, string_1, string_2, a[1],a[2])
print('CHECKING SCORE: {} \n'.format(score))
if score != recent_score:
    print(string_1 + ' and ' + string_2 + ' do not have matching alignments...')

