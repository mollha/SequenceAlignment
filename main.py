# --------------------------------------------- ADDITIONAL FUNCTIONS --------------------------------------------------
def get_score(alpha: str, scoring: list, char_s: str, char_t: str) -> int:
    alpha += '_'
    return scoring[alpha.index(char_s)][alpha.index(char_t)]


def check_score(alpha: str, scoring: list, seq_s: str, seq_t: str, alignment_s: list, alignment_t: list) -> int:
    total_score = 0
    for i in range(alignment_s[0], alignment_s[-1]):
        if i not in alignment_s:
            total_score += get_score(alpha, scoring, seq_s[i], '_')

    for i in range(alignment_t[0], alignment_t[-1]):
        if i not in alignment_t:
            total_score += get_score(alpha, scoring, '_', seq_t[i])

    while alignment_s and alignment_t:
        total_score += get_score(alpha, scoring, seq_s[alignment_s[0]], seq_t[alignment_t[0]])
        alignment_s = alignment_s[1:]
        alignment_t = alignment_t[1:]
    return total_score


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


def banded_sw(alpha, scoring, seq_s, seq_t, st_pair):
    band_width = 30
    best_index, max_score = (0, 0), float('-inf')
    u, v, x, y_corr = st_pair[0], st_pair[1], st_pair[0] + 1, st_pair[1] + 1
    shift = max(0, st_pair[0] - min(st_pair[0], st_pair[1]) - band_width), \
            max(0, st_pair[1] - min(st_pair[0], st_pair[1]) - band_width)
    cell_set = set()

    # ---------------------------- INITIALIZE CELL SET ----------------------------
    def get_cells(element1, element2):
        update_cells = set()
        for i in range(-band_width, band_width + 1):
            for j in range(-band_width, band_width + 1):
                if 0 <= element1 + i < len(seq_s) and 0 <= element2 + j < len(seq_t):
                    update_cells.add((element1 + i, element2 + j))
        return update_cells

    for s_index in range(min(st_pair)):
        u, v = u - s_index, v - s_index
        cell_set.update(get_cells(u, v))

    for t_index in range(min(len(seq_s) - st_pair[0] + 1, len(seq_t) - st_pair[1] + 1)):
        x, y_corr = x + t_index, y_corr + t_index
        cell_set.update(get_cells(u, v))
    # -----------------------------------------------------------------------------

    seq_s, seq_t = seq_s[shift[0]:len(seq_s) - shift[1]], seq_t[shift[1]:len(seq_t) - shift[0]]

    values = [[0 for _ in range(len(seq_s) + 1)] for _ in range(len(seq_t) + 1)]
    len_val, len_val_zero = len(values), len(values[0])
    paths = [['R' for _ in range(len_val_zero)] for _ in range(len_val)]

    for i in range(len(seq_s)):
        values[0][i + 1] = max(0, values[0][i] + get_score(alpha, scoring, seq_s[i], '_'))

    for i in range(len(seq_t)):
        values[i + 1][0] = max(0, values[i][0] + get_score(alpha, scoring, seq_t[i], '_'))

    for y in range(len(seq_t)):
        y_corr, processed = y + 1, 0
        for x in range(1, len(seq_s) + 1):
            if (x + shift[0], y_corr + shift[1]) in cell_set:
                processed += 1
                # ----------------------------- CALCULATE VALUES -------------------------------------
                diag = values[y_corr - 1][x - 1] + get_score(alpha, scoring, seq_t[y_corr - 1], seq_s[x - 1])
                up = values[y_corr - 1][x] + get_score(alpha, scoring, seq_t[y_corr - 1], '_')
                left = values[y_corr][x - 1] + get_score(alpha, scoring, '_', seq_s[x - 1])
                max_val = max([diag, up, left, 0])
                values[y_corr][x] = max_val
                path_val = ['D', 'L', 'U', 'R'][[diag, left, up, 0].index(max_val)]
                # -----------------------------------------------------------------------------------
                tuples = [(x - 1 + shift[0], y_corr - 1 + shift[1]),
                          (x - 1 + shift[0], y_corr + shift[1]),
                          (x + shift[0], y_corr - 1 + shift[1])]
                if path_val != 'R' and tuples[['D', 'L', 'U', 'R'].index(path_val)] in cell_set:
                    paths[y_corr][x] = path_val
                if max_score < max_val:
                    max_score = max_val
                    best_index = (y_corr, x)
            else:
                values[y_corr][x] = 0
                if processed:
                    break
    alignment_t, alignment_s = backtrack(paths, best_index)
    return max_score, alignment_s, alignment_t


# --------------------------------------------------------------------------------------------------------------------


def dynprog(alphabet: str, scoring_matrix: list, seq_s: str, seq_t: str) -> tuple:
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
    return high_score, alignment_s, alignment_t


def dynproglin(alphabet: str, scoring_matrix: list, seq_s: str, seq_t: str) -> tuple:

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
    return (high_score, alignment_s, alignment_t)


def heuralign(alphabet: str, scoring_matrix: list, seq_s: str, seq_t: str) -> tuple:
    ktup = max(3, int(-(-min(len(seq_s), len(seq_t)) // 50)))

    while True:
        seed_dictionary = {}
        seeds = []
        for index in range(len(seq_s) - ktup + 1):
            subword = seq_s[index: index + ktup]
            if subword in seed_dictionary:
                seed_dictionary[subword].append(index)
            else:
                seed_dictionary[subword] = [index]
        for index in range(len(seq_s) - ktup + 1):
            subword = seq_t[index: index + ktup]
            score_list = [get_score(alphabet, scoring_matrix, char, char) for char in subword]
            seed_score = sum(score_list)

            if subword in seed_dictionary:
                seeds += [(s_index, index, seed_score) for s_index in seed_dictionary[subword]]
        if not seeds and ktup > 1:
            ktup -= 1
        else:
            break

    # if there are NO seeds, then just choose the middle diagonal
    if not seeds:
        return banded_sw(alphabet, scoring_matrix, seq_s, seq_t, (0, 0))

    diagonals = {}
    for seed in seeds:
        difference = seed[1] - seed[0]
        if difference in diagonals:
            diagonals[difference].append(seed)
        else:
            diagonals[difference] = [seed]

    def extend_diagonal(diagonal_seeds: list) -> tuple:
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
                    # overlap exists
                    if existing_seed[2] < seed_score:
                        total_score -= existing_seed[2]
                        clone_set.remove(existing_seed)
            extended_seeds = clone_set
            extended_seeds.add((s_start, t_start, seed_score, length))
            total_score += seed_score
        return total_score, list(extended_seeds)

    diagonal_scores = []
    for diagonal_key in diagonals:
        diagonal_scores.append((extend_diagonal(diagonals[diagonal_key]), diagonal_key))
    diagonal_scores.sort(key=lambda x: x[0][0], reverse=True)
    top_3 = diagonal_scores[0: min(3, len(diagonal_scores))]  # get top 3 diagonals
    tuples = [triple[0][1][0] for triple in top_3]
    best_seeds = [(seed_tuple[0], seed_tuple[1]) for seed_tuple in tuples]

    max_score = -float('inf')
    response = None
    for seed in best_seeds:
        score, alignment_s, alignment_t = banded_sw(alphabet, scoring_matrix, seq_s, seq_t, seed)
        if max_score < score:
            max_score = score
            response = score, alignment_s, alignment_t
    return check_score(alphabet, scoring_matrix, seq_s, seq_t, response[1], response[2]), response[1], response[2]
