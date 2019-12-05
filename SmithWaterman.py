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
    return score


def banded_sw(alpha, scoring, seq_s, seq_t, st_pair):
    band_width = 30
    u, v, x, y = st_pair[0], st_pair[1], st_pair[0] + 1, st_pair[1] + 1
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
        x, y = x + t_index, y + t_index
        cell_set.update(get_cells(u, v))
    # -----------------------------------------------------------------------------

    seq_s, seq_t = seq_s[shift[0]:len(seq_s) - shift[1]], seq_t[shift[1]:len(seq_t) - shift[0]]

    values, paths = [], []
    for i in range(len(seq_s) + 1):
        values.append([])
        paths.append([])
        for j in range(len(seq_t) + 1):
            values[i].append(0)
            if not i and not j:
                paths[i].append('R')
            else:
                paths[i].append('U' if not j else 'L' if not i else 'R')

    values[0] = [0 for _ in range(len(values[0]))]

    for i in range(len(values)):
        values[i][0] = 0

    max_index, max_score = (0, 0), float('-inf')

    # Iterate over scoring matrix and generate scoring (start at 1,1 and work from there)
    for y in range(1, len(seq_t) + 1):  # y -> seq_t
        flag = False
        for x in range(1, len(seq_s) + 1):  # x -> seq_s
            # Check if current scoring cell lies in diagonal
            if (x + shift[0], y + shift[1]) in cell_set:
                flag = True
                diag = values[y - 1][x - 1] + get_score(alpha, scoring, seq_t[y - 1], seq_s[x - 1])
                up = values[y - 1][x] + get_score(alpha, scoring, seq_t[y - 1], '_')
                left = values[y][x - 1] + get_score(alpha, scoring, '_', seq_s[x - 1])
                max_val = max([diag, up, left, 0])
                # Update scoring matrix
                values[y][x] = max_val
                # Get index of max
                path_val = ['D', 'L', 'U', 'R'][[diag, left, up, 0].index(max_val)]

                # Update backtrack matrix if score it come from is a valid cell
                if path_val == 'D' and (x - 1 + shift[0], y - 1 + shift[1]) in cell_set:
                    paths[y][x] = 'D'
                elif path_val == 'U' and (x + shift[0], y - 1 + shift[1]) in cell_set:
                    paths[y][x] = 'U'
                elif path_val == 'L' and (x - 1 + shift[0], y + shift[1]) in cell_set:
                    paths[y][x] = 'L'
                # Check if new greatest score seen (score for vals outside diagonals score = 0)
                if max_val > max_score:
                    max_score = max_val
                    max_index = (y, x)
            else:
                # If cell doesn't lie in diagonal -> score = 0 (local alignment still)
                values[y][x] = 0
                # If flag = True, have passed over banded region and back into region outside of band -> can break
                if flag:
                    break
    alignment_s, alignment_t = backtrack(paths, max_index)
    return max_score, alignment_s, alignment_t


def heuralign(alphabet: str, scoring_matrix: list, seq_s: str, seq_t: str):
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
    return response


if __name__ == "__main__":
    alphabet = "ABCD"
    scoring_matrix = [
        [1, -5, -5, -5, -1],
        [-5, 1, -5, -5, -1],
        [-5, -5, 5, -5, -4],
        [-5, -5, -5, 6, -4],
        [-1, -1, -4, -4, -9]]
    sequence1 = "DDCDDCCCDCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCDDDDDDAAAACADDCDADCDCDCDCD"
    sequence2 = "DDCDDCCCDCBCCCCDDDCDBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBDCDCDCDCD"

    # TODO check what happens when no matches

    # Strip to ensure no whitespace
    sequence1, sequence2 = sequence1.strip(), sequence2.strip()
    print("Seq 1 - {0} ".format(sequence1))
    print("Seq 2 - {0}".format(sequence2))
    print("------------")

    # Part 2 - O(n) dynamic prog. (space)
    score, out2_indices, out1_indices = heuralign(alphabet, scoring_matrix, sequence1, sequence2)

    # Output - print results
    print("Score: {0}".format(score))
    print("Indices: {0} | {1}".format(out1_indices, out2_indices))
    score = check_score(alphabet + '_', scoring_matrix, sequence1, sequence2, out1_indices, out2_indices)
    print('CHECKING SCORE: {} \n'.format(score))
