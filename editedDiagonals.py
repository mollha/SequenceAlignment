from random import choice
import helper_functions

"""
1) Basic dynamic programming that runs in quadratic time and space [up to 50 marks].
2) Dynamic programming that runs in linear space [up to 65 marks].
3)A Heuristic procedure that runs in sub-quadratic time (similar to FASTA and BLAST) [up to 85 marks].
They should output three items: the score of the best local alignment found by the algorithm plus two lists of indices,
one for each input sequences, that realise the matches/mismatches in the alignment.
"""

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
    print(sequence_s)
    print(sequence_t)
    return score


def banded_smith_waterman(alphabet, scoring_matrix, seq_s, seq_t, diagonal: int, search_width) -> tuple:
    high_score = -float('inf')
    max_indices = None
    values, paths = [], []
    j_pairs = {}

    for i in range(max(-diagonal, 0), min(len(seq_s), len(seq_t) - diagonal)):
        values.append([])
        paths.append([])
        grid_i = i - max(-diagonal, 0)
        for j in range(i + diagonal - search_width, min(i + diagonal + search_width + 1, len(seq_t))):
            grid_j = len(values[grid_i])
            path_val = ''

            if (not grid_i and not grid_j) or j < 0:
                val = 0
                path_val = 'R'
            elif not grid_i:
                val = max(0, values[grid_i][grid_j - 1] + get_score(alphabet, scoring_matrix, '_', seq_t[j - 1]))
                if not val:
                    path_val = 'R'
                else:
                    path_val = 'L'
            elif not grid_j:
                # # up is now [i - 1][j]
                val = max(0, values[grid_i - 1][grid_j] + get_score(alphabet, scoring_matrix, seq_s[i - 1], '_'))
                if not val:
                    path_val = 'R'
                else:
                    path_val = 'U'
            else:
                # [i - 1][j] is now diagonal
                diag = values[grid_i - 1][grid_j - 1] + get_score(alphabet, scoring_matrix, seq_s[i - 1], seq_t[j - 1])
                up = values[grid_i - 1][grid_j] + get_score(alphabet, scoring_matrix, seq_s[i - 1], '_')
                left = 0 if not values[grid_i] else values[grid_i][grid_j - 1] + get_score(alphabet, scoring_matrix, '_', seq_t[j - 1])
                val = max(0, diag, up, left)

                if val == diag:
                    path_val = 'D'
                    j_pairs[i] = j
                elif val == up:
                    path_val = 'U'
                elif val == left:
                    path_val = 'L'
                elif val == 0:
                    path_val = 'R'
            values[grid_i].append(val)
            paths[grid_i].append(path_val)

            if val > high_score:
                high_score = val
                max_indices = (grid_i, grid_j)

    print('-----')
    for row in values:
        print(row)
    for path in paths:
        print(path)

    alignment_s, alignment_t = backtrack(paths, max_indices)
    return high_score, alignment_s, alignment_t


def heuralign(alphabet: str, scoring_matrix: list, seq_s: str, seq_t: str) -> tuple:
    band_width = 3

    def tune_ktup() -> int:
        min_length = min(len(seq_s), len(seq_t))
        # TODO tune ktup somehow
        return 3

    def get_seeds(ktup_val: int) -> tuple:
        while True:
            seed_dictionary = {}
            seed_list = []
            for index in range(len(seq_s) - ktup_val + 1):
                subword = seq_s[index: index + ktup_val]
                if subword in seed_dictionary:
                    seed_dictionary[subword].append(index)
                else:
                    seed_dictionary[subword] = [index]
            for index in range(len(seq_s) - ktup_val + 1):
                subword = seq_t[index: index + ktup_val]
                score_list = [get_score(alphabet, scoring_matrix, char, char) for char in subword]
                seed_score = sum(score_list)

                if subword in seed_dictionary:
                    seed_list += [(s_index, index, seed_score) for s_index in seed_dictionary[subword]]
            if not seed_list and ktup_val > 1:
                ktup_val -= 1
            else:
                break
        return ktup_val, seed_list

    def get_diagonals(seed_list: list) -> dict:
        seed_diagonals = {}
        for seed in seed_list:
            # seeds are of length ktup, and are in the form (s_index, t_index, score)
            difference = seed[1] - seed[0]    # t_index - s_index
            if difference in seed_diagonals:
                seed_diagonals[difference].append(seed)
            else:
                seed_diagonals[difference] = [seed]
        return seed_diagonals

    def extend_diagonal(diagonal_seeds: list) -> tuple:
        # extend diagonals and merge seeds if necessary
        # seeds start at length ktup, and begin in the form ((s_start, t_start), score)
        # they end as length >= ktup, and are translated to the form ((s_start, t_start), score, length)
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
                    # Overlap exists
                    if existing_seed[2] < seed_score:
                        total_score -= existing_seed[2]
                        clone_set.remove(existing_seed)
            extended_seeds = clone_set
            extended_seeds.add((s_start, t_start, seed_score, length))
            total_score += seed_score
        return total_score, list(extended_seeds)

    #--- 1) Seed sequences - find word sequences the sequences have in common ---
    ktup, seeds = get_seeds(tune_ktup())

    # if there are NO seeds, then just choose the middle diagonal
    if not seeds:
        return banded_smith_waterman(alphabet, scoring_matrix, seq_s, seq_t, 0, band_width)

    diagonals = get_diagonals(seeds)
    diagonal_scores = []
    for diagonal_key in diagonals:
        diagonal_scores.append((extend_diagonal(diagonals[diagonal_key]), diagonal_key))
    diagonal_scores.sort(key=lambda x: x[0][0], reverse=True)
    diagonal_scores = diagonal_scores[0: min(3, len(diagonal_scores))]  # get top 3 diagonals
    diagonals = [x[1] for x in diagonal_scores]
    print(diagonals)
    highest_score = float('-inf')
    best_alignment = None
    for diff in diagonals:
        print('Diff', diff)
        diff_score, alignment_s, alignment_t = banded_smith_waterman(alphabet, scoring_matrix, seq_s, seq_t, diff, band_width)
        if diff_score > highest_score:
            best_alignment = diff_score, alignment_s, alignment_t
    return best_alignment



if __name__ == "__main__":
    # # Debug input 1
    # alphabet = "ABC"
    # scoring_matrix = [[1, -1, -2, -1], [-1, 2, -4, -1], [-2, -4, 3, -2], [-1, -1, -2, 0]]
    # sequence1 = "AABBAACA"
    # sequence2 = "CBACCCBA"
    # # Debug input 2
    # alphabet = "ABCD"
    # scoring_matrix = [
    #         [ 1,-5,-5,-5,-1],
    #         [-5, 1,-5,-5,-1],
    #         [-5,-5, 5,-5,-4],
    #         [-5,-5,-5, 6,-4],
    #         [-1,-1,-4,-4,-9]]
    # sequence1 = "AAAAACCDDCCDDAAAAACC"
    # sequence2 = "CCAAADDAAAACCAAADDCCAAAA"
    # # Debug input 3
    # alphabet = "ABCD"
    # scoring_matrix = [
    #         [ 1,-5,-5,-5,-1],
    #         [-5, 1,-5,-5,-1],
    #         [-5,-5, 5,-5,-4],
    #         [-5,-5,-5, 6,-4],
    #         [-1,-1,-4,-4,-9]]
    # sequence1 = "AACAAADAAAACAADAADAAA"
    # sequence2 = "CDCDDD"
    # Debug input 4
    alphabet = "ABCD"
    scoring_matrix = [
        [1, 1, -5, -5, -1],
        [1, 1, -5, -5, -1],
        [-5, -5, 5, -5, -4],
        [-5, -5, -5, 6, -4],
        [-1, -1, -4, -4, -9]]
    x = 5
    sequence1 = "".join(choice(list(alphabet)) for i in range(x))
    sequence2 = "".join(choice(list(alphabet)) for i in range(x))
    sequence1 = "ACABA"
    sequence2 = "CDABD"

    print("Starting:")
    # Strip to ensure no whitespace
    sequence1, sequence2 = sequence1.strip(), sequence2.strip()
    print("Seq 1 - {0} ".format(sequence1))
    print("Seq 2 - {0}".format(sequence2))
    print("------------")

    # TODO: improve runtime of scoring by using dicts!

    #  Part 3 - < O(n^2) heuristic procedure, similar to FASTA and BLAST (time)
    score, out1_indices, out2_indices = heuralign(alphabet, scoring_matrix, sequence1, sequence2)

    # Output - print results
    print("------------")
    print("Score: {0}".format(score))
    print("Indices: {0} {1}".format(out1_indices, out2_indices))
    score = check_score(alphabet + '_', scoring_matrix, sequence1, sequence2, out1_indices, out2_indices)
    print('CHECKING SCORE: {} \n'.format(score))