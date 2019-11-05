
def smith_waterman(unique_letters: str, scoring_matrix: list, seq_s: str, seq_t: str):
    """
    Applies the smith waterman algorithm to produce an alignment of two sequences
    :param unique_letters: A string of unique letters of length p
    :param scoring_matrix: a (p + 1) × (p + 1) scoring matrix (list of lists)
    :param seq1: sequence 1
    :param seq2: sequence 2
    :return: list containing score of best local alignment, and two lists of indices describing the alignment
    """
    indel = len(unique_letters)
    unique_letters += '-'
    print('Unique Letters', unique_letters)
    print([x for x in unique_letters])
    for index, row in enumerate(scoring_matrix):
        print(row)
    print('Sequence S:', seq_s)
    print('Sequence T:', seq_t)

    high_score = 0
    max_indices = None

    values = []
    paths = []
    for i in range(0, len(seq_s)+1):
        values.append([])
        paths.append([])
        for j in range(0, len(seq_t)+1):
            index_s = unique_letters.index(seq_s[i - 1])
            index_t = unique_letters.index(seq_t[j - 1])
            path_val = ''

            if not i and not j:
                val = 0
                path_val = 'R'
            elif not i:
                val = max(0, values[i][j-1] + scoring_matrix[indel][index_t])
                if val == values[i][j-1] + scoring_matrix[indel][index_t]:
                    path_val = 'L'
                else:
                    path_val = 'R'
            elif not j:
                val = max(0, values[i-1][j] + scoring_matrix[index_s][indel])
                if val == values[i-1][j] + scoring_matrix[index_s][indel]:
                    path_val = 'U'
                else:
                    path_val = 'R'
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
    for row in values:
        print(row)
    print('-----------------')
    for row in paths:
        print(row)
    print('-----------------')

    # ---- BACKTRACK -----
    i, j = max_indices
    print(high_score)
    print(max_indices)
    align_seq_s, align_seq_t = '', ''
    alignment_s, alignment_t = [], []

    while True:       # while we haven't yet had to restart
        path = paths[i][j]
        if path == 'D':
            align_seq_s += seq_s[i - 1]
            align_seq_t += seq_t[j - 1]
            alignment_s.append(i-1)
            alignment_t.append(j-1)
            i, j, = i - 1, j - 1
        elif path == 'U':
            align_seq_s += seq_s[i - 1]
            align_seq_t += '-'
            i = i - 1
        elif path == 'L':
            align_seq_s += '-'
            align_seq_t += seq_t[j - 1]
            j = j - 1
        else:
            break
    align_seq_s = align_seq_s[::-1]
    align_seq_t = align_seq_t[::-1]
    alignment_s.reverse()
    alignment_t.reverse()
    print('SEQ S: ', align_seq_s)
    print('SEQ T: ', align_seq_t)
    return [high_score, alignment_s, alignment_t]



# score = smith_waterman('ABC', [[1,-1,-2,-1],[-1,2,-4,-1],[-2,-4,3,-2],[-1,-1,-2,0]], 'ABCACA', 'BAACB')
# score = smith_waterman('AGC', [[1,-1,-2,-1],[-1,2,-4,-1],[-2,-4,3,-2],[-1,-1,-2,0]], 'AAAC', 'AGC')
# score = smith_waterman('ACT', [[1,-1,-1,-2],[-1,1,-1,-2],[-1,-1,1,-2],[-2,-2,-2,0]], 'TAATA', 'TACTAA')

score = smith_waterman('CTGA', [[1,-1,-1,-1,-5], [-1,1,-1,-1,-5], [-1,-1,1,-1,-5],
                        [-1,-1,-1,1,-5], [-5,-5,-5,-5,-5]], 'CATTCAC', 'CTCGCAGC')

print(score)