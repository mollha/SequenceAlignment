def dynprog(unique_letters: str, scoring_matrix: list, seq_s: str, seq_t: str):
    indel = len(unique_letters)
    unique_letters += '-'
    high_score = -1000
    max_indices = None

    # Compute scores
    values, paths = [], []
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
    # print(high_score)
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
    return [high_score, alignment_s, alignment_t]

# SMITH WATERMAN TEST CASES
# score = smith_waterman('ABC', [[1,-1,-2,-1],[-1,2,-4,-1],[-2,-4,3,-2],[-1,-1,-2,0]], 'ABCACA', 'BAACB')
# score = smith_waterman('AGC', [[1,-1,-2,-1],[-1,2,-4,-1],[-2,-4,3,-2],[-1,-1,-2,0]], 'AAAC', 'AGC')
# score = smith_waterman('ACT', [[1,-1,-1,-2],[-1,1,-1,-2],[-1,-1,1,-2],[-2,-2,-2,0]], 'TAATA', 'TACTAA')
# score = smith_waterman('CTGA', [[1,-1,-1,-1,-5], [-1,1,-1,-1,-5], [-1,-1,1,-1,-5],
#                         [-1,-1,-1,1,-5], [-5,-5,-5,-5,-5]], 'CATTCAC', 'CTCGCAGC')

def dynproglin(unique_letters: str, scoring_matrix: list, seq_s: str, seq_t: str):
    indel = len(unique_letters)
    unique_letters += '-'
    high_score = -1000
    max_indices = None

    # Compute matrix values based on previous columns

    col1, col2 = [], []
    # Compute scores
    for j in range(0, len(seq_t) + 1):
        for i in range(0, len(seq_s)+1):
            index_s = unique_letters.index(seq_s[i - 1])
            index_t = unique_letters.index(seq_t[j - 1])

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
                up = col2[i-1] + scoring_matrix[index_s][indel]
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

    # ---- BACKTRACK -----
    i, j = max_indices
    print(high_score)
    print(max_indices)


a = dynprog ("ABC", [[1,-1,-2,-1],[-1,2,-4,-1],[-2,-4,3,-2],[-1,-1,-2,0]], "AABBAACA", "CBACCCBA")
print("Score:   ", a[0])
print("Indices: ", a[1],a[2])

a = dynproglin ("ABC", [[1,-1,-2,-1],[-1,2,-4,-1],[-2,-4,3,-2],[-1,-1,-2,0]], "AABBAACA", "CBACCCBA")
