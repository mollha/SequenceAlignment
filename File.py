def smith_waterman(unique_letters: str, scoring_matrix: list, seq_s: str, seq_t: str):
    """
    Applies the smith waterman algorithm to produce an alignment of two sequences
    :param unique_letters: A string of unique letters of length p
    :param scoring_matrix: a (p + 1) × (p + 1) scoring matrix (list of lists)
    :param seq1: sequence 1
    :param seq2: sequence 2
    :return:
    """
    indel = len(unique_letters)
    unique_letters += '-'
    print('Unique Letters', unique_letters)
    print([x for x in unique_letters])
    for index, row in enumerate(scoring_matrix):
        print(row)
    print('Sequence S:', seq_s)
    print('Sequence T:', seq_t)

    values = []
    for i in range(0, len(seq_s)+1):
        values.append([])
        for j in range(0, len(seq_t)+1):
            index_s = unique_letters.index(seq_s[i - 1])
            index_t = unique_letters.index(seq_t[j - 1])

            if not i and not j:
                val = 0
                print('corner')
            elif not i:
                val = max(0, values[i][j-1] + scoring_matrix[indel][index_t])
                if val == values[i][j-1] + scoring_matrix[indel][index_t]:
                    print('left')
                else:
                    print('restart')
            elif not j:
                val = max(0, values[i-1][j] + scoring_matrix[index_s][indel])
                if val == values[i-1][j] + scoring_matrix[index_s][indel]:
                    print('up')
                else:
                    print('restart')
            else:
                val = max(0,
                          values[i-1][j-1] + scoring_matrix[index_s][index_t],
                          values[i-1][j] + scoring_matrix[index_s][indel],
                          values[i][j-1] + scoring_matrix[indel][index_t])
                if val == values[i-1][j-1] + scoring_matrix[index_s][index_t]:
                    print('diag')
                elif val == values[i-1][j] + scoring_matrix[index_s][indel]:
                    print('up')
                elif val == values[i][j-1] + scoring_matrix[indel][index_t]:
                    print('left')
                elif val == 0:
                    print('restart')

            values[i].append(val)
            for row in values:
                print(row)
            print('-----------------')




smith_waterman('ABC', [[1,-1,-2,-1],[-1,2,-4,-1],[-2,-4,3,-2],[-1,-1,-2,0]], 'ABCACA', 'BAACB')

# smith_waterman('AGC', [[1,-1,-2,-1],[-1,2,-4,-1],[-2,-4,3,-2],[-1,-1,-2,0]], 'AAAC', 'AGC')
# smith_waterman('ACT', [[1,-1,-1,-2],[-1,1,-1,-2],[-1,-1,1,-2],[-2,-2,-2,0]], 'TAATA', 'TACTAA')

# smith_waterman('CTGA', [[1,-1,-1,-1,-5], [-1,1,-1,-1,-5], [-1,-1,1,-1,-5],
#                         [-1,-1,-1,1,-5], [-5,-5,-5,-5,-5]], 'CATTCAC', 'CTCGCAGC')

