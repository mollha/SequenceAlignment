def smith_waterman(unique_letters: str, scoring_matrix: list, seq_s: str, seq_t: str):
    """
    Applies the smith waterman algorithm to produce an alignment of two sequences
    :param unique_letters: A string of unique letters of length p
    :param scoring_matrix: a (p + 1) × (p + 1) scoring matrix (list of lists)
    :param seq1: sequence 1
    :param seq2: sequence 2
    :return:
    """
    # PRINT SCORING MATRIX
    # Scoring matrix is used to retrieve the static scores for matching certain elements
    my_row = [' '] + [a for a in unique_letters] + ['-']
    print(my_row)
    for index, line in enumerate(scoring_matrix):
        line = [str(x) for x in line]
        char = unique_letters[index]
        line = [char] + line
        print(line)

    unique_dict = {char:index for (index, char) in enumerate(unique_letters)}
    unique_dict['-'] = len(unique_letters)
    print(unique_dict)

    # ---- FILL MATRIX
    m = len(seq_s)
    n = len(seq_t)
    # values_matrix = []
    # values_matrix.append([] * n)
    # for i in range(len(seq_s) - 1):
    #     values_matrix.append([])

    values_matrix = []
    for i in range(m + 1):
        values_matrix.append([])
        for j in range(n + 1):
            print('i = ' + str(i) + ' and j = ' + str(j))
            if not i and not j:
                values_matrix[i].append(0)

            index_s = unique_dict[seq_s[i]]
            index_t = unique_dict[seq_t[j]]
            print('Score of matching ' + seq_s[i] + ' with ' + seq_t[j] + ': ', scoring_matrix[index_s][index_t])

            #
            #
            # if not j:
            #     max(0, values_matrix[][])


    values_matrix = []  # Create values matrix
    values_matrix.append([0] * len(seq_t))
    for i in range(len(seq_s) - 1):
        values_matrix.append([0])

    for row in values_matrix:
        print(row)
    unique_letters += '-'

    # PRINT SCORING MATRIX
    # Scoring matrix is used to retrieve the static scores for matching certain elements

smith_waterman('abc', [[1,-1,-2,-1],[-1,2,-4,-1],[-2,-4,3,-2],[-1,-1,-2,0]], 'ab', 'bc')


def index_scoring(letter_s, letter_t):
    scoring_matrix = [[1, -1, -2, -1], [-1, 2, -4, -1], [-2, -4, 3, -2], [-1, -1, -2, 0]]
    unique_letters = 'abc'
    index_s = None
    index_t = None
    if letter_s == '-':
        index_s = len(unique_letters)
    else:
        index_s = unique_letters.index(letter_s)
    if letter_t == '-':
        index_t = len(unique_letters)
    else:
        index_t = unique_letters.index(letter_t)
    return scoring_matrix[index_s][index_t]




# def smith_waterman1(unique_letters: str, scoring_matrix: list):
#     unique_letters += '-'
#     """
#     Applies the smith waterman algorithm to produce an alignment of two sequences
#     :param unique_letters: A string of unique letters of length p
#     :param scoring_matrix: a (p + 1) × (p + 1) scoring matrix (list of lists)
#     :param seq1: sequence 1
#     :param seq2: sequence 2
#     :return:
#     """
#
#

#
# smith_waterman1('abc', [[1,-1,-2,-1],[-1,2,-4,-1],[-2,-4,3,-2],[-1,-1,-2,0]])
# letter_s = 'c'; letter_t = '-'

