


# def smith_waterman(unique_letters: str, scoring_matrix: list, seq1: str, seq2: str):
#     """
#     Applies the smith waterman algorithm to produce an alignment of two sequences
#     :param unique_letters: A string of unique letters of length p
#     :param scoring_matrix: a (p + 1) × (p + 1) scoring matrix (list of lists)
#     :param seq1: sequence 1
#     :param seq2: sequence 2
#     :return:
#     """
#     # PRINT SCORING MATRIX
#     # Scoring matrix is used to retrieve the static scores for matching certain elements
#     for line in scoring_matrix:
#         for num in line:
#             print(num)


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




def smith_waterman1(unique_letters: str, scoring_matrix: list):
    unique_letters += '-'
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
    print([' '] + [a for a in unique_letters])
    for index, line in enumerate(scoring_matrix):
        line = [str(x) for x in line]
        char = unique_letters[index]
        line = [char] + line
        print(line)

smith_waterman1('abc', [[1,-1,-2,-1],[-1,2,-4,-1],[-2,-4,3,-2],[-1,-1,-2,0]])
letter_s = 'c'; letter_t = '-'
print('Score of matching ' + letter_s + ' with ' + letter_t + ': ', index_scoring(letter_s, letter_t))
