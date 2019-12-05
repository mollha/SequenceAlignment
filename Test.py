from random import choice
from main import heuralign

def get_score(alpha: str, scoring: list, char_s: str, char_t: str) -> int:
    alpha += '_'
    return scoring[alpha.index(char_s)][alpha.index(char_t)]


def check_score(alpha, scoring, seq_s, seq_t, alignment_s, alignment_t):
    total_score = 0
    for i in range(alignment_s[0], alignment_s[-1]):
        if i not in alignment_s:
            total_score += get_score(alphabet, scoring_matrix, seq_s[i], '_')

    for i in range(alignment_t[0], alignment_t[-1]):
        if i not in alignment_t:
            total_score += get_score(alphabet, scoring_matrix, '_', seq_t[i])

    while alignment_s and alignment_t:
        total_score += get_score(alphabet, scoring_matrix, seq_s[alignment_s[0]], seq_t[alignment_t[0]])
        alignment_s = alignment_s[1:]
        alignment_t = alignment_t[1:]
    return score


if __name__ == "__main__":
    alphabet = "ABCD"
    scoring_matrix = [
        [1, 1, -5, -5, -1],
        [-2, 1, -5, -5, -1],
        [-5, -5, 5, -5, -4],
        [-5, -5, -5, 6, -4],
        [-1, -1, -4, -4, -9]]
    x = 1000
    sequence1 = "".join(choice(list(alphabet)) for i in range(45))
    sequence2 = "".join(choice(list(alphabet)) for i in range(x))


    print("Starting:")
    # Strip to ensure no whitespace
    print("Seq 1 - {0} ".format(sequence1))
    print("Seq 2 - {0}".format(sequence2))
    print("------------")

    score, out1_indices, out2_indices = heuralign(alphabet, scoring_matrix, sequence1, sequence2)

    # Output - print results
    print("------------")
    print("Score: {0}".format(score))
    print("Indices: {0} {1}".format(out1_indices, out2_indices))
    score = check_score(alphabet + '_', scoring_matrix, sequence1, sequence2, out1_indices, out2_indices)
    print('CHECKING SCORE: {} \n'.format(score))
