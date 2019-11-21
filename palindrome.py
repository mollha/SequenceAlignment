# Given a string, we want to determine if removing a character will making it a palindrome in O(n)

def palindrome(sequence: str):
    """
    Given a string, we want to determine if removing a character will make it a palindrome in O(n)
    Note: True is returned if the string is already a palindrome
    :param sequence: a string denoting the sequence of characters
    :return: boolean, true or false, new sequence and removed character
    """
    original_sequence = sequence

    # If we wanted to only return palindromes formed after removing a letter, if true but the life is intact,
    # then return False

    length = len(sequence)
    life = True
    for index in range(0, length // 2):
        first_char = sequence[0]
        last_char = sequence[-1]

        if first_char != last_char:
            if sequence[1] == last_char and life:
                print('removed character ' + first_char)
                life = False
                sequence = sequence[1:-1]
            elif sequence[-2] == first_char:
                print('removed character ' + last_char)
                life = False
                sequence = sequence[-2:]
            else:
                print('String: ' + sequence + ' is not palindromeable')
                return False
        sequence = sequence[1:-1]

    print(original_sequence + ' is / can be turned into a palindrome')
    return True



palindrome('ACBBBCA')
palindrome('ACCCBCA')
palindrome('ACBCACBCA')
palindrome('ACBABCA')

palindrome('ACBAABCA')

