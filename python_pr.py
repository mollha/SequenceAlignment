def find_missing(l1: list, l2: list) -> int:
    new_set = set(l1).symmetric_difference(l2)
    assert(type(l1) is list)
    assert (type(l2) is list)
    return new_set.pop()

def find_missing_new(l1, l2):


a = find_missing({1, 2}, [4, 9, 12, 6])
print(a)
