import random

def f(n):
    print(n)

def number_split(len_s, N, len_subsequence):
    beginnings = []
    for i in range(N):
        begin_i = random.uniform(0, len_s - len_subsequence - 1)
        j = 0
        while j < i:
            if abs(beginnings[j] - begin_i) < len_subsequence:
                j = 0
                begin_i = random.uniform(0, len_s - len_subsequence - 1)
            else:
                j += 1
        beginnings[i] = begin_i
    return beginnings

def string_split(S, N, len_subsequence):
    beginnings = number_split(len(S), N, len_subsequence)
    set_of_beginnings = [S[i:i + len_subsequence] for i in number_split]
    return set_of_beginnings

#print(number_split(100, 5, 10))
print(f(5))
