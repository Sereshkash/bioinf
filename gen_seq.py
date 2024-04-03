import seq_split
import random   
import os

def generation_seq(len_s, number_repeat, len_repeat, alphabet):
    repeat_seq = ''.join(random.choice(alphabet) for _ in range(len_repeat))
    beginning_number_repeat = seq_split.number_split(len_s, number_repeat, len_repeat)
    result_seq = ''.join(random.choice(alphabet) for _ in range(beginning_number_repeat[0])) + repeat_seq
    for i in range(1, number_repeat):
        rand_gen = ''.join(random.choice(alphabet) for _ in range(beginning_number_repeat[i] - beginning_number_repeat[i - 1] - len_repeat))
        result_seq +=   rand_gen + repeat_seq
    result_seq +=  ''.join(random.choice(alphabet) for _ in range(len_s - len_repeat - beginning_number_repeat[-1]))
    return result_seq

s = generation_seq(10 ** 7, 10 ** 3, 400, 'ATGC')
file = open('example.txt','w')
file.write(s)
file.close()