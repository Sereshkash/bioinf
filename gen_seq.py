import seq_split
import random   
import subprocess

def insert_del(repeat_seq, number_mut, alphabet = 'GTAC'):
    insert_del_choice = [random.randint(0, 1) for i in range(number_mut)]
    for i in range(number_mut):
        if insert_del_choice[i] == 0:
            index_del = random.randint(0, len(repeat_seq) - 1)
            repeat_seq = repeat_seq[:index_del] + repeat_seq[index_del + 1:] 
        elif insert_del_choice[i] == 1:
            index_ins = random.randint(0, len(repeat_seq) - 1)
            repeat_seq = repeat_seq[:index_ins] + random.choice(alphabet) + repeat_seq[index_ins:]
    return repeat_seq

def change(repeat_seq, number_mut, alphabet = 'GTAC'):
    for i in range(number_mut):
        index_change = random.randint(0, len(repeat_seq) - 1)
        repeat_seq = repeat_seq[:index_change] + random.choice(alphabet) + repeat_seq[index_change + 1:]
    return repeat_seq

def mutation(repeat_seq, number_repeat, x = 1,  alphabet = 'GTAC'):
    set_repeat = [0] * number_repeat
    number_insert_del = len(repeat_seq) // 100
    number_change = int(len(repeat_seq) * x / 2)
    for i in range(number_repeat):
        repeat_i = insert_del(repeat_seq, number_insert_del, alphabet)
        repeat_i = change(repeat_i, number_change, alphabet)
        set_repeat[i] = repeat_i
    return set_repeat



def generation_seq(len_s, number_repeat, len_repeat, x = 0, alphabet = 'GTAC'):
    repeat_seq = ''.join(random.choice(alphabet) for _ in range(len_repeat))
    set_mutation_repeat = mutation(repeat_seq, number_repeat, x, alphabet)

    beginning_number_repeat = seq_split.number_split(len_s, number_repeat, len_repeat, len_repeat // 100)


    f1 = open('num_repeat', 'w')
    s = ' '.join(map(str, beginning_number_repeat))
    f1.write(s)
    f1.close()


    result_seq = ''.join(random.choice(alphabet) for _ in range(beginning_number_repeat[0])) + set_mutation_repeat[0]
   
    for i in range(1, number_repeat):
        rand_gen = ''.join(random.choice(alphabet) for _ in range(beginning_number_repeat[i] - beginning_number_repeat[i - 1] - len(set_mutation_repeat[i - 1])))
        result_seq +=  rand_gen + set_mutation_repeat[i]
    result_seq +=  ''.join(random.choice(alphabet) for _ in range(len_s - len(set_mutation_repeat[-1]) - beginning_number_repeat[-1]))
    return result_seq


#print(len(generation_seq(10**7, 10**3, 400, x = 1, alphabet = 'GTAC')))