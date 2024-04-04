import random
import subprocess

def number_split(len_s, N, len_subsequence):
    beginnings = [0] * N
    for i in range(N):
        begin_i = random.randint(0, len_s - len_subsequence - 1)
        j = 0
        while j < i:
            if abs(beginnings[j] - begin_i) < len_subsequence:
                j = 0
                begin_i = random.randint(0, len_s - len_subsequence - 1)
            else:
                j += 1
        beginnings[i] = begin_i
    return sorted(beginnings)

def string_split(S, N, len_subsequence):
    beginnings = number_split(len(S), N, len_subsequence)
    set_of_subseq = [S[i:i + len_subsequence] for i in beginnings]
    return set_of_subseq

def creating_fasta_set(set_of_subseq, name = 'example.fa'):
    string_set_of_subset = ''
    for i in range(len(set_of_subseq)):
        string_set_of_subset += '>i ' + str(i + 1) + ';\n' +str(set_of_subseq[i]) + '\n'
    file = open(name,'w')
    file.write(string_set_of_subset)
    file.close()

def muscle_msa(name_in, name_out = 'seqs1.afa'):
    process = subprocess.run(['muscle', '-in', name_in, '-out', name_out], capture_output=True, text = True)
    return process



# set = string_split('CACAGGACTGTCGAAAGCAG', 3, 4)
# name = "example.fa"
# creating_fasta_set(set, name)
# muscle_msa(name, "seqs2.afa")



#creating_fasta_set(['ATTGC', 'AATTGC', 'ATGC'])

