import random
import subprocess
from Bio import SeqIO

def number_split(len_s, N, len_subsequence, dist = 0):
    beginnings = [0] * N
    for i in range(N):
        begin_i = random.randint(0, len_s - len_subsequence - 1)
        j = 0
        while j < i:
            if abs(beginnings[j] - begin_i) < (len_subsequence + dist):
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

def string_from_file(name_file):
    f1 = open(name_file, 'r')
    all_string = f1.readlines()
    if len(all_string) > 1:
        s = all_string[1]
    else:
        s = all_string[0]
    f1.close()
    return s

def creating_test_fasta_set(name_S, N, len_subsequence, name_fasta_set = 'example.fa'):
    S = string_from_file(name_S)
    
    f1 = open('num_repeat', 'r')
    s = f1.readline()
    beginnings_all = list(map(int, s.split()))
    f1.close()
    beginnings_1 = random.choices(beginnings_all, k = 1)
    beginnings_1 = beginnings_1 * 4
    #beginnings_1[0] = beginnings_1[0] - 20
    #beginnings_1[1] = beginnings_1[1] + 20
    #beginnings_1 = [beginnings_1[i] + random.randint(-20, 20) for i in range(len(beginnings_1))]
    beginnings_2 = number_split(len(S), N - 1, len_subsequence)
    beginnings_2 = beginnings_2 * 4
    beginnings = beginnings_1 + beginnings_2
    set_of_subseq = [S[i:i + len_subsequence] for i in beginnings]

    string_set_of_subset = ''
    for i in range(len(set_of_subseq)):
        string_set_of_subset +=  '>i' + str(i) + ' iteration;\n' +str(set_of_subseq[i]) + '\n'
        #string_set_of_subset += '>i' + str(i + 1) + ' Name;\n' +str(set_of_subseq[i]) + '\n'
    file = open(name_fasta_set,'w')
    file.write(string_set_of_subset)
    file.close()


def creating_fasta_set(name_S, N, len_subsequence, repeat_set = 1, name_fasta_set = 'example.fa'):
    S = string_from_file(name_S)
    set_of_subseq = string_split(S, N, len_subsequence) * repeat_set
    string_set_of_subset = ''
    for i in range(len(set_of_subseq)):
        string_set_of_subset +=  '>i' + str(i) + ' iteration;\n' +str(set_of_subseq[i]) + '\n'
        #string_set_of_subset += '>i' + str(i + 1) + ' Name;\n' +str(set_of_subseq[i]) + '\n'
    file = open(name_fasta_set,'w')
    file.write(string_set_of_subset)
    file.close()

def creating_fasta_new_nhmmer_set(name_seq, name_out_posision, N = 20, name_fasta_set = 'example.fa'):
    set_of_subseq = creating_new_set_from_out_nhmmer(name_seq, name_out_posision, N)
    string_set_of_subset = ''
    for i in range(len(set_of_subseq)):
        string_set_of_subset +=  '>i' + str(i) + ' iteration;\n' +str(set_of_subseq[i]) + '\n'
        #string_set_of_subset += '>i' + str(i + 1) + ' Name;\n' +str(set_of_subseq[i]) + '\n'
    file = open(name_fasta_set,'w')
    file.write(string_set_of_subset)
    file.close()

def creating_txt_res_nhmmer_set(name_seq, name_out_posision, min_len_subseq = 300, max_len_subseq = 500, N = 20, name_res_txt = 'example.fa'):
    set_of_subseq = creating_new_set_from_out_nhmmer(name_seq, name_out_posision, N)
    string_set_of_subset = ''
    for i in range(len(set_of_subseq)):
        if (len(set_of_subseq[i]) >= min_len_subseq) and (len(set_of_subseq[i]) <= max_len_subseq):
            string_set_of_subset +=  str(i) + ': ' + str(len(str(set_of_subseq[i]))) + '\n' + str(set_of_subseq[i]) + '\n'
        #string_set_of_subset += '>i' + str(i + 1) + ' Name;\n' +str(set_of_subseq[i]) + '\n'
    file = open(name_res_txt,'w')
    file.write(string_set_of_subset)
    file.close()

# def creating_txt_res_nhmmer_set(name_seq, name_out_posision, min_len_subseq = 300, max_len_subseq = 500, N = 20, name_res_txt = 'example.fa'):
#     set_of_subseq = creating_new_set_from_out_nhmmer(name_seq, name_out_posision, N)
#     string_set_of_subset = ''
#     for i in range(len(set_of_subseq)):
#         if (len(set_of_subseq[i]) >= min_len_subseq) and (len(set_of_subseq[i]) <= max_len_subseq):
#             string_set_of_subset +=  str(set_of_subseq[i]) + '\n'
#         #string_set_of_subset += '>i' + str(i + 1) + ' Name;\n' +str(set_of_subseq[i]) + '\n'
#     file = open(name_res_txt,'w')
#     file.write(string_set_of_subset)
#     file.close()



def txt_into_fasta(name_txt, name_fasta):
    f = open(name_txt, 'r')
    string = '>i1 seq;\n' + f.readline()
    f.close()

    f1 = open(name_fasta, 'w')
    f1.write(string)
    f1.close()

def muscle_msa(name_in, name_out = 'seqs1.afa'):
    process = subprocess.run(['muscle', '-in', name_in, '-out', name_out], capture_output=True, text = True)
    return process

def fasta_into_stockholm(name_in = 'example.fa', name_out = 'fasta_file6.sto'):
    records = SeqIO.parse(name_in, 'fasta')
    count = SeqIO.write(records, name_out, 'stockholm')

def nhmmer_on_stock_msa(name_set_msa = '1mult_sto.sto', name_seq = '1ex.fa', name_out_inf = 'result_hmm.fa', name_out_posision = 'www.fa', p_val = 0.001):
    command = 'nhmmer -o ' + name_out_inf + ' --aliscoresout ' + name_out_posision + ' --noali --notextw --singlemx --dna --incE ' + str(p_val) + ' ' + name_set_msa + ' ' + name_seq
    process = subprocess.run(command, shell = True, capture_output = True, text = True)
    return name_out_posision

def creating_new_set_from_out_nhmmer(name_seq, name_out_posision, N = 20):
    s = string_from_file(name_seq)

    f = open(name_out_posision, 'r')
    all_string = f.readlines()
    f.close()
    len_set = len(all_string)
    # if len_set > N:
    #     new_set = [0] * N
    # else:
    #     new_set = [0] * len_set
    new_set = []
    for i in range(len_set):
        string = all_string[i]
        index_sep = string.rfind(':') - 1
        prev_space = string.rfind(' ', 0, index_sep - 1)
        end = int(string[prev_space + 1:index_sep])
        prev_space2 = string.rfind(' ', 0, prev_space - 1)
        begin = int(string[prev_space2 + 1:prev_space])
        if i > (N - 1):
            break
        #if abs(begin - 1 - end) < 500 and abs(begin - 1 - end) > 300:
        new_set.append(s[begin - 1:end])
        
    # new_set = []
    # i = 0
    # if len_set == min(len_set, N):
    #     string = all_string[i]
    #     index_sep = string.rfind(':') - 1
    #     prev_space = string.rfind(' ', 0, index_sep - 1)
    #     end = int(string[prev_space + 1:index_sep])
    #     prev_space2 = string.rfind(' ', 0, prev_space - 1)
    #     begin = int(string[prev_space2 + 1:prev_space])
    #     if (begin - 1 - end) < 500 and (begin - 1 - end) > 300:
    #         new_set.append(s[begin - 1:end])
        
    # while i <= min(len_set, N):
    #     string = all_string[i]
    #     index_sep = string.rfind(':') - 1
    #     prev_space = string.rfind(' ', 0, index_sep - 1)
    #     end = int(string[prev_space + 1:index_sep])
    #     prev_space2 = string.rfind(' ', 0, prev_space - 1)
    #     begin = int(string[prev_space2 + 1:prev_space])
        
    #     new_set[i] = s[begin - 1:end]
    return new_set



#creating_new_set_from_out_nhmmer('1ex.fa','www.fa')
# set = string_split('CACAGGACTAGGATCGAAAGGCAG', 3, 4)
# # name = "example.fa"
# # creating_fasta_set(set, name)
# # muscle_msa('seqs.afa', "seqs2.fa")
# print("WWWWW")
# fasta_into_stockholm('example.fa', 'seqs.sto')

# set = ['GACCGTTATGACCCATGA', 'GACCGTTATGACCCATGA', 'GACCGTTATGACCCATGA']
# name = "1example.fa"
# creating_fasta_set(set, name)
# muscle_msa(name, "1mult.fa")
# fasta_into_stockholm("1mult.fa", "1mult_sto.sto")


# nhmmer -o result_hmm.fa --noali --notextw --singlemx --dna --incT 1 1mult_sto.sto 1ex.fa



#nhmmer_on_stock_msa()

#creating_fasta_set(['ATTGC', 'AATTGC', 'ATGC'])

#nhmmer -o result_hmm.fa --noali --notextw --singlemx --dna --incE 0.1 1mult_sto.sto 1ex.fa
