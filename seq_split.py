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
    beginnings_1 = random.choices(beginnings_all, k = 2)
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

def creating_fasta_new_nhmmer_set(name_seq, name_out_posision, N = 20, name_fasta_set = 'example.fa', left_step = 0, right_step = 0):
    set_of_subseq = creating_new_set_from_out_nhmmer(name_seq, name_out_posision, N, left_step = left_step, right_step = right_step)
    string_set_of_subset = ''
    for i in range(len(set_of_subseq)):
        string_set_of_subset +=  '>i' + str(i) + ' iteration;\n' +str(set_of_subseq[i]) + '\n'
        #string_set_of_subset += '>i' + str(i + 1) + ' Name;\n' +str(set_of_subseq[i]) + '\n'
    file = open(name_fasta_set,'w')
    file.write(string_set_of_subset)
    file.close()

def creating_txt_res_nhmmer_set(name_seq, name_out_posision, min_len_subseq = 300, max_len_subseq = 500, N = 20, name_res_txt = 'example.fa', left_step = 0, right_step = 0):
    set_of_subseq = creating_new_set_from_out_nhmmer(name_seq, name_out_posision, N, left_step = left_step, right_step = right_step)
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

def creating_new_set_from_out_nhmmer(name_seq, name_out_posision, N = 20, left_step = 0, right_step = 0, name_position_subseq = 'position_subseq1.fa'):
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
    f = open(name_position_subseq, 'w')
    f.write('')
    f.close()
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
        if (begin - left_step) <= 0:
            begin = 0
        else:
            begin = begin - left_step
        
        if (end + right_step) >= (len(s) - 1):
            end = len(s) - 1
        else:
            end = end + right_step
        #print(i, begin, end, 'diff',end - begin)
        if end - begin > 0:
            f = open(name_position_subseq, 'a')
            f.write('>' + str(begin) + ' ' + str(end) + ':\n' + str(s[begin:end]) + '\n')
            f.close()
            new_set.append(s[begin:end])
        
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

def increase_len_subseq(s_name_txt, name_result_nhmmer_posision, name_result_nhmmer_posision_new = 'increase2/name_result_pos.fa',left_step = 0, right_step = 0, min_len_subseq = 20, max_len_subseq = 1000, N = 150):
    name_folder = 'increase/'
    name_set_subseq0 = name_folder + 'name_old_set_fasta.fa'
    name_set_subseq0_msa = name_folder + 'name_old_set_fasta_msa.fa'
    name_set_subseq0_msa_stockh = name_folder + 'name_old_set_fasta_msa.sto'
    s_name_fasta = name_folder + 's_name.fa'
    name_result_nhmmer_info = name_folder + 'name_result_nhmmer_info.fa'
    #name_result_nhmmer_posision_new = name_folder + 'name_result_pos.fa'
    txt_into_fasta(s_name_txt, s_name_fasta)
    mean_len_subseq = 0

    creating_txt_res_nhmmer_set(s_name_txt, name_result_nhmmer_posision, min_len_subseq = 0, max_len_subseq = 1000, N =  10 ** 6, name_res_txt = name_folder + 'first_subseq.txt')

    set_subseq = creating_new_set_from_out_nhmmer(s_name_txt, name_result_nhmmer_posision, N = 10 ** 6, left_step = 0, right_step = 0, name_position_subseq = name_folder + 'pos_sub2.fa')
    min_x = 1000
    max_x = 0
    for i in range(len(set_subseq)):
        mean_len_subseq += len(set_subseq[i])
        if min_x > len(set_subseq[i]):
            min_x = len(set_subseq[i])
        if max_x < len(set_subseq[i]):
            max_x = len(set_subseq[i])
    mean_len_subseq = mean_len_subseq // len(set_subseq)
    #print(mean_len_subseq, len(set_subseq), min_x, max_x)
    creating_fasta_new_nhmmer_set(s_name_txt, name_result_nhmmer_posision, N, name_set_subseq0, left_step = left_step, right_step = right_step)

    muscle_msa(name_set_subseq0, name_set_subseq0_msa)
    fasta_into_stockholm(name_set_subseq0_msa, name_set_subseq0_msa_stockh)
    nhmmer_on_stock_msa(name_set_subseq0_msa_stockh, s_name_fasta, name_result_nhmmer_info, name_result_nhmmer_posision_new, 50)

    creating_txt_res_nhmmer_set(s_name_txt, name_result_nhmmer_posision_new, min_len_subseq = 0, max_len_subseq = 1000, N =  10 ** 6, name_res_txt = name_folder + 'second_subseq4.txt')
    set_subseq2 = creating_new_set_from_out_nhmmer(s_name_txt, name_result_nhmmer_posision_new, N = 10 ** 6, left_step = 0, right_step = 0)
    
    mean_len_subseq2 = 0
    min_x = 1000
    max_x = 0
    for i in range(len(set_subseq2)):
        mean_len_subseq2 += len(set_subseq2[i])
        if min_x > len(set_subseq2[i]):
            min_x = len(set_subseq2[i])
        if max_x < len(set_subseq2[i]):
            max_x = len(set_subseq2[i])
    mean_len_subseq2 = mean_len_subseq2 // len(set_subseq2)
    #print(mean_len_subseq2, len(set_subseq2), min_x, max_x)



def increase_len_subseq_msa(s_name_txt, name_result_nhmmer_posision, name_result_nhmmer_posision_new = 'increase2/name_result_pos.fa',left_step = 0, right_step = 0, min_len_subseq = 20, max_len_subseq = 1000, N = 150):
    name_folder = 'increase/'
    name_set_subseq0 = name_folder + 'name_old_set_fasta.fa'
    name_set_subseq0_msa = name_folder + 'name_old_set_fasta_msa.fa'
    name_set_subseq0_msa_stockh = name_folder + 'name_old_set_fasta_msa.sto'
    s_name_fasta = name_folder + 's_name.fa'
    name_result_nhmmer_info = name_folder + 'name_result_nhmmer_info.fa'
    name_position_subseq = name_folder + 'pos_sub2.fa'
    name_position_subseq_msa = name_folder + 'pos_sub2_msa.fa'
    #name_result_nhmmer_posision_new = name_folder + 'name_result_pos.fa'
    txt_into_fasta(s_name_txt, s_name_fasta)
    
    mean_len_subseq = 0

    creating_txt_res_nhmmer_set(s_name_txt, name_result_nhmmer_posision, min_len_subseq = 0, max_len_subseq = 1000, N =  10 ** 6, name_res_txt = name_folder + 'first_subseq.txt')

    set_subseq = creating_new_set_from_out_nhmmer(s_name_txt, name_result_nhmmer_posision, N = N, left_step = 0, right_step = 0, name_position_subseq = name_position_subseq)
    #process = subprocess.run(['muscle', '-in', name_position_subseq, '-out', name_position_subseq_msa], capture_output=True, text = True)
    #print(mean_len_subseq, len(set_subseq), min_x, max_x)
    muscle_msa(name_position_subseq,name_position_subseq_msa)
    f = open(name_position_subseq_msa, 'r')
    msa_pos = f.read()
    f.close()
    s = string_from_file(s_name_txt)
    set_of_subseq = []
    i = 0
    while i != -1:

        
        g = msa_pos.find(':', i + 1)
        begin, end = map(int, msa_pos[i + 1:g].split())

        i = msa_pos.find('>', i + 1)
        seq = msa_pos[g + 2:i - 1]
        seq = seq.replace('\n', '')
        
        j_left = 0
        k = seq[0]
        while k == '-':
            j_left += 1
            k = seq[j_left]
        begin -= j_left

        j_right = -1
        k = seq[-1]
        while k == '-':
            j_right -= 1
            k = seq[j_right]
        end -= j_right + 1
        
        if (begin - left_step) <= 0:
            begin = 0
        else:
            begin = begin - left_step
        
        if (end + right_step) >= (len(s) - 1):
            end = len(s) - 1
        else:
            end = end + right_step
        set_of_subseq.append(s[begin:end])

    string_set_of_subset = ''
    for i in range(len(set_of_subseq)):
        string_set_of_subset +=  '>i' + str(i) + ' iteration;\n' + str(set_of_subseq[i]) + '\n'
        #string_set_of_subset += '>i' + str(i + 1) + ' Name;\n' +str(set_of_subseq[i]) + '\n'
    file = open(name_set_subseq0,'w')
    file.write(string_set_of_subset)
    file.close()

    muscle_msa(name_set_subseq0, name_set_subseq0_msa)
    fasta_into_stockholm(name_set_subseq0_msa, name_set_subseq0_msa_stockh)
    nhmmer_on_stock_msa(name_set_subseq0_msa_stockh, s_name_fasta, name_result_nhmmer_info, name_result_nhmmer_posision_new, 50)

    name_res_txt = name_folder + 'second_subseq6.txt'
    creating_txt_res_nhmmer_set(s_name_txt, name_result_nhmmer_posision_new, min_len_subseq = 0, max_len_subseq = 1000, N =  10 ** 6, name_res_txt = name_res_txt)
    
    f = open(name_res_txt, 'r')
    res = f.readlines()
    f.close()
    
    return len(res) // 2

def increase(s_name_txt, name_result_nhmmer_posision, name_result_nhmmer_posision_new = 'increase2/name_result_pos.fa', mean_predict = 400, N = 150):
    name_folder = 'increase/'
    creating_txt_res_nhmmer_set(s_name_txt, name_result_nhmmer_posision, min_len_subseq = 0, max_len_subseq = 1000, N =  10 ** 6, name_res_txt = name_folder + 'first_subseq.txt')
    f = open(name_folder + 'first_subseq.txt', 'r')
    set_subseq = f.read()
    f.close()
    min_x = 1000
    max_x = 0
    for i in range(1, max(len(set_subseq), N), 2):
        if min_x > len(set_subseq[i]):
            min_x = len(set_subseq[i])
        if max_x < len(set_subseq[i]):
            max_x = len(set_subseq[i])
    diff = mean_predict - max_x
    res_number_all = []
    max_res_number = 0
    max_res_number_index = 0
    for i, left_diff in enumerate([j for j in range(0, diff + 1, diff // 10)]):
        name_result_nhmmer_posision_new_i = name_folder + str(i) +'_name_result_pos.fa'
        res_number = increase_len_subseq_msa(s_name_txt,  name_result_nhmmer_posision, name_result_nhmmer_posision_new = name_result_nhmmer_posision_new_i, left_step = left_diff, right_step = diff - left_diff, min_len_subseq = 20, max_len_subseq = 1000, N = 150)
        res_number_all.append(res_number)
        if res_number >= max_res_number:
            max_res_number_index = i
            max_res_number = res_number
    print(res_number_all)
    
    f = open(name_folder + str(max_res_number_index) +'_name_result_pos.fa', 'r')
    res = f.read()
    f.close()
    
    f = open(name_result_nhmmer_posision_new, 'w')
    f.write(res)
    f.close()

    creating_txt_res_nhmmer_set(s_name_txt, name_result_nhmmer_posision_new, min_len_subseq = 0, max_len_subseq = 1000, N =  10 ** 6, name_res_txt = name_folder + 'sw_subseq.txt')
    #creating_new_set_from_out_nhmmer(s_name_txt, name_result_nhmmer_posision, N = 10 ** 6, left_step = left_step, right_step = right_step)
    
    #process = subprocess.run(['muscle', '-in', name_old_set_fasta, '-out', name_old_set_fasta_msa], capture_output=True, text = True)

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
