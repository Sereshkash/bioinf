import random
import subprocess
from Bio import SeqIO

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
        string_set_of_subset +=  '>i' + str(i) + ' iteration;\n' +str(set_of_subseq[i]) + '\n'
        #string_set_of_subset += '>i' + str(i + 1) + ' Name;\n' +str(set_of_subseq[i]) + '\n'
    file = open(name,'w')
    file.write(string_set_of_subset)
    file.close()

def muscle_msa(name_in, name_out = 'seqs1.afa'):
    process = subprocess.run(['muscle', '-in', name_in, '-out', name_out], capture_output=True, text = True)
    return process

def fasta_into_stockholm(name_in = 'example.fa', name_out = 'fasta_file6.sto'):
    records = SeqIO.parse(name_in, 'fasta')
    count = SeqIO.write(records, name_out, 'stockholm')

def nhmmer_on_stock_msa(name_set_msa = '1mult_sto.sto', name_seq = '1ex.fa', name_out_inf = 'result_hmm.fa', p_val = 0.1):
    command = 'nhmmer -o ' + name_out_inf + ' --noali --notextw --singlemx --dna --incE ' + str(p_val) + ' ' + name_set_msa + ' ' + name_seq
    process = subprocess.run(command, shell = True, capture_output = True, text = True)
    return process

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



nhmmer_on_stock_msa()

#creating_fasta_set(['ATTGC', 'AATTGC', 'ATGC'])

#nhmmer -o result_hmm.fa --noali --notextw --singlemx --dna --incE 0.1 1mult_sto.sto 1ex.fa
