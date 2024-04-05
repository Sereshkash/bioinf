import seq_split as sq
import os

s_name_txt = 'example.txt'
s_name_fasta = 'example.fa'

sq.txt_into_fasta(s_name_txt, s_name_fasta)

N = 20
len_subseq = 410
E_val = 0.1
nubmer_set_of_subseq = 20

#for i in range(nubmer_set_of_subseq):

j = 1
name_folder = 'it' + str(j) + '/'

if not os.path.isdir(name_folder):
     os.mkdir(name_folder)

name_set_subseq0 = name_folder + 'example_subseq.fa'
name_set_subseq0_msa = name_folder + 'example_subseq_msa.fa'
name_set_subseq0_msa_stockh = name_folder + 'example_subseq_msa.sto'
name_result_nhmmer_info_0 = name_folder + 'name_result_nhmmer_info_0.fa'
name_result_nhmmer_info = name_folder + 'name_result_nhmmer_info.fa'
name_result_nhmmer_posision = name_folder + 'name_result_nhmmer_posision.fa'
name_set_res = name_folder + 'name_set_res.fa'

sq.creating_fasta_set(s_name_txt, N, len_subseq, name_set_subseq0)
sq.muscle_msa(name_set_subseq0, name_set_subseq0_msa)
sq.fasta_into_stockholm(name_set_subseq0_msa, name_set_subseq0_msa_stockh)
sq.nhmmer_on_stock_msa(name_set_subseq0_msa_stockh, s_name_fasta, name_result_nhmmer_info_0, name_result_nhmmer_posision, 0.1)
sq.creating_fasta_new_nhmmer_set(s_name_txt, name_result_nhmmer_posision, N, name_set_subseq0)

for i in range(10):
    #sq.creating_fasta_set(s_name_txt, N, len_subseq, name_set_subseq0)
    sq.muscle_msa(name_set_subseq0, name_set_subseq0_msa)
    sq.fasta_into_stockholm(name_set_subseq0_msa, name_set_subseq0_msa_stockh)
    sq.nhmmer_on_stock_msa(name_set_subseq0_msa_stockh, s_name_fasta, name_result_nhmmer_info, name_result_nhmmer_posision, 0.1)
    sq.creating_fasta_new_nhmmer_set(s_name_txt, name_result_nhmmer_posision, N, name_set_subseq0)

sq.nhmmer_on_stock_msa(name_set_subseq0_msa_stockh, s_name_fasta, name_result_nhmmer_info, name_result_nhmmer_posision, 0.01)
sq.creating_fasta_new_nhmmer_set(s_name_txt, name_result_nhmmer_posision, 10 ** 6, name_set_res)