import seq_split as sq
import os
import time

# name_folder = 'increase/'
# name_nhmm_pos = name_folder + 'name_result_nhmmer_posision.fa'
# s_name_txt = name_folder + 'example_new.txt'


name_folder = 'increase/'
name_nhmm_pos = name_folder + 'name_result_nhmmer_posision.fa'
s_name_txt = name_folder + 'example_new.txt'
s_name_fasta = name_folder + 'example_new.fa'
name_set_subseq1 = name_folder + '1example_subseq.fa'
name_set_subseq1_msa = name_folder + '1example_subseq_msa2.fa'
name_set_subseq1_msa_stockh = name_folder + '1example_subseq_msa2.sto'

name_result_nhmmer_posision = name_folder + 'name_result_nhmmer_posision_new.fa'

# for i in range(1, 21):
#     f = open(name_folder + '1example_subseq.fa', 'r')
#     full = f.readlines()
#     f.close()

#     name_set_subseq1 = name_folder + '1example_subseq' + str(i) + '.fa'
#     name_set_subseq1_msa = name_folder + '1example_subseq_msa' + str(i) + '.fa'
#     name_result_nhmmer_posision = name_folder + 'name_result_nhmmer_posision_new' + str(i) + '.fa'

#     f = open(name_set_subseq1, 'w')
#     s = '\n'.join(full[:6*i + 1])
#     f.write(s)
#     f.close()

#     name_result_nhmmer_info_0 = name_folder + 'name_result_nhmmer_info_0_new_f' + str(i) + '.fa'

#     sq.txt_into_fasta(s_name_txt, s_name_fasta)
#     sq.muscle_msa(name_set_subseq1, name_set_subseq1_msa)
#     sq.fasta_into_stockholm(name_set_subseq1_msa, name_set_subseq1_msa_stockh)
#     sq.nhmmer_on_stock_msa(name_set_subseq1_msa_stockh, s_name_fasta, name_result_nhmmer_info_0, name_result_nhmmer_posision, 20)

time_st = time.time()
sq.increase(s_name_txt, name_nhmm_pos, name_result_nhmmer_posision, mean_predict = 400, N = 200)
time_f = time.time()
print('time = {}, N = 200'.format(round((time_f - time_st) / 60, 2)))

#sq.increase_len_subseq_msa(s_name_txt, name_nhmm_pos, name_result_nhmmer_posision, left_step = 200, right_step = 200, min_len_subseq =  20, max_len_subseq = 1000)