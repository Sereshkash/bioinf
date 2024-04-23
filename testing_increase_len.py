import seq_split as sq
import os

name_folder = 'increase/'
name_nhmm_pos = name_folder + 'name_result_nhmmer_posision.fa'
s_name_txt = name_folder + 'example_new.txt'


sq.increase_len_subseq(s_name_txt, name_nhmm_pos, left_step = 200, right_step = 200, min_len_subseq =  200, max_len_subseq = 1000)