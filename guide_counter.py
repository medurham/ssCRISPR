import sys
from Bio import SeqIO
from Bio.Seq import Seq
import csv


def write_guide_representation(filename, guide_counter):
    with open (filename, 'w') as fi:
        writer = csv.writer(fi)
        for key, value in guide_counter.items():
            writer.writerow([key[0],value])


def count_guide_sequences(filename, guide_key_list):
    guide_key_counter = {}
    for guide_key in guide_key_list:
        guide_key_counter[(guide_key, guide_key.split("_")[-1])] = 0
    counter_dic_keys = list(guide_key_counter.keys())
    
    for seq_record in SeqIO.parse(filename, "fastq"):
        my_dna = str(seq_record.seq)
        normal_seq = my_dna[1:21]
        left_shift = my_dna[0:20]       
        for counter_key in counter_dic_keys:
            if normal_seq in counter_key[1]:
                guide_key_counter[counter_key] += 1
            if left_shift in counter_key[1]:
                guide_key_counter[counter_key] += 1
        
    
    return guide_key_counter
    


def get_keys(filename):
    guide_dict = SeqIO.to_dict(SeqIO.parse(filename, "fasta"))
    guide_key_dic = {}
    for key, val in guide_dict.items():
        seq_str = str(val.seq)
        if seq_str not in guide_key_dic:
            guide_key_dic[seq_str] = key + "_" + seq_str
        else:
            guide_key_dic[seq_str] = guide_key_dic[seq_str] + "_" + key + "_" + seq_str
            if "MULT_" not in guide_key_dic[seq_str]:
                guide_key_dic[seq_str] = "MULT_" + guide_key_dic[seq_str]
                
    return list(guide_key_dic.values())


def main():
    fasta_filename = sys.argv[1]
    fastq_filename = sys.argv[2]
    output_file = sys.argv[3]
    #print("GOT THE FASTA FILE PATH: " + fasta_filename)
    guide_keys = get_keys(fasta_filename)
    #print(guide_keys)
    guide_key_counter = count_guide_sequences(fastq_filename, guide_keys)
    #print(guide_key_counter)
    write_guide_representation(output_file, guide_key_counter)    
    
if __name__ == "__main__":
    main()
    
