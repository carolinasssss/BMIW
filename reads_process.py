import sys 
import gzip 
import statistics
import os
import matplotlib.pyplot as plt
import numpy as np

from Bio import SeqIO 

def parsing():
    if len(sys.argv) != 4:
        print("Usage: python seq.py <illumina_file_1> <illumina_file_2> <nanopore_file")
        sys.exit(1)

    for file in sys.argv[1:]:
        if not (file.endswith(".fastq") or file.endswith(".fastq.gz")):
            print("wrong file format")
            sys.exit(1)

    return sys.argv[1], sys.argv[2], sys.argv[3]

def gc_content(seq):
    seq = seq.upper()
    gc = seq.count('G') + seq.count('C')
    gc_content = (gc/len(seq)) * 100 

    if len(seq) == 0:
        return 0

    return gc_content

def illumina(illumina_file):
    if illumina_file.endswith('.gz'):
        handle = gzip.open(illumina_file, 'rt')
    else:
        handle = open(illumina_file, 'r')

    count = sum(1 for _ in SeqIO.parse(handle,"fastq"))
    handle.close()

    print(f" {os.path.basename(illumina_file)}: {count} reads")

    return count

def nanopore(nanopore_file):
    if nanopore_file.endswith('.gz'):
        handle = gzip.open(nanopore_file, 'rt')
    else: 
        handle = open(nanopore_file,'r')

    lengths = []
    nt = 0
    gc_values = []

    for record in SeqIO.parse(handle, "fastq"):
        length = len(record.seq)
        lengths.append(length)

        nt += length

        gc_values.append(gc_content(str(record.seq)))

    handle.close()

    if lengths: 
        mean_len = statistics.mean(lengths)
        median_len = statistics.median(lengths)
        gc_np = statistics.mean(gc_values)

        print(f" nucleotides: {nt}")
        print(f" avr. length: {mean_len:.2f}")
        print(f" median of length: {median_len:.2f}")
        print(f" gc content: {gc_np:.2f}")

        histogram(lengths)

        return nt, mean_len, median_len, gc_np

def histogram(lengths):


    plt.hist(lengths, bins=30, color='blue')
    plt.xlabel('single read lentgh (bp)')
    plt.ylabel('number of reads')
    plt.title('nanopore read length distribution')
    #plt.ylim(0,1000)
    #upper_limit = np.percentile(lengths, 99)
    #plt.xlim(0, upper_limit)
    plt.savefig('histogram.png')
    plt.show()

def main():

    print("-" * 10)
    illumina1_file, illumina2_file, nanopore_file = parsing()

    print("ILLUMINA")
    count1 = illumina(illumina1_file)
    count2 = illumina(illumina2_file)
    print(f" illumina reads: {count1+ count2}")

    print("-" * 10)
    print("NANOPORE")
    nanopore(nanopore_file)


if __name__ == "__main__":
    main()



