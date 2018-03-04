import argparse
import gzip
import sys

def join(f1, f2):
        for line1, line2 in zip(f1, f2):
                i7=line1.decode('utf8').strip().split(sep=" ")[1].split(":")[3]
                seq1=next(f1).decode('utf8').strip()
                seq2=next(f2).decode('utf8').strip()
                print(seq1+i7+seq2)
                next(f1)
                next(f1)
                next(f2)
                next(f2)
        f1.close()
        f2.close()

if __name__ =="__main__":
        parser=argparse.ArgumentParser()
        parser.add_argument("-fastq1",nargs="+",help="read1")
        parser.add_argument("-fastq2",nargs="+",help="read2")
        args=parser.parse_args()
        file_list1 = list(args.fastq1)
        file_list2 = list(args.fastq2)
        for f in zip(file_list1, file_list2):
                f1=gzip.open(f[0])
                f2=gzip.open(f[1])
                join(f1, f2)
