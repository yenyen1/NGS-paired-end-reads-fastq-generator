# usage: python fastq_generator.py fasta -t1 target_name -o output -r 150 -f 450 -d 30
# python version above 3
# author: yenyen.wang
# created by 2020.10.12

import argparse, sys
from os import path
from pysam import FastxFile
from random import getrandbits, normalvariate

def parseFasta(filepath:str, ref_dict:dict) -> dict:
    fasta = FastxFile(filepath)
    for ref in fasta:
        ref_dict[ref.name] = ref.sequence
    return ref_dict
    
def to_string(first:str, second:str) -> str:
    forth = "I"*len(second)
    ss = "@" + first + "\n"
    ss += second + "\n"
    ss += "+\n" + forth + "\n"
    return ss

def getRevComplement(seq:str) -> str:
    ss = ""
    for c in seq:
        if c == "A":
            ss += "T"
        elif c== "T":
            ss += "A"
        elif c == "G":
            ss += "C"
        elif c == "C":
            ss += "G"
        else:
            ss += "N"
    return ss[::-1]

def getGap(coverage: float, read_len: int, frag_len: int, ref_len: int) -> int:
    gap = (ref_len - frag_len + 1) * read_len * 2 / (ref_len * coverage)
    return int(round(gap))

def generateFastq(is_hetero:bool, no:int, target:str, refSeq:str, output:str, read_len:int, frag_len:int, depth:int):
    if len(refSeq) > 0:
        f1 = open(output+"1.fq", 'a+')
        f2 = open(output+"2.fq", 'a+')

        gap = getGap(depth/2, read_len, frag_len, len(refSeq)) if is_hetero else getGap(depth, read_len, frag_len, len(refSeq))
        refpos = int( (len(refSeq) - frag_len + 1) /gap )
        
        first_readname = ""
        second_read1 = ""
        second_read2 = ""
        for i in range(refpos):
            first_readname = "simulate_" + target + "_" + str(no) + "_read" + str(i)
            shift = round(normalvariate(0.,1.))
            pos = i*gap + shift if (i > 0 and i < refpos -1) else  i*gap
            
            second_read1 = refSeq[pos:pos+read_len]
            second_read2 = refSeq[pos+frag_len-read_len: pos+frag_len]
            orient = bool(getrandbits(1))
            if  orient:
                f1.write(to_string(first_readname + "/1", second_read1))
                f2.write(to_string(first_readname + "/2", getRevComplement(second_read2)))
            else:
                f1.write(to_string(first_readname + "/1", getRevComplement(second_read2)))
                f2.write(to_string(first_readname + "/2", second_read1))

        f1.close()
        f2.close()
    else:
        print(target + " is not found in fasta file.")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("FASTA_PATH", help = "A fasta file which includes the reference sequences that you're interested.")
    parser.add_argument("-t1", "--target1", required = True, type = str, help = "Specified your target reference name. (must be identical with the name in fasta file)")
    parser.add_argument("-t2", "--target2", default = "*", type = str, help = "Specified another target reference name if it's heterozigous")
    parser.add_argument("-o", "--OUTPUT", required = True, type = str, help = "Output file name (without extension)")

    parser.add_argument("-r", "--read_len", default = 150, type = int, help = "The average length of reads. (Default: 150)")
    parser.add_argument("-f", "--frag_len", default = 450, type = int, help = "The average length of DNA fragments. (Default: 450)")
    parser.add_argument("-d", "--depth", default = 30, type = int, help = "The average alignment depth. (Default: 30)")

    args = parser.parse_args()
    if path.exists(args.OUTPUT):
        print("Error: output file already exists.")
        sys.exit()

    hla_dict = {}
    hla_dict = parseFasta(args.FASTA_PATH, hla_dict)

    is_hetero = (args.target2 != "*")
    if is_hetero:
        print("generate heterozygous fastq for " + args.target1 + " " + args.target2 + " with read length " 
            + str(args.read_len) + ", fragment length " + str(args.frag_len) + " and depth " + str(args.depth) + "X" )
        generateFastq(is_hetero, 1, args.target1, hla_dict.get(args.target1, ""), args.OUTPUT, args.read_len, args.frag_len, args.depth)
        generateFastq(is_hetero, 2, args.target2, hla_dict.get(args.target2, ""), args.OUTPUT, args.read_len, args.frag_len, args.depth)
    else:
        print("generate homozygous fastq for " + args.target1 + " with read length " 
            + str(args.read_len) + ", fragment length " + str(args.frag_len) + " and depth " + str(args.depth) + "X" )
        generateFastq(is_hetero, 1, args.target1, hla_dict.get(args.target1, ""), args.OUTPUT, args.read_len, args.frag_len, args.depth)         
    print("fin.")

if __name__ == '__main__':
    main()