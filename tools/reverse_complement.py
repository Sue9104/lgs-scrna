import argparse
parser = argparse.ArgumentParser( description='Get Reverse Complement Sequence')
parser.add_argument( 'infile', type=str, help='input  file')
parser.add_argument( 'outfile', type=str, help='output file')
args = parser.parse_args()

def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in dna[::-1]])

with open(args.infile, 'r') as infile:
    lines = infile.readlines()
    outlines = [reverse_complement(line.rstrip()) + '\n' for line in lines]
    with open(args.outfile, 'w') as outfile:
        outfile.writelines(outlines)
