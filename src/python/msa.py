import numpy as np
import subprocess
import io
from subprocess import PIPE
import argparse
import os
def paup_distance(in_fasta, output_phylip, operating_system, model="logDet"):
    
    input_str = ("ToNEXUS format=FASTA fromFile={0} toFile=/tmp/nexus;Quit;")
    paup_input = input_str.format(in_fasta)

    cmd = [os.getcwd()+"/{0}", "-n"]
    if operating_system == "ubuntu":
        cmd[0] = cmd[0].format("paup4a164_ubuntu64")
    elif operating_system == "osx":
        cmd[0] = cmd[0].format("paup4a164_osx")
    else:
        assert False
    p = subprocess.Popen(cmd,stdin=PIPE, stdout=PIPE,stderr=PIPE)
    stdout,stderr = p.communicate(paup_input.encode()) 
    p.terminate()

    input_str = ("exe /tmp/nexus; DSet distance={0}; "
          "SaveDist format=RelPHYLIP file={1} triangle=both "
          "diagonal=yes; Quit;")
    
    paup_input = input_str.format(model, output_phylip)
    p = subprocess.Popen(cmd,stdin=PIPE, stdout=PIPE,stderr=PIPE)
    stdout,stderr = p.communicate(paup_input.encode())
    # print(stderr)  
    p.terminate()

    subprocess.call(["rm","/tmp/nexus"])

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--in-fasta', help='the folder containing the 20 subfolders of replicates')
    parser.add_argument('--out', help='desired output file name of the PHYLIP distance matrix')
    parser.add_argument('--os', choices=['ubuntu', 'osx'], help='your operating system, either ubuntu or osx')
    parser.add_argument('--model', default='logDet', help='the distance model, for input into PAUP')
    args = parser.parse_args()

    print("creating matrix...")
    paup_distance(args.in_fasta, args.out, args.os, args.model)
    print("done.")

##############################################################################
###                                                                        ###
### The code below is a slow python implementation of what PAUP does above ###
###                                                                        ###
##############################################################################

def read_fasta(file_name):
    """Given FASTA file name FILE_NAME, return two lists, 
        one contains the names of the sequences
        one contains the corresponding sequence"""
    names = []
    seqs = []
    with open(file_name) as file:
        line = file.readline()
        # read to the first '>'
        while line[0] != '>':
            line = file.readline()
        names.append(line[1:].strip())

        # start reading the first sequence
        line = file.readline()
        curr_seq = []
        while line:
            if line[0] == ';':
                # ignore this line
                pass
            elif line[0] == '>':
                # finished previous sequence
                seqs.append(''.join(curr_seq))
                curr_seq = []
                # add new name, after the '>'
                names.append(line[1:].strip())
            else:
                # remove newline and continue
                curr_seq.append(line.strip())
            line = file.readline()
        # append the last sequence
        seqs.append(''.join(curr_seq))

    # double check validity
    n_names, n_seqs = len(names), len(seqs)
    assert n_names == n_seqs, "%d names but %d sequences" % (n_names, n_seqs)
    assert n_seqs, "Empty file %s" % file_name
    N = len(seqs[0])
    for seq in seqs:
        # print(len(seq))
        assert len(seq) == N, "Sequences have different length"
    return names, seqs

def distance_func(metric):
    """decorator that takes function METRIC(seq1,seq2) and returns 
        a function distances(sequences)"""
    def distances(sequences):
        """ given alist of N sequences, return an N x N numpy matrix
            of the pairwise distances according to METRIC"""
        num_seqs = len(sequences)
        seq_len = len(sequences[0])
        d = np.zeros((num_seqs,num_seqs))
        for i in range(num_seqs):
            for j in range(i):
                assert len(sequences[i]) == len(sequences[j])
                d[i,j] = metric(sequences[i],sequences[j])
                d[j,i] = d[i,j]
        largest = np.nanmax(d)
        d[np.where(np.isnan(d))] = largest * 2.0
        return d       
    return distances

@distance_func
def jc_distances(s1,s2):
    """ return the jc distance between s1 and s2
        NaN if p >= 3/4 """
    p = 0
    for k in range(len(s1)):
        c1, c2 = s1[k], s2[k]
        p += not (c1 == c2 or c1 == '-' or c2 == '-')
    p /= len(s1)
    
    x = 1.0 - 4.0 / 3.0 * p
    if x <= 0:
        print("NAN")
        return np.nan
    return - 0.75 * np.log(x)

@distance_func
def logdet_distances(s1,s2):
    """Return the logdet distance between s1 and s2
        NaN if det(F) <= 0 """
    F = np.zeros((4,4))
    char2ind = {'A':0,'C':1,'T':2,'G':3}
    for k in range(len(s1)):
        c1,c2 = s1[k],s2[k]
        if c1 == '-' or c2 =='-':
            continue
        F[char2ind[c1],char2ind[c2]] += 1
    print(F)
    F /= len(s1)
    x = np.linalg.det(F)
    if x <= 0:
        print("NAN")
        return np.nan
    return - np.log(x)
    

def is_transition(c1,c2):
    return ((c1 == 'G' and c2 == 'A') or 
            (c1 == 'A' and c2 == 'G') or 
            (c1 == 'C' and c2 == 'T') or 
            (c1 == 'T' and c2 == 'C'))

@distance_func
def k2p_distances(s1,s2):
    """return the k2p distance between s1 and s2
        NaN if too large"""
    p = 0
    q = 0
    for k in range(len(s1)):
        c1, c2 = s1[k], s2[k]
        if c1 == c2 or c1 == '-' or c2 == '-':
            # ignore
            pass
        elif is_transition(c1, c2):
            p += 1
        else: 
            #transversion is the only other possibility
            q += 1
    p /= len(s1)
    q /= len(s1)
    
    x = (1.0 - 2.0*p - q)* np.sqrt(1.0 - 2.0*q)
    if x <= 0:
        print("NAN")
        return np.nan
    return -0.5 * np.log(x)


