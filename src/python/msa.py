import numpy as np
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

def is_transition(c1,c2):
    return ((c1 == 'G' and c2 == 'A') or 
            (c1 == 'A' and c2 == 'G') or 
            (c1 == 'C' and c2 == 'T') or 
            (c1 == 'T' and c2 == 'C'))


def k2p_distances(sequences):
    """Given a list of N sequences, return an N x N numpy matrix
        of the pairwise kimura 2 parameter distances"""
    num_seqs = len(sequences)
    seq_len = len(sequences[0])
    d = np.zeros((num_seqs,num_seqs))
    for i in range(num_seqs):
        for j in range(i):
            s1 = sequences[i]
            s2 = sequences[j]
            p = 0
            q = 0
            for k in range(seq_len):
                c1, c2 = s1[k], s2[k]
                if c1 == c2 or c1 == '-' or c2 == '-':
                    # ignore
                    pass
                elif is_transition(c1, c2):
                    p += 1
                else: 
                    #transversion is the only other possibility
                    q += 1
            p /= seq_len
            q /= seq_len
            d[i,j] = -0.5 * np.log((1 - 2*p - q)* np.sqrt(1 - 2*q))
            d[j,i] = d[i,j]
    return d


