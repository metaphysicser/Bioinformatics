import numpy as np
import time
from PSA_Kband import PSA_AGP_Kband
from FASTA import readfasta
from score import spscore

def findCenterSeq(strs: list):
    """
    to find the center sequence
    Returns:
        index of center sequence
    """
    s_psa = [[-float('Inf')] * len(strs) for _ in range(len(strs))]

    for i in range(len(strs)):
        if i % (len(strs) // 10 + 1) == 0.0:
            print(" => ", end="", flush=True)

        for j in range(len(strs)):
            if j > i:
                tmp, _, _,_,_ = PSA_AGP_Kband(strs[i], strs[j])
                s_psa[i][j] = s_psa[j][i] = tmp
            elif j == i:
                s_psa[i][i] = 0

    idxC = np.argmax(np.sum(s_psa, axis=0))

    return idxC

def gap_count(gap):
    new_gap = [0 for i in range((len(gap)-sum(gap)+1))]
    counter = 0
    index = 0

    for i in gap:
        if i == 1:
            counter += 1

        elif i == 0:
            new_gap[index] = counter
            counter = 0
            index += 1
    new_gap[index] = counter
    return new_gap


def psa(strs: list, idxC: int):
    """
    align center sequence with others
    """
    strsAligned = []
    gap_counter = []
    for i in range(len(strs)):
        if i != idxC:
            _, tmp1, tmp2,gap_A,gap_B = PSA_AGP_Kband(strs[idxC], strs[i])
            strsAligned.append(tmp2)
            print(len(tmp2))
            gap_counter.append(gap_count(gap_A))
        else:
            strsAligned.append(strs[idxC])
            gap_counter.append([0 for i in range(len(strs[idxC])+1)])
            print(len(strs[idxC]))

    return strsAligned, gap_counter

def insert_gap(gap,strs):
    length = len(gap)
    new_str = ""
    for i in range(len(strs)):
        if i<length and gap[i] != 0:
            new_str += "-"*gap[i]
        new_str += strs[i]
    return new_str

def MSA_star(strs):
    sTime = time.time()
    # 1. to find the center sequence
    print("-----------RUN-----------")
    print("Loading", end="")
    idxC = findCenterSeq(strs)
    print("Loaded")
    print("center seq:", ''.join(strs[idxC]))

    # 2. do pairwise alignments
    strsAligned,gap_counter = psa(strs, idxC)

    gap = np.array(gap_counter)

    gap_matrix = np.tile(np.max(gap, axis=0), (gap.shape[0],1)) - gap

    new_strs = []

    for i in range(len(gap_counter)):
        new_str = insert_gap(gap_matrix[i,:],strsAligned[i])
        new_strs.append(new_str)


    Value_SP = spscore(new_strs)
    eTime = time.time()
    print("Run time : %.2f s" % (eTime - sTime))
    print("SP : ", Value_SP)
    # for str in strsAligned:
    #     print(' '.join(str))
    print("-----------END-----------")

    return Value_SP, strsAligned


data = readfasta('data/16srRNA(small).fasta')[1][:50]
for i in range(len(data)):
    data[i] = data[i].upper()
m = MSA_star(data)

