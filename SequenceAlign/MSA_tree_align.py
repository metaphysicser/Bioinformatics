import Levenshtein
import numpy as np
import time
from PSA_Kband import PSA_AGP_Kband
from FASTA import readfasta
from score import spscore
from tqdm import trange

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


class MSA_tree(object):
    def __init__(self, strs):
        self.strs = strs
        self.len = len(strs)
        self.tree = []
        self.build_tree()


    def pair_distance(self, strs: list):
        """
        align center sequence with others
        """
        print("begin to build distance matrix")
        dis_matrix = np.zeros((len(strs), len(strs)))
        for i in trange(len(strs)):
            for j in range(1, len(strs)):
                distance = Levenshtein.distance(strs[i], strs[j])
                dis_matrix[i][j] = dis_matrix[j][i] = distance
        return dis_matrix
    


    def compile_tree(self, root):
        leaf1, leaf2 = root
        if isinstance(leaf1,int) and isinstance(leaf2,int):
            _, self.strs[leaf1], self.strs[leaf2],_,_ = PSA_AGP_Kband(self.strs[leaf1], self.strs[leaf2])
            return [leaf1, leaf2]

        if isinstance(leaf1,int) and isinstance(leaf2,list):
            elements = self.compile_tree(leaf2)
            near1, near2 = self.find_nearest([leaf1],elements)
            _, self.strs[near1], self.strs[near2], A_gap_loc, B_gap_loc = PSA_AGP_Kband(self.strs[near1], self.strs[near2])
            self.alignset(near2, elements, B_gap_loc)
            return elements + [leaf1]

        if isinstance(leaf1,list) and isinstance(leaf2,int):
            elements = self.compile_tree(leaf1)
            near1, near2 = self.find_nearest(elements,[leaf2])
            _, self.strs[near1], self.strs[near2], A_gap_loc, B_gap_loc = PSA_AGP_Kband(self.strs[near1], self.strs[near2])
            self.alignset(near1, elements, A_gap_loc)
            return elements + [leaf2]

        if isinstance(leaf1, list) and isinstance(leaf2, list):
            elements1 = self.compile_tree(leaf1)
            elements2 = self.compile_tree(leaf2)
            near1, near2 = self.find_nearest(elements1,elements2)
            _, self.strs[near1], self.strs[near2],A_gap_loc, B_gap_loc = PSA_AGP_Kband(self.strs[near1], self.strs[near2])
            self.alignset(near1, elements1, A_gap_loc)
            self.alignset(near2, elements2, B_gap_loc)
            return elements1+elements2
    def find_nearest(self, set1,set2):
        min_value = 100000
        near1 = 0
        near2 = 0
        for i in set1:
            for j in set2:
                if self.init_distance[i][j] < min_value:
                    min_value = self.init_distance[i][j]
                    near1 = i
                    near2 = j
        return near1, near2
    def alignset(self, id, set, gap_loc):
        for seq in set:
            if seq != id:
                new = ""
                index = 0
                for gap in gap_loc:
                    if gap == 0:
                        try:
                            new += self.strs[seq][index]
                            index += 1
                        except Exception as e:
                            print("index: {}".format(index))
                            print("gap_length: {}".format(len(gap_loc)-sum(gap_loc)))
                            print("str length: {}".format(len(self.strs[seq])))
                            print(e)
                            new += "-"
                    else:
                        new += "-"
                self.strs[seq] = new




    def build_tree(self):
        """
        ?????????????????????????????? https://blog.csdn.net/m0_49960764/article/details/121495721
        :return:
        """
        sTime = time.time()
        print("begin to generate guide tree")
        self.dis_matrix = self.pair_distance(self.strs)  # ?????????????????????????????????????????? Levenshtein ??????
        self.init_distance = self.dis_matrix # ??????????????????????????????????????????????????????????????????????????????????????????
        self.map = {}  # ??????????????????index??????????????????????????????
        for i in range(len(self.strs)):
            self.map[i] = i

        self.r_dis = {} # ???????????????
        list_len = self.len
        for seq in range(self.len):
            r_distance = 0
            for j in range(self.len):
                if self.dis_matrix[seq][j] is not None:
                    r_distance += self.dis_matrix[seq][j]
            self.r_dis[seq] = r_distance  # ????????????????????????????????????
        if self.len > 2:  # ????????????????????????????????????????????????
            M_matrix = np.zeros((list_len, list_len))
            for i in self.r_dis.keys():
                for j in self.r_dis.keys():
                    if i != j:
                        """
                        ???????????????????????????????????????????????????????????????????????????????????????????????????
                        ????????????????????????????????????????????????????????????????????????????????????sp????????????????????????????????????????????????????????????????????????????????????20?????????????????????????????????
                        ??????????????? time???31.95s   SP???5669004
                        ????????????????????? time???27.44s SP???5667473
                        """
                        M_matrix[i][j] = M_matrix[i][j] = self.dis_matrix[i][j] # - (self.r_dis[i] + self.r_dis[j]) / (list_len - 2)

        min_index= 0
        min_value = 0
        value = 0
        for i in self.map.keys():
            for j in self.map.keys():
                if i != j:
                    value += M_matrix[i][j]
            if value < min_value:
                min_value = value
                min_index = i
        idxC = min_index
        # 2. do pairwise alignments
        strsAligned, gap_counter = psa(self.strs, idxC)

        gap = np.array(gap_counter)

        gap_matrix = np.tile(np.max(gap, axis=0), (gap.shape[0], 1)) - gap

        new_strs = []

        for i in range(len(gap_counter)):
            new_str = insert_gap(gap_matrix[i, :], strsAligned[i])
            new_strs.append(new_str)

        Value_SP = spscore(new_strs)
        eTime = time.time()
        print("Run time : %.2f s" % (eTime - sTime))
        print("SP : ", Value_SP)


                        


        # while self.len > 2:
        #     min_value = 10000
        #
        #     for i in self.map.keys():
        #         for j in self.map.keys():
        #             if i != j and M_matrix[i][j] < min_value:
        #                 min_value = M_matrix[i][j]
        #                 pair = [i, j]
        #     l1, l2 = pair[0], pair[1]  # ????????????????????????????????????
        #
        #     self.map[l1] = [self.map[l1], self.map.pop(l2)] # [1,2,[3,4]] ->[1,[2,[3,4]]] 2???[3???4]????????????
        #
        #     self.len -= 1
        #
        #     for i in self.map.keys():
        #         if i != l1:
        #             # dis = self.dis_matrix[l1][i] if self.dis_matrix[l1][i] > self.dis_matrix[l2][i] else self.dis_matrix[l2][i]
        #             # self.dis_matrix[l1][i] = self.dis_matrix[i][l1] = dis
        #             """
        #             ??????????????????????????????????????????????????????
        #             """
        #             self.dis_matrix[l1][i] = self.dis_matrix[i][l1] = (self.dis_matrix[l1][i] + self.dis_matrix[l2][i] -  self.dis_matrix[l1][l2]) / 2
        #
        #     M_matrix = np.zeros((list_len, list_len))
        #     # ??????M??????
        #     for i in self.map.keys():
        #         for j in  self.map.keys():
        #             if i != j:
        #                 M_matrix[i][j] = M_matrix[i][j] = self.dis_matrix[i][j] # - (self.r_dis[i] + self.r_dis[j]) / (list_len - 2)
        #
        # l1, l2 = self.map.keys()  # ???????????????????????????
        #
        # self.tree = [self.map[l1], self.map[l2]] # ???????????????
        #
        # self.compile_tree(self.tree)  # ??????????????????
        #
        # Value_SP = spscore(self.strs)
        # eTime = time.time()
        #
        # print("Run time : %.2f s" % (eTime - sTime))
        # print("SP : ", Value_SP)







# data = ["AAATTT","TTCAA","FUJHJK","FYUYG","ghvefajih","FGYTVYUIFVG"]
#
data = readfasta('data/16srRNA(small).fasta')[1][:20]
for i in range(len(data)):
    data[i] = data[i].upper()
# data = readfasta('data/dna500.fasta')[1]

m = MSA_tree(data)


