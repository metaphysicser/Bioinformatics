import Levenshtein
import numpy as np
import time
from PSA_Kband import PSA_AGP_Kband
from FASTA import readfasta
from score import spscore
from tqdm import trange

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
        生成指导树，参考博客 https://blog.csdn.net/m0_49960764/article/details/121495721
        :return:
        """
        sTime = time.time()
        print("begin to generate guide tree")
        self.dis_matrix = self.pair_distance(self.strs)  # 计算序列之间的距离矩阵，使用 Levenshtein 距离
        self.init_distance = self.dis_matrix # 之后每对齐序列，距离矩阵都会更改，因此记录下最开始的距离矩阵
        self.map = {}  # 距离矩阵中的index和实际序列之间的映射
        for i in range(len(self.strs)):
            self.map[i] = i

        self.r_dis = {} # 净分化距离
        list_len = self.len
        for seq in range(self.len):
            r_distance = 0
            for j in range(self.len):
                if self.dis_matrix[seq][j] is not None:
                    r_distance += self.dis_matrix[seq][j]
            self.r_dis[seq] = r_distance  # 记录每个序列的净分化距离
        if self.len > 2:  # 如果只有两条序列，就没必要计算了
            M_matrix = np.zeros((list_len, list_len))
            for i in self.r_dis.keys():
                for j in self.r_dis.keys():
                    if i != j:
                        """
                        计算分类群的距离矩阵，在实验中，我发现当我删去距离公式的后半部分。
                        也就是该行代码的后半部分，代码运行速度有了一些提升，但是sp值变化没有特别明显，我认为有趣便保留了下来。我进行了使用20条序列进行了简单实验：
                        博客思路： time：31.95s   SP：5669004
                        我修改的代码： time：27.44s SP：5667473
                        """
                        M_matrix[i][j] = M_matrix[i][j] = self.dis_matrix[i][j] # - (self.r_dis[i] + self.r_dis[j]) / (list_len - 2)

        while self.len > 2:
            min_value = 10000

            for i in self.map.keys():
                for j in self.map.keys():
                    if i != j and M_matrix[i][j] < min_value:
                        min_value = M_matrix[i][j]
                        pair = [i, j]
            l1, l2 = pair[0], pair[1]  # 找到距离最近都两个分类群

            self.map[l1] = [self.map[l1], self.map.pop(l2)] # [1,2,[3,4]] ->[1,[2,[3,4]]] 2和[3，4]距离最短

            self.len -= 1

            for i in self.map.keys():
                if i != l1:
                    # dis = self.dis_matrix[l1][i] if self.dis_matrix[l1][i] > self.dis_matrix[l2][i] else self.dis_matrix[l2][i]
                    # self.dis_matrix[l1][i] = self.dis_matrix[i][l1] = dis
                    """
                    更新距离矩阵，可以尝试不同定义的距离
                    """
                    self.dis_matrix[l1][i] = self.dis_matrix[i][l1] = (self.dis_matrix[l1][i] + self.dis_matrix[l2][i] -  self.dis_matrix[l1][l2]) / 2

            M_matrix = np.zeros((list_len, list_len))
            # 更新M矩阵
            for i in self.map.keys():
                for j in  self.map.keys():
                    if i != j:
                        M_matrix[i][j] = M_matrix[i][j] = self.dis_matrix[i][j] # - (self.r_dis[i] + self.r_dis[j]) / (list_len - 2)

        l1, l2 = self.map.keys()  # 最后只剩两个分类群

        self.tree = [self.map[l1], self.map[l2]] # 生成指导树

        self.compile_tree(self.tree)  # 递归进行比对

        Value_SP = spscore(self.strs)
        eTime = time.time()

        print("Run time : %.2f s" % (eTime - sTime))
        print("SP : ", Value_SP)







# data = ["AAATTT","TTCAA","FUJHJK","FYUYG","ghvefajih","FGYTVYUIFVG"]
#
data = readfasta('data/2019nCoVR_20200301/2019nCoVR_20200301.fasta')[1][:20]
for i in range(len(data)):
    data[i] = data[i].upper()
# data = readfasta('data/dna500.fasta')[1]

m = MSA_tree(data)


