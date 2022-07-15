"""
-*- coding: utf-8 -*-
@Time    : 2022/7/8 20:57
@Author  : 夕照深雨
@File    : HAlign_new.py
@Software: PyCharm

Attention：

"""
import multiprocessing
import pickle

from HAlign.FASTA import readfasta
from HAlign.PSA_Kband import PSA_AGP_Kband
from HAlign.TrieTree import Trie
import numpy as np

from HAlign.score import spscore


def getGapsLoc(strsAligned: list, markInsertion: list, idxC: int):
    """
    compute the gaps location
    """
    for str in strsAligned:
        i = 0
        counter = 0
        for c in str[0]:
            if c == '-':
                counter += 1
            else:
                markInsertion[i] = max(markInsertion[i], counter)
                counter = 0
                i += 1
            markInsertion[i] = max(markInsertion[i], counter)
    return markInsertion

def insertGap(mark: list, seq: str):
    """
    insert gaps to sequences

    Args:
        mark: gaps loc list
        seq: sequence

    Returns:
        sequence
    """
    res = ""
    length = len(mark)
    for i in range(length):
        res += "-" * mark[i]
        if i < length - 1:
            res += seq[i]
    return res

def insertSeqsGap(strsAligned: list, markInsertion: list, strs: list, idxC):
    """
    insert gaps to all sequences
    """
    S_aligned = [""] * (len(strs))
    S_aligned[idxC] = insertGap(markInsertion, strs[idxC])
    idx = 0
    for str2 in strsAligned:
        mark = [0] * (len(str2[0]) + 1)
        total = 0
        pi = 0
        pj = 0
        for c in str2[0]:
            if c == '-':
                total += 1
            else:
                mark[pi] = markInsertion[pj] - total
                pi += 1
                pj += 1
                while total != 0:
                    pi += 1
                    total -= 1
        mark[pi] = markInsertion[pj] - total
        if idx >= idxC:
            S_aligned[idx + 1] = insertGap(mark, str2[1])
        else:
            S_aligned[idx] = insertGap(mark, str2[1])
        idx += 1
    return S_aligned

def psa(strs: list, idxC: int):
        """
        align center sequence with others
        """
        strsAligned = []
        for i in range(len(strs)):
            if i != idxC:
                _, tmp1, tmp2, _, _ = PSA_AGP_Kband(strs[idxC], strs[i])
                strsAligned.append([tmp1, tmp2])
        return strsAligned

class HAlign_new:
    def __init__(self, strs):
        self.strs = strs
        self.homologous_threshold = 2.266  # 同源区段判定阈值
        # self.multi_lines = multiprocessing.cpu_count()
        self.multi_lines = 10
        self.split_strs()  # 分割序列
        self.find_center()  # 寻找中心序列
        self.pair_align()  # 两两对齐

    def merge_clips(self, split_strs):
        """
        合并被分割的序列
        :param split_strs: 被分割点序列
        :return:
        """
        clip_num = len(split_strs)
        total_count = np.sum(split_strs)
        homologous_clips = split_strs
        """
        寻找同源区段，大于阈值的是同源区段，否则是非同源区段，合并方法参考FM-index
        """
        homologous_clips = homologous_clips >= total_count * self.homologous_threshold / clip_num

        clip_type = None
        begin = None
        end = None

        result = []
        for index, i in enumerate(homologous_clips):
            if clip_type is None:
                clip_type = i
                begin = index
                end = index
            if clip_type == i:
                end = index
            if clip_type != i:
                temp = [begin, end]
                result.append(temp)
                begin = index
                end = index
                clip_type = i
        temp = [begin, end]
        result.append(temp)

        merge_strs = []
        merge_count = []
        """
        按照划分的结果合并序列片段，并计算每一个片段的中心序列
        """
        for index, sequence in enumerate(self.new_strs):
            tmp = []
            count = []
            for index2, clip in enumerate(result):
                tmp_strs = ""
                tmp_count = 0
                begin, end = clip
                for i in range(begin, end+1):
                    tmp_strs += sequence[i]
                    tmp_count += self.count_matrix[index][i]
                tmp.append(tmp_strs)
                count.append(tmp_count)

            merge_strs.append(tmp)
            merge_count.append(count)

        self.merge_strs = merge_strs
        self.merge_count = merge_count
        self.center_id = np.array(self.merge_count).argmax(axis=0)

        return result

    def pair_align(self):
        total_str = []
        for i in range(len(self.center_id)):
            strs = [x[i] for x in self.merge_strs]
            print(len(strs[i]))
            # 2. do pairwise alignments
            strsAligned = psa(strs, self.center_id[i])

            # 3. build the multiple alignment
            markInsertion = [0] * (len(strs[self.center_id[i]]) + 1)
            markInsertion = getGapsLoc(strsAligned, markInsertion, self.center_id[i])
            strsAligned = insertSeqsGap(strsAligned, markInsertion, strs, self.center_id[i])
            Value_SP = spscore(strsAligned)
            print(Value_SP)
            total_str.append(strsAligned)

        self.result = []
        for i in range(len(self.strs)):
            tmp = ""
            for j in range(len(self.center_id)):
                tmp += total_str[j][i]
            self.result.append(tmp)

        # for i in self.result:
        #     print(len(i))

        Value_SP = spscore(self.result)
        print(Value_SP)
        # 2. do pairwise alignments
        # strs = self.strs
        # idxC = self.center_id
        # strsAligned = psa(strs, idxC)
        #
        # # 3. build the multiple alignment
        # markInsertion = [0] * (len(strs[idxC]) + 1)
        # markInsertion = getGapsLoc(strsAligned, markInsertion, idxC)
        # strsAligned = insertSeqsGap(strsAligned, markInsertion, strs, idxC)
        #
        # # 4. compute the SP value
        # Value_SP = spscore(strsAligned)
        # print("SP : ", Value_SP)
        # print(Value_SP / len(self.strs))




    def find_center(self):
        count_matrix = np.zeros((len(self.strs), self.multi_lines))
        for i in range(len(self.strs)):
            for j in range(self.multi_lines):
                tree = self.trie_tree[i][j]
                for m in range(len(self.strs)):
                    if m !=i:
                        for k_str in self.str_k_partition[m][j]:
                            if tree.search(k_str):
                                count_matrix[i][j] += 1
        self.count_matrix = count_matrix  # 同源序列计数
        # total_matrix = count_matrix.sum(axis=1)  #
        split_matrix = np.max(count_matrix, axis=0) # 找到每一个区段同源序列最多的片段
        # with open("../data/psa_30.pkl", "rb") as f:
        #     true_matrix = np.array(pickle.load(f)).sum(axis=1)
        # corr = np.corrcoef(total_matrix, true_matrix, "valid")
        # print(corr[0][1])
        # print(total_matrix)
        # print(true_matrix)
        merge_clips = self.merge_clips(split_matrix)
        print(merge_clips)
        # self.center_id = total_matrix.argmax()



    def k_str_partition(self, str, k=20):
        assert len(str) > k
        t = Trie()
        str_partition = []
        for i in range(len(str) - k):
            tmp_str = str[i:i + k]
            str_partition.append(tmp_str)
            t.insert(tmp_str)
        return str_partition, t

    def split_strs(self):
        new_strs = []  # 切分后的序列
        str_k_partition = []  # 划分为k元组的序列
        trie_tree = []  # 前缀树
        for str in self.strs:
            tmp_strs = []
            tmp_k_partition = []
            tmp_tree = []
            str_len = len(str)
            for i in range(0, self.multi_lines):
                begin = i * str_len // self.multi_lines  # 开始
                end = (i + 1) * str_len // self.multi_lines  # 结束
                str_split = str[begin:end]
                tmp_strs.append(str_split)
                k_partition, trie = self.k_str_partition(str_split)  # 生成k元组和前缀树
                tmp_k_partition.append(k_partition)
                tmp_tree.append(trie)

            new_strs.append(tmp_strs)
            str_k_partition.append(tmp_k_partition)
            trie_tree.append(tmp_tree)
        self.new_strs = new_strs
        self.str_k_partition = str_k_partition
        self.trie_tree = trie_tree


if __name__ == '__main__':
    data = readfasta('../SequenceAlign/data/16srRNA(small).fasta')[1][:30]
    for i in range(len(data)):

        data[i] = data[i].upper()
    h = HAlign_new(data)
    print(h.center_id)

