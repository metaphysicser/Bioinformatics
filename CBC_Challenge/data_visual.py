"""
-*- coding: utf-8 -*-
@Time    : 2022/7/17 11:40
@Author  : 夕照深雨
@File    : data_visual.py
@Software: PyCharm

Attention：

"""
import multiprocessing
import pickle
import random
from collections import Counter
import re
from multiprocessing import Process
from multiprocessing import Manager
import time

reference_path = "data/chr1-3_hg38.fasta"
pattern = re.compile('(@S[0-9]+\n)([AaTtGgCc]+)')  # re.I 忽略大小写
reference_pattern = re.compile('(>chr\d)([AaGgCcTtNn]*)')
N_parttern = re.compile('[Nn]+')
ACTG_parttern = re.compile('AaCcGgTt')


def match_fasta(str_dict, str):
    """
    使用正则表达式匹配，统计长read
    :param str_dict: 共享进程变量
    :return:
    """
    # m = pattern.match(str)  # 正则表达式匹配
    results = re.findall(pattern, str)
    for res in results:
        id = res[0][1:-1]
        content = res[1].replace("\n", "").upper()
        str_dict[id] = content


def data_process(index):
    """
    将原始数据处理为pkl文件，分30个小文件，每个小文件用字典统计
    格式为 S001 ：AATT...GGCC
    """
    start_time = time.time()  # 开始时间

    data_path = "F:\CBCdata/split.fastq." + str(index)
    # data_path = "data/sample/simulate_read.sample.align.fastq"
    print("loading data split.fastq." + str(index))
    # data = readfasta(data_path)[1]  # 读取原始数据文件
    with open(data_path, 'r') as f:
        data = f.read()
        str_dict = {}
        print(Counter(data))
        match_fasta(str_dict, data)
        #  写入pkl文件
        with open(f"data/data_process/read{index}.pkl", "wb") as f:
            pickle.dump(str_dict, f)

    end_time = time.time()  # 结束时间
    print("Run time : %.2f s" % (end_time - start_time))


def chr_process():
    data_path = "F:/CBCdata/chr1-3_hg38.fasta"
    with open(data_path, 'r') as f:
        data = f.read().replace('\n', "")
        str_dict = {}
        results = re.findall(reference_pattern, data)
        for res in results:
            id = res[0][1:]
            content = res[1].upper()
            str_dict[id] = content

        with open(f"data/chr.pkl", "wb") as f:
            pickle.dump(str_dict, f)


def split_chr():
    """
    把染色体中的 NNNNNNN 去除
    :return:
    """
    with open(f"data/chr.pkl", "rb") as f:
        data = pickle.load(f)

    result = {}

    for key, value in data.items():
        results = N_parttern.finditer(value)

        """
        经过观察，选取23作为阈值比较合适
        """
        split_begin = 0
        new_chr = []
        for match in results:
            new_value = {}
            first, last = match.span()
            last = last - 1
            length = last - first + 1
            """
            NNNN...NNNAATTTTCCGGGNNNNNNNNNNN
            """
            if length > 23:
                split_end = first
                if split_end > split_begin:
                    split_value = value[split_begin:split_end]
                    if "N" in split_value:
                        print("N" in split_value)
                        split_value.replace("N", random.choice("ACGT"))
                    new_value["data"] = split_value
                    new_value["begin_index"] = split_begin
                    new_chr.append(new_value)
                split_begin = last + 1
            else:
                replace_str = ""
                for i in range(length):
                    replace_str += random.choice('ATCG')
                value = value[:first] + replace_str + value[last + 1:]

        result[key + "+"] = new_chr
        # ----------------------------------------------------------------------
        """
        将染色体倒序进行处理
        """
        reverse_results = N_parttern.finditer(value[::-1])
        value = value[::-1]
        split_begin = 0
        new_chr = []
        for match in reverse_results:
            new_value = {}
            first, last = match.span()
            last = last - 1
            length = last - first + 1
            """
            NNNN...NNNAATTTTCCGGGNNNNNNNNNNN
            """
            if length > 23:
                split_end = first
                if split_end > split_begin:
                    split_value = value[split_begin:split_end]
                    if "N" in split_value:
                        print("N" in split_value)
                        print(value[split_end - 5:split_end + 5])
                        split_value.replace("N", random.choice("ACGT"))
                    new_value["data"] = split_value
                    new_value["begin_index"] = split_begin
                    new_chr.append(new_value)
                split_begin = last + 1
            else:
                replace_str = ""
                for i in range(length):
                    replace_str += random.choice('ATCG')
                value = value[:first] + replace_str + value[last + 1:]

        result[key + "-"] = new_chr

    with open(f"data/new_chr.pkl", "wb") as f:
        pickle.dump(result, f)


def decrease_length():
    """
    将染色体区段最大长度切割， 避免BWT变换时间复杂度过高
    :return:
    """
    with open("data/new_chr.pkl", "rb") as f:
        data = pickle.load(f)

    base = 5000000  # 最大长度为5000000
    split_chr = {}
    for chr, reads in data.items():
        add_reads = []
        for read in reads:
            #  如果该区段的染色体长度大于一千万，需要进行切割
            if len(read["data"]) > 10000000:
                begin_index = read["begin_index"]
                read_data = read["data"]
                split_num = len(read["data"]) // base + 1  # 切割的区段数目
                begin = 0
                end = base
                for i in range(split_num):
                    new_data = read_data[begin:end]
                    new_begin_index = begin_index + begin
                    short_read = {}
                    short_read["data"] = new_data
                    short_read["begin_index"] = new_begin_index
                    add_reads.append(short_read)
                    begin = end
                    end = end + base
            else:
                add_reads.append(read)
        reads = add_reads
        split_chr[chr] = reads

    with open("data/split_chr.pkl", "wb") as f:
        pickle.dump(split_chr, f)


#
# def split_chr(longread):
#     first_index = None
#     last_index = None
#     for index, i in enumerate(longread):
#         if i != "N":
#             first_index = index
#             break
#     for index, i in enumerate(longread[::-1]):
#         if i != "N":
#             last_index = index
#             break
#     return first_index, len(longread) - last_index - 1


if __name__ == '__main__':
    # for i in range(1, 4):  # 使用电脑最大cpu数
    #     p = Process(target=data_process, args=(i,))
    #     p.start()
    #     p.join()

    # split_chr()
    # with open("data/data_process/read7.pkl", "rb") as f:
    #     data = pickle.load(f)
    # for i in data.values():
    #     print(len(i))

    data_process(1)
