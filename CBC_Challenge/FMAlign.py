"""
-*- coding: utf-8 -*-
@Time    : 2022/7/18 20:30
@Author  : 夕照深雨
@File    : FMAlign.py
@Software: PyCharm

Attention：

"""
import os
import pickle

from multiprocessing import Manager, Process
from tqdm import tqdm
import threading
from CBC_Challenge.BWT import BWT

K = 40



def BWT_transform(chr_name,data, BWT_dict, BWT_index_dict):
    values = data[chr_name]
    for value in tqdm(values):
        index = value["begin_index"]
        name = chr_name + "_" + str(index)
        if os.path.exists(f"data/BWT/{name}"+".pkl"):
            continue
        else:
            long_read = value["data"]
            BWT_trie = BWT(long_read)
            with open(f"data/BWT/{name}"+".pkl", "wb") as f:
                pickle.dump(BWT_trie, f)




def FMAlign(long_read, k = K):
    assert len(long_read) > k
    str_partition = []
    for i in range(len(long_read) - k):
        tmp_str = long_read[i:i + k]
        str_partition.append(tmp_str)


    pass

if __name__ == '__main__':
    # with open("data/new_chr.pkl", "rb") as f:
    #     data = pickle.load(f)
    # values = data["chr1+"][0]
    # b = BWT(values["data"])
    # with open("data/test.pkl", "wb") as f:
    #     pickle.dump(b, f)
    # with open("data/test.pkl", "rb") as f:
    #     data = pickle.load(f)
    # res,_ = data._select_CommonStrings("ATTTTAGGGCT")
    # print(res)

    BWT_dict = Manager().dict()
    BWT_index_dict = Manager().dict()
    with open("data/new_chr.pkl", "rb") as f:
        data = pickle.load(f)
    for i in ["chr1+", "chr1-","chr2+", "chr2-","chr3+", "chr3-"]:  # 使用电脑最大cpu数
        p = threading.Thread(target=BWT_transform, args=(i, data, BWT_dict, BWT_index_dict,))
        p.start()
        # p.join()

    # with open(f"data/chr_BWT.pkl", "wb") as f:
    #     pickle.dump(BWT_dict, f)
    #
    # with open(f"data/chr_BWT_index.pkl", "wb") as f:
    #     pickle.dump(BWT_index_dict, f)