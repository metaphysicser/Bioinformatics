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
from multiprocessing.pool import Pool

from tqdm import tqdm
from BWT import BWT

RANGE = 4  # 多进程划分

def BWT_transform(chr_name, begin, data):
    """
    将分割的染色体进行BWT变换
    :param chr_name: 名称
    :param begin: 开始位置
    :param data: 数据
    :return:
    """
    values = data[chr_name]
    length = len(values)
    #  将数据进行划分， 便于多进程运行
    values = values[begin * length // RANGE: (begin + 1) * length // RANGE]
    for value in tqdm(values):
        index = value["begin_index"]
        name = chr_name + "_" + str(index)
        if os.path.exists(f"data/BWT/{name}"+".pkl"):
            continue
        else:
            long_read = value["data"]
            BWT_trie = BWT(long_read)  # BWT变换
            with open(f"data/BWT/{name}"+".pkl", "wb") as f:
                pickle.dump(BWT_trie, f)
    return

def append_res(new_res, results):
    """
    判断匹配的区段能否成为最终结果
    :param new_res: 新的匹配区段
    :param results: 已有的匹配区段
    :return: 结果
    """
    new_res_value = list(new_res.values())[0]
    new_length = sum([x[1] for x in new_res_value])  # 新区段匹配长度
    new_range = [new_res_value[0][0], new_res_value[-1][0] + new_res_value[-1][1]] # 新区段匹配范围
    tmp = []  # 暂存影响区段
    count = 0  # 影响区段长度之和
    for match in results:
        res_value = list(match.values())[0]
        res_range = [res_value[0][0], res_value[-1][0]+res_value[-1][1]]
        length = sum([x[1] for x in res_value])
        # 判断和已有区段是否有冲突
        if new_range[0] > res_range[1] or new_range[1] < res_range[0]:
            continue
        else:
            tmp.append(match)
            count += length

    if len(tmp) == 0:  # 没有任何冲突
        results.append(new_res)
    elif count <= new_length:  # 判断影响区段和该区段哪个覆盖长度更多
        for i in tmp:
            results.remove(i)  # 移除影响区段
        results.append(new_res)  # 加入新区段

    return results


def find_rough_range(chr, res_dict):
    """
    找到区段的大致区间
    :param chr: 染色体
    :return:
    """
    chr_path = "data/BWT/" + chr  # 经过BWT变换后的染色体路径
    read_path = "data/data_process/"  # 长read序列路径
    res_path = "data/rough_res/"

    if not os.path.exists(res_path):
        os.makedirs(res_path)

    chr_list = os.listdir(chr_path)  # 获取文件列表
    read_list = os.listdir(read_path)

    # for read in tqdm(read_list):
    for read in tqdm(["sample.unalign.pkl"]):
        # 读取长序列
        with open(read_path + read, "rb") as f:
            long_read = pickle.load(f)
        for chr in tqdm(chr_list):
            # 读取染色体文件
            chr_name = chr.split('.')[0]
            with open(os.path.join(chr_path, chr), "rb") as f:
                chr_BWT = pickle.load(f)
            # 遍历每一个长序列
            for key, value in long_read.items():
                # 寻找能和长read匹配的染色体区段
                # [[index, length, start]...], index是长read的开始位置， length是匹配的长度， start是染色体开始位置
                res, _ = chr_BWT._select_CommonStrings(value)
                if len(res) == 0:
                    continue
                elif key not in res_dict.keys():
                    res_dict[key] = [{chr_name: res}]
                else:
                    res_dict[key] = append_res({chr_name: res}, res_dict[key])



def align():
    """
    确定比对区间，然后进行序列比对
    :return:
    """
    with open("data/data_process/sample.align.pkl", "rb") as f:
        read_data = pickle.load(f)
    with open("data/rough_res/sample.align.pkl", "rb") as f:
        data = pickle.load(f)

    count = 0

    for read_name, read_range in data.items():
        range_num = len(read_range)
        if range_num > 3:
            continue
        str_length = len(read_data[read_name])
        range_len = 0
        for match in read_range:
            for key, value in match.items():
                length = sum([x[1] for x in value])
                range_len += length
        rate = range_len/str_length
        if rate < 0.01:
            continue




if __name__ == '__main__':
    # res_dict = Manager().dict()
    # pool = Pool(processes=10)
    # for i in ["chr1+", "chr1-", "chr2+", "chr2-", "chr3+", "chr3-"]:
    #     pool.apply_async(find_rough_range, (i,res_dict,))  # 非阻塞
    #
    # print("Started processes")
    # pool.close()  # 需要关闭进程池，防止池其他任务的提交，注意！这里不是关闭进程。简单来说就是关掉了屋外的大门，但是各个房间在运行。
    # pool.join()  # 等待进程池里面的进程运行完
    #
    # res_dict = dict(res_dict)
    #
    # with open("data/rough_res/sample.unalign.pkl" , "wb") as f:
    #     pickle.dump(res_dict, f)



    #
    # process_list = []
    # # find_rough_range("chr2-", res_dict)
    # for i in ["chr1+", "chrtarget=find_rough_range, args=(i,res_dict,)1-", "chr2+", "chr2-", "chr3+", "chr3-"]:  # 开启5个子进程执行fun1函数
    #     p = Process()  # 实例化进程对象
    #     p.start()
    #     process_list.append(p)
    #
    # for i in process_list:
    #     p.join()
    #
    # with open("data/rough_res/sample.align.pkl" , "wb") as f:
    #     pickle.dump(res_dict, f)
    #
    # with open("data/data_process/sample.align.pkl", "rb") as f:
    #     data = pickle.load(f)
    # for i in data.values():
    #     print(len(i))



    with open("data/rough_res/sample.align.pkl", "rb") as f:
        data = pickle.load(f)

    align()






    # with open("data/new_chr.pkl", "rb") as f:
    #     data = pickle.load(f)
    #
    # process_list = []
    # for i in ["chr3+", "chr3-"]:  # 开启5个子进程执行fun1函数
    #     for j in range(RANGE):
    #         p = Process(target=BWT_transform, args=(i,j, data,))  # 实例化进程对象
    #         p.start()
    #         process_list.append(p)
    #
    # for i in process_list:
    #     p.join()




