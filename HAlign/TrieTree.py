"""
-*- coding: utf-8 -*-
@Time    : 2022/7/6 20:49
@Author  : 夕照深雨
@File    : TrieTree.py
@Software: PyCharm

Attention：

"""


class Trie:
    def __init__(self):
        """
        Initialize the structure here.
        """
        self.root = {}
        self.sequence_end = -1

    def insert(self, sequence):
        """
        Inserts a sequence into the trie.
        :type sequence: str
        :rtype: void
        """
        curNode = self.root
        for c in sequence:
            if not c in curNode:
                curNode[c] = {}
            curNode = curNode[c]

        curNode[self.sequence_end] = True

    def search(self, sequence):
        """
        Returns if the sequence is in the trie.
        :type sequence: str
        :rtype: bool
        """
        curNode = self.root
        for c in sequence:
            if not c in curNode:
                return False
            curNode = curNode[c]

        # Doesn't end here
        if self.sequence_end not in curNode:
            return False

        return True

    def startsWith(self, prefix):
        """
        Returns if there is any sequence in the trie that starts with the given prefix.
        :type prefix: str
        :rtype: bool
        """
        curNode = self.root
        for c in prefix:
            if not c in curNode:
                return False
            curNode = curNode[c]

        return True


