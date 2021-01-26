## Adapted from https://www.geeksforgeeks.org/bloom-filters-introduction-and-python-implementation/
import math
import mmh3
from bitarray import bitarray

class BestBloomFilter(object):
    def __init__(self, items_count, fp_prob, get_digest):
        self.fp_prob = fp_prob
        self.size = self.get_size(items_count, fp_prob)
        self.hash_count = self.get_hash_count(self.size, items_count)
        self.bit_array = bitarray(self.size)
        self.bit_array.setall(0)
        self.get_digest = get_digest

    def add(self, item):
        for i in range(self.hash_count):
            digest = self.get_bucket(item, i)
            self.bit_array[digest] = 1

    def check(self, item):
        for i in range(self.hash_count):
            digest = self.get_bucket(item, i)
            if (self.bit_array[digest] == 0):
                return False
        return True

    def get_size(self,n,p):
        m = -(n * math.log(p))/(math.log(2)**2)
        return int(m)

    def get_hash_count(self, m, n):
        k = (m/n) * math.log(2)
        return int(k)

    def get_bucket(self, item, i):
        return self.get_digest(item, i) % self.size

class BloomFilter(object):
    def __init__(self, k, m, get_digest):
        self.size = max(0, int(m))
        self.hash_count = max(0, int(math.ceil(k)))
        self.bit_array = bitarray(self.size)
        self.bit_array.setall(0)
        self.get_digest = get_digest

    def add(self, item):
        for i in range(self.hash_count):
            digest = self.get_bucket(item, i)
            self.bit_array[digest] = 1

    def check(self, item):
        for i in range(self.hash_count):
            digest = self.get_bucket(item, i)
            if (self.bit_array[digest] == 0):
                return False
        return True
    def set_k(self, k):
        self.hash_count = k
    def insert_many(self, key_arr):
        for key in key_arr:
            self.add(key)
        

    def get_bucket(self, item, i):
        return self.get_digest(item, i) % self.size

    def query_many(self, key_arr):
        queries = [False for i in range(len(key_arr))]
        for i in range(len(key_arr)):
            queries[i] = self.check(key_arr[i])
        return queries

