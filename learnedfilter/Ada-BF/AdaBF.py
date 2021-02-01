import sys
sys.path.append("../lib") 
from dataload import *
from GRUModel import GRUModel
from utils import *
import pandas as pd
import numpy as np
import pickle
from keras.models import Sequential, load_model
import xxhash
from sklearn.utils import murmurhash3_32
from random import randint


def hashfunc(m):
    ss = randint(1, 99999999)
    def hash_m(x):
        return xxhash.xxh128(x,ss).intdigest() % m
    return hash_m

class Ada_BloomFilter():
    def __init__(self, n, hash_len, k_max):
        self.n = n
        self.hash_len = int(hash_len)
        self.h = []
        for i in range(int(k_max)):
            self.h.append(hashfunc(self.hash_len))
        self.table = np.zeros(self.hash_len, dtype=int)
    def insert(self, key, k):
        for j in range(int(k)):
            t = self.h[j](key)
            self.table[t] = 1
    def test(self, key, k):
        test_result = 0
        match = 0
        for j in range(int(k)):
            t = self.h[j](key)
            match += 1*(self.table[t] == 1)
        if match == k:
            test_result = 1
        return test_result

def R_size(count_key, count_nonkey, R0):
    R = [0]*len(count_key)
    R[0] = R0
    for k in range(1, len(count_key)):
        R[k] = max(int(count_key[k] * (np.log(count_nonkey[0]/count_nonkey[k])/np.log(0.618) + R[0]/count_key[0])), 1)
    return R

def Find_Optimal_Parameters(c_min, c_max, num_group_min, num_group_max, R_sum, train_negative, positive_sample):
    c_set = np.arange(c_min, c_max+10**(-6), 0.1)
    FP_opt = train_negative.shape[0]

    k_min = 0
    for k_max in range(num_group_min, num_group_max+1):
        for c in c_set:
            tau = sum(c ** np.arange(0, k_max - k_min + 1, 1))
            n = positive_sample.shape[0]
            hash_len = R_sum
            bloom_filter = Ada_BloomFilter(n, hash_len, k_max)
            thresholds = np.zeros(k_max - k_min + 1)
            thresholds[-1] = 1.1
            num_negative = sum(train_negative['score'] <= thresholds[-1])
            num_piece = int(num_negative / tau) + 1
            score = train_negative.loc[(train_negative['score'] <= thresholds[-1]), 'score']
            score = np.sort(score)
            for k in range(k_min, k_max):
                i = k - k_min
                score_1 = score[score < thresholds[-(i + 1)]]
                if int(num_piece * c ** i) < len(score_1):
                    thresholds[-(i + 2)] = score_1[-int(num_piece * c ** i)]

            url = positive_sample['url']
            score = positive_sample['score']

            for score_s, url_s in zip(score, url):
                ix = min(np.where(score_s < thresholds)[0])
                k = k_max - ix
                bloom_filter.insert(url_s, k)
            bloom_filter_opt = bloom_filter
            thresholds_opt = thresholds
            k_max_opt = k_max
    return bloom_filter_opt, thresholds_opt, k_max_opt

class AdaBF(object):
    def __init__(self, model, data, num_group = 16, c=3, R_sum = 22000000, is_train = True):
        self.model = model
        self.is_train = is_train
        self.num_group = num_group
        self.c = c
        self.R_sum = R_sum
        if self.is_train:
            self.fit(data.positives, data.negatives)
        else:
            self.model.loadmodel()
            print("load model success!!!")

        

    def create_Ada_BloomFilter(self, positive_sample, negative_sample):
        train_negative = negative_sample
        self.bloom_filter_opt, self.thresholds_opt, self.k_max_opt = Find_Optimal_Parameters(self.c, self.c, self.num_group, self.num_group, self.R_sum, train_negative, positive_sample)

    def check(self, item):
        score = self.model.predict(item)
        return self.checkbyScore(score, item)

    def checkbyScore(self, score, item):
        if score >= self.thresholds_opt[-2]:
            return True
        ix = min(np.where(score < self.thresholds_opt)[0])
        k = self.k_max_opt - ix
        return self.bloom_filter_opt.test(item, k)

    def gen_store_Scores(self, positives, negatives):
        print("gen_store_Scores")
        scores= np.array(self.model.predicts(negatives))
        # print(scores)
        negatives_dict = {'url': [p for p in negatives], 'score': scores.flatten().tolist()}
        negative_sample = pd.DataFrame(negatives_dict)

        scores=np.array(self.model.predicts(positives))
        # print(scores)
        positives_dict = {'url': [p for p in positives], 'score': scores.flatten().tolist()}
        positive_sample = pd.DataFrame(positives_dict)
        
        return positive_sample, negative_sample

    def fit(self, positives, negatives):
        shuffled = shuffle_for_training(negatives, positives)
        self.model.fit(shuffled[0], shuffled[1])
        print("Done fitting")



