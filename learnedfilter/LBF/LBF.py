import sys
sys.path.append("../lib") 
from BloomFilter import BestBloomFilter, BloomFilter
import math
import random
from utils import *
import pickle


class LBF(object):
    def __init__(self, model, data, using_Fpr = True, fp_rate = 0.01, total_size = 100000, model_size = int(70*1024*8), is_train = True):
        self.model = model
        self.threshold = 0.9
        self.using_Fpr = using_Fpr
        self.is_train = is_train
        (s1, s2) = split_negatives(data)
        if self.is_train:
            self.fit(data.positives, data.negatives)
        else:
            self.model.loadmodel()
            print("load model success!!!")

        if using_Fpr:
            self.fp_rate = float(fp_rate)
            self.create_best_bloom_filter(data, s2)
        else:
            self.m = total_size - model_size
            self.create_bloom_filter(data, s2)


    def check(self, item):
        if self.model.predict(item) > self.threshold:
            return True
        return self.bloom_filter.check(item)


    def create_bloom_filter(self, data, test_negatives):
        print("Creating bloom filter")
        print("predicts positives")
        positives_preds = self.model.predicts(data.positives)
        positives_preds_sort = positives_preds.copy()
        positives_preds_sort.sort()
        tn_preds = self.model.predicts(test_negatives)
        tn_preds.sort()
        thresh_o = 1
        thresh_max = 1
        k_o = 4
        for k in range(4,30):
            n = int((self.m / k) * math.log(2))
            p = math.exp(- (self.m / n) * (math.log(2)**2))
            if n >= len(positives_preds_sort):
                continue
            positive_thresh = positives_preds_sort[n]
            tn_idx = math.ceil((len(tn_preds) * (1 - p)))
            if tn_idx <= 0 or tn_idx>=len(tn_preds):
                continue
            tn_thresh =  tn_preds[tn_idx]
            if abs(tn_thresh - positive_thresh) < thresh_max:
                thresh_max = abs(tn_thresh - positive_thresh)
                thresh_o = positive_thresh
                k_o = k

        self.bloom_filter = BloomFilter(k_o, self.m)
        n_o = int((self.m / k_o) * math.log(2))
        self.threshold = thresh_o
        print("thresh_o",thresh_o)
        print("k_o", k_o)
        print("n", n_o)
        false_negatives = []
        for i in range(len(data.positives)):
            if positives_preds[i] <= self.threshold:
                false_negatives.append(data.positives[i])
        # false_negatives = data.positives[:n_o]
        print("Number of false negatives at bloom time", len(false_negatives))

        for fn in false_negatives:
            self.bloom_filter.add(fn)
        
        print("Created bloom filter")

    def create_best_bloom_filter(self, data, test_negatives):
        print("Creating bloom filter")
        self.get_threshold(test_negatives)
        print("threshold: %f" %self.threshold)
        false_negatives = []
        preds = self.model.predicts(data.positives)
        print(len(preds))
        for i in range(len(data.positives)):
            if preds[i] <= self.threshold:
                false_negatives.append(data.positives[i])
        print("Number of false negatives at bloom time", len(false_negatives))
        self.bloom_filter = BestBloomFilter(
            len(false_negatives),
            self.fp_rate / 2,
            string_digest
        )
        for fn in false_negatives:
            self.bloom_filter.add(fn)
        print("Created bloom filter")
    
    def fit(self, positives, negatives):
        shuffled = shuffle_for_training(negatives, positives)
        self.model.fit(shuffled[0], shuffled[1])
        print("Done fitting")

    def get_threshold(self, test_negatives):
        fp_index = math.ceil((len(test_negatives) * (1 - self.fp_rate/2)))
        predictions = self.model.predicts(test_negatives)
        predictions.sort()

        self.threshold = predictions[fp_index]

     

