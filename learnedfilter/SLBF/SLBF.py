import sys
sys.path.append("../lib") 
from BloomFilter import BloomFilter
import math
import random
from utils import *
import mmh3
import pickle
class SLBF(object):
    def __init__(self, model, data, b, threshold = 0.9, is_train = True):
        self.model = model
        self.threshold = threshold
        self.is_train = is_train
        (s1, s2) = split_negatives(data)
        if self.is_train:
            print("Training model with train, dev, positives", len(s1), len(s2), len(data.positives))   
            self.fit(data.positives, s1)
        else:
            self.model.loadmodel()
            print("load model success!!!")
        
        positives_predictions = self.measure_learned_Fn(data.positives);
        self.measure_learned_Fp(s2);
        self.create_bloom_filter(b, data)
        self.insert_many(data.positives, positives_predictions)

    def check(self, item):
        if not self.initial_filter.check(item):
            return False
        if self.model.predict(item) > self.threshold:
            return True
        if self.backup_filter.check(item):
            return True
        return False

    def getNegativesCache(self, negatives):
        negatives_preds = self.model.predicts(negatives)

    def create_bloom_filter(self, b, data):
        print("Creating bloom filter")
        alpha = 0.6185
        big_expr = (self.learned_Fp) / (1 - self.learned_Fp) / (1 / self.learned_Fn - 1)
        b_2 = min(b, self.learned_Fn * math.log(big_expr, alpha))
        b_1 = b - b_2
        n = len(data.positives)
        print("b1, b2", b_1, b_2)
        initial_k = int(math.ceil(0.6931 * b_1))
        backup_k = int(math.ceil(0.6931 * b_2))
        self.initial_filter = BloomFilter(initial_k, b_1 * n, string_digest)
        self.backup_filter = BloomFilter(backup_k, b_2 * n, string_digest)
        print("Created bloom filter")

    def insert_many(self, positives, positives_predictions):
        print("insert data")
        self.initial_filter.insert_many(positives)
        neg = 0
        neg_keys = []
        for i in range(len(positives_predictions)):
            if positives_predictions[i] <= self.threshold:
                neg += 1
                neg_keys.append(positives[i])
        print("neg_keys num: %d" % neg)
        if len(neg_keys) == 0:
            new_k = 1
        else:
            new_k = int(math.ceil(0.6931 * self.backup_filter.size / len(neg_keys)))
        self.backup_filter.set_k(new_k)
        self.backup_filter.insert_many(neg_keys)

    def test(self, data):
        print("test fpr ...")
        fpr_result = self.query_many(data.negatives)
        fpr = sum(fpr_result) / len(data.negatives)
        print('False positive: %f' % fpr)


    def query_many(self, data):
        initial_results = self.initial_filter.query_many(data)
        learned_result = self.model.predicts(data);
        backup_result = self.backup_filter.query_many(data)

        final_results = [False for i in range(len(data))]

        for i in range(len(final_results)):
            final_results[i] = initial_results[i] and (learned_result[i]>self.threshold or backup_result[i])
        
        return final_results

    def measure_learned_Fp(self, ts2):
        predictions = self.model.predicts(ts2)
        tot = 0
        for val in predictions:
            if val > self.threshold:
                tot += 1
        self.learned_Fp = float(tot) / len(ts2)
        print("learned_Fp: %f", self.learned_Fp)

    def measure_learned_Fn(self, tp):
        predictions = self.model.predicts(tp)
        tot = 0
        for val in predictions:
            if val <= self.threshold:
                tot += 1
        self.learned_Fn = float(tot) / len(tp)
        print("learned_Fn: %f", self.learned_Fn)
        return predictions
    
    def fit(self, positives, negatives):
        shuffled = shuffle_for_training(negatives, positives)
        self.model.fit(shuffled[0], shuffled[1])
        print("Done fitting")

