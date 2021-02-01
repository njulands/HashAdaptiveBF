import sys
sys.path.append("../lib") 
from GRUModel import GRUModel, get_shalla_model, get_ycsb_model
from SLBF import SLBF
import json
from dataload import *
from utils import *
from time import *
import os
import tensorflow as tf
    
def test_shalla_cost(total_size = 200000, model_size = 10000):
    positives, negatives, neg_cost = getShalla()
    begin_time = time()
    b = (total_size - model_size) / len(positives)
    model = get_shalla_model()
    data = Data(positives, negatives)
    slbf = SLBF(model, data, b, is_train = False)
    end_time = time()
    run_time = end_time - begin_time
    print ('SLBF buildtime：',run_time)
    # cost test
    cost = 0
    count = 0
    total_count = 0
    for i in range(len(negatives)):
        if slbf.check(negatives[i]):
            count += 1
            cost += neg_cost[i]
        total_count += neg_cost[i]
    print("weighted FPR:", cost / total_count)

def test_ycsb(total_size = 200000, model_size = 10000):
    positives, negatives, neg_cost = getYCSB()
    begin_time = time()
    b = (total_size - model_size) / len(positives)
    model = get_ycsb_model()
    data = Data(positives, negatives)
    slbf = SLBF(model, data, b, is_train = False)
    end_time = time()
    run_time = end_time - begin_time
    print ('SLBF buildtime：',run_time)
    # cost test
    cost = 0
    count = 0
    total_count = 0
    for i in range(len(negatives)):
        if slbf.check(negatives[i]):
            count += 1
            cost += neg_cost[i]
        total_count += neg_cost[i]
    print("weighted FPR:", cost / total_count)


if __name__=='__main__':
    GPUFLAG = True
    # GPU config
    if GPUFLAG:
        os.environ["CUDA_VISIBLE_DEVICES"]="0" 
        sess = tf.Session(config=tf.ConfigProto(log_device_placement=True))
    '''Shalla test'''
    # b = 10
    # test_shalla(b)
    
    '''Shalla Cost test'''
    total_size = int(1.375*1024*1024*8)
    model_size = int(70*1024*8)
    test_shalla_cost(total_size = total_size, model_size = model_size)

    '''ycsb test'''
    # total_size = int(20*1024*1024*8)
    # model_size = int(3.57*1024*1024*8)
    # test_ycsb(total_size = total_size, model_size = model_size)

