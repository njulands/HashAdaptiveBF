import sys
sys.path.append("../lib") 
from GRUModel import GRUModel, get_shalla_model, get_ycsb_model
from LBF import LBF
import json
from dataload import *
from utils import *
from time import *
import os
import tensorflow as tf

def test_shalla_cost(total_size=100000, model_size = int(70*1024*8)):
    positives, negatives, neg_cost = getShalla()
    begin_time = time()
    data = Data(positives, negatives)
    model = get_shalla_model()
    lbf = LBF(model, data, using_Fpr = False, total_size = total_size, model_size = model_size, is_train = True)
    end_time = time()
    run_time = end_time - begin_time
    print ('LBF buildtime：',run_time)
    # cost test
    print("cost testing")
    cost = 0
    for i in range(len(negatives)):
        cost += lbf.check(negatives[i])*neg_cost[i]
    print("cost:", cost)

    # cost test
    print("cost testing")
    cost = 0
    count = 0
    total_cost = 0
    for i in range(len(negatives)):
        if lbf.check(negatives_cache[i], negatives[i]):
            count+=1
            cost += neg_cost[i]
        total_cost+=neg_cost[i]
    print("Weighted FPR:", 100*cost / total_cost)
    print("fpr:", 100*count/ (len(negatives)))



def test_ycsb(total_size=100000, model_size = int(70*1024*8)):
    positives, negatives, neg_cost = getYCSB()
    begin_time = time()
    data = Data(positives, negatives)
    model = get_ycsb_model()
    lbf = LBF(model, data, using_Fpr = False, total_size = total_size, model_size = model_size, is_train = True)
    end_time = time()
    run_time = end_time - begin_time
    print ('LBF buildtime：',run_time)
   
    # cost test by cache
    print("cost testing")
    cost = 0
    count = 0
    total_cost = 0
    for i in range(len(negatives)):
        if lbf.check(negatives[i]):
            count+=1
            cost += neg_cost[i]
        total_cost+=neg_cost[i]
    print("Weighted FPR:", 100*cost / total_cost)
    print("fpr:", 100*count/ (len(negatives)))


if __name__=='__main__':
    GPUFLAG = True
    # GPU config
    if GPUFLAG:
        os.environ["CUDA_VISIBLE_DEVICES"]="0" 
        sess = tf.Session(config=tf.ConfigProto(log_device_placement=True))

    '''Shalla test'''
    total_size = int(1.5*1024*1024*8)
    model_size = int(70*1024*8)
    test_shalla_cost(total_size = total_size, model_size = model_size)

    '''ycsb test'''
    # total_size = int(15*1024*1024*8)
    # model_size = int(3.57*1024*1024*8)
    # print("total_size:", total_size/(1024*1024*8))
    # test_ycsb(total_size = total_size, model_size = model_size)


