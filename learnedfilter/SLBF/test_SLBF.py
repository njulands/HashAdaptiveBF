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

def test_shalla(b):
    positives, negatives = getShalla()
    begin_time = time()
    model = get_shalla_model()
    data = Data(positives, negatives)
    slbf = SLBF(model, data, b, is_train = True, using_cache = True, cache_path=data_dir+"tmp/SLBF/", cache_name = "shalla")
    end_time = time()
    run_time = end_time - begin_time
    print ('SLBF buildtime：',run_time)
    print("SLBF size(MB): %f" % (len(data.positives)*b/(8*1024*1024)+0.068))

    # slbf.test(data)
    begin_time = time()
    count = 1000
    for i in range(int(count/2)):
        slbf.check(positives[i])
    for i in range(int(count/2)):
        slbf.check(negatives[i])
    end_time = time()
    query_perkey_time = (end_time-begin_time) / count
    print("SLBF query time/key: ", query_perkey_time)

    
def test_shalla_cost(total_size = 200000, model_size = 10000):
    positives, negatives, neg_cost = getShalla()
    begin_time = time()
    b = (total_size - model_size) / len(positives)
    model = get_shalla_model()
    data = Data(positives, negatives)
    slbf = SLBF(model, data, b, is_train = False, using_cache = False, cache_path=data_dir+"tmp/SLBF/", cache_name = "shalla_cost")
    end_time = time()
    run_time = end_time - begin_time
    print ('SLBF buildtime：',run_time)
    negatives_cache = slbf.getNegativesCache(data.negatives)
    # cost test
    cost = 0
    count = 0
    for i in range(len(negatives)):
        if slbf.checkbyScore(negatives_cache[i], negatives[i]):
            count += 1
            cost += neg_cost[i]
    print("cost:", cost)
    print("fpr:", count / (len(negatives)))
    # begin_time = time()
    # count = 1000
    # for i in range(int(count/2)):
    #     slbf.check(positives[i])
    # for i in range(int(count/2)):
    #     slbf.check(negatives[i])
    # end_time = time()
    # query_perkey_time = (end_time-begin_time) / count
    # print("SLBF query time/key: ", query_perkey_time)

def test_ycsb(total_size = 200000, model_size = 10000):
    positives, negatives, neg_cost = getYCSB()
    begin_time = time()
    b = (total_size - model_size) / len(positives)
    model = get_ycsb_model()
    data = Data(positives, negatives)
    slbf = SLBF(model, data, b, is_train = False, using_cache = False, cache_path=data_dir+"/tmp/SLBF/", cache_name = "ycsb")
    end_time = time()
    run_time = end_time - begin_time
    print ('SLBF buildtime：',run_time)
    # negatives_cache = slbf.getNegativesCache(data.negatives)
    # # cost test
    # cost = 0
    # count = 0
    # for i in range(len(negatives)):
    #     if slbf.checkbyScore(negatives_cache[i], negatives[i]):
    #         count += 1
    #         cost += neg_cost[i]
    # print("cost:", cost)
    # print("fpr:", count / (len(negatives)))
    begin_time = time()
    count = 1000
    for i in range(int(count/2)):
        slbf.check(positives[i])
    for i in range(int(count/2)):
        slbf.check(negatives[i])
    end_time = time()
    query_perkey_time = (end_time-begin_time) / count
    print("SLBF query time/key: ", query_perkey_time)


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

