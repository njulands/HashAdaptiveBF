from AdaBF import *
from time import *
from GRUModel import GRUModel, get_shalla_model, get_ycsb_model
import os
import tensorflow as tf


def test_shalla_cost(num_group=16, c=3, R_sum=20000000):
    positives, negatives, neg_cost = getShalla()
    begin_time = time()
    data = Data(positives, negatives)
    model = get_shalla_model()
    ada_bf = AdaBF(model, data, num_group = num_group, c=c, R_sum =R_sum, is_train = False)
    positive_sample, negative_sample = ada_bf.gen_store_Scores(data.positives, data.negatives)
    ada_bf.create_Ada_BloomFilter(positive_sample, negative_sample)
    print("AdaBF size(MB): %f" % (R_sum/(8*1024*1024)+0.068))
    end_time = time()
    run_time = end_time - begin_time
    print ('Ada-BF buildtime：',run_time)
    
    print("cost testing")
    cost = 0
    count = 0
    total_cost = 0
    for i in range(len(negatives)):
        if ada_bf.check(negatives[i]):
            count+=1
            cost += neg_cost[i]
        total_cost+=neg_cost[i]
    print("Weighted FPR:", 100*cost / total_cost)




def test_ycsb(num_group=16, c=3, R_sum=20000000):
    positives, negatives, neg_cost = getYCSB()
    begin_time = time()
    data = Data(positives, negatives)
    model = get_ycsb_model()
    ada_bf = AdaBF(model, data, num_group = num_group, c=c, R_sum =R_sum, is_train = False)
    positive_sample, negative_sample = ada_bf.gen_store_Scores(data.positives, data.negatives)
    ada_bf.create_Ada_BloomFilter(positive_sample, negative_sample)
    end_time = time()
    run_time = end_time - begin_time
    print ('Ada-BF buildtime：',run_time)
    # cost test
    cost = 0
    count = 0
    total_cost = 0
    for i in range(len(negatives)):
        if ada_bf.check(negatives[i]):
            count+=1
            cost += neg_cost[i]
        total_cost+=neg_cost[i]
    print("Weighted FPR:", 100*cost / total_cost)

if __name__=='__main__':
    GPUFLAG = True
    # GPU config
    if GPUFLAG:
        os.environ["CUDA_VISIBLE_DEVICES"]="0" 
        sess = tf.Session(config=tf.ConfigProto(log_device_placement=True))
    
    '''shalla test'''
    total_size = int(1.5*1024*1024*8)
    model_size = int(70*1024*8)
    num_group = 16
    c=8
    R_sum = total_size - model_size
    test_shalla_cost(num_group=num_group, c=c, R_sum=R_sum)

    '''ycsb test'''
    # total_size = int(1.25*1024*1024*8)
    # model_size = int(3.57*1024*1024*8)
    # num_group = 8
    # c=3
    # R_sum = total_size - model_size
    # test_ycsb(num_group=num_group, c=c, R_sum=R_sum)



