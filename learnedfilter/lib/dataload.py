import pandas as pd
import json
"""
load shalla.json
"""
data_dir = "../../data/"
'''
load shalla_cost.txt
'''
def getShalla():
    print("Shalla_Cost reading...")
    data=pd.read_csv(data_dir+"shalla_cost.txt", sep=' ', names=['type','ID', 'cost'], header=None)
    positives = data.loc[data['type']==1].iloc[:,1].values
    negatives = data.loc[data['type']==0].iloc[:,1].values
    neg_cost = data.loc[data['type']==0].iloc[:,2].values
    print("positives size: ", len(positives))
    print("negatives size: ", len(negatives))
    return positives, negatives, neg_cost

'''
load ycsb
'''
def getYCSB():
    print("YCSB reading...")
    data=pd.read_csv(data_dir+"ycsbt.txt", sep=' ', names=['type','ID', 'cost'], header=None)
    positives = data.loc[data['type']=="FILTERKEY"].iloc[:,1].values
    negatives = data.loc[data['type']=="OTHERKEY"].iloc[:,1].values
    neg_cost = data.loc[data['type']=="OTHERKEY"].iloc[:,2].values
    print("positives size: ", len(positives))
    print("negatives size: ", len(negatives))
    return positives, negatives, neg_cost

# getYCSB()