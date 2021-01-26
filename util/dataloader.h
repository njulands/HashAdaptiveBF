#ifndef DATALOADER
#define DATALOADER

#include <fstream>
#include <iostream>
#include <random>
#include <chrono>
#include <assert.h>
#include <algorithm>
#include "key.h"

// Dataset directory
#define SHALLA_PATH  "../data/shalla_cost"
#define YCSB_PATH  "../data/ycsbt"

#define RANDOM_KEYSTR_PATH "../util/randomKeyStr.txt"
#define RANDOM_COST_TYPE zipf

enum {uniform, hotcost, normal, zipf};

// data loader for datasets and random data
class dataloader{
    public:
        std::vector<Slice *> pos_keys_; // positive keys
        std::vector<Slice *> neg_keys_; // negative keys

        ~dataloader();
        bool load(std::string data_name_, bool using_cost_);
        bool load(std::string data_name_, bool using_cost_, std::string epcho);
        bool loadRandomKey(int positives_number, int negatives_number, bool using_cost_);

    private:
        bool loadShalla(bool using_cost_, std::string epcho);
        bool loadYCSB(bool using_cost_, std::string epcho);
};

// generate random keys and costs
class KeyBuilder{
    public:
        KeyBuilder();
        std::string GetKeyStr();
        bool ReadKeys(std::vector<Slice *> &v, int start_position_);
        void GenKeyStrAndToFile();
        void GenKeysUniformCosts(std::vector<Slice *> &keys, int interval);
        void GenKeysHotCosts(std::vector<Slice *> &keys, double hotNumberpro, int hotcost, int coldcost);
        void GenKeysNormalCosts(std::vector<Slice *> &keys, int u, int d);
        void GenKeysZipfCosts(std::vector<Slice *> &keys, double a, double c);
    private:
        std::vector<std::string> key_strs;
};

bool dataloader::load(std::string data_name_, bool using_cost_){
    if(data_name_ == "shalla") return loadShalla(using_cost_, "");
    else if(data_name_ == "ycsb") return loadYCSB(using_cost_, "");
    else return false;
}

bool dataloader::load(std::string data_name_, bool using_cost_, std::string epcho){
    if(data_name_ == "shalla") return loadShalla(using_cost_, epcho);
    else if(data_name_ == "ycsb") return loadYCSB(using_cost_, epcho);
    else return false;
}
bool dataloader::loadShalla(bool using_cost_, std::string epcho){
    std::cout << "shalla reading..."  << std::endl;
    std::ifstream is(SHALLA_PATH+epcho+".txt", std::ios::binary);
    if(is){
        std::string optype, keystr;
        double cost;
        while (is >> optype >> keystr >> cost) 
        { 
            Slice * key = new Slice();    
            key->str = keystr;   
            if(optype == "1") pos_keys_.push_back(key);
            else if(optype == "0"){    
                key->cost = using_cost_ ? cost : 1;
                neg_keys_.push_back(key); 
            }
        }
        is.close();
        return true;
    }
    is.close();
    return false;
}

bool dataloader::loadYCSB(bool using_cost_, std::string epcho){
    std::cout << "ycsb reading..."  << std::endl;
    std::ifstream is(YCSB_PATH+epcho+".txt");
    if(is){
        std::string optype, keystr;
        double cost;
        while (is >> optype >> keystr >> cost) 
        { 
            Slice * key = new Slice();    
            key->str = keystr;   
            if(optype == "FILTERKEY" || optype == "1") pos_keys_.push_back(key);
            else if(optype == "OTHERKEY" || optype == "0"){          
                key->cost = using_cost_ ? cost : 1;
                neg_keys_.push_back(key); 
            }
        }
        is.close();
        return true;
    }
    is.close();
    return false;
}

bool dataloader::loadRandomKey(int positives_number, int negatives_number, bool using_cost_){
    KeyBuilder kb;
    for(int i=0; i<positives_number; i++){
        Slice *key = new Slice();
        pos_keys_.push_back(key);
    }
    for(int i=0; i<negatives_number; i++){
        Slice *key = new Slice();
        neg_keys_.push_back(key);
    }
    if(!kb.ReadKeys(pos_keys_, 0)) return false; 
    if(!kb.ReadKeys(pos_keys_, 0)) return false; 
    if(using_cost_)
        switch (RANDOM_COST_TYPE)
        {
        case uniform:
            kb.GenKeysUniformCosts(neg_keys_, 5);
            break;
        case hotcost:
            kb.GenKeysHotCosts(neg_keys_, 0.01, 100, 1);
            break;
        case normal:
            kb.GenKeysNormalCosts(neg_keys_, 20, 10);
            break;
        case zipf:
            kb.GenKeysZipfCosts(neg_keys_, 1.25, 1.0);
            break;
        default:
            break;
        }
    return true;
}

dataloader::~dataloader(){
    for(Slice *key : pos_keys_)
        delete key;
    for(Slice *key : neg_keys_)
        delete key;
}

KeyBuilder::KeyBuilder(){
    std::fstream ifs(RANDOM_KEYSTR_PATH);
    if(ifs.is_open()){
        std::string s;
        while(std::getline(ifs,s))
            key_strs.push_back(s);
    }else{
        std::cout << "Keystr file not exists, generate again..." << std::endl;
        GenKeyStrAndToFile();
    }
    ifs.close();
}
std::string KeyBuilder::GetKeyStr(){
    int k=rand()%10+1;
    char arr[10];
    for(int i=1;i<=k;i++){
        int x,s;                         
        s=rand()%2;                     
        if(s==1) x=rand()%('Z'-'A'+1)+'A';        
        else x=rand()%('z'-'a'+1)+'a';      
        arr[i-1] = x;                  
    }
    return std::string(arr,k);
}

bool KeyBuilder::ReadKeys(std::vector<Slice *> &v, int start_position_){
    int size_ = v.size();
    if(start_position_ + size_ >= key_strs.size()) return false;
    for(int j=0; j<size_; j++)
        v[j]->str = key_strs[start_position_+j];
    return true;
}

void KeyBuilder::GenKeyStrAndToFile(){
    int gen_key_size_ = 200000;
    int i = key_strs.size();
    std::ofstream ofs(RANDOM_KEYSTR_PATH);
    while(i < gen_key_size_){
        if(0 == i%10000) std::cout << i << "keys have been created..." << std::endl;
        std::string str = GetKeyStr();
        if(std::find(key_strs.begin(), key_strs.end(), str) == key_strs.end()){
            key_strs.push_back(str);
            ofs << str << std::endl;
            i++;
        }
    }
    ofs.close();
}

void KeyBuilder::GenKeysUniformCosts(std::vector<Slice *> &keys, int interval){
    for(int i=0; i<keys.size(); i++)
        keys[i]->cost = 1+(i+1)*interval;
}
void KeyBuilder::GenKeysHotCosts(std::vector<Slice *> &keys, double hotNumberpro, int hotcost, int coldcost){
    int hotNumSize = hotNumberpro * keys.size();
    for(int i=0; i<keys.size(); i++){
        if(i <= hotNumSize) 
            keys[i]->cost = hotcost;
        else 
            keys[i]->cost = coldcost;
    }
}
void KeyBuilder::GenKeysNormalCosts(std::vector<Slice *> &keys, int u, int d){
    for(int i=0; i<keys.size(); i++){
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::default_random_engine generator(seed);
        std::normal_distribution<double> distribution(u, d);
        keys[i]->cost = distribution(generator);
    }
}

void KeyBuilder::GenKeysZipfCosts(std::vector<Slice *> &keys, double a, double c){
    int r = 10000;
    double pf[10000];
    double sum = 0.0;
    for (int i = 0; i < r; i++)        
        sum += c/pow((double)(i+2), a);  
    for (int i = 0; i < r; i++){ 
        if (i == 0)
            pf[i] = c/pow((double)(i+2), a)/sum;
        else
            pf[i] = pf[i-1] + c/pow((double)(i+2), a)/sum;
    }
     for (int i = 0; i < keys.size(); i++){
        int index = 0;
        double data = (double)rand()/RAND_MAX;  
        while (data > pf[index])  
            index++;
        keys[i]->cost = index;
    }
}

#endif