#ifndef HABF
#define HABF

#include <iostream>
#include <cmath>
#include <assert.h>
#include <algorithm>
#include <set>
#include <exception> 
#include <string.h>
#include <chrono> 
#include "../util/key.h"
#include "../util/hashutil.h"

#define HASH_NUM 3
#define ALLOCATION_RATIO 4 // Allocation ratio between the size of bloom and hashexpressor
#define HASHEXPRESSOR_ENTRY_HASH xxhash64// entry_hash_ of hashexpressor

namespace habf{

const std::vector<HashFamily> hashset_ = {BOB, xxhash128, cityhash64, SuperFastHash, BKDRHash, murmurhash, APHash};


struct HABFilterKey{
    Slice * data_;
    uint8_t hash_map_[HASH_NUM];
};

struct Tuple_V{
    bool singleflag_;
    HABFilterKey * keyid;
    Tuple_V(): singleflag_(true), keyid(NULL){};
};

struct HashExpressorCell{
    bool endbit;
    uint8_t hash_idx_;
};

class BloomFilter{
    public:
        void CreateNewFilter(uint32_t &bit_size_);
        void AddAll(std::vector<HABFilterKey> &keys);
        bool Contain(HABFilterKey &key);

        void initByV_(Tuple_V * V_);

        uint32_t getSize(){return 8 * data_.size();};

    private:
        std::string data_;
        HashSetUtil hsu_;  
};

class HashExpressor{
    public:
        HashExpressor(uint8_t k_, HashFamily entry_hash_);
        void CreateNewFilter(uint32_t bit_size_);
        bool Add(HABFilterKey &key);
        void Add(HABFilterKey &key, std::vector<uint8_t> &insert_sequence_);
        bool GetInsertSequence(HABFilterKey &key, std::vector<uint8_t> &insert_sequence_);
        void Query(HABFilterKey &key);

        bool AddWithMaxPattarn(HABFilterKey &key);
        bool GetInsertSequenceWithMaxPattarn(HABFilterKey &key, std::vector<uint8_t> &insert_sequence_);
        void QueryWithMaxPattarn(HABFilterKey &key);
        
        int getSize(){return 8 * data_.size();};
        
        void AssignDefaultHashset_(HABFilterKey *keys_, uint32_t &size_);
        void AssignDefaultHashset_(HABFilterKey &key);

    private:
        void setbit(const int &bitinx);
        bool getbit(const int &bitinx);
        void getCell(const int &cellinx, HashExpressorCell &cell);
        void setCell(const int &cellinx, const HashExpressorCell &cell);

        void FindInsertSequence(HABFilterKey &key, std::vector<uint8_t> &hash_map_, std::vector<uint8_t> &insert_sequence_, int left_, int right_);
        void clearHashSet(std::vector<uint8_t> &hash_set);


        uint32_t cell_num_;
        uint8_t cell_size_;
        uint8_t k_;
        HashFamily entry_hash_;
        HashSetUtil hsu_;
        std::string data_;
        std::vector<uint8_t> default_hashset_;


};

class HABFilter{
    public:
        HABFilter(float bits_per_key_, int pos_count_); // In order to fine-tune the size of the filter space for comparing, we allow bits_per_key_ to be float
        void AddAndOptimize(const std::vector<Slice *> &pos_keys_, const std::vector<Slice *> &neg_keys_);
        bool Contain(Slice &key);

    private:
        void Optimize(Tuple_V * V_, std::vector<std::vector<HABFilterKey *>> &T_, std::vector<HABFilterKey *> &CQ_);
        void DetectConflict(Tuple_V * V_, std::vector<std::vector<HABFilterKey *>> &T_, int bitpos_, 
        std::vector<HABFilterKey *> &cost_keys_, double &cost);

        void UpdateStructure(Tuple_V * V_, std::vector<std::vector<HABFilterKey *>> &T_,  HABFilterKey *ck,  HABFilterKey *optimal_pk_, 
        std::vector<uint8_t> &optimal_insert_sequence_, std::pair<uint32_t, uint32_t> &pre_next_,
        std::vector<HABFilterKey *> &optimal_cost_keys_);

        uint8_t k_;
        BloomFilter bloom_;
        HashExpressor hash_expressor_;
        HashSetUtil hsu_;        
};
uint32_t getHashByidx(const std::string &str_,const uint8_t &hash_idx_, HashSetUtil &hsu_){
    return hsu_.GetHash(hashset_[hash_idx_ - 1], str_);
}

bool compare_func(HABFilterKey * key1, HABFilterKey * key2){
    return key1->data_->cost > key2->data_->cost;
}

void BloomFilter::CreateNewFilter(uint32_t &bit_size_){
    if (bit_size_ < 64) bit_size_ = 64;
    uint32_t bytes = (bit_size_ + 7) / 8;
    bit_size_ = bytes * 8;
    data_.resize(bytes, 0);
}

void BloomFilter::AddAll(std::vector<HABFilterKey> &keys){
    uint32_t bits = data_.size() * 8;
    char *array = &(data_[0]);
    for(int i=0; i<keys.size(); i++){
        for(uint8_t j=0; j<HASH_NUM; j++){          
            uint32_t h = getHashByidx(keys[i].data_->str, keys[i].hash_map_[j], hsu_);
            uint32_t bitpos = h % bits;
            array[bitpos / 8] |= (1 << (bitpos % 8));
        }
    }
}

bool BloomFilter::Contain(HABFilterKey &key){
    uint32_t bytes = data_.size();
    uint32_t bits = bytes * 8;
    for(int i=0 ; i<HASH_NUM; i++){
        uint32_t h = getHashByidx(key.data_->str, key.hash_map_[i], hsu_);
        uint32_t bitpos = h % bits;
        if ((data_[bitpos / 8] & (1 << (bitpos % 8))) == 0) return 0;
    }
    return 1;
}

void BloomFilter::initByV_(Tuple_V * V_){
    char *array = &(data_[0]);
    int bits = data_.size() * 8;
    for(int i=0; i<bits; i++){    
        if(V_[i].keyid != NULL){
            array[i / 8] |= (1 << (i % 8));
        }
    }
}

HashExpressor::HashExpressor(uint8_t k_, HashFamily entry_hash_) : k_(k_), entry_hash_(entry_hash_),
cell_size_(log(hashset_.size())/log(2) + 2), default_hashset_(std::vector<uint8_t>(k_)){
    for(int i=0; i<k_; i++) 
        default_hashset_[i] = i + 1;
}

void HashExpressor::CreateNewFilter(uint32_t bit_size_){
    cell_num_ = bit_size_ / cell_size_ + 1;
    data_.resize(cell_num_ * cell_size_ / 8 + 1, 0);
}

bool HashExpressor::Add(HABFilterKey &key){
    std::vector<uint8_t> insert_sequence_;

    if(GetInsertSequence(key, insert_sequence_))
        Add(key, insert_sequence_);
    else
        return false;
    return true;
}

void HashExpressor::Query(HABFilterKey &key){
    HashExpressorCell cell;

    //get first cell by entry_hash_
    uint32_t h =  hsu_.GetHash(entry_hash_, key.data_->str) % cell_num_; 

    int queryed_n_ = 0;
    for(int i=0; i<k_; i++){
        getCell(h, cell); 
        // cell is empty or enbit of last cell is false, return default hashset
        if(cell.hash_idx_ == 0 || (!cell.endbit && i == k_ - 1))
            break;
        key.hash_map_[queryed_n_++] = cell.hash_idx_;
        h = getHashByidx(key.data_->str, cell.hash_idx_, hsu_) % cell_num_; // mapped to next cell
    }
    if(queryed_n_ != k_)
        AssignDefaultHashset_(key);
}

void HashExpressor::Add(HABFilterKey &key, std::vector<uint8_t> &insert_sequence_){
    HashExpressorCell cell;
    uint32_t h =  hsu_.GetHash(entry_hash_, key.data_->str) % cell_num_; 

    for(int i=0; i<insert_sequence_.size()-1; i++){
        cell.endbit = 0;
        cell.hash_idx_ = insert_sequence_[i];
        setCell(h, cell);
        h = getHashByidx(key.data_->str, insert_sequence_[i], hsu_) % cell_num_;
    }
    cell.endbit = 1; 
    cell.hash_idx_ = insert_sequence_.back();
    setCell(h, cell);
}

bool HashExpressor::GetInsertSequence(HABFilterKey &key, std::vector<uint8_t> &insert_sequence_){

    clearHashSet(insert_sequence_);

    //get first cell by entry_hash_
    HashExpressorCell cell;
    uint32_t h =  hsu_.GetHash(entry_hash_, key.data_->str) % cell_num_;
    getCell(h, cell);

    std::set<uint32_t> hashpos_; // hash position for key
    hashpos_.insert(h);
    for(int i=0; i<k_; i++){
        h = getHashByidx(key.data_->str, key.hash_map_[i], hsu_) % cell_num_;
        hashpos_.insert(h);
    }
    if(hashpos_.size() < k_ + 1) return false; // existing duplicate hashing

    for(uint8_t i=0; i<k_; i++){
        if(cell.hash_idx_ == 0 || cell.hash_idx_ == key.hash_map_[i]){ // if can be inserted into first cell

            std::vector<uint8_t> hash_map_(key.hash_map_, key.hash_map_+k_);
            std::swap(hash_map_[i], hash_map_[0]); // i-th element to be inserted firstly 

            FindInsertSequence(key, hash_map_, insert_sequence_, 1, k_ - 1); // Find the remainder k_-1 insert sequence

            if(insert_sequence_.size() > 0) return true;
        }
    }
    return false;
}

bool HashExpressor::AddWithMaxPattarn(HABFilterKey &key){
    std::vector<uint8_t> insert_sequence_;

    if(GetInsertSequenceWithMaxPattarn(key, insert_sequence_))
        Add(key, insert_sequence_);
    else
        return false;
    return true;
}

bool HashExpressor::GetInsertSequenceWithMaxPattarn(HABFilterKey &key, std::vector<uint8_t> &insert_sequence_){
    assert(k_>=2);

    clearHashSet(insert_sequence_);

    //get first cell by entry_hash_
    HashExpressorCell cell;
    uint32_t h =  hsu_.GetHash(entry_hash_, key.data_->str) % cell_num_;
    getCell(h, cell);

    std::set<uint32_t> hashpos_; // hash position for key
    hashpos_.insert(h);
    for(int i=0; i<k_; i++){
        h = getHashByidx(key.data_->str, key.hash_map_[i], hsu_) % cell_num_;
        hashpos_.insert(h);
    }
    if(hashpos_.size() < k_ + 1) return false; // existing duplicate hashing

    uint8_t maxidx = std::distance(key.hash_map_, std::max_element(key.hash_map_, key.hash_map_+k_));
    std::swap(key.hash_map_[maxidx], key.hash_map_[k_-1]); // max element is always to be inserted lastly 
    for(uint8_t i=0; i<k_-1; i++){
        if(i == maxidx) continue;
        if(cell.hash_idx_==0 || cell.hash_idx_ == key.hash_map_[i]){ // if can be inserted into first cell
            std::vector<uint8_t> hash_map_(key.hash_map_, key.hash_map_ + k_);
            std::swap(hash_map_[i], hash_map_[0]); // i-th element to be inserted firstly 

            FindInsertSequence(key, hash_map_, insert_sequence_, 1, k_ - 2); // Find the remainder k_-2 insert sequence

            if(insert_sequence_.size() > 0){
                // judge whether the last cell can be inserted 
                HashExpressorCell c_;
                h = getHashByidx(key.data_->str, insert_sequence_[k_ - 2], hsu_) % cell_num_;
                getCell(h, c_);
                if(c_.hash_idx_==0 || c_.hash_idx_ == insert_sequence_.back()) return true;
            }
        }
    }
    return false;
}

void HashExpressor::QueryWithMaxPattarn(HABFilterKey &key){
    HashExpressorCell cell;

    //get first cell by entry_hash_
    uint32_t h =  hsu_.GetHash(entry_hash_, key.data_->str) % cell_num_; 

    uint8_t queryed_n_ = 0;
    for(int i=0; i<k_; i++){
        getCell(h, cell); 
        // cell is empty or enbit of last cell is false, return default hashset
        if(cell.hash_idx_ == 0 || (!cell.endbit && i == k_ - 1))
            break;
        key.hash_map_[queryed_n_++] = cell.hash_idx_;
        h = getHashByidx(key.data_->str, cell.hash_idx_, hsu_) % cell_num_; // mapped to next cell
    }

    if(queryed_n_ != k_)
        AssignDefaultHashset_(key);
    else{
        uint8_t maxidx = *(std::max_element(key.hash_map_, key.hash_map_ + k_));
        // max element need to be the last item
        if(key.hash_map_[k_-1] != maxidx)
            AssignDefaultHashset_(key);
    }
}

void HashExpressor::FindInsertSequence(HABFilterKey &key, std::vector<uint8_t> &hash_map_, std::vector<uint8_t> &insert_sequence_, int left_, int right_){
    if(insert_sequence_.size() > 0) return; // assign already
    if(left_ == right_){ // one possible insert sequence
        uint32_t h;
        HashExpressorCell cell;
        for(int i=0; i<right_; ++i){
            h = getHashByidx(key.data_->str, hash_map_[i], hsu_) % cell_num_;
            getCell(h, cell); //get the cell mapped by hash of i-th element

            if(cell.hash_idx_!=0 && cell.hash_idx_ != hash_map_[i+1]) return; 
            // insert only when cell is empty or its hash_idx_ equals the next element
        }
        insert_sequence_.assign(hash_map_.begin(), hash_map_.end());

    }else if(left_ < right_){
        for(int i=left_;i<=right_;++i){
            std::swap(hash_map_[i], hash_map_[left_]);
            FindInsertSequence(key, hash_map_, insert_sequence_, left_ + 1, right_);
            std::swap(hash_map_[i], hash_map_[left_]);
        }
    }else return;
}

void HashExpressor::setbit(const int &bitinx){
    assert(bitinx<=data_.size()*8);
    char *array = &(data_[0]);
    int charIdx = (bitinx) / 8; 
    int b = (bitinx) % 8;
    array[charIdx] |= (1 << b);
}

bool HashExpressor::getbit(const int &bitinx){
    assert(bitinx<=data_.size()*8);
    char *array = &(data_[0]);
    int charIdx = (bitinx) / 8; 
    int b = (bitinx) % 8;
    return (array[charIdx] & (1 << b)) > 0;
}

void HashExpressor::getCell(const int &cellinx, HashExpressorCell &cell){
    cell.endbit = getbit(cellinx * cell_size_);
    uint8_t w = 1;
    cell.hash_idx_ = 0;
    for(uint8_t i=cell_size_-1; i>0; i--){
        cell.hash_idx_ += w*getbit(cellinx*cell_size_ + i);
        w*=2; 
    }
}

void HashExpressor::setCell(const int &cellinx, const HashExpressorCell &cell){
    if(cell.endbit) setbit(cellinx * cell_size_);
    uint8_t i = cell_size_ - 1;
    uint8_t hash_index_ = cell.hash_idx_;
    while(hash_index_ > 0){
        assert(i>0);
        if(hash_index_ % 2 > 0) setbit(cellinx * cell_size_ + i);
        i--;
        hash_index_ = hash_index_ / 2;
    }
}

void HashExpressor::clearHashSet(std::vector<uint8_t> &hash_set){
    for(std::vector<uint8_t>::iterator it=hash_set.begin(); it!= hash_set.end();)
        hash_set.erase(it);
}

void HashExpressor::AssignDefaultHashset_(HABFilterKey *keys_, uint32_t &size_){
    for(int i=0; i<size_; i++)
        for(int j=0; j<k_; j++)
            keys_[i].hash_map_[j] = j + 1;
}

void HashExpressor::AssignDefaultHashset_(HABFilterKey &key){
    for(int j=0; j<k_; j++)
            key.hash_map_[j] = j + 1;
}

HABFilter::HABFilter(float bits_per_key_, int pos_count_) : k_(HASH_NUM), hash_expressor_(HashExpressor(HASH_NUM, HASHEXPRESSOR_ENTRY_HASH)){
    uint32_t bits = bits_per_key_ * pos_count_;
    uint32_t bloom_size_ = bits * (ALLOCATION_RATIO) / (ALLOCATION_RATIO + 1);
    uint32_t hashexpressor_size_ = bits / (ALLOCATION_RATIO + 1);
    bloom_.CreateNewFilter(bloom_size_);
    hash_expressor_.CreateNewFilter(hashexpressor_size_);
    std::cout << "HABF Size:" << (bloom_size_ + hashexpressor_size_) / (double)(8*1024*1024)  << std::endl;
}

void HABFilter::AddAndOptimize(const std::vector<Slice *> &pos_keys_, const std::vector<Slice *> &neg_keys_){

    // Keys with hash_idx_
    // std::vector<HABFilterKey> habf_pos_keys_(pos_keys_.size());
    // std::vector<HABFilterKey> habf_neg_keys_(neg_keys_.size());

    uint32_t pos_n = pos_keys_.size();
    uint32_t neg_n = neg_keys_.size();
    HABFilterKey* habf_pos_keys_ = new HABFilterKey[pos_n];
    HABFilterKey* habf_neg_keys_ = new HABFilterKey[neg_n];

    
    for(int i=0; i<pos_n; ++i)
        habf_pos_keys_[i].data_ = pos_keys_[i];
    for(int i=0; i<neg_n; ++i)
        habf_neg_keys_[i].data_ = neg_keys_[i];
    
    // Assign default hashset H_0 for each key 
    hash_expressor_.AssignDefaultHashset_(habf_pos_keys_, pos_n);
    hash_expressor_.AssignDefaultHashset_(habf_neg_keys_, neg_n);

    std::cout << "Assign all keys..." << std::endl;
    // Init V_
    // auto t1 = std::chrono::steady_clock::now();

    // define struct V_ and T_ (Gamma in paper)
    uint32_t bloom_size_ = bloom_.getSize();
    Tuple_V * V_ = new Tuple_V[bloom_size_];
    std::vector<std::vector<HABFilterKey *>> T_(bloom_size_);

    // hash position cache for init
    uint32_t (* pos_cache)[HASH_NUM] = new uint32_t[pos_keys_.size()][HASH_NUM]; 
    uint32_t (* neg_cache)[HASH_NUM] = new uint32_t[neg_keys_.size()][HASH_NUM]; 
    for(int i=0; i<pos_keys_.size(); i++){
        for(int j=0; j<k_; j++){
            pos_cache[i][j] = getHashByidx(habf_pos_keys_[i].data_->str, habf_pos_keys_[i].hash_map_[j], hsu_) % bloom_size_;
        }
    }
    for(int i=0; i<neg_keys_.size(); i++){
        for(int j=0; j<k_; j++){
            neg_cache[i][j] = getHashByidx(habf_neg_keys_[i].data_->str, habf_neg_keys_[i].hash_map_[j], hsu_) % bloom_size_;
        }
    }

    std::cout << "Init V_... " << std::flush;

    uint32_t h;
    for(int i=0; i<pos_n; i++){
        for(uint32_t &h : pos_cache[i]){
            if(V_[h].singleflag_){
                if(NULL == V_[h].keyid)
                    V_[h].keyid = &habf_pos_keys_[i];
                else
                    V_[h].singleflag_ = false;
            }
        }
    }

    std::cout << "\r                       \r" << std::flush;
    std::cout << "Init T_... " << std::flush;

    // Init T_ and collision queue CQ_
    std::vector<HABFilterKey *> CQ_;
    bool inBloom;
    for(int i=0; i<neg_n; i++){
        inBloom = true;
        // membership testing by V_
        for(uint32_t &h : neg_cache[i]){
            if(V_[h].keyid == NULL){
                inBloom = false;
                break;
            }
        }
        if(inBloom)
            CQ_.emplace_back(&habf_neg_keys_[i]); // misjudge and pushed into CQ
        else{
            for(uint32_t &h : neg_cache[i]){ //judge correctly and pushed into T_
                T_[h].emplace_back(&habf_neg_keys_[i]);
            }
        }
    }
    delete [] pos_cache;
    delete [] neg_cache;

    std::cout << "\r                       \r" << std::flush;
    std::cout << "Init CQ_... " << std::flush;

    // init CQ_ sorted by the cost of collision keys
    std::sort(CQ_.begin(), CQ_.end(), compare_func);

    std::cout << "\r                       \r" << std::flush;
    // std::cout << "init all struct" << std::endl;

    // std::cout << "CQ_ size: " << CQ_.size() << std::endl;
    
    // auto t2 = std::chrono::steady_clock::now();

    // double construct_time_ = std::chrono::duration<double>(t2 - t1).count();
    // std::cout<<"init_time_:"<< construct_time_ << std::endl;

    // Start Optimization
    Optimize(V_, T_, CQ_);

    // Insert optimized positive keys into bloom_
    bloom_.initByV_(V_);
    // bloom_.AddAll(habf_pos_keys_);

    delete [] habf_pos_keys_;
    delete [] habf_neg_keys_;
    delete [] V_;
}

void HABFilter::Optimize(Tuple_V * V_, std::vector<std::vector<HABFilterKey *>> &T_, std::vector<HABFilterKey *> &CQ_){
    
    HABFilterKey *ck; // current optimizing collision key

    uint32_t bloom_size_ = bloom_.getSize();

    int count_ = 0;
    while(true){
        if(count_ == CQ_.size()){
            break;
        }
        ck = CQ_[count_++];
        double min_cost_ = __DBL_MAX__;
        bool not_empty_flag = false;
        HABFilterKey * optimal_pk_; // optimal positive key to adjust hash func
        std::vector<uint8_t> optimal_insert_sequence_; // optimal hash idx insertion sequence of optimal_pk_ into hashexpressor
        std::vector<HABFilterKey *> optimal_cost_keys_; // negative keys in T_ which will caused FPR as the adjustment of optimal_pk_
        std::pair<uint32_t, uint32_t> pre_next_; // previous hash position of optimal_pk_ and next hash idx which optimal_pk_ will adjust to

        for(uint8_t &hashidx_ck_ : ck->hash_map_){
            if(not_empty_flag) break; // the minimum adjustment has been got already
            uint32_t h1 = getHashByidx(ck->data_->str, hashidx_ck_, hsu_) % bloom_size_;
            if(V_[h1].keyid == NULL) break; // ck has been optimized by front optimization
            if(!V_[h1].singleflag_) continue; // V_[h1] is mapped more than once positive key,

            HABFilterKey *pk = V_[h1].keyid;
            
            // hash func mapped to h is changed to the hash func which is not in hashset of fk with minimum cost
            double cost;
            for(uint8_t hashidx_ = 1; hashidx_ < hashset_.size(); hashidx_++){ // Traverse all candidate hash func
                if(std::find(pk->hash_map_, pk->hash_map_ + k_, hashidx_) == pk->hash_map_ + k_){ 
                    uint32_t h2 = getHashByidx(pk->data_->str, hashidx_, hsu_) % bloom_size_;

                    HABFilterKey * tmp_key_ = new HABFilterKey(); // temporary key to judge if new hashset of fk can be inserted into HashExpressor
                    std::vector<uint8_t> insert_sequence_;

                    tmp_key_->data_ = pk->data_;
                    uint8_t i = 0;
                    tmp_key_->hash_map_[i++] = hashidx_;
                    for(uint8_t &hashidx_pk_ : pk->hash_map_){
                        uint32_t h3 = getHashByidx(pk->data_->str, hashidx_pk_, hsu_) % bloom_size_;
                        if(h3 != h1) 
                            tmp_key_->hash_map_[i++] = hashidx_pk_;
                    }

                    if(hash_expressor_.GetInsertSequence(*tmp_key_, insert_sequence_)){
                        if(V_[h2].keyid != NULL){  // new hash is mapped into a empty cell
                            min_cost_ = 0;
                            optimal_pk_ = pk;
                            pre_next_.first = h1;
                            pre_next_.second = h2;
                            optimal_insert_sequence_.assign(insert_sequence_.begin(), insert_sequence_.end());
                            not_empty_flag = true;
                        }else{ // new hash is mapped into a cell which is not empty, so no cost will be caused
                            std::vector<HABFilterKey *> cost_keys_; 
                            DetectConflict(V_, T_, h2, cost_keys_, cost);
                            // if cost less than ck and get the minimum cost
                            if(cost < ck->data_->cost && cost < min_cost_){
                                min_cost_ = cost;
                                optimal_pk_ = pk;
                                pre_next_.first = h1;
                                pre_next_.second = h2;
                                optimal_insert_sequence_.assign(insert_sequence_.begin(), insert_sequence_.end());
                                optimal_cost_keys_.assign(cost_keys_.begin(), cost_keys_.end());
                            }
                        }
                    }
                    delete tmp_key_;
                    if(not_empty_flag) break; // the minimum adjustment has been got alread
                }
            }
        }
        if(optimal_pk_ != NULL && optimal_insert_sequence_.size() > 0){ // optimize successfully, update structure
            UpdateStructure(V_, T_, ck, optimal_pk_, optimal_insert_sequence_, pre_next_, optimal_cost_keys_);
        }
    }
}

// detect the bitpos_ -th bucket of T_
void HABFilter::DetectConflict(Tuple_V * V_, std::vector<std::vector<HABFilterKey *>> &T_, int bitpos_, 
std::vector<HABFilterKey *> &cost_keys_, double &cost){
    cost = 0.0;
    for(HABFilterKey * key : T_[bitpos_]){
        uint8_t col_count = 0;
        uint8_t bitpos_func_count = 0; // number of hash functions map into bitpos
        for(uint8_t &hash_idx_ : key->hash_map_){
            uint32_t h = getHashByidx(key->data_->str, hash_idx_, hsu_) % bloom_.getSize();
            if(h != bitpos_ && V_[h].keyid != NULL)
                col_count++;
            if(h == bitpos_)
                bitpos_func_count++;
        }
        //all hash func except hash functions map into bitpos conflict with positive keys
        if(col_count == k_ - bitpos_func_count){
            cost += key->data_->cost;
            cost_keys_.emplace_back(key);
        }
    }
}

void HABFilter::UpdateStructure(Tuple_V * V_, std::vector<std::vector<HABFilterKey *>> &T_, HABFilterKey *ck, HABFilterKey *optimal_pk_, 
std::vector<uint8_t> &optimal_insert_sequence_, std::pair<uint32_t, uint32_t> &pre_next_,
std::vector<HABFilterKey *> &optimal_cost_keys_){
    // push cost_keys_ into CQ_ and delete from T_ 
    for(HABFilterKey *cost_key_ : optimal_cost_keys_){
        for(uint8_t &hashidx_ : cost_key_->hash_map_){
            uint32_t h = getHashByidx(cost_key_->data_->str, hashidx_, hsu_) % bloom_.getSize();
            auto it = T_[h].begin();
            while(it != T_[h].end()){
                if(*it == cost_key_)
                    it = T_[h].erase(it);
                else
                    it++;
            }
        }
    }
    // push optimized ck into T_
    for(uint8_t &hashidx_ : ck->hash_map_){
        uint32_t h = getHashByidx(ck->data_->str, hashidx_, hsu_) % bloom_.getSize();
        T_[h].emplace_back(ck);
    }
    // insert changed hash functions of optimal_pk_ into hashexpressor with insert_sequence_
    hash_expressor_.Add(*optimal_pk_, optimal_insert_sequence_);

    // change hash of pk
    for(int j=0; j<k_; j++)
        optimal_pk_->hash_map_[j] = optimal_insert_sequence_[j];


    // update V_

    V_[pre_next_.first].keyid = NULL;
    if(V_[pre_next_.second].singleflag_){
        if(V_[pre_next_.second].keyid == NULL)
            V_[pre_next_.second].keyid = optimal_pk_;
        else
            V_[pre_next_.second].singleflag_ = false;
    }
}

bool HABFilter::Contain(Slice &key){
    HABFilterKey habf_key_;
    habf_key_.data_ = &key;
    hash_expressor_.AssignDefaultHashset_(habf_key_);
    if(!bloom_.Contain(habf_key_)){
        hash_expressor_.Query(habf_key_);
        return bloom_.Contain(habf_key_);
    }
    return true;
}


}
#endif