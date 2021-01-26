#ifndef FASTHABF
#define FASTHABF

#include <iostream>
#include <cmath>
#include <assert.h>
#include <algorithm>
#include <set>
#include <limits.h>
#include <exception> 
#include <string.h>
#include <chrono> 

#include "../util/key.h"
#include "../util/hashutil.h"
#define HASH_NUM 3
#define ALLOCATION_RATIO 4 // Allocation ratio between the size of bloom and hashexpressor

namespace fasthabf{

const int hashset_size_ = 15;

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
};

class HashExpressor{
    public:
        HashExpressor(uint8_t k_);
        void CreateNewFilter(uint32_t bit_size_);
        bool Add(HABFilterKey &key);
        void Add(HABFilterKey &key, std::vector<uint8_t> &insert_sequence_);
        bool GetInsertSequence(HABFilterKey &key, std::vector<uint8_t> &insert_sequence_);
        void Query(HABFilterKey &key);
        
        int getSize(){return 8 * data_.size();};
        
        void AssignDefaultHashset_(HABFilterKey *keys_, uint32_t &size_);
        void AssignDefaultHashset_(HABFilterKey &key);

    private:
        void setbit(const int &bitinx);
        bool getbit(const int &bitinx);
        void getCell(const int &cellinx, HashExpressorCell &cell);
        void setCell(const int &cellinx, const HashExpressorCell &cell);

        void FindInsertSequence(uint64_t &hash_key_, std::vector<uint8_t> &hash_map_, std::vector<uint8_t> &insert_sequence_, int left_, int right_);
        void clearHashSet(std::vector<uint8_t> &hash_set);


        uint32_t cell_num_;
        uint8_t cell_size_;
        uint8_t k_;
        std::string data_;
        std::vector<uint8_t> default_hashset_;


};

class FastHABFilter{
    public:
        FastHABFilter(float bits_per_key_, int pos_count_); // In order to fine-tune the size of the filter space for comparing, we allow bits_per_key_ to be float
        void AddAndOptimize(const std::vector<Slice *> &pos_keys_, const std::vector<Slice *> &neg_keys_);
        bool Contain(Slice &key);

    private:
        void Optimize(Tuple_V * V_, std::vector<HABFilterKey *> &CQ_);

        void UpdateStructure(Tuple_V * V_,  HABFilterKey *ck,  HABFilterKey *optimal_pk_, 
        std::vector<uint8_t> &optimal_insert_sequence_, std::pair<uint32_t, uint32_t> &pre_next_);

        uint8_t k_;
        BloomFilter bloom_;
        HashExpressor hash_expressor_;     
};

uint64_t rotl64(uint64_t n, unsigned int c) {
    // assumes width is a power of 2
    const unsigned int mask = (CHAR_BIT * sizeof(n) - 1);
    // assert ( (c<=mask) &&"rotate by type width or more");
    c &= mask;
    return (n << c) | ( n >> ((-c) & mask));
}

uint32_t reduce(uint32_t hash, uint32_t n) {
    // http://lemire.me/blog/2016/06/27/a-fast-alternative-to-the-modulo-reduction/
    return (uint32_t) (((uint64_t) hash * n) >> 32);
}

uint32_t GetHashFromHash(uint64_t hash, uint32_t index, uint32_t blockLength) {
    uint32_t r = rotl64(hash, index * 21);
    return (uint32_t) reduce(r, blockLength);
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
        uint64_t hash_key_ = XXH3_128bits(&keys[i].data_->str[0], keys[i].data_->str.size()).low64; 
        for(uint8_t j=0; j<HASH_NUM; j++){         
            uint32_t bitpos = GetHashFromHash(hash_key_, keys[i].hash_map_[j], bits);
            array[bitpos / 8] |= (1 << (bitpos % 8));
        }
    }
}

bool BloomFilter::Contain(HABFilterKey &key){
    uint32_t bytes = data_.size();
    uint32_t bits = bytes * 8;
    uint64_t hash_key_ = XXH3_128bits(&key.data_->str[0], key.data_->str.size()).low64; 
    for(int i=0 ; i<HASH_NUM; i++){
        uint32_t bitpos = GetHashFromHash(hash_key_, key.hash_map_[i], bits);
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

HashExpressor::HashExpressor(uint8_t k_) : k_(k_), cell_size_(log(hashset_size_)/log(2) + 2), default_hashset_(std::vector<uint8_t>(k_)){
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
    uint64_t hash_key_ = XXH3_128bits(&key.data_->str[0], key.data_->str.size()).low64;
    uint32_t h =  GetHashFromHash(hash_key_, 0, cell_num_);

    int queryed_n_ = 0;
    for(int i=0; i<k_; i++){
        getCell(h, cell); 
        // cell is empty or enbit of last cell is false, return default hashset
        if(cell.hash_idx_ == 0 || (!cell.endbit && i == k_ - 1))
            break;
        key.hash_map_[queryed_n_++] = cell.hash_idx_;
        h = GetHashFromHash(hash_key_, cell.hash_idx_, cell_num_); // mapped to next cell
    }
    if(queryed_n_ != k_)
        AssignDefaultHashset_(key);
}

void HashExpressor::Add(HABFilterKey &key, std::vector<uint8_t> &insert_sequence_){
    HashExpressorCell cell;
    uint64_t hash_key_ = XXH3_128bits(&key.data_->str[0], key.data_->str.size()).low64;
    uint32_t h =  GetHashFromHash(hash_key_, 0, cell_num_);

    for(int i=0; i<insert_sequence_.size()-1; i++){
        cell.endbit = 0;
        cell.hash_idx_ = insert_sequence_[i];
        setCell(h, cell);
        h = GetHashFromHash(hash_key_, insert_sequence_[i], cell_num_);
    }
    cell.endbit = 1; 
    cell.hash_idx_ = insert_sequence_.back();
    setCell(h, cell);
}

bool HashExpressor::GetInsertSequence(HABFilterKey &key, std::vector<uint8_t> &insert_sequence_){

    clearHashSet(insert_sequence_);

    //get first cell by entry_hash_
    HashExpressorCell cell;
    uint64_t hash_key_ = XXH3_128bits(&key.data_->str[0], key.data_->str.size()).low64;
    uint32_t h =  GetHashFromHash(hash_key_, 0, cell_num_);
    getCell(h, cell);

    std::set<uint32_t> hashpos_; // hash position for key
    hashpos_.insert(h);
    for(int i=0; i<k_; i++){
        h = GetHashFromHash(hash_key_, key.hash_map_[i], cell_num_); 
        hashpos_.insert(h);
    }
    if(hashpos_.size() < k_ + 1) return false; // existing duplicate hashing

    for(uint8_t i=0; i<k_; i++){
        if(cell.hash_idx_ == 0 || cell.hash_idx_ == key.hash_map_[i]){ // if can be inserted into first cell

            std::vector<uint8_t> hash_map_(key.hash_map_, key.hash_map_+k_);
            std::swap(hash_map_[i], hash_map_[0]); // i-th element to be inserted firstly 

            FindInsertSequence(hash_key_, hash_map_, insert_sequence_, 1, k_ - 1); // Find the remainder k_-1 insert sequence

            if(insert_sequence_.size() > 0) return true;
        }
    }
    return false;
}


void HashExpressor::FindInsertSequence(uint64_t &hash_key_, std::vector<uint8_t> &hash_map_, std::vector<uint8_t> &insert_sequence_, int left_, int right_){
    if(insert_sequence_.size() > 0) return; // assign already
    if(left_ == right_){ // one possible insert sequence
        uint32_t h;
        HashExpressorCell cell;
        for(int i=0; i<right_; ++i){
            h = GetHashFromHash(hash_key_, hash_map_[i], cell_num_);
            getCell(h, cell); //get the cell mapped by hash of i-th element

            if(cell.hash_idx_!=0 && cell.hash_idx_ != hash_map_[i+1]) return; 
            // insert only when cell is empty or its hash_idx_ equals the next element
        }
        insert_sequence_.assign(hash_map_.begin(), hash_map_.end());

    }else if(left_ < right_){
        for(int i=left_;i<=right_;++i){
            std::swap(hash_map_[i], hash_map_[left_]);
            FindInsertSequence(hash_key_, hash_map_, insert_sequence_, left_ + 1, right_);
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


FastHABFilter::FastHABFilter(float bits_per_key_, int pos_count_) : k_(HASH_NUM), hash_expressor_(HashExpressor(HASH_NUM)){
    uint32_t bits = bits_per_key_ * pos_count_;
    uint32_t bloom_size_ = bits * (ALLOCATION_RATIO) / (ALLOCATION_RATIO + 1);
    uint32_t hashexpressor_size_ = bits / (ALLOCATION_RATIO + 1);
    bloom_.CreateNewFilter(bloom_size_);
    hash_expressor_.CreateNewFilter(hashexpressor_size_);
    std::cout << "FastHABF Size:" << (bloom_size_ + hashexpressor_size_) / (double)(8*1024*1024)  << std::endl;
}

void FastHABFilter::AddAndOptimize(const std::vector<Slice *> &pos_keys_, const std::vector<Slice *> &neg_keys_){

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

    uint64_t * pos_cache = new uint64_t[pos_n]; 
    uint64_t * neg_cache = new uint64_t[neg_n]; 

    for(int i=0; i<pos_n; i++)
        pos_cache[i] = XXH3_128bits(&habf_pos_keys_[i].data_->str[0], habf_pos_keys_[i].data_->str.size()).low64;
    for(int i=0; i<neg_n; i++)
        neg_cache[i] = XXH3_128bits(&habf_neg_keys_[i].data_->str[0], habf_neg_keys_[i].data_->str.size()).low64;

        
    // std::cout << "Assign all keys..." << std::endl;
    // Init V_
    auto t1 = std::chrono::steady_clock::now();

    // define struct V_ 
    uint32_t bloom_size_ = bloom_.getSize();
    Tuple_V * V_ = new Tuple_V[bloom_size_];

    std::cout << "Init V_... " << std::flush;

    uint32_t h;
    for(int i=0; i<pos_n; i++){
        for(int j=1; j<=k_; j++){
            h = GetHashFromHash(pos_cache[i], j, bloom_size_);
            if(V_[h].singleflag_){
                if(NULL == V_[h].keyid)
                    V_[h].keyid = &habf_pos_keys_[i];
                else
                    V_[h].singleflag_ = false;
            }
        }
    }

    std::cout << "\r                       \r" << std::flush;
    std::cout << "Init CQ_... " << std::flush;

    // Init collision queue CQ_
    std::vector<HABFilterKey *> CQ_;
    bool inBloom;
    for(int i=0; i<neg_n; i++){
        inBloom = true;
        // membership testing by V_
        for(int j=1; j<=k_; j++){
            h = GetHashFromHash(neg_cache[i], j, bloom_size_);
            if(V_[h].keyid == NULL){
                inBloom = false;
                break;
            }
        }
        if(inBloom)
            CQ_.emplace_back(&habf_neg_keys_[i]); // misjudge and pushed into CQ
    }

    delete [] pos_cache;
    delete [] neg_cache;

    // init CQ_ sorted by the cost of collision keys
    std::sort(CQ_.begin(), CQ_.end(), compare_func);

    // for(int i=0;i<CQ_.size();i++){
    //     std::cout << CQ_[i]->data_->cost << ",";
    // }

    std::cout << "\r                       \r" << std::flush;
    // std::cout << "init all struct" << std::flush;

    // std::cout << "CQ_ size: " << CQ_.size() << std::endl;
    
    auto t2 = std::chrono::steady_clock::now();

    double construct_time_ = std::chrono::duration<double>(t2 - t1).count();
    // std::cout<<"init_time_:"<< construct_time_ << std::endl;

    // Start Optimization
    Optimize(V_, CQ_);

    // Insert optimized positive keys into bloom_
    bloom_.initByV_(V_);
    // bloom_.AddAll(habf_pos_keys_);

    delete [] habf_pos_keys_;
    delete [] habf_neg_keys_;
    delete [] V_;
}

void FastHABFilter::Optimize(Tuple_V * V_, std::vector<HABFilterKey *> &CQ_){
    
    HABFilterKey *ck; // current optimizing collision key

    uint32_t bloom_size_ = bloom_.getSize();

    int count_ = 0;
    while(true){
        if(count_ == CQ_.size()){
            break;
        }
        ck = CQ_[count_++];
        bool min_cost_ = false;
        HABFilterKey * optimal_pk_; // optimal positive key to adjust hash func
        std::vector<uint8_t> optimal_insert_sequence_; // optimal hash idx insertion sequence of optimal_pk_ into hashexpressor
        std::pair<uint32_t, uint32_t> pre_next_; // previous hash position of optimal_pk_ and next hash idx which optimal_pk_ will adjust to

        uint64_t hash_ck_ = XXH3_128bits(&ck->data_->str[0], ck->data_->str.size()).low64;
        for(uint8_t &hashidx_ck_ : ck->hash_map_){
            if(min_cost_) break; // the minimum adjustment has been got already

            uint32_t h1 = GetHashFromHash(hash_ck_, hashidx_ck_, bloom_size_);
            if(V_[h1].keyid == NULL) break; // ck has been optimized by front optimization
            if(!V_[h1].singleflag_) continue; // V_[h1] is mapped more than once positive key,
            // std::cout << (V_[h1].keyid == NULL) << std::endl;
            HABFilterKey *pk = V_[h1].keyid;
            uint64_t hash_pk_ = XXH3_128bits(&pk->data_->str[0], pk->data_->str.size()).low64;
            // hash func mapped to h is changed to the hash func which is not in hashset of fk with minimum cost
            for(uint8_t hashidx_ = 1; hashidx_ < hashset_size_; hashidx_++){ // Traverse all candidate hash func
                if(std::find(pk->hash_map_, pk->hash_map_ + k_, hashidx_) == pk->hash_map_ + k_){ 
                    uint32_t h2 = GetHashFromHash(hash_pk_, hashidx_, bloom_size_);
                    if(V_[h2].keyid == NULL) continue; // new hash is mapped into a cell which is empty, not consider
                    HABFilterKey * tmp_key_ = new HABFilterKey(); // temporary key to judge if new hashset of fk can be inserted into HashExpressor
                    std::vector<uint8_t> insert_sequence_;
                    tmp_key_->data_ = pk->data_;
                    int i = 0;
                    tmp_key_->hash_map_[i++] = hashidx_;
                    for(uint8_t &hashidx_pk_ : pk->hash_map_){
                        uint32_t h3 = GetHashFromHash(hash_pk_, hashidx_pk_, bloom_size_);
                        if(h3 != h1) 
                            tmp_key_->hash_map_[i++] = hashidx_pk_;
                    }
                    if(hash_expressor_.GetInsertSequence(*tmp_key_, insert_sequence_)){
                        min_cost_ = true;
                        optimal_pk_ = pk;
                        pre_next_.first = h1;
                        pre_next_.second = h2;
                        optimal_insert_sequence_.assign(insert_sequence_.begin(), insert_sequence_.end());
                    }
                    delete tmp_key_;
                    if(min_cost_ ) break; // the minimum adjustment has been got alread
                }
            }
        }
        if(optimal_pk_ != NULL && optimal_insert_sequence_.size() > 0){ // optimize successfully, update structure
            UpdateStructure(V_, ck, optimal_pk_, optimal_insert_sequence_, pre_next_);
        }
    }
}


void FastHABFilter::UpdateStructure(Tuple_V * V_, HABFilterKey *ck, HABFilterKey *optimal_pk_, 
std::vector<uint8_t> &optimal_insert_sequence_, std::pair<uint32_t, uint32_t> &pre_next_){

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

bool FastHABFilter::Contain(Slice &key){
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