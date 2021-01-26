#ifndef DYNAMIXHABF
#define DYNAMIXHABF

#include <iostream>
#include <cmath>
#include <assert.h>
#include <algorithm>
#include <set>
#include <unordered_map>
#include <limits.h>
#include <exception> 
#include <string.h>
#include <chrono> 
#include <fstream>

#include "../util/key.h"
#include "../util/hashutil.h"
#define HASH_NUM 3
#define ALLOCATION_RATIO 4 // Allocation ratio between the size of bloom and hashexpressor

namespace dynamichabf{

const int hashset_size_ = 7;

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
        void Add(HABFilterKey &key);
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
        bool GetInsertSequence(HABFilterKey &key, std::vector<uint8_t> &insert_sequence_, int &min_cell_cout_);
        void Query(HABFilterKey &key);
        
        int getSize(){return 8 * data_.size();};
        int getCellNum(){return cell_num_;};


        void AssignDefaultHashset_(HABFilterKey *keys_, uint32_t &size_);
        void AssignDefaultHashset_(HABFilterKey &key);

    private:
        void setbit(const int &bitinx);
        bool getbit(const int &bitinx);
        void getCell(const int &cellinx, HashExpressorCell &cell);
        void setCell(const int &cellinx, const HashExpressorCell &cell);

        void FindInsertSequence(uint64_t &hash_key_, std::vector<uint8_t> &hash_map_, std::vector<uint8_t> &insert_sequence_, int left_, int right_, int init_cell_count_, int &min_cell_count_);
        void clearHashSet(std::vector<uint8_t> &hash_set);


        uint32_t cell_num_;
        uint8_t cell_size_;
        uint8_t k_;
        uint32_t t;
        std::string data_;
        std::vector<uint8_t> default_hashset_;


};

class DynamicHABFilter{
    public:
        DynamicHABFilter(float bits_per_key_, int neg_count_, std::vector<Slice *> &neg_keys_); // In order to fine-tune the size of the filter space for comparing, we allow bits_per_key_ to be float
        void Add(Slice &key);
        bool Contain(Slice &key);
        bool ContainTest(Slice &key);
        int m=0;
        int getHashExpressorCellNum(){return hash_expressor_.getCellNum();}
        void writeToFile();

    private:
        void combine_inner(std::vector<uint8_t> &data, int start, int n, int m, int depth, std::vector<uint8_t> temp,std::vector<std::vector<uint8_t>> &result);
        uint8_t k_;
        BloomFilter bloom_;
        HashExpressor hash_expressor_;
        double total_cost = 0.0;
        std::vector<double> bit_positions; 
        std::vector<double> bit_positions_tmp;
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

// for scene of few good hash functions or caculate quickly 
uint32_t GetHashFromHash(uint64_t hash, uint32_t index, uint32_t blockLength) {
    uint32_t r = rotl64(hash, index * 21);
    return (uint32_t) reduce(r, blockLength);
}

bool compare_func(Slice * key1, Slice * key2){
    return key1->cost > key2->cost;
}

void BloomFilter::CreateNewFilter(uint32_t &bit_size_){
    if (bit_size_ < 64) bit_size_ = 64;
    uint32_t bytes = (bit_size_ + 7) / 8;
    bit_size_ = bytes * 8;
    data_.resize(bytes, 0);
}


void BloomFilter::Add(HABFilterKey &key){
    uint32_t bits = data_.size() * 8;
    char *array = &(data_[0]);
    for(uint8_t j=0; j<HASH_NUM; j++){         
        uint32_t bitpos = XXH3_128bits_withSeed(&key.data_->str[0], key.data_->str.size(), key.hash_map_[j]).low64 % bits;
        array[bitpos / 8] |= (1 << (bitpos % 8));
    }
}

bool BloomFilter::Contain(HABFilterKey &key){
    uint32_t bytes = data_.size();
    uint32_t bits = bytes * 8;
    for(int i=0 ; i<HASH_NUM; i++){
        uint32_t bitpos = XXH3_128bits_withSeed(&key.data_->str[0], key.data_->str.size(), key.hash_map_[i]).low64 % bits;
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
    t = 0;
    for(int i=0; i<k_; i++) 
        default_hashset_[i] = i + 1;
}

void HashExpressor::CreateNewFilter(uint32_t bit_size_){
    cell_num_ = bit_size_ / cell_size_ + 1;
    // std::cout<< "cell_num" << cell_num_ << std::endl;
    data_.resize(cell_num_ * cell_size_ / 8 + 1, 0);
}

bool HashExpressor::Add(HABFilterKey &key){
    std::vector<uint8_t> insert_sequence_;
    int min_cell_count_;
    if(GetInsertSequence(key, insert_sequence_, min_cell_count_)){
        Add(key, insert_sequence_);
    }
    else
        return false;
    return true;
}

void HashExpressor::Query(HABFilterKey &key){
    HashExpressorCell cell;

    //get first cell by entry_hash_
    uint64_t hash_key_ = XXH3_64bits(&key.data_->str[0], key.data_->str.size());
    uint32_t h =  GetHashFromHash(hash_key_, 0, cell_num_);

    int queryed_n_ = 0;
    for(int i=0; i<k_; i++){
        getCell(h, cell); 
        // cell is empty or enbit of last cell is false, return default hashset
        if(cell.hash_idx_ == 0 || (!cell.endbit && i == k_ - 1))
            break;
        if(std::find(key.hash_map_, key.hash_map_+queryed_n_, cell.hash_idx_) != key.hash_map_+queryed_n_)
            break;
        key.hash_map_[queryed_n_++] = cell.hash_idx_;
        h = GetHashFromHash(hash_key_, cell.hash_idx_, cell_num_); // mapped to next cell
    }
    if(queryed_n_ != k_)
        AssignDefaultHashset_(key);
}

void HashExpressor::Add(HABFilterKey &key, std::vector<uint8_t> &insert_sequence_){
    HashExpressorCell cell;
    uint64_t hash_key_ = XXH3_64bits(&key.data_->str[0], key.data_->str.size());
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

bool HashExpressor::GetInsertSequence(HABFilterKey &key, std::vector<uint8_t> &insert_sequence_, int &min_cell_cout_){

    clearHashSet(insert_sequence_);

    //get first cell by entry_hash_
    HashExpressorCell cell;
    uint64_t hash_key_ = XXH3_64bits(&key.data_->str[0], key.data_->str.size());
    uint32_t h =  GetHashFromHash(hash_key_, 0, cell_num_);
    getCell(h, cell);

    std::set<uint32_t> hashpos_; // hash position for key
    hashpos_.insert(h);
    for(int i=0; i<k_; i++){
        h = GetHashFromHash(hash_key_, key.hash_map_[i], cell_num_); 
        hashpos_.insert(h);
    }
    if(hashpos_.size() < k_ + 1) return false; // existing duplicate hashing

    min_cell_cout_ = HASH_NUM+1;
    for(uint8_t i=0; i<k_; i++){
        if(cell.hash_idx_ == key.hash_map_[i]){ // if can be inserted into first cell
            std::vector<uint8_t> hash_map_(key.hash_map_, key.hash_map_+k_);
            std::swap(hash_map_[i], hash_map_[0]); // i-th element to be inserted firstly 

            FindInsertSequence(hash_key_, hash_map_, insert_sequence_, 1, k_ - 1, 0,  min_cell_cout_); // Find the remainder k_-1 insert sequence
            if(insert_sequence_.size() > 0) return true;
        }
    }
    for(uint8_t i=0; i<k_; i++){
        if(cell.hash_idx_ == 0){ // if can be inserted into first cell
            std::vector<uint8_t> hash_map_(key.hash_map_, key.hash_map_+k_);
            std::swap(hash_map_[i], hash_map_[0]); // i-th element to be inserted firstly 

            FindInsertSequence(hash_key_, hash_map_, insert_sequence_, 1, k_ - 1, 1, min_cell_cout_); // Find the remainder k_-1 insert sequence
            if(insert_sequence_.size() > 0) return true;
        }
    }
    return false;
}


void HashExpressor::FindInsertSequence(uint64_t &hash_key_, std::vector<uint8_t> &hash_map_, std::vector<uint8_t> &insert_sequence_, int left_, int right_, int init_cell_count_, int &min_cell_count_){
    if(insert_sequence_.size() > 0 && min_cell_count_ == init_cell_count_) return; // assign already
    if(left_ == right_){ // one possible insert sequence
        uint32_t h;
        HashExpressorCell cell;
        int cell_count = init_cell_count_;
        for(int i=0; i<right_; ++i){
            h = GetHashFromHash(hash_key_, hash_map_[i], cell_num_);
            getCell(h, cell); //get the cell mapped by hash of i-th element
            if(cell.hash_idx_!=0 && cell.hash_idx_ != hash_map_[i+1]){
                return; 
            }
            if(cell.hash_idx_ == 0){
                cell_count += 1;
            }
        }
        if(cell_count < min_cell_count_){
            min_cell_count_ = cell_count;
            insert_sequence_.assign(hash_map_.begin(), hash_map_.end());
        }
        

    }else if(left_ < right_){
        for(int i=left_;i<=right_;++i){
            std::swap(hash_map_[i], hash_map_[left_]);
            FindInsertSequence(hash_key_, hash_map_, insert_sequence_, left_ + 1, right_, init_cell_count_, min_cell_count_);
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
    HashExpressorCell tmp_cell;
    getCell(cellinx, tmp_cell);
    if(tmp_cell.hash_idx_ == 0){
        t++;
        // std::cout << t << std::endl;
    } 

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


DynamicHABFilter::DynamicHABFilter(float bits_per_key_, int neg_count_, std::vector<Slice *> &neg_keys_) : k_(HASH_NUM), hash_expressor_(HashExpressor(HASH_NUM)){
    uint32_t bloom_size_ = 8 * neg_count_;

    bloom_.CreateNewFilter(bloom_size_);

    hash_expressor_.CreateNewFilter(2 * neg_count_);

    bit_positions.resize(bloom_size_, 0.0);
    bit_positions_tmp.resize(bloom_size_, 0.0);
    for(int i=0; i<neg_keys_.size(); i++){
        for(int j=1; j<=k_; j++){
            uint32_t h = XXH3_128bits_withSeed(&neg_keys_[i]->str[0], neg_keys_[i]->str.size(), j).low64 % bloom_size_;
            bit_positions[h]+=neg_keys_[i]->cost;
            bit_positions_tmp[h]+=neg_keys_[i]->cost;
            total_cost += neg_keys_[i]->cost;
        }
    }

}

void DynamicHABFilter::combine_inner(std::vector<uint8_t> &data, int start, int n, int m, int depth, std::vector<uint8_t> temp,std::vector<std::vector<uint8_t>> &result)
{
    if (depth == m - 1)
    {
        for (int i = start; i < n - (m - depth - 1); ++i)
        {
            temp[depth] = data[i];
            result.push_back(temp);
        }
    }
    else
        for (int i = start; i < n - (m - depth - 1);++i)
    {
        temp[depth] = data[i];
        combine_inner(data,i + 1, n, m, depth+1,temp,result);
    }
}


void DynamicHABFilter::Add(Slice &key){
    HABFilterKey * pk  = new HABFilterKey;
    pk->data_ = &key;
    std::vector<uint8_t> candate_hash(hashset_size_);

    //is default hash functions conflicting
    double cost_0 = 0.0;
    double max_cell_cost = 0.0;
    int new_cell_num = 0;
    std::vector<uint8_t> optimal_insert_sequence_;
    for(int j=1; j<=k_; j++){
        uint32_t h = XXH3_128bits_withSeed(&key.str[0], key.str.size(), j).low64 % bloom_.getSize();
        if(bit_positions[h] != __DBL_MAX__){
            cost_0+=bit_positions[h];
        }
    } 

    for(int hashidx_ = 1; hashidx_ <= hashset_size_; hashidx_++){
        candate_hash[hashidx_] = hashidx_;
    }
    std::vector<std::vector<uint8_t>> vec_res;
    std::vector<uint8_t> temp(k_,0);
    combine_inner(candate_hash,0, candate_hash.size(), k_, 0, temp, vec_res);
    for(std::vector<uint8_t> &v:vec_res){
        for(int i=0; i<v.size(); i++){
            pk->hash_map_[i] = v[i];
        }
        std::vector<uint8_t> insert_sequence;
        int min_cell_cout_;
        hash_expressor_.GetInsertSequence(*pk, insert_sequence, min_cell_cout_);
        if(insert_sequence.size()>0){
            double cost_1 = 0.0;
            for(int j=0;j<HASH_NUM;j++){
                uint32_t h = XXH3_128bits_withSeed(&key.str[0], key.str.size(), insert_sequence[j]).low64 % bloom_.getSize();
                if(bit_positions[h] != __DBL_MAX__){
                    cost_1+=bit_positions[h];
                }
            }
            if(min_cell_cout_ == 0){
                max_cell_cost = __DBL_MAX__;
                optimal_insert_sequence_.assign(insert_sequence.begin(),insert_sequence.end());
                new_cell_num = 0;
                break;
            }else{
                double cell_cost = (cost_0-cost_1) / (double)min_cell_cout_;
                if(cell_cost>max_cell_cost){
                    optimal_insert_sequence_.assign(insert_sequence.begin(),insert_sequence.end());
                    max_cell_cost = cell_cost;
                    new_cell_num = min_cell_cout_;
                }
            }
        }
    }
    if(max_cell_cost > (total_cost / ((hash_expressor_.getCellNum()/1.3)))){
        hash_expressor_.Add(*pk, optimal_insert_sequence_);
        hash_expressor_.Query(*pk);
        m+=new_cell_num;
    }else {
        hash_expressor_.AssignDefaultHashset_(*pk); 
    }
    for(int i=0; i<HASH_NUM; i++){
        uint32_t h = XXH3_128bits_withSeed(&key.str[0], key.str.size(), pk->hash_map_[i]).low64 % bloom_.getSize();
        bit_positions[h] = __DBL_MAX__;
    }
    bloom_.Add(*pk);
    delete pk;
}

void DynamicHABFilter::writeToFile(){
    std::ofstream ofs("bitposition.txt");
    for(int i=0; i<bit_positions.size(); i++){
        ofs << bit_positions_tmp[i] << ",";
        if(bit_positions[i] == __DBL_MAX__){
            ofs << "0" << ",";
        }else{
            ofs << "1" << ",";
        }
        ofs << std::endl;
    }
    ofs.close();
}

bool DynamicHABFilter::Contain(Slice &key){
    HABFilterKey habf_key_;
    habf_key_.data_ = &key;
    hash_expressor_.AssignDefaultHashset_(habf_key_);
    if(bloom_.Contain(habf_key_)) return true;
    hash_expressor_.Query(habf_key_);
    if(bloom_.Contain(habf_key_)) return true;
    return false;
}

bool DynamicHABFilter::ContainTest(Slice &key){
    HABFilterKey habf_key_;
    habf_key_.data_ = &key;
    hash_expressor_.AssignDefaultHashset_(habf_key_);
    if(bloom_.Contain(habf_key_)) return true;
    hash_expressor_.Query(habf_key_);
    for(int i=0; i<k_; i++){
        std::cout << (int)habf_key_.hash_map_[i] << ",";
    }
    std::cout << std::endl;
    if(bloom_.Contain(habf_key_)) return true;
    return false;
}

}
#endif