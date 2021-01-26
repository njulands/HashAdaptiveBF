#ifndef BLOOM
#define BLOOM

#include <algorithm>
#include <math.h>
#include "../util/key.h"
#include "../util/hashutil.h"

#define WBF_COST_RADIO 0.001

namespace wbf{

struct WBFKey{
    Slice * data_;
    uint8_t k_;
};

double getEr(double cost_e, double xe, double s){
    return (1-xe)*cost_e / (xe*cost_e + s);
}

bool compare_func(WBFKey &key1, WBFKey &key2){
    return key1.data_->cost > key2.data_->cost;
}

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

class BloomFilter{
    public:
        void CreateNewFilter(uint32_t &bit_size_);
        void AddAll(std::vector<Slice *> &keys, uint32_t k_);
        bool Contain(Slice &key, uint32_t k_);

        uint32_t getSize(){return 8 * data_.size();};

    private:
        std::string data_;
};

class WBF{
    public:
        WBF(float bits_per_key_) : bits_per_key_(bits_per_key_){};
        void AddAll(std::vector<Slice *> &pos_keys_, std::vector<Slice *> &neg_keys_); 
        bool Contain(Slice &key);
        double testFPR();
    private:
        float bits_per_key_;
        BloomFilter bloom_;
        std::vector<Slice *> cost_list_;
        std::vector<Slice *> no_cost_list_;
        double s;
        double s2;
        uint32_t ke;
        uint32_t pos_n;
        uint32_t N;
};

void BloomFilter::CreateNewFilter(uint32_t &bit_size_){
    if (bit_size_ < 64) bit_size_ = 64;
    uint32_t bytes = (bit_size_ + 7) / 8;
    bit_size_ = bytes * 8;
    data_.resize(bytes, 0);
}

void BloomFilter::AddAll(std::vector<Slice *> &keys, uint32_t k_){
    uint32_t bits = data_.size() * 8;
    char *array = &(data_[0]);
    for(int i=0; i<keys.size(); i++){
        uint64_t hash = XXH3_64bits(&keys[i]->str[0], keys[i]->str.size()); 
        uint64_t a = (hash >> 32) | (hash << 32);
        uint64_t b = hash;
        for(uint8_t j=0; j<k_; j++){         
            uint32_t bitpos = a % bits;
            array[bitpos / 8] |= (1 << (bitpos % 8));
            a+=b;
        }
    }
}

bool BloomFilter::Contain(Slice &key, uint32_t k_){
    uint32_t bytes = data_.size();
    uint32_t bits = bytes * 8;
    uint64_t hash = XXH3_64bits(&key.str[0], key.str.size()); 
    uint64_t a = (hash >> 32) | (hash << 32);
    uint64_t b = hash;
    for(int i=0 ; i<k_; i++){
        uint32_t bitpos = a % bits;
        if ((data_[bitpos / 8] & (1 << (bitpos % 8))) == 0) return 0;
        a+=b;
    }
    return 1;
}

void WBF::AddAll(std::vector<Slice *> &pos_keys_, std::vector<Slice *> &neg_keys_){
    pos_n = pos_keys_.size();
    uint32_t neg_n = neg_keys_.size();
    WBFKey* wbf_pos_keys_ = new WBFKey[pos_n];
    WBFKey* wbf_neg_keys_ = new WBFKey[neg_n];
    
    uint32_t bits = pos_n * bits_per_key_;
    uint32_t cost_storing_num = neg_n*WBF_COST_RADIO;
    uint32_t cost_bits = cost_storing_num * 8 * sizeof(Slice); 
    uint32_t bloom_bits = bits - cost_bits; 
    
    bloom_.CreateNewFilter(bloom_bits);

    std::cout << "bloom size:" << (double)bloom_bits/(double)(8*1024*1024) << std::endl;

    for(int i=0; i<pos_n; ++i)
        wbf_pos_keys_[i].data_ = pos_keys_[i];
    for(int i=0; i<neg_n; ++i)
        wbf_neg_keys_[i].data_ = neg_keys_[i];

    std::sort(wbf_neg_keys_, wbf_neg_keys_+neg_n, compare_func);

    cost_list_.resize(cost_storing_num);
    for(int i=0; i<cost_storing_num; ++i)
        cost_list_[i] = wbf_neg_keys_[i].data_;

    for(int i=cost_storing_num; i<neg_n; ++i)
        no_cost_list_.push_back(wbf_neg_keys_[i].data_);

    N = pos_n + neg_n;
    s = 0;

    for(int i=0; i< neg_n; i++)
        s += (1.0-(double)pos_n / (double)N) * (neg_keys_[i]->cost + 1.0);
    for(int i=0; i< pos_n; i++)
        s += (1.0-(double)pos_n / (double)N) * 1.0;

    // std::cout << s << std::endl;

    s2 = 0;
    for(int i=0; i< neg_n; i++){
        double eri = getEr(neg_keys_[i]->cost + 1.0, (double)pos_n / (double)N, s);
        s2 += log(eri) / (double)N;
    }
    for(int i=0; i< pos_n; i++){
        double eri = getEr(1.0, (double)pos_n / (double)N, s);
        s2 += log(eri) / (double)N;
    }
    // std::cout << s2 << std::endl;

    double Ere = getEr(1, (double)pos_n / (double)N, s);
    ke = ceil((double)bloom_.getSize() / (double)pos_n *log(2) + (log(Ere) - s2)/log(2));
    
    uint32_t op_ke = uint32_t((double)bloom_.getSize() / (double)pos_n * log(2));

    if(ke < op_ke) ke = op_ke;
    std::cout << ke << std::endl;

    bloom_.AddAll(pos_keys_, ke); 

    delete [] wbf_pos_keys_;
    delete [] wbf_neg_keys_;
}

bool WBF::Contain(Slice &key){
    uint32_t k = ke;
    for(Slice * key2 : cost_list_){
        if(key.str == key2->str){
            double Ere = getEr(key2->cost+1.0, (double)pos_n / (double)N, s);
            k = ceil((double)bloom_.getSize() / (double)pos_n *log(2) + (log(Ere) - s2)/log(2));
            if (k<ke)
                k = int((double)bloom_.getSize() / (double)pos_n * log(2));
            // std::cout << k << std::endl;
            break;
        }
    }
    return bloom_.Contain(key, k);
}

double WBF::testFPR(){
    uint32_t k = ke;
    double fpr_c = 0;
    double total_cost_ = 0.0;
    for(Slice * key2 : cost_list_){
        double Ere = getEr(key2->cost+1.0, (double)pos_n / (double)N, s);
        k = ceil((double)bloom_.getSize() / (double)pos_n *log(2) + (log(Ere) - s2)/log(2));
        if (k<ke)
            k = int((double)bloom_.getSize() / (double)pos_n * log(2));
        if(bloom_.Contain(*key2, k))
            fpr_c += key2->cost;
        total_cost_ += key2->cost;
    }
    for(Slice * key2 : no_cost_list_){
        if(bloom_.Contain(*key2, ke))
            fpr_c += key2->cost;
        total_cost_ += key2->cost;
    }
    return (double)(fpr_c*100) / total_cost_;
}
}

#endif