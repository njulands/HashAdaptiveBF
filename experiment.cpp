
#include <time.h>

#include "util/dataloader.h"
#include "habf/habf.h"
#include "habf/fasthabf.h"
#include "nonlearnedfilter/bloom.h"
#include "nonlearnedfilter/wbf.h"
#include "nonlearnedfilter/xorfilter/xorfilter_2.h"


using namespace std;
using namespace std::chrono;

void testBloom_XXH128(dataloader &dl, double space_size_);
void testHABF(dataloader &dl, double space_size_);
void testFastHABF(dataloader &dl, double bits_per_key);
void testXorFilter(dataloader &dl, double bits_per_key);
void testWBF(dataloader &dl, double bits_per_key);
int main(){

    const string dataset_ = "shalla"; // shalla or ycsb
    const double space_size_ = 1.5; // total space size of data structure (MB)
    dataloader dl;
    if(dl.load(dataset_, true)){
        cout << "keys size (positives, negatives): (" << dl.pos_keys_.size() << "," << dl.neg_keys_.size() << ")" << endl;
        std::cout << "Space size allocated for each data structure: (" << space_size_ << " MB)" << endl;
        double bits_per_key = (double)(space_size_ * 1024 * 1024 * 8) / (double) dl.pos_keys_.size();
        std::cout << "bits_per_key: (" << bits_per_key << ")" << endl;

        if(bits_per_key < 5){
            std::cout << "bits_per_key is too small" << endl;
            return 0;
        }

        cout <<"Bloom (XXH128) testing......" << endl;
        testBloom_XXH128(dl, bits_per_key);

        std::cout << "---------------------------------------------------" << endl;

        cout <<"Xorfilter testing......" << endl;
        testXorFilter(dl, bits_per_key);

        std::cout << "---------------------------------------------------" << endl;

        cout <<"f-HABFfilter testing......" << endl;
        testFastHABF(dl, bits_per_key);

        std::cout << "---------------------------------------------------" << endl;

        cout <<"HABFfilter testing......" << endl;
        testHABF(dl, bits_per_key);

        // std::cout << "---------------------------------------------------" << endl;

        // cout <<"WBF testing......" << endl;
        // testWBF(dl, bits_per_key);

    }else{
        cout << "dataset reading failed"  << endl;
    }
}

void testWBF(dataloader &dl, double bits_per_key){
    wbf::WBF wbf(bits_per_key);
    auto t1 = steady_clock::now();

    wbf.AddAll(dl.pos_keys_, dl.neg_keys_);

    auto t2 = steady_clock::now();
    double construct_time_ = duration<double>(t2 - t1).count();
    std::cout<<"construct_time_:"<< construct_time_ <<endl;

    // query testing
    double fpr_c = 0;
    int fnr_c = 0;
    double total_cost_ = 0.0;
    t1 = steady_clock::now();
    for(int i=0; i<dl.pos_keys_.size(); i++)
        if(!wbf.Contain(*dl.pos_keys_[i])) fnr_c++;
    for(int i=0; i<dl.neg_keys_.size(); i++){
        if(wbf.Contain(*dl.neg_keys_[i]))
            fpr_c += dl.neg_keys_[i]->cost;
        total_cost_ += dl.neg_keys_[i]->cost;
    }
        
    t2 = steady_clock::now();
    double query_perkey_time =  1000000000.0 * duration<double>(t2 - t1).count() / (double)(dl.pos_keys_.size()+dl.neg_keys_.size());
    std::cout<<"query time/key (ns):"<< query_perkey_time <<endl;
    std::cout<<"Weighted FPR:"<< (double)(fpr_c*100) / total_cost_ << "%" <<endl;
    std::cout<<"FNR:"<< (double)(fnr_c*100) / (double)dl.pos_keys_.size() << "%" <<endl;
}

void testBloom_XXH128(dataloader &dl, double bits_per_key){
    bloomfilter::BloomFilter<Slice, false> bf(dl.pos_keys_.size(), bits_per_key);

    // construction testing
    auto t1 = steady_clock::now();

    bf.AddAll(dl.pos_keys_);

    auto t2 = steady_clock::now();

    double construct_time_ = duration<double>(t2 - t1).count();
    std::cout<<"construct time/key (ns):"<< 1000000000.0 * construct_time_ / (double)(dl.pos_keys_.size()) <<endl;

    // query testing
    double fpr_c = 0;
    int fnr_c = 0;
    double total_cost_ = 0.0;
    t1 = steady_clock::now();
    for(int i=0; i<dl.pos_keys_.size(); i++)
        if(bf.Contain(*dl.pos_keys_[i]) != bloomfilter::Ok) fnr_c++;
    for(int i=0; i<dl.neg_keys_.size(); i++){
        if(bf.Contain(*dl.neg_keys_[i]) == bloomfilter::Ok)
            fpr_c += dl.neg_keys_[i]->cost;
        total_cost_ += dl.neg_keys_[i]->cost;
    }
        
    t2 = steady_clock::now();
    std::cout << total_cost_ << std::endl;
    double query_perkey_time =  1000000000.0 * duration<double>(t2 - t1).count() / (double)(dl.pos_keys_.size()+dl.neg_keys_.size());
    std::cout<<"query time/key (ns):"<< query_perkey_time <<endl;
    std::cout<<"Weighted FPR:"<< (double)(fpr_c*100) / total_cost_ << "%" <<endl;
    std::cout<<"FNR:"<< (double)(fnr_c*100) / (double)dl.pos_keys_.size() << "%" <<endl;
}

//test HABF
void testHABF(dataloader &dl, double bits_per_key){

    // construction testing
    auto t1 = steady_clock::now();

    habf::HABFilter habf(bits_per_key, dl.pos_keys_.size());

    habf.AddAndOptimize(dl.pos_keys_, dl.neg_keys_);

    auto t2 = steady_clock::now();

    double construct_time_ = duration<double>(t2 - t1).count();
    std::cout<<"construct time/key (ns):"<< 1000000000.0 * construct_time_ / (double)(dl.pos_keys_.size()) <<endl;

    // query testing
    double fpr_c = 0;
    int fnr_c = 0;
    double total_cost_ = 0.0;
    t1 = steady_clock::now();
    for(int i=0; i<dl.pos_keys_.size(); i++)
        if(!habf.Contain(*dl.pos_keys_[i])) fnr_c++;
    for(int i=0; i<dl.neg_keys_.size(); i++){
        if(habf.Contain(*dl.neg_keys_[i]))
            fpr_c += dl.neg_keys_[i]->cost;
        total_cost_ += dl.neg_keys_[i]->cost;
    }
        
    t2 = steady_clock::now();
    double query_perkey_time =  1000000000.0 * duration<double>(t2 - t1).count() / (double)(dl.pos_keys_.size()+dl.neg_keys_.size());
    std::cout<<"query time/key (ns):"<< query_perkey_time <<endl;
    std::cout<<"Weighted FPR:"<< (double)(fpr_c*100) / total_cost_ << "%" <<endl;
    std::cout<<"FNR:"<< (double)(fnr_c*100) / (double)dl.pos_keys_.size() << "%" <<endl;
}

//test f-HABF
void testFastHABF(dataloader &dl, double bits_per_key){

    // construction testing
    auto t1 = steady_clock::now();

    fasthabf::FastHABFilter fhabf(bits_per_key, dl.pos_keys_.size());

    fhabf.AddAndOptimize(dl.pos_keys_, dl.neg_keys_);

    auto t2 = steady_clock::now();

    double construct_time_ = duration<double>(t2 - t1).count();
    std::cout<<"construct time/key (ns):"<< 1000000000.0 * construct_time_ / (double)(dl.pos_keys_.size()) <<endl;

    // query testing
    double fpr_c = 0;
    int fnr_c = 0;
    double total_cost_ = 0.0;
    t1 = steady_clock::now();
    for(int i=0; i<dl.pos_keys_.size(); i++)
        if(!fhabf.Contain(*dl.pos_keys_[i])) fnr_c++;
    for(int i=0; i<dl.neg_keys_.size(); i++){
        if(fhabf.Contain(*dl.neg_keys_[i]))
            fpr_c += dl.neg_keys_[i]->cost;
        total_cost_ += dl.neg_keys_[i]->cost;
    }
        
    t2 = steady_clock::now();
    double query_perkey_time =  1000000000.0 * duration<double>(t2 - t1).count() / (double)(dl.pos_keys_.size()+dl.neg_keys_.size());
    std::cout<<"query time/key (ns):"<< query_perkey_time <<endl;
    std::cout<<"Weighted FPR:"<< (double)(fpr_c*100) / total_cost_ << "%" <<endl;
    std::cout<<"FNR:"<< (double)(fnr_c*100) / (double)dl.pos_keys_.size() << "%" <<endl;
    
}

// test xorfilter
void testXorFilter(dataloader &dl, double bits_per_key){
    int fingerprint_size_ = bits_per_key / 1.23;
    double mul = bits_per_key / (double) fingerprint_size_;
    xorfilter2::XorFilter2<uint64_t, uint16_t, NBitArray<uint16_t>, SimpleMixSplit> xorf(dl.pos_keys_.size(), mul, fingerprint_size_);
    
    // construction testing
    auto t1 = steady_clock::now();

    xorf.AddAll(dl.pos_keys_, 0, dl.pos_keys_.size());

    auto t2 = steady_clock::now();

    double construct_time_ = duration<double>(t2 - t1).count();
    std::cout<<"construct time/key (ns):"<< 1000000000.0 * construct_time_ / (double)(dl.pos_keys_.size()) <<endl;


    // query testing
    double fpr_c = 0;
    int fnr_c = 0;
    double total_cost_ = 0.0;
    t1 = steady_clock::now();
    for(int i=0; i<dl.pos_keys_.size(); i++)
        if(xorf.Contain(*dl.pos_keys_[i]) != xorfilter2::Ok) fnr_c++;
    for(int i=0; i<dl.neg_keys_.size(); i++){
        if(xorf.Contain(*dl.neg_keys_[i]) == xorfilter2::Ok)
            fpr_c += dl.neg_keys_[i]->cost;
        total_cost_ += dl.neg_keys_[i]->cost;
    }

    t2 = steady_clock::now();
    double query_perkey_time =  1000000000.0 * duration<double>(t2 - t1).count() / (double)(dl.pos_keys_.size()+dl.neg_keys_.size());
    std::cout<<"query time/key (ns):"<< query_perkey_time <<endl;
    std::cout<<"Weighted FPR:"<< (double)(fpr_c*100) / total_cost_ << "%" <<endl;
    std::cout<<"FNR:"<< (double)(fnr_c*100) / (double)dl.pos_keys_.size() << "%" <<endl;
}
