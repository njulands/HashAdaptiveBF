#ifndef HASH_SET
#define HASH_SET
#include <iostream>
#include <vector>
#include <string>
#include <stdint.h>
#include <assert.h>
#include "cityhash/city.h"
#include "murmurhash/MurMurHash.h"

// #define XXH_STATIC_LINKING_ONLY
// #include "xxhash/xxhash.h"

#define XXH_INLINE_ALL
#include "xxhash/xxhash.h"

class HashSetUtil;

typedef uint32_t (HashSetUtil::* HashFunction)(const char *str);
typedef uint32_t (HashSetUtil::* HashSeedFunction)(const char *str, uint32_t seed);

enum HashFamily {SDBMHash, RSHash, JSHash, PJWHash, BKDRHash, DJBHash, APHash, ELFHash, BOB, 
JenkinsOAAT, sBoxHash, TWMXHash, HsiehHash, NDJBHash, DEKHash, crc32Hash, BRPHash, PYHash, FNVHash, murmurhash, 
SuperFastHash, cityhash32, cityhash64, cityhash128, xxhash64, xxhash128};

// HashSetUtil
class HashSetUtil{
    public: 
        uint32_t GetHash(HashFamily hf, const std::string &str);
        uint32_t GetHashWithSeed(HashFamily hf, const std::string &str, uint32_t seed);

    private:
        uint32_t SDBMHash(const char *str);
        uint32_t RSHash(const char *str);
        uint32_t JSHash(const char *str);
        uint32_t PJWHash(const char *str);
        uint32_t BKDRHash(const char *str);
        uint32_t DJBHash(const char *str);
        uint32_t APHash(const char *str);
        uint32_t ELFHash(const char *str);
        uint32_t BOB(const char * str);
        uint32_t JenkinsOAAT(const char *str);
        uint32_t sBoxHash(const char *str);
        uint32_t TWMXHash(const char *str);
        uint32_t HsiehHash(const char *str);
        uint32_t NDJBHash(const char *str);
        uint32_t DEKHash(const char *str);
        uint32_t crc32Hash(const char * str);
        uint32_t BRPHash(const char * str);
        uint32_t PYHash(const char * str);
        uint32_t FNVHash(const char * str);
        uint32_t SuperFastHash(const char * str);

        uint32_t SDBMHash(const char *str, uint32_t seed);
        uint32_t RSHash(const char *str, uint32_t seed);
        uint32_t JSHash(const char *str, uint32_t seed);
        uint32_t PJWHash(const char *str, uint32_t seed);
        uint32_t BKDRHash(const char *str, uint32_t seed);
        uint32_t DJBHash(const char *str, uint32_t seed);
        uint32_t APHash(const char *str, uint32_t seed);
        uint32_t ELFHash(const char *str, uint32_t seed);
        uint32_t BOB(const char * str, uint32_t seed);
        uint32_t JenkinsOAAT(const char *str, uint32_t seed);
        uint32_t sBoxHash(const char *str, uint32_t seed);
        uint32_t TWMXHash(const char *str, uint32_t seed);
        uint32_t HsiehHash(const char *str, uint32_t seed);
        uint32_t NDJBHash(const char *str, uint32_t seed);
        uint32_t DEKHash(const char *str, uint32_t seed);
        uint32_t crc32Hash(const char * str, uint32_t seed);
        uint32_t BRPHash(const char * str, uint32_t seed);
        uint32_t PYHash(const char * str, uint32_t seed);
        uint32_t FNVHash(const char * str, uint32_t seed);
        uint32_t SuperFastHash(const char * str, uint32_t seed);

};


#endif