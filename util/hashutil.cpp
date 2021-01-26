#include "hashutil.h"
#include <string.h>

uint32_t HashSetUtil::GetHash(HashFamily hf, const std::string &str){
    const char * p = &str[0];
    switch (hf)
    {
        case HashFamily::SDBMHash:
            return SDBMHash(p);
        case HashFamily::RSHash:
            return RSHash(p);
        case HashFamily::JSHash:
            return JSHash(p);
        case HashFamily::PJWHash:
            return PJWHash(p);
        case HashFamily::BKDRHash:
            return BKDRHash(p);
        case HashFamily::DJBHash:
            return DJBHash(p);
        case HashFamily::APHash:
            return APHash(p);
        case HashFamily::ELFHash:
            return ELFHash(p);
        case HashFamily::BOB:
            return BOB(p);
        case HashFamily::JenkinsOAAT:
            return JenkinsOAAT(p);
        case HashFamily::TWMXHash:
            return TWMXHash(p);
        case HashFamily::HsiehHash:
            return HsiehHash(p);
        case HashFamily::NDJBHash:
            return NDJBHash(p);
        case HashFamily::DEKHash:
            return DEKHash(p);
        case HashFamily::crc32Hash:
            return crc32Hash(p);
        case HashFamily::BRPHash:
            return BRPHash(p);
        case HashFamily::PYHash:
            return PYHash(p);
        case HashFamily::FNVHash:
            return FNVHash(p);
        case HashFamily::murmurhash:
            return MurmurHash64A(p, str.size(), 0);
        case HashFamily::SuperFastHash:
            return SuperFastHash(p);
        case HashFamily::cityhash32:
            return CityHash32(p, str.size());
        case HashFamily::cityhash64:
            return CityHash64(p, str.size());
        case HashFamily::cityhash128:
            return Hash128to64(CityHash128(p, str.size()));
        case HashFamily::xxhash64:
            return XXH3_64bits(p, str.size());
        case HashFamily::xxhash128:
            return XXH3_128bits(p, str.size()).low64;
        default:
            return CityHash64(p, str.size());
    }
}

uint32_t HashSetUtil::GetHashWithSeed(HashFamily hf, const std::string &str, uint32_t seed){
    const char * p = &str[0];
    switch (hf)
    {
        case HashFamily::SDBMHash:
            return SDBMHash(p, seed);
        case HashFamily::RSHash:
            return RSHash(p, seed);
        case HashFamily::JSHash:
            return JSHash(p, seed);
        case HashFamily::PJWHash:
            return PJWHash(p, seed);
        case HashFamily::BKDRHash:
            return BKDRHash(p, seed);
        case HashFamily::DJBHash:
            return DJBHash(p, seed);
        case HashFamily::APHash:
            return APHash(p, seed);
        case HashFamily::ELFHash:
            return ELFHash(p, seed);
        case HashFamily::BOB:
            return BOB(p, seed);
        case HashFamily::JenkinsOAAT:
            return JenkinsOAAT(p, seed);
        case HashFamily::TWMXHash:
            return TWMXHash(p, seed);
        case HashFamily::HsiehHash:
            return HsiehHash(p, seed);
        case HashFamily::NDJBHash:
            return NDJBHash(p, seed);
        case HashFamily::DEKHash:
            return DEKHash(p, seed);
        case HashFamily::crc32Hash:
            return crc32Hash(p, seed);
        case HashFamily::BRPHash:
            return BRPHash(p, seed);
        case HashFamily::PYHash:
            return PYHash(p, seed);
        case HashFamily::FNVHash:
            return FNVHash(p, seed);
        case HashFamily::SuperFastHash:
            return SuperFastHash(p, seed);
        case HashFamily::cityhash32:
            return CityHash32(p, seed);
        case HashFamily::cityhash64:
            return CityHash64WithSeed(p, str.size(), seed);
        case HashFamily::cityhash128:
            return Hash128to64(CityHash128WithSeed(p, str.size(), std::make_pair((uint64_t)seed, (uint64_t)0)));
        case HashFamily::xxhash128:
            return XXH3_128bits_withSeed(p, str.size(), seed).low64;
        default:
            return XXH3_128bits_withSeed(p, str.size(), seed).low64;
    }
}

uint32_t HashSetUtil::SDBMHash(const char *str)
{
    uint32_t hash = 0;
 
    while (*str)
    {
        // equivalent to: hash = 65599*hash + (*str++);
        hash = (*str++) + (hash << 6) + (hash << 16) - hash;
    }
 
    return (hash & 0x7FFFFFFF);
}

uint32_t HashSetUtil::SDBMHash(const char *str, uint32_t seed)
{
    uint32_t hash = seed;
 
    while (*str)
    {
        // equivalent to: hash = 65599*hash + (*str++);
        hash = (*str++) + (hash << 6) + (hash << 16) - hash;
    }
 
    return (hash & 0x7FFFFFFF);
}

// RS Hash 
uint32_t HashSetUtil::RSHash(const char *str)
{
    uint32_t b = 378551;
    uint32_t a = 63689;
    uint32_t hash = 0;
 
    while (*str)
    {
        hash = hash * a + (*str++);
        a *= b;
    }
 
    return (hash & 0x7FFFFFFF);
}

uint32_t HashSetUtil::RSHash(const char *str, uint32_t seed)
{
    uint32_t b = 378551;
    uint32_t a = 63689;
    uint32_t hash = seed;
 
    while (*str)
    {
        hash = hash * a + (*str++);
        a *= b;
    }
 
    return (hash & 0x7FFFFFFF);
}

// JS Hash 
uint32_t HashSetUtil::JSHash(const char *str)
{
    uint32_t hash = 1315423911;
 
    while (*str)
    {
        hash ^= ((hash << 5) + (*str++) + (hash >> 2));
    }
 
    return (hash & 0x7FFFFFFF);
}

uint32_t HashSetUtil::JSHash(const char *str, uint32_t seed)
{
    uint32_t hash = seed;
 
    while (*str)
    {
        hash ^= ((hash << 5) + (*str++) + (hash >> 2));
    }
 
    return (hash & 0x7FFFFFFF);
}

// P. J. Weinberger Hash 
uint32_t HashSetUtil::PJWHash(const char *str)
{
    uint32_t BitsInUnignedInt = (uint32_t)(sizeof(uint32_t) * 8);
    uint32_t ThreeQuarters    = (uint32_t)((BitsInUnignedInt  * 3) / 4);
    uint32_t OneEighth = (uint32_t)(BitsInUnignedInt / 8);
    uint32_t HighBits = (uint32_t)(0xFFFFFFFF) << (BitsInUnignedInt 
                                               - OneEighth);
    uint32_t hash    = 0;
    uint32_t test    = 0;
 
    while (*str)
    {
        hash = (hash << OneEighth) + (*str++);
        if ((test = hash & HighBits) != 0)
        {
            hash = ((hash ^ (test >> ThreeQuarters)) & (~HighBits));
        }
    }
 
    return (hash & 0x7FFFFFFF);
}

 uint32_t HashSetUtil::PJWHash(const char *str, uint32_t seed)
{
    uint32_t BitsInUnignedInt = (uint32_t)(sizeof(uint32_t) * 8);
    uint32_t ThreeQuarters    = (uint32_t)((BitsInUnignedInt  * 3) / 4);
    uint32_t OneEighth = (uint32_t)(BitsInUnignedInt / 8);
    uint32_t HighBits = (uint32_t)(0xFFFFFFFF) << (BitsInUnignedInt 
                                               - OneEighth);
    uint32_t hash    = seed;
    uint32_t test    = 0;
 
    while (*str)
    {
        hash = (hash << OneEighth) + (*str++);
        if ((test = hash & HighBits) != 0)
        {
            hash = ((hash ^ (test >> ThreeQuarters)) & (~HighBits));
        }
    }
 
    return (hash & 0x7FFFFFFF);
}

// ELF Hash 
uint32_t HashSetUtil::ELFHash(const char *str)
{
    uint32_t hash = 0;
    uint32_t x    = 0;
 
    while (*str)
    {
        hash = (hash << 4) + (*str++);
        if ((x = hash & 0xF0000000L) != 0)
        {
            hash ^= (x >> 24);
            hash &= ~x;
        }
    }
 
    return (hash & 0x7FFFFFFF);
}

 uint32_t HashSetUtil::ELFHash(const char *str, uint32_t seed)
{
    uint32_t hash = seed;
    uint32_t x    = 0;
 
    while (*str)
    {
        hash = (hash << 4) + (*str++);
        if ((x = hash & 0xF0000000L) != 0)
        {
            hash ^= (x >> 24);
            hash &= ~x;
        }
    }
 
    return (hash & 0x7FFFFFFF);
}

// BKDR Hash 
uint32_t HashSetUtil::BKDRHash(const char *str)
{
    uint32_t seed = 131; // 31 131 1313 13131 131313 etc..
    uint32_t hash = 0;
 
    while (*str)
    {
        hash = hash * seed + (*str++);
    }
 
    return (hash & 0x7FFFFFFF);
}

 uint32_t HashSetUtil::BKDRHash(const char *str,  uint32_t seed)
{
    uint32_t hash = 0;
 
    while (*str)
    {
        hash = hash * seed + (*str++);
    }
 
    return (hash & 0x7FFFFFFF);
}

// DJB Hash 
uint32_t HashSetUtil::DJBHash(const char *str)
{
    uint32_t hash = 5381;
 
    while (*str)
    {
        hash += (hash << 5) + (*str++);
    }
 
    return (hash & 0x7FFFFFFF);
}

uint32_t HashSetUtil::DJBHash(const char *str,  uint32_t seed)
{
    uint32_t hash = seed;
 
    while (*str)
    {
        hash += (hash << 5) + (*str++);
    }
 
    return (hash & 0x7FFFFFFF);
}

// AP Hash 
uint32_t HashSetUtil::APHash(const char *str)
{
    uint32_t hash = 0;
    int i;
 
    for (i=0; *str; i++)
    {
        if ((i & 1) == 0)
        {
            hash ^= ((hash << 7) ^ (*str++) ^ (hash >> 3));
        }
        else
        {
            hash ^= (~((hash << 11) ^ (*str++) ^ (hash >> 5)));
        }
    }
 
    return (hash & 0x7FFFFFFF);
}

uint32_t HashSetUtil::APHash(const char *str,  uint32_t seed)
{
    uint32_t hash = seed;
    int i;
 
    for (i=0; *str; i++)
    {
        if ((i & 1) == 0)
        {
            hash ^= ((hash << 7) ^ (*str++) ^ (hash >> 3));
        }
        else
        {
            hash ^= (~((hash << 11) ^ (*str++) ^ (hash >> 5)));
        }
    }
 
    return (hash & 0x7FFFFFFF);
}

#define bobmix(a,b,c) \
{ \
  a -= b; a -= c; a ^= (c>>13); \
  b -= c; b -= a; b ^= (a<<8); \
  c -= a; c -= b; c ^= (b>>13); \
  a -= b; a -= c; a ^= (c>>12);  \
  b -= c; b -= a; b ^= (a<<16); \
  c -= a; c -= b; c ^= (b>>5); \
  a -= b; a -= c; a ^= (c>>3);  \
  b -= c; b -= a; b ^= (a<<10); \
  c -= a; c -= b; c ^= (b>>15); \
}

uint32_t HashSetUtil::BOB(const char * str)
{
	//register ub4 a,b,c,len;
	uint32_t a,b,c;
	uint32_t initval = 0;
	/* Set up the internal state */
	//len = length;
	a = b = 0x9e3779b9;  /* the golden ratio; an arbitrary value */
	c = initval;         /* the previous hash value */
    uint32_t len = strlen(str);
	/*---------------------------------------- handle most of the key */
	while (len >= 12)
	{
		a += (str[0] +((uint32_t)str[1]<<8) +((uint32_t)str[2]<<16) +((uint32_t)str[3]<<24));
		b += (str[4] +((uint32_t)str[5]<<8) +((uint32_t)str[6]<<16) +((uint32_t)str[7]<<24));
		c += (str[8] +((uint32_t)str[9]<<8) +((uint32_t)str[10]<<16)+((uint32_t)str[11]<<24));
		bobmix(a,b,c);
		str += 12; len -= 12;
	}

	/*------------------------------------- handle the last 11 bytes */
	c += len;
	switch(len)              /* all the case statements fall through */
	{
		case 11: c+=((uint32_t)str[10]<<24);
		case 10: c+=((uint32_t)str[9]<<16);
		case 9 : c+=((uint32_t)str[8]<<8);
		/* the first byte of c is reserved for the length */
		case 8 : b+=((uint32_t)str[7]<<24);
		case 7 : b+=((uint32_t)str[6]<<16);
		case 6 : b+=((uint32_t)str[5]<<8);
		case 5 : b+=str[4];
		case 4 : a+=((uint32_t)str[3]<<24);
		case 3 : a+=((uint32_t)str[2]<<16);
		case 2 : a+=((uint32_t)str[1]<<8);
		case 1 : a+=str[0];
		/* case 0: nothing left to add */
	}
	bobmix(a,b,c);
	/*-------------------------------------------- report the result */
	return c;
}

uint32_t HashSetUtil::BOB(const char * str,  uint32_t seed)
{
	//register ub4 a,b,c,len;
	uint32_t a,b,c;
	uint32_t initval = seed;
	/* Set up the internal state */
	//len = length;
	a = b = 0x9e3779b9;  /* the golden ratio; an arbitrary value */
	c = initval;         /* the previous hash value */
    uint32_t len = strlen(str);
	/*---------------------------------------- handle most of the key */
	while (len >= 12)
	{
		a += (str[0] +((uint32_t)str[1]<<8) +((uint32_t)str[2]<<16) +((uint32_t)str[3]<<24));
		b += (str[4] +((uint32_t)str[5]<<8) +((uint32_t)str[6]<<16) +((uint32_t)str[7]<<24));
		c += (str[8] +((uint32_t)str[9]<<8) +((uint32_t)str[10]<<16)+((uint32_t)str[11]<<24));
		bobmix(a,b,c);
		str += 12; len -= 12;
	}

	/*------------------------------------- handle the last 11 bytes */
	c += len;
	switch(len)              /* all the case statements fall through */
	{
		case 11: c+=((uint32_t)str[10]<<24);
		case 10: c+=((uint32_t)str[9]<<16);
		case 9 : c+=((uint32_t)str[8]<<8);
		/* the first byte of c is reserved for the length */
		case 8 : b+=((uint32_t)str[7]<<24);
		case 7 : b+=((uint32_t)str[6]<<16);
		case 6 : b+=((uint32_t)str[5]<<8);
		case 5 : b+=str[4];
		case 4 : a+=((uint32_t)str[3]<<24);
		case 3 : a+=((uint32_t)str[2]<<16);
		case 2 : a+=((uint32_t)str[1]<<8);
		case 1 : a+=str[0];
		/* case 0: nothing left to add */
	}
	bobmix(a,b,c);
	/*-------------------------------------------- report the result */
	return c;
}

//https://github.com/demerphq/smhasher/blob/d621653eea179344de4c75bc3cdb316eb96b1642/Hashes.cpp
uint32_t HashSetUtil::JenkinsOAAT(const char *str)
{
  uint32_t hash = 0;

  while (*str) {
    hash += *str++;
    hash += (hash << 10);
    hash ^= (hash >> 6);
  }
  hash += (hash << 3);
  hash ^= (hash >> 11);
  hash = hash + (hash << 15);
  return hash;
}

uint32_t HashSetUtil::JenkinsOAAT(const char *str,uint32_t seed)
{
  uint32_t hash = seed;

  while (*str) {
    hash += *str++;
    hash += (hash << 10);
    hash ^= (hash >> 6);
  }
  hash += (hash << 3);
  hash ^= (hash >> 11);
  hash = hash + (hash << 15);
  return hash;
}

uint32_t sBoxTable[256] = {
    0x4660c395, 0x3baba6c5, 0x27ec605b, 0xdfc1d81a, 0xaaac4406, 0x3783e9b8, 0xa4e87c68, 0x62dc1b2a,
    0xa8830d34, 0x10a56307, 0x4ba469e3, 0x54836450, 0x1b0223d4, 0x23312e32, 0xc04e13fe, 0x3b3d61fa,
    0xdab2d0ea, 0x297286b1, 0x73dbf93f, 0x6bb1158b, 0x46867fe2, 0xb7fb5313, 0x3146f063, 0x4fd4c7cb,
    0xa59780fa, 0x9fa38c24, 0x38c63986, 0xa0bac49f, 0xd47d3386, 0x49f44707, 0xa28dea30, 0xd0f30e6d,
    0xd5ca7704, 0x934698e3, 0x1a1ddd6d, 0xfa026c39, 0xd72f0fe6, 0x4d52eb70, 0xe99126df, 0xdfdaed86,
    0x4f649da8, 0x427212bb, 0xc728b983, 0x7ca5d563, 0x5e6164e5, 0xe41d3a24, 0x10018a23, 0x5a12e111,
    0x999ebc05, 0xf1383400, 0x50b92a7c, 0xa37f7577, 0x2c126291, 0x9daf79b2, 0xdea086b1, 0x85b1f03d,
    0x598ce687, 0xf3f5f6b9, 0xe55c5c74, 0x791733af, 0x39954ea8, 0xafcff761, 0x5fea64f1, 0x216d43b4,
    0xd039f8c1, 0xa6cf1125, 0xc14b7939, 0xb6ac7001, 0x138a2eff, 0x2f7875d6, 0xfe298e40, 0x4a3fad3b,
    0x066207fd, 0x8d4dd630, 0x96998973, 0xe656ac56, 0xbb2df109, 0x0ee1ec32, 0x03673d6c, 0xd20fb97d,
    0x2c09423c, 0x093eb555, 0xab77c1e2, 0x64607bf2, 0x945204bd, 0xe8819613, 0xb59de0e3, 0x5df7fc9a,
    0x82542258, 0xfb0ee357, 0xda2a4356, 0x5c97ab61, 0x8076e10d, 0x48e4b3cc, 0x7c28ec12, 0xb17986e1,
    0x01735836, 0x1b826322, 0x6602a990, 0x7c1cef68, 0xe102458e, 0xa5564a67, 0x1136b393, 0x98dc0ea1,
    0x3b6f59e5, 0x9efe981d, 0x35fafbe0, 0xc9949ec2, 0x62c765f9, 0x510cab26, 0xbe071300, 0x7ee1d449,
    0xcc71beef, 0xfbb4284e, 0xbfc02ce7, 0xdf734c93, 0x2f8cebcd, 0xfeedc6ab, 0x5476ee54, 0xbd2b5ff9,
    0xf4fd0352, 0x67f9d6ea, 0x7b70db05, 0x5a5f5310, 0x482dd7aa, 0xa0a66735, 0x321ae71f, 0x8e8ad56c,
    0x27a509c3, 0x1690b261, 0x4494b132, 0xc43a42a7, 0x3f60a7a6, 0xd63779ff, 0xe69c1659, 0xd15972c8,
    0x5f6cdb0c, 0xb9415af2, 0x1261ad8d, 0xb70a6135, 0x52ceda5e, 0xd4591dc3, 0x442b793c, 0xe50e2dee,
    0x6f90fc79, 0xd9ecc8f9, 0x063dd233, 0x6cf2e985, 0xe62cfbe9, 0x3466e821, 0x2c8377a2, 0x00b9f14e,
    0x237c4751, 0x40d4a33b, 0x919df7e8, 0xa16991a4, 0xc5295033, 0x5c507944, 0x89510e2b, 0xb5f7d902,
    0xd2d439a6, 0xc23e5216, 0xd52d9de3, 0x534a5e05, 0x762e73d4, 0x3c147760, 0x2d189706, 0x20aa0564,
    0xb07bbc3b, 0x8183e2de, 0xebc28889, 0xf839ed29, 0x532278f7, 0x41f8b31b, 0x762e89c1, 0xa1e71830,
    0xac049bfc, 0x9b7f839c, 0x8fd9208d, 0x2d2402ed, 0xf1f06670, 0x2711d695, 0x5b9e8fe4, 0xdc935762,
    0xa56b794f, 0xd8666b88, 0x6872c274, 0xbc603be2, 0x2196689b, 0x5b2b5f7a, 0x00c77076, 0x16bfa292,
    0xc2f86524, 0xdd92e83e, 0xab60a3d4, 0x92daf8bd, 0x1fe14c62, 0xf0ff82cc, 0xc0ed8d0a, 0x64356e4d,
    0x7e996b28, 0x81aad3e8, 0x05a22d56, 0xc4b25d4f, 0x5e3683e5, 0x811c2881, 0x124b1041, 0xdb1b4f02,
    0x5a72b5cc, 0x07f8d94e, 0xe5740463, 0x498632ad, 0x7357ffb1, 0x0dddd380, 0x3d095486, 0x2569b0a9,
    0xd6e054ae, 0x14a47e22, 0x73ec8dcc, 0x004968cf, 0xe0c3a853, 0xc9b50a03, 0xe1b0eb17, 0x57c6f281,
    0xc9f9377d, 0x43e03612, 0x9a0c4554, 0xbb2d83ff, 0xa818ffee, 0xf407db87, 0x175e3847, 0x5597168f,
    0xd3d547a7, 0x78f3157c, 0xfc750f20, 0x9880a1c6, 0x1af41571, 0x95d01dfc, 0xa3968d62, 0xeae03cf8,
    0x02ee4662, 0x5f1943ff, 0x252d9d1c, 0x6b718887, 0xe052f724, 0x4cefa30b, 0xdcc31a00, 0xe4d0024d,
    0xdbb4534a, 0xce01f5c8, 0x0c072b61, 0x5d59736a, 0x60291da4, 0x1fbe2c71, 0x2f11d09c, 0x9dce266a
};

//https://github.com/f-xyz/sbox-hash
uint32_t HashSetUtil::sBoxHash(const char *str) {
    uint32_t length = strlen(str);
    uint32_t len = length;
    uint32_t seed = 0;
    uint32_t hash = len + seed;
    for (; (len + 1)!=0 && (len + 1)!=1; len -= 2) {
        
        uint32_t iInput = length - len;
        unsigned char c1 = str[iInput] % 256;
        unsigned char c2 = str[iInput+1] % 256;
        hash = 3 * (((hash ^ sBoxTable[c1]) * 3) ^ sBoxTable[c2]);
    }

    if (len != 0) {
        hash = (hash ^ sBoxTable[str[0]]) * 3;
    }

    hash += (hash >> 22) ^ (hash << 4);

    return hash;
}

uint32_t HashSetUtil::sBoxHash(const char *str, uint32_t seed) {
    uint32_t length = strlen(str);
    uint32_t len = length;
    uint32_t hash = len + seed;
    for (; (len + 1)!=0 && (len + 1)!=1; len -= 2) {
        
        uint32_t iInput = length - len;
        unsigned char c1 = str[iInput] % 256;
        unsigned char c2 = str[iInput+1] % 256;
        hash = 3 * (((hash ^ sBoxTable[c1]) * 3) ^ sBoxTable[c2]);
    }

    if (len != 0) {
        hash = (hash ^ sBoxTable[str[0]]) * 3;
    }

    hash += (hash >> 22) ^ (hash << 4);

    return hash;
}

#define TWMXmix(key) \
{ \
  key += ~(key << 15); \
  key ^=  ((key & 0x7FFFFFFF) >> 10); \
  key +=  (key << 3); \
  key ^=  ((key & 0x7FFFFFFF) >> 6); \
  key += ~(key << 11); \
  key ^=  ((key & 0x7FFFFFFF) >> 16); \
}

// https://github.com/tubav/impd4e/blob/ddcec28414bbafb5cfb225688fe7c9d7b33de481/src/twmx.c
uint32_t HashSetUtil::TWMXHash(const char *str)
{
   const char *k = str;
   uint32_t len = strlen(str);
   uint32_t hash = 0;

  if (str == NULL) return 0;

   while (len >= 4)
   {
      hash += (k[0] +((uint32_t)k[1]<<8) +((uint32_t)k[2]<<16) +((uint32_t)k[3]<<24));
      TWMXmix(hash);
      k += 4; len -= 4;
   }

   switch(len)              /* all the case statements fall through */
   {
   case 3 : hash+=((uint32_t)k[2]<<16);
   case 2 : hash+=((uint32_t)k[1]<<8);
   case 1 : hash+=k[0];
     /* case 0: nothing left to add */
   }
   TWMXmix(hash);

   return hash;
}

uint32_t HashSetUtil::TWMXHash(const char *str, uint32_t seed)
{
   const char *k = str;
   uint32_t len = strlen(str);
   uint32_t hash = seed;

  if (str == NULL) return 0;

   while (len >= 4)
   {
      hash += (k[0] +((uint32_t)k[1]<<8) +((uint32_t)k[2]<<16) +((uint32_t)k[3]<<24));
      TWMXmix(hash);
      k += 4; len -= 4;
   }

   switch(len)              /* all the case statements fall through */
   {
   case 3 : hash+=((uint32_t)k[2]<<16);
   case 2 : hash+=((uint32_t)k[1]<<8);
   case 1 : hash+=k[0];
     /* case 0: nothing left to add */
   }
   TWMXmix(hash);

   return hash;
}

//https://github.com/tubav/impd4e/blob/ddcec28414bbafb5cfb225688fe7c9d7b33de481/src/hsieh.c
#define get16bits(d) ((((uint32_t)(((const unsigned char *)(d))[1])) << 8)\
                       +(uint32_t)(((const unsigned char *)(d))[0]) )
uint32_t HashSetUtil::HsiehHash(const char * str) {
    uint32_t len = strlen(str);
    uint32_t hash = len, tmp;
    int rem;

    if (len <= 0 || str == NULL) return 0;

    rem = len & 3;
    len >>= 2;

    /* Main loop */
    for (;len > 0; len--) {
        hash  += get16bits (str);
        tmp    = (get16bits (str+2) << 11) ^ hash;
        hash   = (hash << 16) ^ tmp;
        str  += 2*sizeof (uint32_t);
        hash  += hash >> 11;
    }

    /* Handle end cases */
    switch (rem) {
        case 3: hash += get16bits (str);
                hash ^= hash << 16;
                hash ^= str[sizeof (unsigned short)] << 18;
                hash += hash >> 11;
                break;
        case 2: hash += get16bits (str);
                hash ^= hash << 11;
                hash += hash >> 17;
                break;
        case 1: hash += *str;
                hash ^= hash << 10;
                hash += hash >> 1;
    }

    /* Force "avalanching" of final 127 bits */
    hash ^= hash << 3;
    hash += hash >> 5;
    hash ^= hash << 4;
    hash += hash >> 17;
    hash ^= hash << 25;
    hash += hash >> 6;

    return hash;
}

uint32_t HashSetUtil::HsiehHash(const char * str, uint32_t seed) {
    uint32_t len = strlen(str);
    uint32_t hash = seed, tmp;
    int rem;

    if (len <= 0 || str == NULL) return 0;

    rem = len & 3;
    len >>= 2;

    /* Main loop */
    for (;len > 0; len--) {
        hash  += get16bits (str);
        tmp    = (get16bits (str+2) << 11) ^ hash;
        hash   = (hash << 16) ^ tmp;
        str  += 2*sizeof (uint32_t);
        hash  += hash >> 11;
    }

    /* Handle end cases */
    switch (rem) {
        case 3: hash += get16bits (str);
                hash ^= hash << 16;
                hash ^= str[sizeof (unsigned short)] << 18;
                hash += hash >> 11;
                break;
        case 2: hash += get16bits (str);
                hash ^= hash << 11;
                hash += hash >> 17;
                break;
        case 1: hash += *str;
                hash ^= hash << 10;
                hash += hash >> 1;
    }

    /* Force "avalanching" of final 127 bits */
    hash ^= hash << 3;
    hash += hash >> 5;
    hash ^= hash << 4;
    hash += hash >> 17;
    hash ^= hash << 25;
    hash += hash >> 6;

    return hash;
}

//https://github.com/alvatar/snippets/blob/89ce0c71ff1554a9307ed0726b84dfa5782babe6/scheme/hashes/.svn/text-base/NDJBHash.scm.svn-base
uint32_t HashSetUtil::NDJBHash(const char *str)
{
    uint32_t hash = 5381;
    uint32_t len = strlen(str);
    for ( ; len; str++, len--) {
        hash = (hash * 33) ^ ((uint32_t) *str);
    }

    return hash;
}

uint32_t HashSetUtil::NDJBHash(const char *str, uint32_t seed)
{
    uint32_t hash = seed;
    uint32_t len = strlen(str);
    for ( ; len; str++, len--) {
        hash = (hash * 33) ^ ((uint32_t) *str);
    }

    return hash;
}

//https://github.com/alvatar/snippets/blob/89ce0c71ff1554a9307ed0726b84dfa5782babe6/scheme/hashes/.svn/text-base/DEKHash.scm.svn-base
uint32_t HashSetUtil::DEKHash(const char *str)
{
    uint32_t len = strlen(str);
    uint32_t hash = len;

    for ( ; len; str++, len--) {
        hash = ((hash << 5) ^ (hash >> 27)) ^ ((uint32_t) *str);
    }


    return hash;
}

uint32_t HashSetUtil::DEKHash(const char *str, uint32_t seed)
{
    uint32_t len = strlen(str);
    uint32_t hash = seed;

    for ( ; len; str++, len--) {
        hash = ((hash << 5) ^ (hash >> 27)) ^ ((uint32_t) *str);
    }


    return hash;
}

static const uint32_t crc_table[256] = {
  0x00000000L, 0x77073096L, 0xee0e612cL, 0x990951baL, 0x076dc419L,
  0x706af48fL, 0xe963a535L, 0x9e6495a3L, 0x0edb8832L, 0x79dcb8a4L,
  0xe0d5e91eL, 0x97d2d988L, 0x09b64c2bL, 0x7eb17cbdL, 0xe7b82d07L,
  0x90bf1d91L, 0x1db71064L, 0x6ab020f2L, 0xf3b97148L, 0x84be41deL,
  0x1adad47dL, 0x6ddde4ebL, 0xf4d4b551L, 0x83d385c7L, 0x136c9856L,
  0x646ba8c0L, 0xfd62f97aL, 0x8a65c9ecL, 0x14015c4fL, 0x63066cd9L,
  0xfa0f3d63L, 0x8d080df5L, 0x3b6e20c8L, 0x4c69105eL, 0xd56041e4L,
  0xa2677172L, 0x3c03e4d1L, 0x4b04d447L, 0xd20d85fdL, 0xa50ab56bL,
  0x35b5a8faL, 0x42b2986cL, 0xdbbbc9d6L, 0xacbcf940L, 0x32d86ce3L,
  0x45df5c75L, 0xdcd60dcfL, 0xabd13d59L, 0x26d930acL, 0x51de003aL,
  0xc8d75180L, 0xbfd06116L, 0x21b4f4b5L, 0x56b3c423L, 0xcfba9599L,
  0xb8bda50fL, 0x2802b89eL, 0x5f058808L, 0xc60cd9b2L, 0xb10be924L,
  0x2f6f7c87L, 0x58684c11L, 0xc1611dabL, 0xb6662d3dL, 0x76dc4190L,
  0x01db7106L, 0x98d220bcL, 0xefd5102aL, 0x71b18589L, 0x06b6b51fL,
  0x9fbfe4a5L, 0xe8b8d433L, 0x7807c9a2L, 0x0f00f934L, 0x9609a88eL,
  0xe10e9818L, 0x7f6a0dbbL, 0x086d3d2dL, 0x91646c97L, 0xe6635c01L,
  0x6b6b51f4L, 0x1c6c6162L, 0x856530d8L, 0xf262004eL, 0x6c0695edL,
  0x1b01a57bL, 0x8208f4c1L, 0xf50fc457L, 0x65b0d9c6L, 0x12b7e950L,
  0x8bbeb8eaL, 0xfcb9887cL, 0x62dd1ddfL, 0x15da2d49L, 0x8cd37cf3L,
  0xfbd44c65L, 0x4db26158L, 0x3ab551ceL, 0xa3bc0074L, 0xd4bb30e2L,
  0x4adfa541L, 0x3dd895d7L, 0xa4d1c46dL, 0xd3d6f4fbL, 0x4369e96aL,
  0x346ed9fcL, 0xad678846L, 0xda60b8d0L, 0x44042d73L, 0x33031de5L,
  0xaa0a4c5fL, 0xdd0d7cc9L, 0x5005713cL, 0x270241aaL, 0xbe0b1010L,
  0xc90c2086L, 0x5768b525L, 0x206f85b3L, 0xb966d409L, 0xce61e49fL,
  0x5edef90eL, 0x29d9c998L, 0xb0d09822L, 0xc7d7a8b4L, 0x59b33d17L,
  0x2eb40d81L, 0xb7bd5c3bL, 0xc0ba6cadL, 0xedb88320L, 0x9abfb3b6L,
  0x03b6e20cL, 0x74b1d29aL, 0xead54739L, 0x9dd277afL, 0x04db2615L,
  0x73dc1683L, 0xe3630b12L, 0x94643b84L, 0x0d6d6a3eL, 0x7a6a5aa8L,
  0xe40ecf0bL, 0x9309ff9dL, 0x0a00ae27L, 0x7d079eb1L, 0xf00f9344L,
  0x8708a3d2L, 0x1e01f268L, 0x6906c2feL, 0xf762575dL, 0x806567cbL,
  0x196c3671L, 0x6e6b06e7L, 0xfed41b76L, 0x89d32be0L, 0x10da7a5aL,
  0x67dd4accL, 0xf9b9df6fL, 0x8ebeeff9L, 0x17b7be43L, 0x60b08ed5L,
  0xd6d6a3e8L, 0xa1d1937eL, 0x38d8c2c4L, 0x4fdff252L, 0xd1bb67f1L,
  0xa6bc5767L, 0x3fb506ddL, 0x48b2364bL, 0xd80d2bdaL, 0xaf0a1b4cL,
  0x36034af6L, 0x41047a60L, 0xdf60efc3L, 0xa867df55L, 0x316e8eefL,
  0x4669be79L, 0xcb61b38cL, 0xbc66831aL, 0x256fd2a0L, 0x5268e236L,
  0xcc0c7795L, 0xbb0b4703L, 0x220216b9L, 0x5505262fL, 0xc5ba3bbeL,
  0xb2bd0b28L, 0x2bb45a92L, 0x5cb36a04L, 0xc2d7ffa7L, 0xb5d0cf31L,
  0x2cd99e8bL, 0x5bdeae1dL, 0x9b64c2b0L, 0xec63f226L, 0x756aa39cL,
  0x026d930aL, 0x9c0906a9L, 0xeb0e363fL, 0x72076785L, 0x05005713L,
  0x95bf4a82L, 0xe2b87a14L, 0x7bb12baeL, 0x0cb61b38L, 0x92d28e9bL,
  0xe5d5be0dL, 0x7cdcefb7L, 0x0bdbdf21L, 0x86d3d2d4L, 0xf1d4e242L,
  0x68ddb3f8L, 0x1fda836eL, 0x81be16cdL, 0xf6b9265bL, 0x6fb077e1L,
  0x18b74777L, 0x88085ae6L, 0xff0f6a70L, 0x66063bcaL, 0x11010b5cL,
  0x8f659effL, 0xf862ae69L, 0x616bffd3L, 0x166ccf45L, 0xa00ae278L,
  0xd70dd2eeL, 0x4e048354L, 0x3903b3c2L, 0xa7672661L, 0xd06016f7L,
  0x4969474dL, 0x3e6e77dbL, 0xaed16a4aL, 0xd9d65adcL, 0x40df0b66L,
  0x37d83bf0L, 0xa9bcae53L, 0xdebb9ec5L, 0x47b2cf7fL, 0x30b5ffe9L,
  0xbdbdf21cL, 0xcabac28aL, 0x53b39330L, 0x24b4a3a6L, 0xbad03605L,
  0xcdd70693L, 0x54de5729L, 0x23d967bfL, 0xb3667a2eL, 0xc4614ab8L,
  0x5d681b02L, 0x2a6f2b94L, 0xb40bbe37L, 0xc30c8ea1L, 0x5a05df1bL,
  0x2d02ef8dL
};

/* ========================================================================= */

#define DO1(buf) crc = crc_table[((int)crc ^ (*buf++)) & 0xff] ^ (crc >> 8);
#define DO2(buf)  DO1(buf); DO1(buf);
#define DO4(buf)  DO2(buf); DO2(buf);
#define DO8(buf)  DO4(buf); DO4(buf);

/* ========================================================================= */
// https://github.com/aappleby/smhasher/blob/master/src/crc.cpp
uint32_t HashSetUtil::crc32Hash(const char * str)
{
  uint32_t seed = 0;
  uint32_t len = strlen(str);
  uint32_t crc = seed ^ 0xffffffffL;

  while (len >= 8)
  {
    DO8(str);
    len -= 8;
  }

  while(len--)
  {
    DO1(str);
  } 

  crc ^= 0xffffffffL;

  return crc;
}

uint32_t HashSetUtil::crc32Hash(const char * str, uint32_t seed)
{
  uint32_t len = strlen(str);
  uint32_t crc = seed ^ 0xffffffffL;

  while (len >= 8)
  {
    DO8(str);
    len -= 8;
  }

  while(len--)
  {
    DO1(str);
  } 

  crc ^= 0xffffffffL;

  return crc;
}
// https://github.com/alvatar/snippets/blob/89ce0c71ff1554a9307ed0726b84dfa5782babe6/scheme/hashes/.svn/text-base/BRPHash.scm.svn-base

#define MASK (((uint32_t) (~0)) >> 2)
uint32_t HashSetUtil::BRPHash(const char * str)
{
        uint32_t len = strlen(str);
        uint32_t hash = 0;

        for ( ; len; str++, len--) {
            hash = (hash & MASK) ^ (hash << 6) ^ ((uint32_t ) *str);
        }
    return hash;

}

uint32_t HashSetUtil::BRPHash(const char * str, uint32_t seed)
{
        uint32_t len = strlen(str);
        uint32_t hash = seed;

        for ( ; len; str++, len--) {
            hash = (hash & MASK) ^ (hash << 6) ^ ((uint32_t ) *str);
        }
    return hash;

}
// https://github.com/alvatar/snippets/blob/89ce0c71ff1554a9307ed0726b84dfa5782babe6/scheme/hashes/PYHash.scm
uint32_t HashSetUtil::PYHash(const char * str)
{
    uint32_t hash = ((uint32_t) *str) << 7;
    uint32_t len = strlen(str);
    for (int i = 0; i < len; str++, i++) {
        hash = (1000003 * hash) ^ ((uint32_t) *str);
    }

    hash ^= len;

    return (((uint32_t) -1) == hash) ? ((uint32_t) -2) : hash;
}
uint32_t HashSetUtil::PYHash(const char * str, uint32_t seed)
{
    uint32_t hash = seed;
    uint32_t len = strlen(str);
    for (int i = 0; i < len; str++, i++) {
        hash = (1000003 * hash) ^ ((uint32_t) *str);
    }

    hash ^= len;

    return (((uint32_t) -1) == hash) ? ((uint32_t) -2) : hash;
}


// https://github.com/alvatar/snippets/blob/89ce0c71ff1554a9307ed0726b84dfa5782babe6/scheme/hashes/FNVHash.scm
uint32_t HashSetUtil::FNVHash(const char * str)
{
    uint32_t hash = ((uint32_t) 0x811c9dc5);
    uint32_t len = strlen(str);

    /* FNV-1 hash each octet in the buffer */
    for (int i = 0; i < len; i++) {

        /* multiply by the 32 bit FNV magic prime mod 2^32 */
        hash += (hash << 1) + (hash << 4) + (hash << 7) + (hash << 8) + (hash << 24);

        /* xor the bottom with the current octet */
        hash ^= ((uint32_t) *str++);
    }


    /* return our new hash value */
    return hash;

}
uint32_t HashSetUtil::FNVHash(const char * str, uint32_t seed)
{
    uint32_t hash = seed;
    uint32_t len = strlen(str);

    /* FNV-1 hash each octet in the buffer */
    for (int i = 0; i < len; i++) {

        /* multiply by the 32 bit FNV magic prime mod 2^32 */
        hash += (hash << 1) + (hash << 4) + (hash << 7) + (hash << 8) + (hash << 24);

        /* xor the bottom with the current octet */
        hash ^= ((uint32_t) *str++);
    }


    /* return our new hash value */
    return hash;

}

unsigned short SFHGet16bits ( const void * p )
{
  return *(const unsigned short*)p;
}
// https://github.com/demerphq/smhasher/blob/d621653eea179344de4c75bc3cdb316eb96b1642/Hashes.cpp
uint32_t HashSetUtil::SuperFastHash(const char * str) {
    uint32_t len = strlen(str);
    uint32_t hash = len;
    uint32_t tmp;
    int rem;

  if (len <= 0 || str == NULL) return 0;

  rem = len & 3;
  len >>= 2;

  /* Main loop */
  for (;len > 0; len--) {
    hash  += SFHGet16bits (str);
    tmp    = (SFHGet16bits (str+2) << 11) ^ hash;
    hash   = (hash << 16) ^ tmp;
    str  += 2*sizeof (unsigned short);
    hash  += hash >> 11;
  }

  /* Handle end cases */
  switch (rem) {
    case 3:	hash += SFHGet16bits (str);
        hash ^= hash << 16;
        hash ^= str[sizeof (unsigned short)] << 18;
        hash += hash >> 11;
        break;
    case 2:	hash += SFHGet16bits (str);
        hash ^= hash << 11;
        hash += hash >> 17;
        break;
    case 1: hash += *str;
        hash ^= hash << 10;
        hash += hash >> 1;
  }

  /* Force "avalanching" of final 127 bits */
  hash ^= hash << 3;
  hash += hash >> 5;
  hash ^= hash << 4;
  hash += hash >> 17;
  hash ^= hash << 25;
  hash += hash >> 6;

  return hash;
}
uint32_t HashSetUtil::SuperFastHash(const char * str, uint32_t seed) {
    uint32_t len = strlen(str);
    uint32_t hash = seed;
    uint32_t tmp;
    int rem;

  if (len <= 0 || str == NULL) return 0;

  rem = len & 3;
  len >>= 2;

  /* Main loop */
  for (;len > 0; len--) {
    hash  += SFHGet16bits (str);
    tmp    = (SFHGet16bits (str+2) << 11) ^ hash;
    hash   = (hash << 16) ^ tmp;
    str  += 2*sizeof (unsigned short);
    hash  += hash >> 11;
  }

  /* Handle end cases */
  switch (rem) {
    case 3:	hash += SFHGet16bits (str);
        hash ^= hash << 16;
        hash ^= str[sizeof (unsigned short)] << 18;
        hash += hash >> 11;
        break;
    case 2:	hash += SFHGet16bits (str);
        hash ^= hash << 11;
        hash += hash >> 17;
        break;
    case 1: hash += *str;
        hash ^= hash << 10;
        hash += hash >> 1;
  }

  /* Force "avalanching" of final 127 bits */
  hash ^= hash << 3;
  hash += hash >> 5;
  hash ^= hash << 4;
  hash += hash >> 17;
  hash ^= hash << 25;
  hash += hash >> 6;

  return hash;
}

