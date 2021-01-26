# ICDE-2021 paper #455 - Hash-Adaptive-Bloom-Filter

# About this repo
This repo contains the source code of HABF and all comparison algorithms in our experiments, which are as shown in the following table.

|Algorithm| Description|
|:----:|----|
|HABF| Implementation - [/habf/habf.h](https://github.com/AnonymousAuthor455/HashAdaptiveBF/blob/master/habf/habf.h)|
|f-HABF| Implementation - [/habf/fasthabf.h](https://github.com/AnonymousAuthor455/HashAdaptiveBF/blob/master/habf/fasthabf.h)|
|Bloom filter|B. H. Bloom, “Space/time trade-offs in hash coding with allowable errors,” Communications of ACM, 1970. Implementation: https://github.com/FastFilter/fastfilter_cpp|
|Xor|T. M. Graf and D. Lemire, “Xor filters: Faster and smaller than bloom and cuckoo filters,” Journal of Experimental Algorithmics, 2020. Implementation: https://github.com/FastFilter/fastfilter_cpp|
|WBF|J. Bruck, J. Gao, and A. Jiang, “Weighted Bloom filter,” in Proceedings of International Symposium on Information Theory. IEEE, 2006. Implementation - [/nonlearnedfilter/wbf.h](https://github.com/AnonymousAuthor455/HashAdaptiveBF/blob/master/nonlearnedfilter/wbf.h)|
|LBF|T. Kraska, A. Beutel, E. H. Chi, J. Dean, and N. Polyzotis, “The case for learned index structures,” in Proceedings of the International Conference on Management of Data. ACM, 2018. Implementation: https://github.com/karan1149/DeepBloom/.|
|SLBF|M. Mitzenmacher, “A model for learned Bloom filters and optimizing by sandwiching,” in Advances in Neural Information Processing Systems. Curran Associates, Inc., 2018.Implementation - [/learnedfilter/SLBF](https://github.com/AnonymousAuthor455/HashAdaptiveBF/blob/master/learnedfilter/SLBF)|
|AdaBF|Z. Dai and A. Shrivastava, “Adaptive learned Bloom filter (Ada-BF): Efficient utilization of the classifier,” arXiv preprint, 2019. Implementation: https://github.com/DAIZHENWEI/Ada-BF|

# Requirement 
   1. cmake@3+
   2. make
# Build

Build benchmarking executable file
```bash
mkdir -p build && cd build && cmake .. && make
```
# BenchMarking
Running benchmark executable file
```Bash
./experiment
```
![image](/data/result.png)

