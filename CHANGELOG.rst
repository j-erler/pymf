Release 1.2 9 Jan 2020

* Version 1.2 adds vectorization to the CMMF and removes all non-vectorized code.


Release 1.1 4 Dez 2018

* Version 1.1 provides vastly improved performance for the MF, CMF and MMF filters that is achieved through additional vectorization techniques. Those improvemnts will be applied to the CMMF soon.
* A bug that led to an error massage when the no_cross variable had been set to True in filter_map_mmf() and filter_map_cmmf() has been fixed.
