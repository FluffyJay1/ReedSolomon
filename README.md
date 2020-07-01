# ReedSolomon 

Basic Reed-Solomon coding algorithim implementation for C++ using [this tutorial](https://en.wikiversity.org/wiki/Reed–Solomon_codes_for_coders), only meant for small pieces of data.

This was made as a programming exercise, I take no responsiblity for damage caused by using this code.

To test, run `ReedSolomonTest.cpp`.

## Basic Usage

Whatever the use, first instantiate a ReedSolomon object, constructing it with the number of bits in the field (a.k.a. set p in GF(2^p) for p < 32). 
Keep in mind that as p approaches 30 and higher, the memory taken by the program increases exponentially, hitting 30+ Gb.

For example:

`ReedSolomon rs(12); // set up RS for GF(2^12)`

Data input/output to the algorithm will be handled in arrays of the RS_WORD type, defined in ReedSolomon.h. This can be customized to your needs, but it is set to be an `unsigned long` by default.

### Encoding

Arrange data in the form of an unsigned int array, then call `ReedSolomon::encode(RS_WORD* out, RS_WORD* data, int k, int nsym)`

`out` is the output buffer for the encoded message (with length `k` + `nsym`), `data` is the message array, `k` is the length of the message, 
and `nsym` is the number of error correction symbols to use.

The first `k` elements of the output buffer will be the message, unchanged, followed by `nsym` error correction elements.

### Decoding

Take the incoming data, corrupted or not, and call `ReedSolomon::decode(RS_WORD* wholeOut, RS_WORD* out, RS_WORD* data, int k, int nsym, vector<unsigned int>* erasePos, bool debug)`

`wholeOut` is the output buffer for the whole message (error correction elements included), `out` is the output buffer for just the `k` elements of the message, 
`data` is the incoming data, `k` is the length of the message, `nsym` is the number of error correction elements, `erasePos` is a vector of known indices of corruption 
(determined by user), and `debug` is whether to print the decoding process.

The method will return `true` if the decoding was successful, `false` otherwise. Keep in mind that a successful decoding doesn't necessarily mean that 
the original data has been recovered.

If only one of `wholeOut` and `out` are needed, replace the uneeded one with `nullptr`.