//
// Created by li002252 on 2/9/22.
//

#ifndef SCANNLSEXT_BAM_H
#define SCANNLSEXT_BAM_H
#include "htslib/sam.h"
#include <vector>
#include <utility>
#include <string>

struct parseCigarResult_t {
    int lt_soft_len{0};
    int rt_soft_len{0};
    int read_match{0};
    int ref_match{0};
    int indel_len{0};
    int query_len{0};
    std::vector<uint> cigartuples_without_soft{};
    std::vector<uint> cigartuples{};
};

struct soft_clipResult_t{
    int lt_soft_len{0};
    std::string seq{};
    int ref_pos{0};
    int mode{0};
};

void printChrome(const bam_hdr_t *har);
void readBam(const char *bamFile);
parseCigarResult_t parseCigar(const char *cigar);

#endif //SCANNLSEXT_BAM_H
