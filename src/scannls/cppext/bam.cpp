//
// Created by li002252 on 2/9/22.
//
#include <iostream>
#include <vector>
#include "bam.h"
#include <pybind11/pybind11.h>



void printChrome(const bam_hdr_t *har) {
    for (int i = 0; i < har->n_targets; i++) {
        std::cout << har->target_name[i] << "\n";
    }
}

void readBam(const char *bamFile) {
//    read bam file
    samFile *sam = sam_open(bamFile, "r");
//    header
    sam_hdr_t *header = sam_hdr_read(sam);
//    one alignment
    bam1_t *b = bam_init1();

    while (sam_read1(sam, header, b) >= 0) {
        const uint32_t *cigar = bam_get_cigar(b);
        const auto cigarLen = bam_cigar2rlen(b->core.n_cigar, cigar);
        std::cout << "cigarLen: " << cigarLen << "\n";

        for (unsigned int i{0}; i < b->core.n_cigar; i++) {
            const auto op{bam_cigar_opchr(cigar[i])};
            const auto len{bam_cigar_oplen(cigar[i])};
            std::cout << op << " " << len << "\n";
        }
    }
    sam_close(sam);
    bam_destroy1(b);
    sam_hdr_destroy(header);
}


//
//#define BAM_CMATCH      0
//#define BAM_CINS        1
//#define BAM_CDEL        2
//#define BAM_CREF_SKIP   3
//#define BAM_CSOFT_CLIP  4
//#define BAM_CHARD_CLIP  5
//#define BAM_CPAD        6
//#define BAM_CEQUAL      7
//#define BAM_CDIFF       8
//#define BAM_CBACK       9
//
//#define BAM_CIGAR_STR   "MIDNSHP=XB"
parseCigarResult_t parseCigar(const char *cigar) {
    parseCigarResult_t result{};
    uint32_t *buf{nullptr};
    size_t m{0};
    if (sam_parse_cigar(cigar, nullptr, &buf, &m) == -1) {
        std::cerr << "Error: Cannot parse Cigar  " << cigar << "\n";
    }

    for (size_t i{0}; i < m; i++) {
        uint op{bam_cigar_op(buf[i])};
        uint32_t len{bam_cigar_oplen(buf[i])};
        result.cigartuples.insert(result.cigartuples.end(), {op, len});
        switch (op) {
            case 0:
                result.ref_match += len;
                result.read_match += len;
                result.query_len += len;
                result.cigartuples_without_soft.insert(result.cigartuples_without_soft.end(), {op, len});
                break;
            case 1:
                result.indel_len -= len;
                result.read_match += len;
                result.query_len += len;
                result.cigartuples_without_soft.insert(result.cigartuples_without_soft.end(), {op, len});
                break;
            case 2:
            case 3:
                result.indel_len += len;
                result.ref_match += len;
                break;
            case 4:
                result.query_len += len;
        }
    }
    if (result.cigartuples[0] == 4)
        result.lt_soft_len = result.cigartuples[1];
    if (result.cigartuples[2 * m - 2] == 4)
        result.rt_soft_len = result.cigartuples[2 * m - 1];
    return result;
}
