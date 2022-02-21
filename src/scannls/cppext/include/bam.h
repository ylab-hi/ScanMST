//
// Created by li002252 on 2/9/22.
//

#ifndef SCANNLSEXT_BAM_H
#define SCANNLSEXT_BAM_H
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "htslib/sam.h"
namespace bam_parser {

  struct parseCigarResult_t {
    uint lt_soft_len{0};
    uint rt_soft_len{0};
    uint read_match{0};
    uint ref_match{0};
    int indel_len{0};  // may be negative
    uint query_len{0};
    std::vector<uint> cigartuples_without_soft{};
    std::vector<uint> cigartuples{};
  };

  [[maybe_unused]] void printChrome(const bam_hdr_t *har);
  [[maybe_unused]] void readBam(const char *bamFile);
  parseCigarResult_t parseCigar(const char *cigar);

  struct read_t {
    std::string query_name{};
    std::string chrom{};
    long ref_start{};
    long ref_end{};
    int mapping_quality{};
    const uint32_t *cigar{};
    uint32_t n_cigar{};
    std::string query_seq{};
    bool is_reverse{};
    read_t(std::string t_query_name, std::string t_chrom, long t_ref_start, long t_ref_end,
           int t_mapping, const uint32_t *t_cigar, uint32_t t_n_cigar, std::string t_query_seq,
           bool t_is_reverse);
  };

  void print_reads(const std::vector<read_t> &reads);

  class bam_handler {
  public:
    samFile *sam_file{nullptr};
    sam_hdr_t *sam_header{nullptr};
    hts_idx_t *sam_index{nullptr};
    bam1_t *sam_record{bam_init1()};

    bam_handler() = default;
    bam_handler(const char *bamFile);
    ~bam_handler();

    [[maybe_unused]] std::vector<read_t> fetch(const char *t_chrom, long t_start, long t_end);
  };

}  // namespace bam_parser
#endif  // SCANNLSEXT_BAM_H
