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

  /**
   * @brief parse cigar string
   * @param cigar string of cigar
   * @return  parseCigarResult_t
   */
  [[maybe_unused]] parseCigarResult_t parseCigar(const char *cigar);

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

  [[maybe_unused]] void print_reads(const std::vector<read_t> &reads);

  class bam_handler {
  public:
    samFile *sam_file{nullptr};
    sam_hdr_t *sam_header{nullptr};
    hts_idx_t *sam_index{nullptr};
    bam1_t *sam_record{bam_init1()};

    bam_handler() = default;
    explicit bam_handler(const char *bamFile);

    /**
     * @breif: release memory of bam_handler
     */
    ~bam_handler();

    /** rule of 5
     * @brief: copy constructor and copy assignment
     * @return: this
     */
    bam_handler(const bam_handler &) = delete;
    bam_handler &operator=(const bam_handler &) = delete;

    /** rule of 5
     * @brief: move constructor and move assignment
     * @return: this
     */
    bam_handler(bam_handler &&src) noexcept
        : sam_file{src.sam_file},
          sam_header{src.sam_header},
          sam_index{src.sam_index},
          sam_record{src.sam_record} {
      src.sam_file = nullptr;
      src.sam_header = nullptr;
      src.sam_index = nullptr;
      src.sam_record = nullptr;
    }
    bam_handler &operator=(bam_handler &&src) noexcept {
      sam_file = src.sam_file;
      sam_header = src.sam_header;
      sam_index = src.sam_index;
      sam_record = src.sam_record;
      src.sam_file = nullptr;
      src.sam_header = nullptr;
      src.sam_index = nullptr;
      src.sam_record = nullptr;
      return *this;
    }

    /**
     * @breif fetch reads from bam file in terms of chromosome and start and end
     * @param t_chrom  chromosome name
     * @param t_start  start position
     * @param t_end  end position
     * @return  vector of read_t
     */
    [[maybe_unused]] std::vector<read_t> fetch(const char *t_chrom, long t_start, long t_end);

    /**
     * @brief count reads from bam file in terms of the position
     * @param t_chrom  chromosome name
     * @param t_start  start position
     * @param t_end  end position
     * @return  number of reads
     */
    int count(const char *t_chrom, long t_start, long t_end) const;
  };

}  // namespace bam_parser
#endif  // SCANNLSEXT_BAM_H
