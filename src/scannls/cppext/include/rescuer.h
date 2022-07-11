//
// Created by li002252 on 2/11/22.
//

#ifndef SCANNLSEXT_RESCUER_H
#define SCANNLSEXT_RESCUER_H
#include <algorithm>
#include <iostream>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include "bam.h"
#include "htslib/sam.h"
#include "output_container.h"
#include "ssw_cpp.h"

namespace rescuer {

  using bam_parser::bam_handler;
  using bam_parser::parseCigarResult_t;
  constexpr int max_seq_len = 100;
  constexpr int max_sr_pos_diff = 10;
  constexpr int min_seq_align_len = 10;
  constexpr uint align_match_score = 2;
  constexpr uint align_mismatch_penalty = 5;
  constexpr uint align_gap_open_penalty = 8;
  constexpr uint align_gap_extend_penalty = 6;

  struct get_softclip_result_t {
    int soft_len{};
    std::string read_seq{};
    long pos{-1};
    int mode{0};
  };

  inline int get_read_max_length(std::string_view t_seq);

  /**
   * @brief Check if cigar string is soft-clipped
   * @param t_cigar cigar string
   * @param t_cigar_len cigar string length
   * @return true if soft-clipped
   */
  bool is_soft_clipped(const uint32_t *t_cigar, size_t t_cigar_len);

  /**
   * @brief Get read sequence
   * @param t_seqs : the read sequence
   * @param t_seq_len : the length of the read sequence
   * @return : read sequence
   */
  std::string get_read_seq(const uint8_t *t_seqs, int t_seq_len);

  /**
   * @brief Parse cigar string from original uint32_t array
   * @param t_cigar_str: uint32_t array of cigar string
   * @param t_cigar_len: length of cigar string array
   * @return: parseCigarResult_t
   */
  parseCigarResult_t parser_cigar(const uint32_t *t_cigar_str, size_t t_cigar_len);

  class Rescuer {
  private:
    const char *m_file_path{};
    bam_handler m_bam_handler{};
    int min_mapq{};
    int min_soft_len{};
    int min_mismatch{};
    double min_identity{};
    int min_seq_align_len{min_seq_align_len};
    int average_read_depth{};
    StripedSmithWaterman::Aligner m_aligner{
        StripedSmithWaterman::Aligner{align_match_score, align_mismatch_penalty,
                                      align_gap_open_penalty, align_gap_extend_penalty}};
    StripedSmithWaterman::Filter m_filter{StripedSmithWaterman::Filter{}};
    mutable StripedSmithWaterman::Alignment m_alignment{};
    mutable std::vector<std::string> m_names_list{};

  public:
    Rescuer() = default;
    Rescuer(const char *t_file, int t_mapq, int t_soft_len, int t_mismatch, double t_identity,
            int t_min_seq_align_len, int t_average_read_depth);

    /**
     * @brief calculate number of sr for every node
     * @param t_chrom chromosome
     * @param t_start start position 0-based
     * @param t_end   end position
     * @param t_mode  mode of node
     * @param t_current_query_name current query names
     * @param t_query_name_list query names list
     * @return number of sr
     */
    int calculate_sr(std::string_view t_chrom, long t_start, long t_end, int t_mode,
                     std::string_view t_strand, long t_read_start,
                     std::vector<std::string> &t_current_query_name,
                     const std::vector<uint> &cigartuples_without_soft) const;

    /**
     * @brief calculate number of sr according to alignment between srlist and
     * svlist
     * @param t_sr sr list
     * @param t_sv sv list
     * @return number of sr
     */
    int determine_num_increment_sr(std::vector<std::string> const &t_sr,
                                   std::vector<std::string> const &t_sv,
                                   std::vector<std::string> const &t_sr_list_names) const;

    /**
     * @brief calculate number of read for a position
     * @param t_chrom chromosome
     * @param t_start start position 0-based
     * @param t_end  end position
     * @return number of read
     */

    [[maybe_unused]] int count_reads(const std::string &t_chrom, long t_start, long t_end) const;

    /**
     * @brief check if two read sequence need to be aligned or not
     * @param t_query query sequence
     * @param t_target target sequence
     * @return true if need to be aligned
     */
    bool check_if_align(std::string_view t_query, std::string_view t_target) const;

    /**
     * @brief reset name list
     * @param t_names_list
     */
    void reset_names_list(std::vector<std::string> const &t_names_list) const;

    /**
     * @brief Add seqs for sr and sv list
     * @param t_sr seq list of sr
     * @param t_sv  seq list of sv
     * @param t_bam  bam handler
     * @param tt_chrom  chrom name
     * @param tt_start  start position 0-based
     * @param tt_end  end position 0-based
     * @param tt_mode  mode
     * @param t_min_mapq  min mapq
     * @param t_min_soft  min soft clip length
     * @param t_current_names  current names list
     * @param t_name_list  names in graph
     */
    std::vector<std::string> add_sr_sv_list(
        std::vector<std::string> &t_sr, std::vector<std::string> &t_sv, std::string_view tt_chrom,
        long tt_start, const long tt_end, const int tt_mode, std::string_view tt_strand,
        long read_start, std::vector<std::string> &t_current_names,
        const std::vector<uint> &tt_cigartuples_without_soft) const;

    std::tuple<get_softclip_result_t, bool> get_softclip(
        const std::string &t_read_seq, int t_mode, const uint32_t *t_cigar_buffer,
        const std::vector<uint> &cigartuples_without_soft) const;

    /**
     * @brief check if the rescued read need to be skipped
     * @param read_seq read sequence of rescued read
     * @param tt_mode mode of node's read
     * @param tt_start start position of node's read for rescuing
     * @param cigar cigar string buff of rescued read
     * @param tt_cigartuples_without_soft cigar tuples without soft clip of node's read
     * @return tuple<get_softclip_result_t, seq_len, is_skip>
     */
    std::tuple<get_softclip_result_t, int, bool> check_if_skip(
        std::string const &read_seq, int tt_mode, long tt_start, uint32_t const *cigar,
        const std::vector<uint> &tt_cigartuples_without_soft) const;

    /**
     * @brief check if rescued reads and node's read have same strand
     * @param strand  strand of node's read
     * @return true if not same strand
     */
    bool check_strand_if_skip(std::string_view strand) const;

    /**
     * @brief check start position for rescued reads and node
     * @param read_start node's read reference start
     * @return true if rescued reads need to be skipped
     */
    bool check_read_start_if_skip(long read_start) const;

    /**
     * @brief check if has same isoform
     * @param left_cigartuples_without_soft cigar tuples without soft clip for one read
     * @param right_cigartuples_without_soft cigar tuples without soft clip for another read
     * @return true if has same isoform
     */
    bool check_if_same_isform(const std::vector<uint> &left_cigartuples_without_soft,
                              const std::vector<uint> &right_cigartuples_without_soft) const;
  };

}  // namespace rescuer

#endif  // SCANNLSEXT_RESCUER_H
