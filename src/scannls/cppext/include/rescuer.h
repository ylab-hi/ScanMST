//
// Created by li002252 on 2/11/22.
//

#ifndef SCANNLSEXT_RESCUER_H
#define SCANNLSEXT_RESCUER_H
#include <algorithm>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include "bam.h"
#include "htslib/sam.h"
#include "ssw_cpp.h"

namespace rescuer {

  using bam_parser::bam_handler;
  using bam_parser::parseCigarResult_t;

  struct get_softclip_result_t {
    int soft_len{};
    std::string read_seq{};
    long pos{-1};
    int mode{0};
  };

  get_softclip_result_t get_softclip(const std::string &t_read_seq, const bam1_t *t_alignment,
                                     int t_mode, const uint32_t *t_cigar_str, size_t t_cigar_len);

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
  void add_sr_sv_list(std::vector<std::string> &t_sr, std::vector<std::string> &t_sv,
                      const bam_handler &t_bam, const std::string &tt_chrom, long tt_start,
                      long tt_end, int tt_mode, int t_min_mapq, int t_min_soft,
                      const std::vector<std::string> &t_current_names,
                      const std::vector<std::string> &t_name_list);
  /**
   * @brief calculate number of sr according to alignment between srlist and
   * svlist
   * @param t_sr sr list
   * @param t_sv sv list
   * @param t_min_frac min fraction of alignment
   * @param t_mode mode
   * @param t_min_mismatch min mismatch of alignment
   * @return number of sr
   */
  int determine_num_increment_sr(const std::vector<std::string> &t_sr,
                                 const std::vector<std::string> &t_sv, double t_min_frac,
                                 int t_min_mismatch);

  class Rescuer {
  private:
    const char * m_file_path{};
    bam_handler m_bam_handler{};
    int min_mapq{};
    int min_soft_len{};
    int min_mismatch{};
    double min_identity{};

  public:
    Rescuer();
    Rescuer(const char *t_file, int t_mapq, int t_soft_len, int t_mismatch, double t_identity);
    ~Rescuer();

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
    int calculate_sr(const std::string &t_chrom, long t_start, long t_end, int t_mode,
                     std::vector<std::string> &t_current_query_name,
                     std::vector<std::string> &t_query_name_list) const;
  };

}  // namespace rescuer

#endif  // SCANNLSEXT_RESCUER_H
