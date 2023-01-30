//     scannls  Copyright (C) 2022  Yangyang Li
//     This program comes with ABSOLUTELY NO WARRANTY; for details type `show
//     w'. This is free software, and you are welcome to redistribute it under
//     certain conditions; type `show c' for details.
//
// The hypothetical commands `show w' and `show c' should show the appropriate
// parts of the General Public License.  Of course, your program's commands
// might be different; for a GUI interface, you would use an "about box".
//
//   You should also get your employer (if you work as a programmer) or school,
// if any, to sign a "copyright disclaimer" for the program, if necessary.
// For more information on this, and how to apply and follow the GNU GPL, see
// <https://www.gnu.org/licenses/>.
//
//   The GNU General Public License does not permit incorporating your program
// into proprietary programs.  If your program is a subroutine library, you
// may consider it more useful to permit linking proprietary applications with
// the library.  If this is what you want to do, use the GNU Lesser General
// Public License instead of this License.  But first, please read
// <https://www.gnu.org/licenses/why-not-lgpl.html>.

#ifndef SCANNLSEXT_RESCUER_H
#define SCANNLSEXT_RESCUER_H
#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <optional>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include "bam.h"
#include "htslib/sam.h"
#include "ssw_cpp.h"
#include "utils.hpp"

namespace cppext {
constexpr int max_align_seq_len = 100;
constexpr int max_sr_pos_diff = 10;
constexpr uint8_t align_match_score = 2;
constexpr uint8_t align_mismatch_penalty = 5;
constexpr uint8_t align_gap_open_penalty = 8;
constexpr uint8_t align_gap_extend_penalty = 6;

/**
 * @brief Parse cigar string from original uint32_t array
 * @param t_cigar_str: uint32_t array of cigar string
 * @param t_cigar_len: length of cigar string array
 * @return: parseCigarResult_t
 */
CigarResult parser_cigar(const uint32_t *t_cigar_str, size_t t_cigar_len);

/**
 * @brief Parse cigar string from record iterator
 * @param iterator
 * @return parseCigarResult_t
 */
CigarResult parser_cigar(BamReader::Iterator const &iterator);

struct Options;
struct BreakPoint;
struct Seqs;

class Rescuer {
private:
  BamReader bam_reader{};
  int min_mapq{};
  int min_soft_len{}; // may be not used
  int min_mismatch{};
  double min_identity{};
  int min_seq_align_len{};
  std::optional<int> average_read_depth{};
  StripedSmithWaterman::Aligner m_aligner{
      align_match_score, align_mismatch_penalty, align_gap_open_penalty,
      align_gap_extend_penalty};
  StripedSmithWaterman::Filter m_filter{StripedSmithWaterman::Filter{}};
  mutable StripedSmithWaterman::Alignment m_alignment{};
  mutable std::vector<std::string> m_names_list{};

public:
  explicit Rescuer(Options const &t_options);

  [[maybe_unused]] std::string to_string() const;

  /**
   * @brief calculate number of sr for every node
   * @param region info including chrom, start, end
   * @param t_mode  mode of node
   * @param current_query_name current query names
   * @param t_query_name_list query names list
   * @return number of sr
   */
  int calculate_sr(const Region &region, BreakPoint const &break_point,
                   std::vector<std::string> &current_query_name,
                   const std::vector<uint> &cigartuples_without_soft) const;

  /**
   * @brief calculate number of sr according to alignment between srlist and
   * @param candidate_list sr list
   * @param reference_list sv list
   * @return number of sr
   */
  int calculate_incremented_sr(
      const std::vector<Seqs> &candidate_list,
      const std::vector<Seqs> &reference_list,
      std::vector<std::string> const &candidate_list_names) const;

  /**
   * @brief calculate number of read for a position
   * @param t_chrom chromosome
   * @param t_start start position 0-based
   * @param t_end  end position
   * @return number of read
   */

  [[maybe_unused]] int count_reads(const std::string &t_chrom, long t_start,
                                   long t_end) const;

  /**
   * @brief check if two read sequence need to be aligned or not
   * @param query query sequence
   * @param target target sequence
   * @return true if need to be aligned
   */
  bool check_rescue(const Seqs &query, const Seqs &target) const;

  bool check_identity(std::string_view query, std::string_view target,
                      double threshold) const;

  bool check_rescue_condition(std::string_view query,
                              std::string_view target) const;

  /**
   * @param query  query sequence
   * @param target  target sequence
   * @return the identity of query
   */
  std::optional<double> calculate_identity(std::string_view query,
                                           std::string_view target) const;

  /**
   * @brief reset name list
   * @param t_names_list
   */
  void reset_names_list(std::vector<std::string> const &t_names_list) const;

  /**
   * @brief Add seqs for sr and sv list
   * @param candidate_seqs seq list of sr
   * @param reference_seqs  seq list of sv
   * @param t_bam  bam handler
   * @param tt_chrom  chrom name
   * @param tt_start  start position 0-based
   * @param tt_end  end position 0-based
   * @param tt_mode  mode
   * @param t_min_mapq  min mapq
   * @param t_min_soft  min soft clip length
   * @param current_names  current names list
   * @param t_name_list  names in graph
   */
  std::vector<std::string>
  add_align_seqs(std::vector<Seqs> &candidate_seqs,
                 std::vector<Seqs> &reference_seqs, Region const &region,
                 BreakPoint const &break_point,
                 std::vector<std::string> &current_names,
                 const std::vector<uint> &cigartuples_without_soft) const;

  std::optional<Seqs>
  get_align_sequences(BamReader::Iterator const &iterator,
                      BreakPoint const &break_point,
                      CigarResult const &cigar_result) const;

  /**
   * @brief check if the rescued read need to be skipped
   * @param read_seq read sequence of rescued read
   * @param tt_mode mode of node's read
   * @param tt_start start position of node's read for rescuing
   * @param  iterator cigar string buff of rescued read
   * @param tt_cigartuples_without_soft cigar tuples without soft clip of node's
   * read
   * @return tuple<get_softclip_result_t, seq_len, is_skip>
   */
  static std::optional<int> get_align_seq_len(uint seq_size,
                                              uint min_align_len);

  /**
   * @brief check start position for rescued reads and node
   * @param break_point node's read reference start
   * @return true if rescued reads need to be skipped
   */
  static bool check_read_pos(BamReader::Iterator const &iterator,
                             const BreakPoint &break_point);

  /**
   * @brief check if has same isoform
   * @param left_cigartuples_without_soft cigar tuples without soft clip for one
   * read
   * @param right_cigartuples_without_soft cigar tuples without soft clip for
   * another read
   * @return true if has same isoform
   */
  static bool
  check_if_same_isform(const std::vector<uint> &left_cigartuples_without_soft,
                       const std::vector<uint> &right_cigartuples_without_soft);
};

struct Options {
  using self = Options;

  Options() = default;

  self &file(std::string_view file);

  self &mapq(int mapq);

  self &soft_len(int soft_len);

  self &mismatch(int mismatch);

  self &identity(double identity);

  self &min_seq_align_len(int min_seq_align_len);

  [[maybe_unused]] self &average_read_depth(int average_read_depth);

  [[maybe_unused]] [[nodiscard]] std::string to_string() const;

  std::string file_{};
  int mapq_{};
  int soft_len_{};
  int mismatch_{};
  double identity_{};
  int min_seq_align_len_{};
  std::optional<int> average_read_depth_{};
};

struct BreakPoint {
  BreakPoint() = default;
  BreakPoint(long read_start_, long read_end_, int mode_, bool is_reverse_,
             bool is_middle_);

  [[maybe_unused]] [[nodiscard]] std::string to_string() const;

  long read_start{}; // read match start position
  long read_end{};   // read match end position
  int mode{};
  bool is_reverse{};
  bool is_middle{}; // is middle node
};

struct Seqs {
  Seqs() = default;
  explicit Seqs(std::string_view seq1_);
  Seqs(std::string_view seq1_, std::string_view seq2_);

  [[maybe_unused]] [[nodiscard]] std::string to_string() const;

  friend std::ostream &operator<<(std::ostream &os, const Seqs &seqs) {
    return os << seqs.to_string();
  }

  std::string seq1{};
  std::optional<std::string> seq2{};
};

} // namespace cppext

#endif // SCANNLSEXT_RESCUER_H
>>>>>>> c5cd6df (fix complict)
