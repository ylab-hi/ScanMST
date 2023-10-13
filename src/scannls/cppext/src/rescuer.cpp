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

#include "rescuer.h"

#include <sstream>

namespace cppext {

Options::self &Options::file(std::string_view file) {
  file_ = file;
  return *this;
}
Options::self &Options::mapq(int mapq) {
  mapq_ = mapq;
  return *this;
}
Options::self &Options::soft_len(int soft_len) {
  soft_len_ = soft_len;
  return *this;
}
Options::self &Options::mismatch(int mismatch) {
  mismatch_ = mismatch;
  return *this;
}
Options::self &Options::identity(double identity) {
  identity_ = identity;
  return *this;
}
Options::self &Options::min_seq_align_len(int min_seq_align_len) {
  min_seq_align_len_ = min_seq_align_len;
  return *this;
}
[[maybe_unused]] Options::self &Options::average_read_depth(
    int average_read_depth) {
  average_read_depth_ = average_read_depth;
  return *this;
}

[[maybe_unused]] std::string Options::to_string() const {
  std::stringstream ss{};
  ss << "Options("
     << "\n";
  ss << "file: " << file_ << "\n";
  ss << "mapq: " << mapq_ << "\n";
  ss << "soft_len: " << soft_len_ << "\n";
  ss << "mismatch: " << mismatch_ << "\n";
  ss << "identity: " << identity_ << "\n";
  ss << "min_seq_align_len: " << min_seq_align_len_ << "\n";
  ss << "average_read_depth: " << average_read_depth_.value_or(-1) << "\n";
  ss << ")"
     << "\n";
  return ss.str();
}

Rescuer::Rescuer(const Options &options)
    : bam_reader(options.file_),
      min_mapq(options.mapq_),
      min_soft_len(options.soft_len_),
      min_mismatch(options.mismatch_),
      min_identity(options.identity_),
      min_seq_align_len(options.min_seq_align_len_),
      average_read_depth(options.average_read_depth_) {}

int Rescuer::calculate_sr(
    const Region &region, BreakPoint const &break_point,
    std::vector<std::string> &current_query_name,
    const std::vector<uint> &cigartuples_without_soft) const {
  std::vector<Seqs> candidate_list{};
  std::vector<Seqs> reference_list{};

  auto sr_list_names =
      add_align_seqs(candidate_list, reference_list, region, break_point,
                     current_query_name, cigartuples_without_soft);

#ifdef SCDEBUG
  std::cout << "candidate list" << '\n';
  for (auto const &i : candidate_list) {
    std::cout << i << '\n';
  }

  std::cout << "reference list" << '\n';
  for (auto const &i : reference_list) {
    std::cout << i << '\n';
  }
#endif

  if (candidate_list.empty() || reference_list.empty()) {
    return 0;
  }

  return calculate_incremented_sr(candidate_list, reference_list,
                                  sr_list_names);
}

int Rescuer::calculate_incremented_sr(
    const std::vector<Seqs> &candidate_list,
    const std::vector<Seqs> &reference_list,
    std::vector<std::string> const &candidate_list_names) const {
  int num_increment_sr{0};
  int sr_index{-1};

  for (const auto &candidate : candidate_list) {
    ++sr_index;
    bool is_sr_increment{false};

    for (const auto &reference : reference_list) {
      if (!check_rescue(candidate, reference)) {
        continue;
      }

      // candidate can be rescued
      ++num_increment_sr;
      is_sr_increment = true;
      break;
    }
    //       add query name of sr to all names list if the query name of sr
    //       update successfully to ensure that the query name of sr  do not
    //       update repeatedly
    if (is_sr_increment) {
      m_names_list.push_back(candidate_list_names[sr_index]);
      // -1 means no limit when running in loose mode
      if (average_read_depth.has_value() &&
          num_increment_sr >= average_read_depth) {
        std::cout << "early stopping " << '\n';
        return num_increment_sr;
      }
    }
  }

  return num_increment_sr;
}

[[maybe_unused]] int Rescuer::count_reads(const std::string &t_chrom,
                                          long t_start, long t_end) const {
  return bam_reader.count(t_chrom.c_str(), t_start, t_end);
}

std::optional<double> Rescuer::calculate_identity(
    std::string_view query, std::string_view target) const {
  auto const target_len{static_cast<int>(target.length())};
  int const masklen{std::max(target_len / 2, 15)};

  bool return_value{m_aligner.Align(query.data(), target.data(), target_len,
                                    m_filter, &m_alignment, masklen)};

  if (double identity{static_cast<double>(m_alignment.query_end -
                                          m_alignment.query_begin + 1) /
                      static_cast<double>(query.length())};
      return_value) {
    return identity;  // do not align
  } else {
    return {};
  }
}

bool Rescuer::check_identity(std::string_view query, std::string_view target,
                             double threshold) const {
  if (auto const identity{calculate_identity(query, target)};
      identity.has_value() && identity.value() >= threshold) {
    return true;
  }

  return false;
}

bool Rescuer::check_rescue(const Seqs &query, const Seqs &target) const {
  bool flag_1{check_identity(query.seq1.substr(0, min_seq_align_len),
                             target.seq1.substr(0, min_seq_align_len),
                             min_identity)};

  if (query.seq2.has_value() && target.seq2.has_value()) {
    if (bool flag_2{check_identity(
            query.seq2.value().substr(0, min_seq_align_len),
            target.seq2.value().substr(0, min_seq_align_len), min_identity)};
        flag_1 && flag_2) {
      // middle node
      bool rescue_flag_1{check_rescue_condition(query.seq1, target.seq1)};
      bool rescue_flag_2{
          check_rescue_condition(query.seq2.value(), target.seq2.value())};
      return rescue_flag_1 && rescue_flag_2;
    }

    return false;
  }

  if (flag_1) {
    // end node
    return check_rescue_condition(query.seq1, target.seq1);
  }

  return false;
}

bool Rescuer::check_rescue_condition(std::string_view query,
                                     std::string_view target) const {
  if (auto const identity{calculate_identity(query, target)};
      identity.has_value() && identity.value() >= min_identity &&
      m_alignment.mismatches <= min_mismatch &&
      m_alignment.query_begin + m_alignment.ref_begin <= 2) {
    return true;
  }
  return false;
}

void Rescuer::reset_names_list(
    std::vector<std::string> const &t_names_list) const {
  m_names_list = t_names_list;
}

bool Rescuer::check_read_pos(BamReader::Iterator const &iterator,
                             const BreakPoint &break_point) {
  if (break_point.is_middle) {
    return std::abs(iterator.pos() - break_point.read_start) <=
               max_sr_pos_diff &&
           std::abs(iterator.end_pos() - break_point.read_end) <=
               max_sr_pos_diff;
  }

  if (break_point.mode == 1) {
    // right soft clipped
    return std::abs(iterator.pos() - break_point.read_start) <= max_sr_pos_diff;
  }

  if (break_point.mode == 2) {
    // left soft clipped
    return std::abs(iterator.end_pos() - break_point.read_end) <=
           max_sr_pos_diff;
  }
  return false;
}

std::vector<std::string> Rescuer::add_align_seqs(
    std::vector<Seqs> &candidate_seqs, std::vector<Seqs> &reference_seqs,
    Region const &region, BreakPoint const &break_point,
    std::vector<std::string> &current_names,
    const std::vector<uint> &cigartuples_without_soft) const {
  std::vector<std::string> candidate_list_names{};

  for (auto iterator = bam_reader.query(region);
       !iterator.is_end() && iterator.same_strand_with(break_point.is_reverse);
       iterator.next()) {
    //    get read seq

    if (iterator.quality_eq_than(min_mapq)) {
      CigarResult cigar_result{parser_cigar(iterator)};

      if (!check_read_pos(iterator, break_point) ||
          !check_if_same_isform(cigartuples_without_soft,
                                cigar_result.cigartuples_without_soft)) {
        continue;
      }

      auto align_sequence =
          get_align_sequences(iterator, break_point, cigar_result);
      if (!align_sequence.has_value()) {
        continue;
      }

#ifdef SCDEBUG
      std::cout << '\n' << "read name: " << iterator.read_name() << '\n';
      std::cout << "current names" << '\n';
      for (auto const &current_name : current_names) {
        std::cout << current_name << '\n';
      }

      std::cout << "m_names_list" << '\n';
      for (auto const &m_names_list_name : m_names_list) {
        std::cout << m_names_list_name << '\n';
      }
#endif

      if (auto read_name{iterator.read_name()};
          find(current_names.begin(), current_names.end(), read_name) !=
          current_names.end()) {
        reference_seqs.push_back(std::move(align_sequence.value()));
      } else if (find(m_names_list.begin(), m_names_list.end(), read_name) ==
                 m_names_list.end()) {
        candidate_list_names.push_back(std::move(read_name));
        candidate_seqs.push_back(std::move(align_sequence.value()));
      }
    }
  }
  return candidate_list_names;
}

std::optional<int> Rescuer::get_align_seq_len(uint seq_size,
                                              uint min_align_len) {
  if (seq_size >= max_align_seq_len) {
    return max_align_seq_len;
  }

  if (seq_size < min_align_len) {
    return {};  // read length is too short skip
  }

  return {seq_size};  // has value not skip
}

// may need to change
std::optional<Seqs> Rescuer::get_align_sequences(
    BamReader::Iterator const &iterator, BreakPoint const &break_point,
    CigarResult const &cigar_result) const {
  Seqs seqs{};

  std::string const read_seq = iterator.sequence();

  if (break_point.is_middle) {
    auto const lt_seq_len = get_align_seq_len(
        cigar_result.lt_soft_len, static_cast<uint>(min_seq_align_len));
    auto const rt_seq_len = get_align_seq_len(
        cigar_result.rt_soft_len, static_cast<uint>(min_seq_align_len));

    // check if seq length is too short
    if (!lt_seq_len.has_value() || !rt_seq_len.has_value()) {
      return {};
    }

    auto &&temp = read_seq.substr(cigar_result.lt_soft_len - lt_seq_len.value(),
                                  lt_seq_len.value());
    seqs.seq1.assign(temp.rbegin(), temp.rend());

    seqs.seq2 = read_seq.substr(
        cigar_result.query_len - cigar_result.rt_soft_len, rt_seq_len.value());

#ifdef SCDEBUG
    std::cout << '\n' << "get seqs:" << '\n';
    std::cout << "read name: " << iterator.read_name() << '\n';
    std::cout << "read seq len: " << read_seq.length() << '\n';
    std::cout << "left seq len: " << lt_seq_len.value() << '\n';
    std::cout << "right seq len: " << rt_seq_len.value() << '\n';
    std::cout << "read seq len: " << read_seq.length() << '\n';
    std::cout << "read seq: " << read_seq << '\n';
    std::cout << "seq: " << seqs.to_string() << '\n';
#endif

    return seqs;
  }

  if (break_point.mode == 1) {
    // right soft clipped
    auto const rt_seq_len = get_align_seq_len(
        cigar_result.rt_soft_len, static_cast<uint>(min_seq_align_len));

    if (!rt_seq_len.has_value()) {
      return {};
    }
    seqs.seq1.assign(read_seq.substr(
        cigar_result.query_len - cigar_result.rt_soft_len, rt_seq_len.value()));

    return seqs;
  }

  if (break_point.mode == 2) {
    // left soft clipped
    auto const lt_seq_len = get_align_seq_len(
        cigar_result.lt_soft_len, static_cast<uint>(min_seq_align_len));
    if (!lt_seq_len.has_value()) {
      return {};
    }

    auto &&temp = read_seq.substr(cigar_result.lt_soft_len - lt_seq_len.value(),
                                  lt_seq_len.value());
    seqs.seq1.assign(temp.rbegin(), temp.rend());

    return seqs;
  }

  return {};
}

bool Rescuer::check_if_same_isform(
    const std::vector<uint> &left_cigartuples_without_soft,
    const std::vector<uint> &right_cigartuples_without_soft) {
  std::vector<uint> left_result{}, right_result{};

  auto const left_size = left_cigartuples_without_soft.size();
  auto const right_size = right_cigartuples_without_soft.size();

  for (std::vector<int>::size_type i = 0; i < left_size; i += 2) {
    if (left_cigartuples_without_soft[i] == BAM_CREF_SKIP) {
      left_result.push_back(left_cigartuples_without_soft[i + 1]);
    }
  }

  for (std::vector<int>::size_type i = 0; i < right_size; i += 2) {
    if (right_cigartuples_without_soft[i] == BAM_CREF_SKIP) {
      right_result.push_back(right_cigartuples_without_soft[i + 1]);
    }
  }

  if (left_result.size() != right_result.size()) {
    return false;
  }
  if (left_result.empty()) {
    return true;
  }

  return std::equal(left_result.begin(), left_result.end(),
                    right_result.begin());
}

[[maybe_unused]] std::string Rescuer::to_string() const {
  std::stringstream ss{};
  ss << "Rescuer("
     << "\n";
  ss << "  min_mapq: " << min_mapq << "\n";
  ss << "  min_soft_len: " << min_soft_len << "\n";
  ss << "  min_mismatch: " << min_mismatch << "\n";
  ss << "  min_identity: " << min_identity << "\n";
  ss << "  min_seq_align_len: " << min_seq_align_len << "\n";
  ss << "  average_read_depth: " << average_read_depth.value_or(-1) << "\n";
  ss << ")"
     << "\n";
  return ss.str();
}

CigarResult parser_cigar(const uint32_t *t_cigar_str, size_t t_cigar_len) {
  CigarResult result{};

  result.cigartuples.reserve(2 * t_cigar_len);
  result.cigartuples_without_soft.reserve(2 * t_cigar_len);

  for (size_t i{0}; i < t_cigar_len; i++) {
    uint op{bam_cigar_op(t_cigar_str[i])};
    uint32_t len{bam_cigar_oplen(t_cigar_str[i])};
    result.cigartuples.insert(result.cigartuples.end(), {op, len});

    switch (op) {
      case BAM_CMATCH:
        result.ref_match += len;
        result.read_match += len;
        result.query_len += len;
        result.cigartuples_without_soft.insert(
            result.cigartuples_without_soft.end(), {op, len});
        break;
      case BAM_CINS:
        result.indel_len -= static_cast<int>(len);
        result.read_match += len;
        result.query_len += len;
        result.cigartuples_without_soft.insert(
            result.cigartuples_without_soft.end(), {op, len});
        break;
      case BAM_CDEL:
      case BAM_CREF_SKIP:
        result.indel_len += static_cast<int>(len);
        result.ref_match += len;
        result.cigartuples_without_soft.insert(
            result.cigartuples_without_soft.end(), {op, len});
        break;
      case BAM_CSOFT_CLIP:
        result.query_len += len;
        break;
    }
  }

  if (result.cigartuples[0] == BAM_CSOFT_CLIP) {
    result.lt_soft_len = result.cigartuples[1];
  }
  if (result.cigartuples[2 * t_cigar_len - 2] == BAM_CSOFT_CLIP) {
    result.rt_soft_len = result.cigartuples[2 * t_cigar_len - 1];
  }
  return result;
}

CigarResult parser_cigar(BamReader::Iterator const &iterator) {
  return parser_cigar(iterator.cigar_buffer(), iterator.cigar_length());
}

BreakPoint::BreakPoint(long read_start_, long read_end_, int mode_,
                       bool is_reverse_, bool is_middle_)
    : read_start{read_start_},
      read_end{read_end_},
      mode{mode_},
      is_reverse{is_reverse_},
      is_middle{is_middle_} {}

[[maybe_unused]] std::string BreakPoint::to_string() const {
  std::stringstream ss{};

  ss << "BreakPoint("
     << "\n"
     << "  read_start: " << read_start << "\n"
     << "  mode: " << mode << "\n"
     << "  is_reverse: " << is_reverse << "\n"
     << "  is_middle: " << is_middle << "\n"
     << ")"
     << "\n";

  return ss.str();
}

Seqs::Seqs(std::string_view seq1_, std::string_view seq2_)
    : seq1{seq1_}, seq2{seq2_} {}
Seqs::Seqs(std::string_view seq1_) : seq1{seq1_} {}
[[maybe_unused]] std::string Seqs::to_string() const {
  std::stringstream ss{};

  ss << "Seqs("
     << "\n"
     << "  seq1: " << seq1 << "\n"
     << "  seq2: " << seq2.value_or(" ") << "\n"
     << ")"
     << "\n";

  return ss.str();
}
}  // namespace cppext
