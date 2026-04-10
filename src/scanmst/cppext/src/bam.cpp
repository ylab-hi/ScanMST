//     ScanMST  Copyright (C) 2025  YangLab
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

#include "bam.h"

#include <sstream>

#include "utils.hpp"

namespace cppext {

//
// #define BAM_CMATCH      0
// #define BAM_CINS        1
// #define BAM_CDEL        2
// #define BAM_CREF_SKIP   3
// #define BAM_CSOFT_CLIP  4
// #define BAM_CHARD_CLIP  5
// #define BAM_CPAD        6
// #define BAM_CEQUAL      7
// #define BAM_CDIFF       8
// #define BAM_CBACK       9
//
// #define BAM_CIGAR_STR   "MIDNSHP=XB"
[[maybe_unused]] CigarResult parseCigar(const char *cigar) {
  CigarResult result{};
  uint32_t *buf{nullptr};
  size_t m{0};
  if (sam_parse_cigar(cigar, nullptr, &buf, &m) == -1) {
    std::cerr << "Error: Cannot parse Cigar  " << cigar << "\n";
  }

  result.cigartuples.reserve(2 * m);
  result.cigartuples_without_soft.reserve(2 * m);

  for (size_t i{0}; i < m; i++) {
    uint op{bam_cigar_op(buf[i])};
    uint32_t len{bam_cigar_oplen(buf[i])};
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
      result.indel_len -= len;
      result.read_match += len;
      result.query_len += len;
      result.cigartuples_without_soft.insert(
          result.cigartuples_without_soft.end(), {op, len});
      break;
    case BAM_CDEL:
    case BAM_CREF_SKIP:
      result.indel_len += len;
      result.ref_match += len;
      result.cigartuples_without_soft.insert(
          result.cigartuples_without_soft.end(), {op, len});
      break;
    case BAM_CSOFT_CLIP:
      result.query_len += len;
      break;
    default:;
    }
  }
  if (result.cigartuples[0] == 4)
    result.lt_soft_len = result.cigartuples[1];
  if (result.cigartuples[2 * m - 2] == 4)
    result.rt_soft_len = result.cigartuples[2 * m - 1];
  free(buf);
  return result;
}

std::ostream &operator<<(std::ostream &os, const CigarResult &cigar_result) {
  return os << cigar_result.to_string();
}

BamReader::BamReader(std::string_view file_path)
    : file(file_path),
      sam_file(sam_open(file_path.data(), "r"), &deleter<samFile>),
      sam_header(sam_hdr_read(sam_file.get()), &deleter<sam_hdr_t>),
      sam_index(sam_index_load(sam_file.get(), file_path.data()),
                &deleter<hts_idx_t>) {}

BamReader::Iterator BamReader::query(std::string_view chrom, long start,
                                     long end) const {
  const int tid = bam_name2id(sam_header.get(), chrom.data());

  return {sam_itr_queryi(sam_index.get(), tid, start - 1, end), sam_file.get(),
          sam_record.get(), sam_header.get()};
}

[[maybe_unused]] BamReader::Iterator
BamReader::query(const Region &region) const {
  return BamReader::query(region.chrom, region.start, region.end);
}

int BamReader::count(std::string_view t_chrom, long start_t, long end_t) const {
  int num_reads{0};
  auto iterator = query(t_chrom, start_t, end_t);

  while (!iterator.is_end()) {
    ++num_reads;
    iterator.next();
  }
  return num_reads;
}

[[maybe_unused]] int BamReader::count(const Region &region) const {
  return count(region.chrom, region.start, region.end);
}

[[maybe_unused]] void BamReader::print_chroms() const {
  if (sam_header != nullptr) {
    for (int i = 0; i < sam_header->n_targets; i++) {
      std::cout << sam_header->target_name[i] << "\n";
    }
  }
}

bool BamReader::Iterator::is_end() const { return is_end_; }
bool BamReader::Iterator::is_reverse() const { return bam_is_rev(sam_record); }

long BamReader::Iterator::pos() const { return sam_record->core.pos; }
uint8_t BamReader::Iterator::quality() const { return sam_record->core.qual; }
uint BamReader::Iterator::cigar_length() const {
  return sam_record->core.n_cigar;
}
uint *BamReader::Iterator::cigar_buffer() const {
  return bam_get_cigar(sam_record);
}
long BamReader::Iterator::end_pos() const { return bam_endpos(sam_record); }
std::string BamReader::Iterator::read_name() const {
  return bam_get_qname(sam_record);
}
std::string BamReader::Iterator::chrom() const {
  return sam_header->target_name[sam_record->core.tid];
}

bool BamReader::Iterator::same_strand_with(bool is_reversed) const {
  return is_reversed == is_reverse();
}

bool BamReader::Iterator::quality_eq_than(int quality_) const {
  return quality() >= quality_;
}

std::string BamReader::Iterator::cigar_string() const {
  std::string cigar_string_;

  const uint32_t *cigar = bam_get_cigar(sam_record);
  auto n_cigar = sam_record->core.n_cigar;

  for (size_t i{0}; i < n_cigar; ++i) {
    const auto op{bam_cigar_opchr(cigar[i])};
    const auto len{bam_cigar_oplen(cigar[i])};
    cigar_string_ += std::to_string(len) + op;
  }

  return cigar_string_;
}

[[maybe_unused]] std::string BamReader::Iterator::to_string() const {
  std::stringstream ss{};
  ss << "chrom " << chrom() << " read name: " << read_name()
     << " is reverse :" << is_reverse() << " cigar:" << cigar_string()
     << " pos:" << pos() << " end pos:" << end_pos() << "\n";
  return ss.str();
}

[[maybe_unused]] void BamReader::Iterator::print() const {
  std::cout << to_string();
}

std::string BamReader::Iterator::sequence() const {
  std::string read_sequence_{};

  const uint8_t *seq = bam_get_seq(sam_record);
  const int l_seq = sam_record->core.l_qseq;
  read_sequence_.resize(l_seq);
  for (int i = 0; i < l_seq; ++i) {
    read_sequence_[i] = seq_nt16_str[bam_seqi(seq, i)];
  }
  return read_sequence_;
}

void BamReader::Iterator::next() {
  if (auto res = sam_itr_next(sam_file, iter.get(), sam_record); res < 0) {
    is_end_ = true;
  }
}

[[maybe_unused]] std::string CigarResult::to_string() const {
  std::stringstream ss{};
  ss << "parseCigarResult_t("
     << "\n";
  ss << "lt_soft_len: " << lt_soft_len << "\n";
  ss << "rt_soft_len: " << rt_soft_len << "\n";
  ss << "read_match: " << read_match << "\n";
  ss << "ref_match: " << ref_match << "\n";
  ss << "indel_len: " << indel_len << "\n";
  ss << "query_len: " << query_len << "\n";
  ss << "cigartuples_without_soft: ";
  for (auto const &c : cigartuples_without_soft) {
    ss << c << " ";
  }
  ss << "\n";
  ss << "cigartuples: ";
  for (auto const &c : cigartuples) {
    ss << c << " ";
  }
  ss << ")"
     << "\n";
  return ss.str();
}

[[maybe_unused]] std::string Region::to_string() const {
  std::stringstream ss{};
  ss << "Region(" << chrom << ":" << start << "-" << end << ")"
     << "\n";
  return ss.str();
}
} // namespace cppext
