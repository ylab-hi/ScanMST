//     ScanMST  Copyright (C) 2025  Yang Lab
//     This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
//     This is free software, and you are welcome to redistribute it
//     under certain conditions; type `show c' for details.
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

#ifndef SCANMSTEXT_BAM_H
#define SCANMSTEXT_BAM_H
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "htslib/sam.h"

namespace cppext {
  using uint = uint32_t;

  struct CigarResult {
    uint lt_soft_len{0};
    uint rt_soft_len{0};
    uint read_match{0};
    uint ref_match{0};
    int indel_len{0};  // may be negative
    uint query_len{0};
    std::vector<uint> cigartuples_without_soft{};
    std::vector<uint> cigartuples{};

    [[maybe_unused]] [[nodiscard]] std::string to_string() const;
    [[maybe_unused]] [[nodiscard]] bool both_soft_clipped() const {
      return lt_soft_len > 0 && rt_soft_len > 0;
    }
    [[maybe_unused]] [[nodiscard]] bool is_right_soft_clipped() const { return rt_soft_len > 0; }
    [[maybe_unused]] [[nodiscard]] bool is_left_soft_clipped() const { return lt_soft_len > 0; }
  };

  std::ostream &operator<<(std::ostream &os, CigarResult const &cigar_result);

  /**
   * @brief parse cigar string
   * @param cigar string of cigar
   * @return  parseCigarResult_t
   */
  [[maybe_unused]] CigarResult parseCigar(const char *cigar);

  template <typename T> constexpr void deleter(T *p) noexcept {
    if constexpr (std::is_same_v<T, sam_hdr_t>) {
      sam_hdr_destroy(p);
    } else if constexpr (std::is_same_v<T, bam1_t>) {
      bam_destroy1(p);
    } else if constexpr (std::is_same_v<T, samFile>) {
      sam_close(p);
    } else if constexpr (std::is_same_v<T, hts_idx_t>) {
      hts_idx_destroy(p);
    } else if constexpr (std::is_same_v<T, hts_itr_t>) {
      hts_itr_destroy(p);
    } else {
      std::cerr << "Unknown type\n";
    }
  }

  template <typename T> using hts_unique_ptr = std::unique_ptr<T, decltype(&deleter<T>)>;

  template <typename T> hts_unique_ptr<T> make_hts_unique_ptr(T *p) {
    return hts_unique_ptr<T>(p, &deleter<T>);
  }

  struct Region {
    std::string chrom;
    long start{};
    long end{};
    Region() = default;
    Region(std::string_view chrom_, long start_, long end_)
        : chrom{chrom_}, start{start_}, end{end_} {}

    [[maybe_unused]] [[nodiscard]] std::string to_string() const;
  };

  class BamReader {
  private:
    std::string file{};
    hts_unique_ptr<samFile> sam_file = make_hts_unique_ptr<samFile>(nullptr);
    hts_unique_ptr<bam_hdr_t> sam_header = make_hts_unique_ptr<sam_hdr_t>(nullptr);
    hts_unique_ptr<hts_idx_t> sam_index = make_hts_unique_ptr<hts_idx_t>(nullptr);
    hts_unique_ptr<bam1_t> sam_record = make_hts_unique_ptr<bam1_t>(bam_init1());

  public:
    struct Iterator;

    BamReader() = default;
    explicit BamReader(std::string_view bamFile);

    /** rule of 5
     * @brief: copy constructor and copy assignment
     * @return: this
     */
    BamReader(const BamReader &) = delete;
    BamReader &operator=(const BamReader &) = delete;

    /** rule of 5
     * @brief: move constructor and move assignment
     * @return: this
     */
    BamReader(BamReader &&src) noexcept = default;
    BamReader &operator=(BamReader &&src) noexcept = default;

    /**
     * @brief count reads from bam file in terms of the position
     * @param t_chrom  chromosome name
     * @param start_t  start position
     * @param end_t  end position
     * @return  number of reads
     */
    [[nodiscard]] int count(std::string_view t_chrom, long start_t, long end_t) const;

    [[maybe_unused]] [[nodiscard]] int count(const Region &region) const;

    [[nodiscard]] Iterator query(std::string_view chrom, long start, long end) const;
    [[maybe_unused]] [[nodiscard]] Iterator query(const Region &region) const;

    [[maybe_unused]] void print_chroms() const;
  };

  struct BamReader::Iterator {
    hts_unique_ptr<hts_itr_t> iter = make_hts_unique_ptr<hts_itr_t>(nullptr);
    samFile *sam_file = nullptr;
    bam1_t *sam_record = nullptr;
    bam_hdr_t *sam_header = nullptr;
    bool is_end_{false};

    Iterator() = default;

    Iterator(hts_itr_t *iter, samFile *sam_file_, bam1_t *sam_record_, bam_hdr_t *sam_header_)
        : iter(iter, &deleter<hts_itr_t>),
          sam_file(sam_file_),
          sam_record(sam_record_),
          sam_header(sam_header_) {
      if (iter == nullptr) {
        is_end_ = true;
      } else {
        next();
      }
    }

    void next();

    [[nodiscard]] bool is_end() const;
    [[nodiscard]] bool is_reverse() const;

    [[maybe_unused]] [[nodiscard]] std::string read_name() const;
    [[nodiscard]] std::string chrom() const;
    [[nodiscard]] long pos() const;
    [[nodiscard]] long end_pos() const;
    [[nodiscard]] uint8_t quality() const;
    [[nodiscard]] std::string sequence() const;
    [[nodiscard]] std::string cigar_string() const;

    [[nodiscard]] uint cigar_length() const;
    [[nodiscard]] uint *cigar_buffer() const;

    [[maybe_unused]] [[nodiscard]] std::string to_string() const;
    [[maybe_unused]] void print() const;

    [[nodiscard]] bool same_strand_with(bool is_reversed) const;
    [[nodiscard]] bool quality_eq_than(int quality_) const;
  };

}  // namespace cppext
#endif  // SCANMSTEXT_BAM_H
