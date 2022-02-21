//
// Created by li002252 on 2/9/22.
//
#include "bam.h"

namespace bam_parser {

  read_t::read_t(std::string t_query_name, std::string t_chrom, long t_ref_start, long t_ref_end,
                 int t_mapping, const uint32_t *t_cigar, uint32_t t_n_cigar,
                 std::string t_query_seq, bool t_is_reverse)
      : query_name{std::move(t_query_name)},
        chrom{std::move(t_chrom)},
        ref_start{t_ref_start},
        ref_end{t_ref_end},
        mapping_quality{t_mapping},
        cigar{t_cigar},
        n_cigar{t_n_cigar},
        query_seq{std::move(t_query_seq)},
        is_reverse{t_is_reverse} {}

  [[maybe_unused]] void printChrome(const bam_hdr_t *har) {
    for (int i = 0; i < har->n_targets; i++) {
      std::cout << har->target_name[i] << "\n";
    }
  }

  [[maybe_unused]] void readBam(const char *bamFile) {
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

  void print_reads(const std::vector<read_t> &reads) {
    for (auto &read : reads) {
      std::cout << read.query_name << "\t" << read.query_seq << "\t" << read.query_seq.size()
                << "\t" << read.ref_start << "\t" << read.ref_end << "\t" << read.chrom << "\t"
                << read.mapping_quality << "\t" << read.is_reverse << "\n";
    }
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
          result.cigartuples_without_soft.insert(result.cigartuples_without_soft.end(), {op, len});
          break;
        case BAM_CINS:
          result.indel_len -= len;
          result.read_match += len;
          result.query_len += len;
          result.cigartuples_without_soft.insert(result.cigartuples_without_soft.end(), {op, len});
          break;
        case BAM_CDEL:
        case BAM_CREF_SKIP:
          result.indel_len += len;
          result.ref_match += len;
          result.cigartuples_without_soft.insert(result.cigartuples_without_soft.end(), {op, len});
          break;
        case BAM_CSOFT_CLIP:
          result.query_len += len;
          break;
        default:;
      }
    }
    if (result.cigartuples[0] == 4) result.lt_soft_len = result.cigartuples[1];
    if (result.cigartuples[2 * m - 2] == 4) result.rt_soft_len = result.cigartuples[2 * m - 1];
    return result;
  }

  bam_handler::bam_handler(const char *file_path) {
    sam_file = sam_open(file_path, "r");
    if (sam_file == nullptr) {
      std::cerr << "Error: Cannot open file " << file_path << "\n";
      exit(EXIT_FAILURE);
    }

    sam_index = sam_index_load(sam_file, file_path);
    if (sam_index == nullptr) {
      std::cerr << "Error: Cannot open index file " << file_path << "\n";
      exit(EXIT_FAILURE);
    }

    sam_header = sam_hdr_read(sam_file);
    if (sam_header == 0) {
      std::cerr << "Error: Cannot read header " << file_path << "\n";
      exit(EXIT_FAILURE);
    }
  }

  bam_handler::~bam_handler() {
    hts_idx_destroy(sam_index);
    sam_close(sam_file);
    sam_hdr_destroy(sam_header);
    bam_destroy1(sam_record);
  }

  // may be more efficient to use bam_read1_core
  [[maybe_unused]] std::vector<read_t> bam_handler::fetch(const char *t_chrom, long t_start,
                                                          long t_end) {
    const int tid = bam_name2id(sam_header, t_chrom);
    hts_itr_t *iter = sam_itr_queryi(sam_index, tid, t_start, t_end);
    std::vector<read_t> reads;

    while (sam_itr_next(sam_file, iter, sam_record) >= 0) {
      const uint32_t *cigar = bam_get_cigar(sam_record);
      const uint8_t *seq = bam_get_seq(sam_record);
      std::string read_seq{};
      read_seq.reserve(sam_record->core.l_qseq);

      for (int i{0}; i < sam_record->core.l_qseq; i++) {
        read_seq += seq_nt16_str[bam_seqi(seq, i)];
      }

      reads.emplace_back(bam_get_qname(sam_record), t_chrom, sam_record->core.pos,
                         bam_endpos(sam_record), sam_record->core.qual, cigar,
                         sam_record->core.n_cigar, read_seq, bam_is_rev(sam_record));
    }
    return reads;
  }

}  // namespace bam_parser
