//
// Created by li002252 on 2/11/22.
//
#include "rescuer.h"

namespace rescuer {

  // Constructor for rescuer
  Rescuer::Rescuer() = default;

  // Destructor for rescuer
  Rescuer::~Rescuer() = default;

  // Constructor for rescuer
  Rescuer::Rescuer(const char *t_file, int t_mapq, int t_soft_len, int t_mismatch,
                   double t_identity)
      : bam_file{t_file},
        min_mapq{t_mapq},
        min_soft_len{t_soft_len},
        min_mismatch{t_mismatch},
        min_identity{t_identity} {}

  int Rescuer::calculate_sr(const std::string &t_chrom, long t_start, long t_end, int t_mode,
                            std::vector<std::string> &t_current_query_name,
                            std::vector<std::string> &t_query_name_list) const {
    std::vector<std::string> sr_list{};
    std::vector<std::string> sv_list{};

    bam_handler bam_reader{bam_file};

    add_sr_sv_list(sr_list, sv_list, bam_reader, t_chrom, t_start, t_end, t_mode, min_mapq,
                   min_soft_len, t_current_query_name, t_query_name_list);

    if (sr_list.empty() || sv_list.empty()) return 0;

    return determine_num_increment_sr(sr_list, sv_list, min_identity, min_mismatch);
  }

  // non-member function
  void add_sr_sv_list(std::vector<std::string> &t_sr, std::vector<std::string> &t_sv,
                      const bam_handler &t_bam, const std::string &tt_chrom, long tt_start,
                      long tt_end, int tt_mode, int t_min_mapq, int t_min_soft,
                      const std::vector<std::string> &t_current_names,
                      const std::vector<std::string> &t_name_list) {
    --tt_start;
    const int tid = bam_name2id(t_bam.sam_header, tt_chrom.c_str());

    hts_itr_t *iter = sam_itr_queryi(t_bam.sam_index, tid, tt_start, tt_end);

    while (sam_itr_next(t_bam.sam_file, iter, t_bam.sam_record) >= 0) {
      const uint32_t *cigar{bam_get_cigar(t_bam.sam_record)};
      const uint8_t *seq{bam_get_seq(t_bam.sam_record)};

      //    get read seq
      if (t_bam.sam_record->core.qual >= t_min_mapq
          && is_soft_clipped(cigar, t_bam.sam_record->core.n_cigar)) {
        std::string read_seq{get_read_seq(seq, t_bam.sam_record->core.l_qseq)};

        auto softclip_result{get_softclip(read_seq, t_bam.sam_record, tt_mode, cigar,
                                          t_bam.sam_record->core.n_cigar)};

        long reference_pos
            = (tt_mode == 2) ? t_bam.sam_record->core.pos : bam_endpos(t_bam.sam_record);

        if (reference_pos == softclip_result.pos
            && find(t_current_names.begin(), t_current_names.end(), bam_get_qname(t_bam.sam_record))
                   != t_current_names.end()) {
          if (tt_mode == 2)
            t_sv.emplace_back(softclip_result.read_seq.rbegin(), softclip_result.read_seq.rend());
          else
            t_sv.push_back(softclip_result.read_seq);
        } else if (reference_pos == softclip_result.pos && softclip_result.soft_len >= t_min_soft
                   && find(t_name_list.begin(), t_name_list.end(), bam_get_qname(t_bam.sam_record))
                          == t_name_list.end()) {
          if (tt_mode == 2)
            t_sr.emplace_back(softclip_result.read_seq.rbegin(), softclip_result.read_seq.rend());
          else
            t_sr.push_back(softclip_result.read_seq);
        }
      }
    }
  }

  bool is_soft_clipped(const uint32_t *t_cigar, size_t t_cigar_len) {
    for (size_t i{0}; i < t_cigar_len; i++) {
      if (bam_cigar_op(t_cigar[i]) == BAM_CSOFT_CLIP) {
        return true;
      }
    }
    return false;
  }

  std::string get_read_seq(const uint8_t *t_seqs, int t_seq_len) {
    std::string result_seq{};
    result_seq.reserve(t_seq_len);
    for (int i{0}; i < t_seq_len; i++) result_seq += seq_nt16_str[bam_seqi(t_seqs, i)];
    return result_seq;
  }

  parseCigarResult_t parser_cigar(const uint32_t *t_cigar_str, size_t t_cigar_len) {
    parseCigarResult_t result{};

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
    if (result.cigartuples[2 * t_cigar_len - 2] == 4)
      result.rt_soft_len = result.cigartuples[2 * t_cigar_len - 1];
    return result;
  }

  get_softclip_result_t get_softclip(const std::string &t_read_seq, const bam1_t *t_alignment,
                                     int t_mode, const uint32_t *t_cigar_str, size_t t_cigar_len) {
    parseCigarResult_t cigar_result{parser_cigar(t_cigar_str, t_cigar_len)};

    long ref_end{t_alignment->core.pos + cigar_result.ref_match};

    if (t_mode == 0) {
      if (cigar_result.lt_soft_len > cigar_result.rt_soft_len)

        return {static_cast<int>(cigar_result.lt_soft_len),
                t_read_seq.substr(0, cigar_result.lt_soft_len), t_alignment->core.pos, 2};

      else if (cigar_result.lt_soft_len < cigar_result.rt_soft_len)
        return {static_cast<int>(cigar_result.rt_soft_len),
                t_read_seq.substr(cigar_result.query_len - cigar_result.rt_soft_len), ref_end, 1};
      else
        return {};

    } else if (t_mode == 1) {
      return {static_cast<int>(cigar_result.rt_soft_len),
              t_read_seq.substr(cigar_result.query_len - cigar_result.rt_soft_len), ref_end, 1};
    }

    else if (t_mode == 2)
      return {static_cast<int>(cigar_result.lt_soft_len),
              t_read_seq.substr(0, cigar_result.lt_soft_len), t_alignment->core.pos, 2};

    return {};
  }

  int determine_num_increment_sr(const std::vector<std::string> &t_sr,
                                 const std::vector<std::string> &t_sv, double t_min_frac,
                                 int t_min_mismatch) {
    static StripedSmithWaterman::Aligner aligner{StripedSmithWaterman::Aligner{}};
    static StripedSmithWaterman::Filter filter{StripedSmithWaterman::Filter{}};
    static StripedSmithWaterman::Alignment alignment{};
    alignment.Clear();

    int num_increment_sr{0};

    for (const auto &r : t_sr) {
      for (const auto &v : t_sv) {
        if (bool return_value = aligner.Align(r.c_str(), v.c_str(), static_cast<int>(r.length()),
                                              filter, &alignment, 16);
            !return_value) {
          std::cout << "No alignment found"
                    << "\n";
          continue;
        }

        //        StripedSmithWaterman::print_alignment(r, v, alignment);
        if (alignment.query_begin + alignment.ref_begin == 0
            && alignment.mismatches <= t_min_mismatch
            && (alignment.query_end - alignment.query_begin + 1) / (int)r.length() >= t_min_frac) {
          ++num_increment_sr;
        } else
          break;
      }
    }
    return num_increment_sr;
  }

}  // namespace rescuer
