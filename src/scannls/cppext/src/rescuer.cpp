//
// Created by li002252 on 2/11/22.
//
#include "rescuer.h"

namespace rescuer {

  // Constructor for rescuer
  Rescuer::Rescuer(const char *t_file, int t_mapq, int t_soft_len, int t_mismatch,
                   double t_identity, int t_min_seq_align_len)
      : m_file_path{t_file},
        m_bam_handler{t_file},
        min_mapq{t_mapq},
        min_soft_len{t_soft_len},
        min_mismatch{t_mismatch},
        min_identity{t_identity},
        min_seq_align_len{t_min_seq_align_len} {}

  int Rescuer::calculate_sr(std::string_view t_chrom, long t_start, long t_end, int t_mode,
                            std::string_view t_strand,
                            std::vector<std::string> &t_current_query_name) const {

    std::vector<std::string> sr_list{};
    std::vector<std::string> sv_list{};

    auto sr_list_names = add_sr_sv_list(sr_list, sv_list, t_chrom, t_start, t_end, t_mode, t_strand, t_current_query_name);

    if (sr_list.empty() || sv_list.empty()) return 0;

    return determine_num_increment_sr(sr_list, sv_list, sr_list_names);
  }

  int Rescuer::determine_num_increment_sr(std::vector<std::string> &t_sr,
                                          std::vector<std::string> &t_sv,
                                          std::vector<std::string> &t_sr_list_names) const {
    int num_increment_sr{0};
    int sr_index{-1};

    for (const auto &r : t_sr) {
      ++sr_index;
      bool is_sr_increment{false};

      for (const auto &v : t_sv) {
        if (!check_if_align(r.substr(0, min_seq_align_len), v.substr(0, min_seq_align_len)))
          continue;
        // reference length
        int v_len{static_cast<int>(v.length())};
        int mask_len = v_len >= 30 ? v_len / 2 : 15;

        // check if conduct alignment
        if (bool return_value
            = m_aligner.Align(r.c_str(), v.c_str(), v_len, m_filter, &m_alignment, mask_len);
            !return_value) {
          std::cout << "No alignment found"
                    << "\n";
          continue;
        }
        double identity{static_cast<double>(m_alignment.query_end - m_alignment.query_begin + 1)
                        / static_cast<double>(r.length())};

        if (identity >= min_identity && (m_alignment.query_begin + m_alignment.ref_begin) <= 2
            && m_alignment.mismatches <= min_mismatch) {
          ++num_increment_sr;
          is_sr_increment = true;
          break;
        }
      }
//       add query name of sr to all names list if the query name of sr update successfully
//       to ensure that the query name of sr  do not update repeatedly
      if (is_sr_increment) {
        m_names_list.push_back(t_sr_list_names[sr_index]);
      }

    }
    return num_increment_sr;
  }
  [[maybe_unused]] int Rescuer::count_reads(const std::string &t_chrom, long t_start,
                                            long t_end) const {
    return m_bam_handler.count(t_chrom.c_str(), t_start, t_end);
  }

  bool Rescuer::check_if_align(std::string_view t_query, std::string_view t_target) const {
    int target_len{static_cast<int>(t_target.length())};
    int masklen{target_len / 2 > 15 ? target_len / 2 : 15};
    bool return_value{m_aligner.Align(std::string(t_query).c_str(), std::string(t_target).c_str(),
                                      target_len, m_filter, &m_alignment, masklen)};

    if (double identity{static_cast<double>(m_alignment.query_end - m_alignment.query_begin + 1)
                        / static_cast<double>(t_query.length())};
        !return_value || identity < 0.8) {
      return false;  // do not align
    }
    return true;
  }

  void Rescuer::reset_names_list(std::vector<std::string> &t_names_list) const {
    m_names_list.swap(t_names_list);
  }

  std::vector<std::string> Rescuer::add_sr_sv_list(std::vector<std::string> &t_sr, std::vector<std::string> &t_sv,
                               std::string_view tt_chrom, long tt_start,
                      long tt_end, int tt_mode, std::string_view tt_strand,
                      std::vector<std::string> &t_current_names) const {

    std::vector<std::string> sr_list_names;

    --tt_start;
    const int tid = bam_name2id(m_bam_handler.sam_header, std::string(tt_chrom).c_str());

    hts_itr_t *iter = sam_itr_queryi(m_bam_handler.sam_index, tid, tt_start, tt_end);

    while (sam_itr_next(m_bam_handler.sam_file, iter, m_bam_handler.sam_record) >= 0) {
      if (m_bam_handler.sam_record->core.flag & BAM_FREVERSE) {
        if (tt_strand == "+") continue;
      } else {
        if (tt_strand == "-") continue;
      }

      const uint32_t *cigar{bam_get_cigar(m_bam_handler.sam_record)};
      const uint8_t *seq{bam_get_seq(m_bam_handler.sam_record)};

      //    get read seq
      if (m_bam_handler.sam_record->core.qual >= min_mapq
          && is_soft_clipped(cigar, m_bam_handler.sam_record->core.n_cigar)) {
        std::string read_seq{get_read_seq(seq, m_bam_handler.sam_record->core.l_qseq)};

        auto softclip_result{get_softclip(read_seq, m_bam_handler.sam_record, tt_mode, cigar,
                                          m_bam_handler.sam_record->core.n_cigar)};

        long reference_pos
            = (tt_mode == 2) ? m_bam_handler.sam_record->core.pos : bam_endpos(m_bam_handler.sam_record);

        int seq_len{get_read_max_length(softclip_result.read_seq)};
        if (seq_len < min_seq_align_len) continue;  // read length is too short

        auto read_name{bam_get_qname(m_bam_handler.sam_record)};

        if (reference_pos == softclip_result.pos
            && find(t_current_names.begin(), t_current_names.end(), read_name)
                   != t_current_names.end()) {


          if (tt_mode == 2)
            t_sv.emplace_back(softclip_result.read_seq.rbegin(),
                              softclip_result.read_seq.rbegin() + seq_len);
          else
            t_sv.emplace_back(softclip_result.read_seq.begin(),
                              softclip_result.read_seq.begin() + seq_len);

        } else if (reference_pos == softclip_result.pos && softclip_result.soft_len >= min_soft_len
                   && find(m_names_list.begin(), m_names_list.end(),read_name )
                          == m_names_list.end()) {


          sr_list_names.emplace_back(read_name);

          if (tt_mode == 2)
            t_sr.emplace_back(softclip_result.read_seq.rbegin(),
                              softclip_result.read_seq.rbegin() + seq_len);
          else
            t_sr.emplace_back(softclip_result.read_seq.begin(),
                              softclip_result.read_seq.begin() + seq_len);
        }
      }
    }
    sam_itr_destroy(iter);

    return sr_list_names;
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
          result.indel_len -= static_cast<int>(len);
          result.read_match += len;
          result.query_len += len;
          result.cigartuples_without_soft.insert(result.cigartuples_without_soft.end(), {op, len});
          break;
        case BAM_CDEL:
        case BAM_CREF_SKIP:
          result.indel_len += static_cast<int>(len);
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

  inline int get_read_max_length(std::string_view t_seq) {
    if (auto seq_len = static_cast<int>(t_seq.size()); seq_len < max_seq_len) return seq_len;
    return max_seq_len;
  }

}  // namespace rescuer
