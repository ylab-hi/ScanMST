#include "rescuer.h"
#include "utils.hpp"

int main(int argc, char* argv[]) {
  using namespace cppext;

  std::cout << argc << argv[0] << '\n';

  constexpr std::string_view file{"../test/data/TDUP_case.bam"};

  constexpr int mapq_threshold = 15;
  constexpr int min_softclip_length = 5;
  constexpr int min_mismatch_count = 3;
  constexpr double min_align_ratio = 0.8;

  auto options = cppext::Options()
                     .file(file)
                     .mapq(mapq_threshold)
                     .soft_len(min_softclip_length)
                     .mismatch(min_mismatch_count)
                     .identity(min_align_ratio)
                     .min_seq_align_len(10);

  Rescuer rescuer(options);
  std::vector<std::string> current_names{"423_1"};
  std::vector<std::string> names_list{"423_3", "423_7", "423_1",  "423_9", "423_4",
                                      "423_8", "423_2", "423_10", "423_5", "423_6"};

  constexpr long start{49892642};
  constexpr long read_start{49892344};
  constexpr long read_end{49892642};
  constexpr std::string_view strand{"-"};
  constexpr int mode{1};

  std::vector<uint> const cigar_tuple_list{0, 83, 2, 1, 0, 197, 2, 1, 0, 16};

  //  chrom,
  //      start,
  //      start,
  //      mode1,
  //      current_node.strand,
  //      current_node.ref_start,
  //      query_name_current,
  //      current_node.cigartuples_without_soft,
  //
  rescuer.reset_names_list(names_list);

  Region const region{"chr20", start - 1, start};
  BreakPoint const breakpoint{read_start, read_end, mode, true, false};

  log(region.to_string());

  auto sr = rescuer.calculate_sr(region, breakpoint, current_names, cigar_tuple_list);
  std::cout << "sr: " << sr << "\n";

  return 0;
}
