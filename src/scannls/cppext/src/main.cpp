#include "rescuer.h"
#include "utils.hpp"

int main(int argc, char* argv[]) {
  using namespace cppext;

  std::cout << argc << argv[0] << '\n';

  //  constexpr std::string_view file{
  //      "/projects/b1171/ylk4626/project/scannls/src/scannls/cppext/test/chr2_214987647.test.bam"};
  //      "test.bam"};
    constexpr std::string_view file{"../test/data/chr2_214987647.test.bam"};
//  constexpr std::string_view file{"../test/data/test.bam"};

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
  std::vector<std::string> current_names{"559_4", "559_2",  "559_3", "559_1", "559_6",
                                         "559_5", "559_10", "559_9", "559_7", "559_8"};
  std::vector<std::string> names_list{"559_2", "559_4", "559_10", "559_8", "559_7",
                                      "559_5", "559_9", "559_1",  "559_6", "559_3"};

  constexpr long const start{216206704};
  constexpr long const read_start{216206403};
  constexpr std::string_view const strand{"-"};
  constexpr int const mode{1};
  std::vector<uint> const cigar_tuple_list{0, 301};

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

  Region const region{"chr2", start - 1, start};

  //  long read_start{};  // read match start position
  //  long read_end{};    // read match end position
  //  int mode{};
  //  bool is_reverse{};
  //  std::string breakpoint_chrom{};
  //  long breakpoint_start{};
  //  std::optional<long> breakpoint_end{};
  BreakPoint const breakpoint{read_start, 0, mode, true, "chr2", start - 1};

  log(region.to_string());

  auto sr = rescuer.calculate_sr(region, breakpoint, current_names, cigar_tuple_list);
  std::cout << "sr: " << sr << "\n";

  return 0;
}
