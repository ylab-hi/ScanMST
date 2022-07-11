#include "rescuer.h"

int main(int argc, char* argv[]) {
  std::cout << argc << argv[0] << '\n';

  std::string file{"/panfs/home/yang4414/li002252/project/scan_data/resue_test.bam"};

  constexpr int mapq_threshold = 15;
  constexpr int min_softclip_length = 5;
  constexpr int min_mismatch_count = 3;
  constexpr double min_align_ratio = 0.8;

  rescuer::Rescuer p_rescuer{file.c_str(),
                             mapq_threshold,
                             min_softclip_length,
                             min_mismatch_count,
                             min_align_ratio,
                             10,
                             -1};

  std::vector<std::string> current_names{"two"};
  std::vector<std::string> names_list{"one", "two"};

  long const start{49798558};
  long const read_start{49791885};
  std::string const strand{"+"};
  int const mode{1};
  std::vector<uint> const cigar_tuple_list{0, 148, 3, 4716, 0, 1809};

  p_rescuer.reset_names_list(names_list);
  auto sr = p_rescuer.calculate_sr("chr17", start, start, mode, strand, read_start, current_names,
                                   cigar_tuple_list);

  std::cout << "sr: " << sr << "\n";

  return 0;
}
