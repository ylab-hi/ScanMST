#include "rescuer.h"

int main(int argc, char* argv[]) {
  std::cout << argc << argv[0] << '\n';

  std::string file{"/panfs/home/yang4414/li002252/project/scan_data/DU145.extract.bam"};

  constexpr int mapq_threshold = 15;
  constexpr int min_softclip_length = 5;
  constexpr int min_mismatch_count = 3;
  constexpr double min_align_ratio = 0.8;

  rescuer::Rescuer p_rescuer{file.c_str(),       mapq_threshold,  min_softclip_length,
                             min_mismatch_count, min_align_ratio, 10};

  std::vector<std::string> current_names{"m64135_201204_204719/171246914/ccs"};
  std::vector<std::string> names_list{"m64135_201204_204719/171246914/ccs"};

  long const start{65499397};
  std::string const strand{"+"};

  p_rescuer.reset_names_list(names_list);
  auto sr = p_rescuer.calculate_sr("chr11", start, start, 2, strand,
                                   current_names);  // chr17:7708250-7708250

  std::cout << "sr: " << sr << "\n";

  return 0;
}
