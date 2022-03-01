#include "rescuer.h"

int main(int argc, char* argv[]) {
  std::cout << argc << argv[0] << std::endl;

  std::string file{"/panfs/home/yang4414/li002252/project/scan_data/resue2.bam"};

  int mapq_threshold = 15;
  int min_softclip_length = 5;
  int min_mismatch_count = 3;
  double min_align_ratio = 0.8;

  std::vector<std::string> names_list{"one_1", "one_3", "one_2"};
  std::vector<std::string> current_names{"one_3", "one_1", "one_2"};

  rescuer::Rescuer p_rescuer{file.c_str(), mapq_threshold, min_softclip_length, min_mismatch_count,
                             min_align_ratio};
  long start{190670461};
  int sr = p_rescuer.calculate_sr("chr2", start, start, 1, current_names,
                                  names_list);  // chr17:7708250-7708250
  p_rescuer.check_if_align("ATAATTGGCC", "TTCCGACGTT");

  std::cout << "sr: " << sr << "\n";

  return 0;
}
