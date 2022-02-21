#include "rescuer.h"

int main(int argc, char* argv[]) {
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " <input_file>" << std::endl;
    return 1;
  }

  int mapq_threshold = 15;
  int min_softclip_length = 5;
  int min_mismatch_count = 3;
  double min_align_ratio = 0.8;
  std::vector<std::string> names_list{"one_hop_back"};
  std::vector<std::string> current_names{"one_hop_back"};

  rescuer::Rescuer p_rescuer{argv[1], mapq_threshold, min_softclip_length, min_mismatch_count,
                             min_align_ratio};

  int sr = p_rescuer.calculate_sr("chr17", 7708250, 7708250, 1, current_names,
                                  names_list);  // chr17:7708250-7708250

  std::cout << "sr: " << sr << "\n";

  return 0;
}
