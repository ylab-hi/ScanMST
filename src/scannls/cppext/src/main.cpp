#include "rescuer.h"

int main(int argc, char* argv[]) {
  std::cout << argc << argv[0] << std::endl;

  std::string file{"/panfs/home/yang4414/li002252/project/scan_data/sr_count.bam"};

  int mapq_threshold = 15;
  int min_softclip_length = 5;
  int min_mismatch_count = 3;
  double min_align_ratio = 0.8;

  rescuer::Rescuer p_rescuer{file.c_str(),       mapq_threshold,  min_softclip_length,
                             min_mismatch_count, min_align_ratio, 10};

  std::vector<std::string> current_names{"m64135_201202_230957/107415515/ccs"};
//  std::vector<std::string> names_list {"m64135_201202_230957/107415515/ccs"};
  std::vector<std::string> names_list { "m64135_201202_230957/107415515/ccs", "m64135_201202_230957/2097763/ccs", "m64135_201202_230957/53543856/ccs", "m64135_201202_230957/79825441/ccs", "m64135_201202_230957/99419441/ccs", "m64135_201202_230957/164889621/ccs", "m64135_201202_230957/170002391/ccs", "m64135_201202_230957/172361197/ccs" };


  int sr{-1};

//  long start{150804249};
  long start{84472987};
  std::string strand{"-"};

  p_rescuer.reset_names_list(names_list);

  sr = p_rescuer.calculate_sr("chr5", start, start, 2, strand,
                              current_names);  // chr17:7708250-7708250

  std::cout << "sr: " << sr << "\n";

  return 0;
}
