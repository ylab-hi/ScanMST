#include <iostream>
#include <string>
#include <vector>
#include "bam.h"

int main(int argc, char *argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <bam file>" << std::endl;
        return 1;
    }
//    readBam(argv[1]);
    const char *cigar = argv[1];
    parseCigarResult_t res{parseCigar(cigar)};


    for (auto &r : res.cigartuples) {
        std::cout << r << " ";
    }
    std::cout << "\n";
    std::cout << "result: " << "\n"
              << "lt len: " << res.lt_soft_len << " tr len: " << res.rt_soft_len << "\n";
    return 0;
}
