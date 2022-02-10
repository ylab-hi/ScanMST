#include "bam.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

PYBIND11_MODULE(cppext, m){
    m.doc() = "BAM file parser";
    m.def("parseCigar", &parseCigar, "parse cigar string");
    py::class_<parseCigarResult_t>(m, "parseCigarResult")
            .def_readonly("cigartuples", &parseCigarResult_t::cigartuples)
            .def_readonly("cigartuples_without_soft", &parseCigarResult_t::cigartuples_without_soft)
            .def_readonly("lt_soft_len", &parseCigarResult_t::lt_soft_len)
            .def_readonly("rt_soft_len", &parseCigarResult_t::rt_soft_len)
            .def_readonly("ref_match", &parseCigarResult_t::ref_match)
            .def_readonly("read_match", &parseCigarResult_t::read_match)
            .def_readonly("query_len", &parseCigarResult_t::query_len)
            .def_readonly("indel_len", &parseCigarResult_t::indel_len);

}
