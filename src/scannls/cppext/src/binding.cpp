#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "bam.h"
#include "rescuer.h"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;
using bam_parser::bam_handler;
using bam_parser::parseCigarResult_t;
using rescuer::Rescuer;

PYBIND11_MODULE(_cppext, m) {
  m.doc() = "Cpp extension for BAM file parser";
  m.def("parseCigar", &bam_parser::parseCigar, "parse cigar string");
  py::class_<parseCigarResult_t>(m, "parseCigarResult")
      .def_readonly("cigartuples", &parseCigarResult_t::cigartuples)
      .def_readonly("cigartuples_without_soft", &parseCigarResult_t::cigartuples_without_soft)
      .def_readonly("lt_soft_len", &parseCigarResult_t::lt_soft_len)
      .def_readonly("rt_soft_len", &parseCigarResult_t::rt_soft_len)
      .def_readonly("ref_match", &parseCigarResult_t::ref_match)
      .def_readonly("read_match", &parseCigarResult_t::read_match)
      .def_readonly("query_len", &parseCigarResult_t::query_len)
      .def_readonly("indel_len", &parseCigarResult_t::indel_len)
      .def("__repr__", [](const parseCigarResult_t &r) { return "parseCigarResult()"; });

  py::class_<StripedSmithWaterman::Alignment>(m, "Alignment")
      .def(py::init<>())
      .def_readwrite("best_score", &StripedSmithWaterman::Alignment::sw_score)
      .def_readwrite("best_score2", &StripedSmithWaterman::Alignment::sw_score_next_best)
      .def_readwrite("reference_begin", &StripedSmithWaterman::Alignment::ref_begin)
      .def_readwrite("reference_end", &StripedSmithWaterman::Alignment::ref_end)
      .def_readwrite("query_begin", &StripedSmithWaterman::Alignment::query_begin)
      .def_readwrite("query_end", &StripedSmithWaterman::Alignment::query_end)
      .def_readwrite("ref_end_next_best", &StripedSmithWaterman::Alignment::ref_end_next_best)
      .def_readwrite("mismatches", &StripedSmithWaterman::Alignment::mismatches)
      .def_readwrite("cigar_string", &StripedSmithWaterman::Alignment::cigar_string)
      .def_readwrite("cigar", &StripedSmithWaterman::Alignment::cigar)
      .def("Clear", &StripedSmithWaterman::Alignment::Clear);

  // Filter Class
  py::class_<StripedSmithWaterman::Filter>(m, "Filter")
      .def_readwrite("report_begin_position", &StripedSmithWaterman::Filter::report_begin_position)
      .def_readwrite("report_cigar", &StripedSmithWaterman::Filter::report_cigar)
      .def_readwrite("score_filter", &StripedSmithWaterman::Filter::score_filter)
      .def_readwrite("distance_filter", &StripedSmithWaterman::Filter::distance_filter)
      .def(py::init<>())
      .def(py::init<const bool &, const bool &, const uint16_t &, const uint16_t &>());

  // Aligner Class
  py::class_<StripedSmithWaterman::Aligner>(m, "Aligner")
      .def(py::init<>())
      .def(py::init<const uint8_t &, const uint8_t &, const uint8_t &, const uint8_t &>())
      .def("SetReferenceSequence", &StripedSmithWaterman::Aligner::SetReferenceSequence);

  py::class_<bam_handler>(m, "bam_handler")
      .def("count", &bam_handler::count, "count(chrom, star, end) -> int");

  py::class_<Rescuer>(m, "Rescuer",
                      "Rescuer(bam_file, min_mapq, min_soft_len, min_mis, min_frac, "
                      "min_seq_align_len, average_read_depth)")
      .def(py::init<const char *, int, int, int, double, int, int>())
      .def("calculate_sr", &Rescuer::calculate_sr,
           "calculate_sr(chrom, start, end, mode, strand, read_start, current_names, "
           "cigartuples_without_soft) -> int")
      .def("count_reads", &Rescuer::count_reads, "count_reads(chrom, start, end) -> int")
      .def("__repr__", [](const Rescuer &r) { return "Rescuer()"; })
      .def("reset_names_list", &Rescuer::reset_names_list, "reset_names_list(names_list) -> None");
}
