#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "bam.h"
#include "rescuer.h"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;
using cppext::BamReader;
using cppext::Options;
using cppext::parseCigarResult_t;
using cppext::Region;
using cppext::Rescuer;

PYBIND11_MODULE(_cppext, m) {
  m.doc() = "Cpp extension for BAM file parser";
  m.def("parseCigar", &cppext::parseCigar, "parse cigar string");
  py::class_<parseCigarResult_t>(m, "parseCigarResult")
      .def_readonly("cigartuples", &parseCigarResult_t::cigartuples)
      .def_readonly("cigartuples_without_soft", &parseCigarResult_t::cigartuples_without_soft)
      .def_readonly("lt_soft_len", &parseCigarResult_t::lt_soft_len)
      .def_readonly("rt_soft_len", &parseCigarResult_t::rt_soft_len)
      .def_readonly("ref_match", &parseCigarResult_t::ref_match)
      .def_readonly("read_match", &parseCigarResult_t::read_match)
      .def_readonly("query_len", &parseCigarResult_t::query_len)
      .def_readonly("indel_len", &parseCigarResult_t::indel_len)
      .def("__repr__", [](const parseCigarResult_t &r) { return r.to_string(); });

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
      .def(py::init<bool, bool, uint16_t, uint16_t>());

  // Aligner Class
  py::class_<StripedSmithWaterman::Aligner>(m, "Aligner")
      .def(py::init<>())
      .def(py::init<uint8_t, uint8_t, uint8_t, uint8_t>())
      .def("SetReferenceSequence", &StripedSmithWaterman::Aligner::SetReferenceSequence);

  py::class_<BamReader>(m, "BamReader")
      .def("count", py::overload_cast<const Region &>(&BamReader::count, py::const_),
           "count reads in a region")
      .def("count", py::overload_cast<std::string_view, long, long>(&BamReader::count, py::const_),
           "count reads in a region");

  py::class_<Region>(m, "Region")
      .def(py::init<>())
      .def(py::init<const std::string &, uint32_t, const int &>())
      .def_readwrite("chrom", &Region::chrom)
      .def_readwrite("start", &Region::start)
      .def_readwrite("end", &Region::end);

  py::class_<Options>(m, "Options")
      .def(py::init<>())
      .def_readwrite("mapq_", &Options::mapq_)
      .def_readwrite("soft_len_", &Options::soft_len_)
      .def_readwrite("mismatch_", &Options::mismatch_)
      .def_readwrite("identity_", &Options::identity_)
      .def_readwrite("min_seq_align_len_", &Options::min_seq_align_len_)
      .def_readwrite("average_read_depth_", &Options::average_read_depth_)
      .def("file", &Options::file)
      .def("mapq", &Options::mapq)
      .def("soft_len", &Options::soft_len)
      .def("mismatch", &Options::mismatch)
      .def("identity", &Options::identity)
      .def("min_seq_align_len", &Options::min_seq_align_len)
      .def("average_read_depth", &Options::average_read_depth)
      .def("__repr__", [](const Options &o) { return o.to_string(); });


  py::class_<Rescuer>(m, "Rescuer")
      .def(py::init<Options const &>())
      .def("calculate_sr", &Rescuer::calculate_sr,
           "calculate_sr(Region, mode, strand, read_start, current_names, "
           "cigartuples_without_soft) -> int")
      .def("count_reads", &Rescuer::count_reads, "count_reads(chrom, start, end) -> int")
      .def("__repr__", [](const Rescuer &r) { return r.to_string(); })
      .def("reset_names_list", &Rescuer::reset_names_list, "reset_names_list(names_list) -> None");
}
