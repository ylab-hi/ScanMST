#include <bam.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <rescuer.h>

#include <functional>
#include <ios>
#include <iterator>
#include <locale>
#include <memory>
#include <optional>
#include <ostream>
#include <sstream>  // __str__
#include <streambuf>
#include <string>
#include <string_view>
#include <vector>

#ifndef BINDER_PYBIND11_TYPE_CASTER
#  define BINDER_PYBIND11_TYPE_CASTER
PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)
PYBIND11_DECLARE_HOLDER_TYPE(T, T *)
PYBIND11_MAKE_OPAQUE(std::shared_ptr<void>)
#endif

// clang-format off

void bind_rescuer(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	// cppext::parser_cigar(const unsigned int *, unsigned long) file:rescuer.h line:55
	M("cppext").def("parser_cigar", (struct cppext::CigarResult (*)(const unsigned int *, unsigned long)) &cppext::parser_cigar, "Parse cigar string from original uint32_t array\n \n\n uint32_t array of cigar string\n \n\n length of cigar string array\n \n\n: parseCigarResult_t\n\nC++: cppext::parser_cigar(const unsigned int *, unsigned long) --> struct cppext::CigarResult", pybind11::arg("t_cigar_str"), pybind11::arg("t_cigar_len"));

	// cppext::parser_cigar(const struct cppext::BamReader::Iterator &) file:rescuer.h line:62
	M("cppext").def("parser_cigar", (struct cppext::CigarResult (*)(const struct cppext::BamReader::Iterator &)) &cppext::parser_cigar, "Parse cigar string from record iterator\n \n\n\n \n\n parseCigarResult_t\n\nC++: cppext::parser_cigar(const struct cppext::BamReader::Iterator &) --> struct cppext::CigarResult", pybind11::arg("iterator"));

	{ // cppext::Rescuer file:rescuer.h line:68
		pybind11::class_<cppext::Rescuer, std::shared_ptr<cppext::Rescuer>> cl(M("cppext"), "Rescuer", "" );
		cl.def( pybind11::init<const struct cppext::Options &>(), pybind11::arg("t_options") );

		cl.def("to_string", (std::string (cppext::Rescuer::*)() const) &cppext::Rescuer::to_string, "C++: cppext::Rescuer::to_string() const --> std::string");
		cl.def("calculate_sr", (int (cppext::Rescuer::*)(const struct cppext::Region &, const struct cppext::BreakPoint &, class std::vector<std::string > &, const class std::vector<unsigned int> &) const) &cppext::Rescuer::calculate_sr, "calculate number of sr for every node\n \n\n info including chrom, start, end\n \n\n  mode of node\n \n\n current query names\n \n\n query names list\n \n\n number of sr\n\nC++: cppext::Rescuer::calculate_sr(const struct cppext::Region &, const struct cppext::BreakPoint &, class std::vector<std::string > &, const class std::vector<unsigned int> &) const --> int", pybind11::arg("region"), pybind11::arg("break_point"), pybind11::arg("current_query_name"), pybind11::arg("cigartuples_without_soft"));
		cl.def("count_reads", (int (cppext::Rescuer::*)(const std::string &, long, long) const) &cppext::Rescuer::count_reads, "calculate number of read for a position\n \n\n chromosome\n \n\n start position 0-based\n \n\n  end position\n \n\n number of read\n\nC++: cppext::Rescuer::count_reads(const std::string &, long, long) const --> int", pybind11::arg("t_chrom"), pybind11::arg("t_start"), pybind11::arg("t_end"));
		cl.def("check_rescue", (bool (cppext::Rescuer::*)(const struct cppext::Seqs &, const struct cppext::Seqs &) const) &cppext::Rescuer::check_rescue, "check if two read sequence need to be aligned or not\n \n\n query sequence\n \n\n target sequence\n \n\n true if need to be aligned\n\nC++: cppext::Rescuer::check_rescue(const struct cppext::Seqs &, const struct cppext::Seqs &) const --> bool", pybind11::arg("query"), pybind11::arg("target"));
		cl.def("check_identity", (bool (cppext::Rescuer::*)(class std::basic_string_view<char>, class std::basic_string_view<char>, double) const) &cppext::Rescuer::check_identity, "C++: cppext::Rescuer::check_identity(class std::basic_string_view<char>, class std::basic_string_view<char>, double) const --> bool", pybind11::arg("query"), pybind11::arg("target"), pybind11::arg("threshold"));
		cl.def("check_rescue_condition", (bool (cppext::Rescuer::*)(class std::basic_string_view<char>, class std::basic_string_view<char>) const) &cppext::Rescuer::check_rescue_condition, "C++: cppext::Rescuer::check_rescue_condition(class std::basic_string_view<char>, class std::basic_string_view<char>) const --> bool", pybind11::arg("query"), pybind11::arg("target"));
		cl.def("calculate_identity", (class std::optional<double> (cppext::Rescuer::*)(class std::basic_string_view<char>, class std::basic_string_view<char>) const) &cppext::Rescuer::calculate_identity, "query sequence\n \n\n  target sequence\n \n\n the identity of query\n\nC++: cppext::Rescuer::calculate_identity(class std::basic_string_view<char>, class std::basic_string_view<char>) const --> class std::optional<double>", pybind11::arg("query"), pybind11::arg("target"));
		cl.def("reset_names_list", (void (cppext::Rescuer::*)(const class std::vector<std::string > &) const) &cppext::Rescuer::reset_names_list, "reset name list\n \n\n\n     \n\nC++: cppext::Rescuer::reset_names_list(const class std::vector<std::string > &) const --> void", pybind11::arg("t_names_list"));
		cl.def("get_align_sequences", (class std::optional<struct cppext::Seqs> (cppext::Rescuer::*)(const struct cppext::BamReader::Iterator &, const struct cppext::BreakPoint &, const struct cppext::CigarResult &) const) &cppext::Rescuer::get_align_sequences, "C++: cppext::Rescuer::get_align_sequences(const struct cppext::BamReader::Iterator &, const struct cppext::BreakPoint &, const struct cppext::CigarResult &) const --> class std::optional<struct cppext::Seqs>", pybind11::arg("iterator"), pybind11::arg("break_point"), pybind11::arg("cigar_result"));
		cl.def_static("get_align_seq_len", (class std::optional<int> (*)(unsigned int, unsigned int)) &cppext::Rescuer::get_align_seq_len, "check if the rescued read need to be skipped\n \n\n read sequence of rescued read\n \n\n mode of node's read\n \n\n start position of node's read for rescuing\n \n\n cigar string buff of rescued read\n \n\n cigar tuples without soft clip of node's read\n \n\n tuple<get_softclip_result_t, seq_len, is_skip>\n\nC++: cppext::Rescuer::get_align_seq_len(unsigned int, unsigned int) --> class std::optional<int>", pybind11::arg("seq_size"), pybind11::arg("min_align_len"));
		cl.def_static("check_read_pos", (bool (*)(const struct cppext::BamReader::Iterator &, const struct cppext::CigarResult &, const struct cppext::BreakPoint &)) &cppext::Rescuer::check_read_pos, "check start position for rescued reads and node\n \n\n node's read reference start\n \n\n true if rescued reads need to be skipped\n\nC++: cppext::Rescuer::check_read_pos(const struct cppext::BamReader::Iterator &, const struct cppext::CigarResult &, const struct cppext::BreakPoint &) --> bool", pybind11::arg("iterator"), pybind11::arg("cigar_result"), pybind11::arg("break_point"));
		cl.def_static("check_if_same_isform", (bool (*)(const class std::vector<unsigned int> &, const class std::vector<unsigned int> &)) &cppext::Rescuer::check_if_same_isform, "check if has same isoform\n \n\n cigar tuples without soft clip for one read\n \n\n cigar tuples without soft clip for another read\n \n\n true if has same isoform\n\nC++: cppext::Rescuer::check_if_same_isform(const class std::vector<unsigned int> &, const class std::vector<unsigned int> &) --> bool", pybind11::arg("left_cigartuples_without_soft"), pybind11::arg("right_cigartuples_without_soft"));
	}
	{ // cppext::Options file:rescuer.h line:197
		pybind11::class_<cppext::Options, std::shared_ptr<cppext::Options>> cl(M("cppext"), "Options", "" );
		cl.def( pybind11::init( [](){ return new cppext::Options(); } ) );
		cl.def( pybind11::init( [](cppext::Options const &o){ return new cppext::Options(o); } ) );
		cl.def_readwrite("file_", &cppext::Options::file_);
		cl.def_readwrite("mapq_", &cppext::Options::mapq_);
		cl.def_readwrite("soft_len_", &cppext::Options::soft_len_);
		cl.def_readwrite("mismatch_", &cppext::Options::mismatch_);
		cl.def_readwrite("identity_", &cppext::Options::identity_);
		cl.def_readwrite("min_seq_align_len_", &cppext::Options::min_seq_align_len_);
		cl.def_readwrite("average_read_depth_", &cppext::Options::average_read_depth_);
		cl.def("file", (struct cppext::Options & (cppext::Options::*)(class std::basic_string_view<char>)) &cppext::Options::file, "C++: cppext::Options::file(class std::basic_string_view<char>) --> struct cppext::Options &", pybind11::return_value_policy::automatic, pybind11::arg("file"));
		cl.def("mapq", (struct cppext::Options & (cppext::Options::*)(int)) &cppext::Options::mapq, "C++: cppext::Options::mapq(int) --> struct cppext::Options &", pybind11::return_value_policy::automatic, pybind11::arg("mapq"));
		cl.def("soft_len", (struct cppext::Options & (cppext::Options::*)(int)) &cppext::Options::soft_len, "C++: cppext::Options::soft_len(int) --> struct cppext::Options &", pybind11::return_value_policy::automatic, pybind11::arg("soft_len"));
		cl.def("mismatch", (struct cppext::Options & (cppext::Options::*)(int)) &cppext::Options::mismatch, "C++: cppext::Options::mismatch(int) --> struct cppext::Options &", pybind11::return_value_policy::automatic, pybind11::arg("mismatch"));
		cl.def("identity", (struct cppext::Options & (cppext::Options::*)(double)) &cppext::Options::identity, "C++: cppext::Options::identity(double) --> struct cppext::Options &", pybind11::return_value_policy::automatic, pybind11::arg("identity"));
		cl.def("min_seq_align_len", (struct cppext::Options & (cppext::Options::*)(int)) &cppext::Options::min_seq_align_len, "C++: cppext::Options::min_seq_align_len(int) --> struct cppext::Options &", pybind11::return_value_policy::automatic, pybind11::arg("min_seq_align_len"));
		cl.def("average_read_depth", (struct cppext::Options & (cppext::Options::*)(int)) &cppext::Options::average_read_depth, "C++: cppext::Options::average_read_depth(int) --> struct cppext::Options &", pybind11::return_value_policy::automatic, pybind11::arg("average_read_depth"));
		cl.def("to_string", (std::string (cppext::Options::*)() const) &cppext::Options::to_string, "C++: cppext::Options::to_string() const --> std::string");
		cl.def("assign", (struct cppext::Options & (cppext::Options::*)(const struct cppext::Options &)) &cppext::Options::operator=, "C++: cppext::Options::operator=(const struct cppext::Options &) --> struct cppext::Options &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // cppext::BreakPoint file:rescuer.h line:227
		pybind11::class_<cppext::BreakPoint, std::shared_ptr<cppext::BreakPoint>> cl(M("cppext"), "BreakPoint", "" );
		cl.def( pybind11::init( [](){ return new cppext::BreakPoint(); } ) );
		cl.def( pybind11::init<long, long, int, bool, bool>(), pybind11::arg("read_start_"), pybind11::arg("read_end_"), pybind11::arg("mode_"), pybind11::arg("is_reverse_"), pybind11::arg("is_middle_"));

		cl.def( pybind11::init( [](cppext::BreakPoint const &o){ return new cppext::BreakPoint(o); } ) );
		cl.def_readwrite("read_start", &cppext::BreakPoint::read_start);
		cl.def_readwrite("read_end", &cppext::BreakPoint::read_end);
		cl.def_readwrite("mode", &cppext::BreakPoint::mode);
		cl.def_readwrite("is_reverse", &cppext::BreakPoint::is_reverse);
                cl.def_readwrite("is_middle", &cppext::BreakPoint::is_middle);
		cl.def("to_string", (std::string (cppext::BreakPoint::*)() const) &cppext::BreakPoint::to_string, "C++: cppext::BreakPoint::to_string() const --> std::string");
		cl.def("assign", (struct cppext::BreakPoint & (cppext::BreakPoint::*)(const struct cppext::BreakPoint &)) &cppext::BreakPoint::operator=, "C++: cppext::BreakPoint::operator=(const struct cppext::BreakPoint &) --> struct cppext::BreakPoint &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // cppext::Seqs file:rescuer.h line:247
		pybind11::class_<cppext::Seqs, std::shared_ptr<cppext::Seqs>> cl(M("cppext"), "Seqs", "" );
		cl.def( pybind11::init( [](){ return new cppext::Seqs(); } ) );
		cl.def( pybind11::init<class std::basic_string_view<char>>(), pybind11::arg("seq1_") );

		cl.def( pybind11::init<class std::basic_string_view<char>, class std::basic_string_view<char>>(), pybind11::arg("seq1_"), pybind11::arg("seq2_") );

		cl.def( pybind11::init( [](cppext::Seqs const &o){ return new cppext::Seqs(o); } ) );
		cl.def_readwrite("seq1", &cppext::Seqs::seq1);
		cl.def_readwrite("seq2", &cppext::Seqs::seq2);
		cl.def("to_string", (std::string (cppext::Seqs::*)() const) &cppext::Seqs::to_string, "C++: cppext::Seqs::to_string() const --> std::string");
		cl.def("assign", (struct cppext::Seqs & (cppext::Seqs::*)(const struct cppext::Seqs &)) &cppext::Seqs::operator=, "C++: cppext::Seqs::operator=(const struct cppext::Seqs &) --> struct cppext::Seqs &", pybind11::return_value_policy::automatic, pybind11::arg(""));

		cl.def("__str__", [](cppext::Seqs const &o) -> std::string { std::ostringstream s; s << o; return s.str(); } );
	}
}

// clang-format on
