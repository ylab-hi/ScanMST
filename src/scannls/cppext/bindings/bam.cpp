#include <bam.h>
#include <htslib/sam.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <functional>
#include <ios>
#include <iterator>
#include <locale>
#include <memory>
#include <ostream>
#include <sstream> // __str__
#include <streambuf>
#include <string>
#include <string_view>

#ifndef BINDER_PYBIND11_TYPE_CASTER
#define BINDER_PYBIND11_TYPE_CASTER
PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)
PYBIND11_DECLARE_HOLDER_TYPE(T, T *)
PYBIND11_MAKE_OPAQUE(std::shared_ptr<void>)
#endif

// clang-format off
void bind_bam(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // cppext::CigarResult file:bam.h line:35
		pybind11::class_<cppext::CigarResult, std::shared_ptr<cppext::CigarResult>> cl(M("cppext"), "CigarResult", "" );
		cl.def( pybind11::init( [](cppext::CigarResult const &o){ return new cppext::CigarResult(o); } ) );
		cl.def( pybind11::init( [](){ return new cppext::CigarResult(); } ) );
		cl.def_readwrite("lt_soft_len", &cppext::CigarResult::lt_soft_len);
		cl.def_readwrite("rt_soft_len", &cppext::CigarResult::rt_soft_len);
		cl.def_readwrite("read_match", &cppext::CigarResult::read_match);
		cl.def_readwrite("ref_match", &cppext::CigarResult::ref_match);
		cl.def_readwrite("indel_len", &cppext::CigarResult::indel_len);
		cl.def_readwrite("query_len", &cppext::CigarResult::query_len);
		cl.def_readwrite("cigartuples_without_soft", &cppext::CigarResult::cigartuples_without_soft);
		cl.def_readwrite("cigartuples", &cppext::CigarResult::cigartuples);
		cl.def("to_string", (std::string (cppext::CigarResult::*)() const) &cppext::CigarResult::to_string, "C++: cppext::CigarResult::to_string() const --> std::string");
		cl.def("both_soft_clipped", (bool (cppext::CigarResult::*)() const) &cppext::CigarResult::both_soft_clipped, "C++: cppext::CigarResult::both_soft_clipped() const --> bool");
		cl.def("is_right_soft_clipped", (bool (cppext::CigarResult::*)() const) &cppext::CigarResult::is_right_soft_clipped, "C++: cppext::CigarResult::is_right_soft_clipped() const --> bool");
		cl.def("is_left_soft_clipped", (bool (cppext::CigarResult::*)() const) &cppext::CigarResult::is_left_soft_clipped, "C++: cppext::CigarResult::is_left_soft_clipped() const --> bool");
		cl.def("assign", (struct cppext::CigarResult & (cppext::CigarResult::*)(const struct cppext::CigarResult &)) &cppext::CigarResult::operator=, "C++: cppext::CigarResult::operator=(const struct cppext::CigarResult &) --> struct cppext::CigarResult &", pybind11::return_value_policy::automatic, pybind11::arg(""));

		cl.def("__str__", [](cppext::CigarResult const &o) -> std::string { std::ostringstream s; s << o; return s.str(); } );
	}
	// cppext::parseCigar(const char *) file:bam.h line:59
	M("cppext").def("parseCigar", (struct cppext::CigarResult (*)(const char *)) &cppext::parseCigar, "parse cigar string\n \n\n string of cigar\n \n\n  parseCigarResult_t\n\nC++: cppext::parseCigar(const char *) --> struct cppext::CigarResult", pybind11::arg("cigar"));

	{ // cppext::Region file:bam.h line:83
		pybind11::class_<cppext::Region, std::shared_ptr<cppext::Region>> cl(M("cppext"), "Region", "" );
		cl.def( pybind11::init( [](){ return new cppext::Region(); } ) );
		cl.def( pybind11::init<class std::basic_string_view<char>, long, long>(), pybind11::arg("chrom_"), pybind11::arg("start_"), pybind11::arg("end_") );

		cl.def( pybind11::init( [](cppext::Region const &o){ return new cppext::Region(o); } ) );
		cl.def_readwrite("chrom", &cppext::Region::chrom);
		cl.def_readwrite("start", &cppext::Region::start);
		cl.def_readwrite("end", &cppext::Region::end);
		cl.def("to_string", (std::string (cppext::Region::*)() const) &cppext::Region::to_string, "C++: cppext::Region::to_string() const --> std::string");
		cl.def("assign", (struct cppext::Region & (cppext::Region::*)(const struct cppext::Region &)) &cppext::Region::operator=, "C++: cppext::Region::operator=(const struct cppext::Region &) --> struct cppext::Region &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // cppext::BamReader file:bam.h line:94
		pybind11::class_<cppext::BamReader, std::shared_ptr<cppext::BamReader>> cl(M("cppext"), "BamReader", "" );
		cl.def( pybind11::init( [](){ return new cppext::BamReader(); } ) );
		cl.def( pybind11::init<class std::basic_string_view<char>>(), pybind11::arg("bamFile") );

		cl.def("count", (int (cppext::BamReader::*)(class std::basic_string_view<char>, long, long) const) &cppext::BamReader::count, "count reads from bam file in terms of the position\n \n\n  chromosome name\n \n\n  start position\n \n\n  end position\n \n\n  number of reads\n\nC++: cppext::BamReader::count(class std::basic_string_view<char>, long, long) const --> int", pybind11::arg("t_chrom"), pybind11::arg("start_t"), pybind11::arg("end_t"));
		cl.def("count", (int (cppext::BamReader::*)(const struct cppext::Region &) const) &cppext::BamReader::count, "C++: cppext::BamReader::count(const struct cppext::Region &) const --> int", pybind11::arg("region"));
		cl.def("query", (struct cppext::BamReader::Iterator (cppext::BamReader::*)(class std::basic_string_view<char>, long, long) const) &cppext::BamReader::query, "C++: cppext::BamReader::query(class std::basic_string_view<char>, long, long) const --> struct cppext::BamReader::Iterator", pybind11::arg("chrom"), pybind11::arg("start"), pybind11::arg("end"));
		cl.def("query", (struct cppext::BamReader::Iterator (cppext::BamReader::*)(const struct cppext::Region &) const) &cppext::BamReader::query, "C++: cppext::BamReader::query(const struct cppext::Region &) const --> struct cppext::BamReader::Iterator", pybind11::arg("region"));
		cl.def("print_chroms", (void (cppext::BamReader::*)() const) &cppext::BamReader::print_chroms, "C++: cppext::BamReader::print_chroms() const --> void");

		{ // cppext::BamReader::Iterator file:bam.h line:139
			auto & enclosing_class = cl;
			pybind11::class_<cppext::BamReader::Iterator, std::shared_ptr<cppext::BamReader::Iterator>> cl(enclosing_class, "Iterator", "" );
			cl.def( pybind11::init( [](){ return new cppext::BamReader::Iterator(); } ) );
			cl.def_readwrite("is_end_", &cppext::BamReader::Iterator::is_end_);
			cl.def("next", (void (cppext::BamReader::Iterator::*)()) &cppext::BamReader::Iterator::next, "C++: cppext::BamReader::Iterator::next() --> void");
			cl.def("is_end", (bool (cppext::BamReader::Iterator::*)() const) &cppext::BamReader::Iterator::is_end, "C++: cppext::BamReader::Iterator::is_end() const --> bool");
			cl.def("is_reverse", (bool (cppext::BamReader::Iterator::*)() const) &cppext::BamReader::Iterator::is_reverse, "C++: cppext::BamReader::Iterator::is_reverse() const --> bool");
			cl.def("read_name", (std::string (cppext::BamReader::Iterator::*)() const) &cppext::BamReader::Iterator::read_name, "C++: cppext::BamReader::Iterator::read_name() const --> std::string");
			cl.def("chrom", (std::string (cppext::BamReader::Iterator::*)() const) &cppext::BamReader::Iterator::chrom, "C++: cppext::BamReader::Iterator::chrom() const --> std::string");
			cl.def("pos", (long (cppext::BamReader::Iterator::*)() const) &cppext::BamReader::Iterator::pos, "C++: cppext::BamReader::Iterator::pos() const --> long");
			cl.def("end_pos", (long (cppext::BamReader::Iterator::*)() const) &cppext::BamReader::Iterator::end_pos, "C++: cppext::BamReader::Iterator::end_pos() const --> long");
			cl.def("quality", (unsigned char (cppext::BamReader::Iterator::*)() const) &cppext::BamReader::Iterator::quality, "C++: cppext::BamReader::Iterator::quality() const --> unsigned char");
			cl.def("sequence", (std::string (cppext::BamReader::Iterator::*)() const) &cppext::BamReader::Iterator::sequence, "C++: cppext::BamReader::Iterator::sequence() const --> std::string");
			cl.def("cigar_string", (std::string (cppext::BamReader::Iterator::*)() const) &cppext::BamReader::Iterator::cigar_string, "C++: cppext::BamReader::Iterator::cigar_string() const --> std::string");
			cl.def("cigar_length", (unsigned int (cppext::BamReader::Iterator::*)() const) &cppext::BamReader::Iterator::cigar_length, "C++: cppext::BamReader::Iterator::cigar_length() const --> unsigned int");
			cl.def("cigar_buffer", (unsigned int * (cppext::BamReader::Iterator::*)() const) &cppext::BamReader::Iterator::cigar_buffer, "C++: cppext::BamReader::Iterator::cigar_buffer() const --> unsigned int *", pybind11::return_value_policy::automatic);
			cl.def("to_string", (std::string (cppext::BamReader::Iterator::*)() const) &cppext::BamReader::Iterator::to_string, "C++: cppext::BamReader::Iterator::to_string() const --> std::string");
			cl.def("print", (void (cppext::BamReader::Iterator::*)() const) &cppext::BamReader::Iterator::print, "C++: cppext::BamReader::Iterator::print() const --> void");
			cl.def("same_strand_with", (bool (cppext::BamReader::Iterator::*)(bool) const) &cppext::BamReader::Iterator::same_strand_with, "C++: cppext::BamReader::Iterator::same_strand_with(bool) const --> bool", pybind11::arg("is_reversed"));
			cl.def("quality_eq_than", (bool (cppext::BamReader::Iterator::*)(int) const) &cppext::BamReader::Iterator::quality_eq_than, "C++: cppext::BamReader::Iterator::quality_eq_than(int) const --> bool", pybind11::arg("quality_"));
		}

	}
}

// clang-format on
