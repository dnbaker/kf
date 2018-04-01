#include "pybind11/pybind11.h"
#include "include/kfreq.h"

namespace py = pybind11;
using namespace kf;
using namespace kf::freq;

using kfs_t = KFreqArray<size_t>;
PYBIND11_MODULE(kf, m) {
    m.doc() = "kmer histogram calculator";
    py::class_<kfs_t> (m, "kf")
        .def(py::init<unsigned>())
        .def(py::init<const char *>())
        .def("clear", &kfs_t::clear, "Clear all entries.")
        .def("count", (size_t (kfs_t::*)(const std::string &)) &kfs_t::count, "Count stuff.")
        .def("count", (size_t (kfs_t::*)(unsigned, std::size_t)) &kfs_t::count, "Count stuff.");
    m.def("str2kmer", &str2kmer<std::size_t>, "Convert a string into an index. Throws std::runtime_error if it contains illegal characters.");
}
