#include "Matrix.hpp"
#include "LowRank.hpp"
#include "HODLR.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

namespace py = pybind11;

class PyMatrix : public Matrix
{
public:
    /* Inherit the constructors */
    using Matrix::Matrix;

    /* Trampoline (need one for each virtual function) */
    dtype getMatrixEntry(int i, int j) override 
    {
        PYBIND11_OVERLOAD_PURE(dtype,          /* Return type */
                               Matrix,         /* Parent class */
                               getMatrixEntry, /* Name of function in C++ (must match Python name) */
                               i, j            /* Arguments */
                              );
    }
};

PYBIND11_MODULE(pyhodlrlib, m) 
{
    m.doc() = "This is the Python Wrapper to HODLRlib";

    py::class_<Matrix, PyMatrix> matrix(m, "Matrix");
    matrix
        .def(py::init<int>())
        .def("getMatrixEntry", &Matrix::getMatrixEntry)
        .def("getMatrix", &Matrix::getMatrix);

    py::class_<LowRank> lowrank(m, "LowRank");
    lowrank
        .def(py::init<Matrix*, std::string>())
        .def("factorize", &LowRank::factorize)
        .def("getL", &LowRank::getL)
        .def("getR", &LowRank::getR);

    py::class_<HODLR> hodlr(m, "HODLR");
    hodlr
        .def(py::init<int, int, double>())
        .def("assemble", &HODLR::assemble)
        .def("matmatProduct", &HODLR::matmatProduct)
        .def("factorize", &HODLR::factorize)
        .def("solve", &HODLR::solve)
        .def("symmetricFactorProduct", &HODLR::symmetricFactorProduct)
        .def("symmetricFactorTransposeProduct", &HODLR::symmetricFactorTransposeProduct)
        .def("getSymmetricFactor", &HODLR::getSymmetricFactor)
        .def("logDeterminant", &HODLR::logDeterminant);
}
