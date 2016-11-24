#include "acmacs-base/pybind11.hh"

// ----------------------------------------------------------------------

PYBIND11_PLUGIN(seqdb_backend)
{
    py::module m("seqdb_backend", "Seqdb access plugin");

      // ----------------------------------------------------------------------

    return m.ptr();
}

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
