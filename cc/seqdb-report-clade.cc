#include <iostream>

#include "seqdb.hh"
using namespace seqdb;

// ----------------------------------------------------------------------

int main(int argc, char* const argv[])
{
    try {
        if (argc != 3)
            throw std::runtime_error(std::string{"Usage: "} + argv[0] + " seqdb.json.xz <clade-name, e.g. DEL2017>");
        Seqdb seqdb;
        seqdb.load(argv[1]);
        for (const auto entry: seqdb) {
            if (entry.seq().has_clade(argv[2]))
                std::cout << entry.make_name() << '\n';
        }
        return 0;
    }
    catch (std::exception& err) {
        std::cerr << "ERROR: " << err.what() << '\n';
        return 1;
    }
}

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
