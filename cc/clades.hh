#pragma once

#include <string>
#include <vector>

#include "sequence-shift.hh"

// ----------------------------------------------------------------------

std::vector<std::string> clades_b_yamagata(std::string aSequence, Shift aShift);
std::vector<std::string> clades_h1pdm(std::string aSequence, Shift aShift);
std::vector<std::string> clades_h3n2(std::string aSequence, Shift aShift);

// ----------------------------------------------------------------------
