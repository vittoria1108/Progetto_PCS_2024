#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include "DFN.hpp"

using namespace std;

namespace FractureLibrary{

bool ImportFracture(const string& filename, DFN& dfn);

}
