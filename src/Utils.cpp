#include "Utils.hpp"

#include <iostream>
#include <filesystem>
namespace fs = std::filesystem;

namespace PolyhedronLibrary{

bool Import_solid(string &solid_name, Polyhedron &P)
{
    if(!import_cell0Ds(solid_name, P)) return false;

    if(!import_cell1Ds(solid_name, P)) return false;

    if(!import_cell2Ds(solid_name, P)) return false;

    return true;

}

}