#ifndef POLYHEDRON_INTEGRATION_HEADER
#define POLYHEDRON_INTEGRATION_HEADER

#include <array>
#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <medusa/Medusa_fwd.hpp>
#include "../../medusa/include/medusa/bits/approximations/WLS.hpp"

/*****************************************************************************\
|* Authors: E. B. Chin and N. Sukumar                                        *|
|* Date   : May 2020                                                         *|
|* This code is provided as Supplementary Material to our paper contribution *|
|* in Computer Aided Geometric Design, 2020                                  *|
\*****************************************************************************/
/*****************************************************************************\
|* Polyhedron: struct defining a polyhedron in a graph                       *|
|*  - coords: (# of coords) vector holding vertex coordinates as an array    *|
|*  - edges: (# of edges) vector holding index vertex values defining        *|
|*           ordered edge connectivity                                       *|
|*  - faces: (# of faces) vector holding index edge values defining          *|
|*                the face connectivity in counterclockwise orientation      *|
|*  - face_dir: (# of faces) vector defining the direction of edges; +1 for  *|
|*              edges in the correct orientation and -1 for edges which must *|
|*              be reversed                                                  *|
\*****************************************************************************/
struct Polyhedron
{
  std::vector<std::vector<unsigned int>> faces;
  std::vector<std::vector<short>> faces_dir;
  std::vector<std::array<unsigned int, 2>> edges;
  std::vector<std::array<double, 3>> coords;
};

/*****************************************************************************\
|* computeIntegrals: integrates monomials over a polyhedron                  *|
|*  Input:                                                                   *|
|*   - polyhedron: Polyhedron defining the domain of integration             *|
|*   - max_order: all monomials <= max_order will be integrated              *|
|*  Output:                                                                  *|
|*   - integrals: Integral values stored in an Integrals struct              *|
|*  Notes:                                                                   *|
|*   - the closest point to the origin on edges and faces are treated as x_0 *|
\*****************************************************************************/
std::vector<std::vector<double>> computeIntegrals(const Polyhedron& polyhedron, size_t max_order);

// auxiliary function to convert OFF to format needed for computeIntegrals
void edge_check(int a, int b, short &sign, int &idx, std::vector<std::array<unsigned int, 2>> &edges);

// convert OFF format to the format which is necessary for computeIntegrals function
Polyhedron OFF2polyInt(std::string file);

// maxOrder = 12 works well in experiments
Eigen::VectorXd compute_pressure_constraint(std::string OFF_file, const mm::DomainDiscretization<mm::Vec3d>& d_p, int max_order);

#endif
