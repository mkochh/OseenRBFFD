#include "polyhedron_integration.h"

std::vector<std::vector<double>> computeIntegrals(
  const Polyhedron& polyhedron, size_t max_order)
{
  // Create data structures to store integrals ////////////////////////////////
  std::vector<std::vector<double>> integrals(max_order + 1);
  auto number_of_verts = polyhedron.coords.size();
  auto number_of_edges = polyhedron.edges.size();
  auto number_of_faces = polyhedron.faces.size();
  std::vector<std::vector<std::vector<double>>> edge_ints(max_order + 1);
  std::vector<std::vector<std::vector<double>>> face_ints(max_order + 1);
  for (size_t o{}; o < max_order + 1; ++o)
  {
    auto number_of_integrals = (o + 2) * (o + 1) / 2;
    integrals[o].resize(number_of_integrals);
    edge_ints[o].resize(number_of_integrals);
    face_ints[o].resize(number_of_integrals);
    for (size_t i{}; i < number_of_integrals; ++i)
    {
      edge_ints[o][i].resize(number_of_edges);
      face_ints[o][i].resize(number_of_faces);
    }
  }
  // Edge quantities //////////////////////////////////////////////////////////
  // vector tangent to edge
  std::vector<std::array<double, 3>> edges_t(number_of_edges);
  // signed distance from 1st vert to x0
  std::vector<double> edges_d1(number_of_edges);
  // signed distance from 2nd vert to x0
  std::vector<double> edges_d2(number_of_edges);
  // vector from origin to closest point on edge
  std::vector<std::array<double, 3>> edges_x0(number_of_edges);
  // Edge quantity computations ///////////////////////////////////////////////
  for (size_t e{}; e < number_of_edges; ++e)
  {
    auto edge_t_norm = 0.0;
    for (size_t i{}; i < 3; ++i)
    {
      edges_t[e][i] = polyhedron.coords[polyhedron.edges[e][0]][i]
        - polyhedron.coords[polyhedron.edges[e][1]][i];
      edge_t_norm = edge_t_norm + edges_t[e][i] * edges_t[e][i];
    }
    edge_t_norm = std::sqrt(edge_t_norm);
    for (size_t i{}; i < 3; ++i)
    {
      edges_t[e][i] = edges_t[e][i] / edge_t_norm;
      edges_d1[e] = edges_d1[e]
        + edges_t[e][i] * polyhedron.coords[polyhedron.edges[e][0]][i];
      edges_d2[e] = edges_d2[e]
        - edges_t[e][i] * polyhedron.coords[polyhedron.edges[e][1]][i];
    }
    for (size_t i{}; i < 3; ++i)
    {
      edges_x0[e][i] = polyhedron.coords[polyhedron.edges[e][0]][i]
        - edges_d1[e] * edges_t[e][i];
    }
  }
  // Face quantities //////////////////////////////////////////////////////////
  // signed distance from x0 on face to edge line
  std::vector<std::vector<double>> faces_edges_d(number_of_faces);
  // vector from origin to closest point on face
  std::vector<std::array<double, 3>> faces_x0(number_of_faces);
  // signed distance from origin to face
  std::vector<double> faces_d(number_of_faces);
  // Face quantity computations ///////////////////////////////////////////////
  for (size_t f{}; f < number_of_faces; ++f)
  {
    auto number_of_face_edges = polyhedron.faces[f].size();
    std::vector<std::array<double, 3>> face_edges_t(number_of_face_edges);
    for (size_t e{}; e < number_of_face_edges; ++e)
    {
      auto v_curr = polyhedron.edges[polyhedron.faces[f][e]][0];
      auto v_next = polyhedron.edges[polyhedron.faces[f][e]][1];
      faces_x0[f][0] = faces_x0[f][0] +
        static_cast<double>(polyhedron.faces_dir[f][e])
        * (polyhedron.coords[v_curr][1] - polyhedron.coords[v_next][1])
        * (polyhedron.coords[v_curr][2] + polyhedron.coords[v_next][2]);
      faces_x0[f][1] = faces_x0[f][1] +
        static_cast<double>(polyhedron.faces_dir[f][e])
        * (polyhedron.coords[v_curr][2] - polyhedron.coords[v_next][2])
        * (polyhedron.coords[v_curr][0] + polyhedron.coords[v_next][0]);
      faces_x0[f][2] = faces_x0[f][2] +
        static_cast<double>(polyhedron.faces_dir[f][e])
        * (polyhedron.coords[v_curr][0] - polyhedron.coords[v_next][0])
        * (polyhedron.coords[v_curr][1] + polyhedron.coords[v_next][1]);
      face_edges_t[e][0] = static_cast<double>(polyhedron.faces_dir[f][e])
        * (polyhedron.coords[v_next][0] - polyhedron.coords[v_curr][0]);
      face_edges_t[e][1] = static_cast<double>(polyhedron.faces_dir[f][e])
        * (polyhedron.coords[v_next][1] - polyhedron.coords[v_curr][1]);
      face_edges_t[e][2] = static_cast<double>(polyhedron.faces_dir[f][e])
        * (polyhedron.coords[v_next][2] - polyhedron.coords[v_curr][2]);
      auto face_edge_t_norm = std::sqrt(face_edges_t[e][0] * face_edges_t[e][0]
        + face_edges_t[e][1] * face_edges_t[e][1]
        + face_edges_t[e][2] * face_edges_t[e][2]);
      face_edges_t[e][0] /= face_edge_t_norm;
      face_edges_t[e][1] /= face_edge_t_norm;
      face_edges_t[e][2] /= face_edge_t_norm;
    }
    auto face_x0_norm = std::sqrt(faces_x0[f][0] * faces_x0[f][0]
      + faces_x0[f][1] * faces_x0[f][1] + faces_x0[f][2] * faces_x0[f][2]);
    faces_x0[f][0] /= face_x0_norm;
    faces_x0[f][1] /= face_x0_norm;
    faces_x0[f][2] /= face_x0_norm;
    std::vector<std::array<double, 3>> face_edges_n(number_of_face_edges);
    for (size_t e{}; e < number_of_face_edges; ++e)
    {
      face_edges_n[e][0] = face_edges_t[e][1] * faces_x0[f][2]
        - face_edges_t[e][2] * faces_x0[f][1];
      face_edges_n[e][1] = face_edges_t[e][2] * faces_x0[f][0]
        - face_edges_t[e][0] * faces_x0[f][2];
      face_edges_n[e][2] = face_edges_t[e][0] * faces_x0[f][1]
        - face_edges_t[e][1] * faces_x0[f][0];
    }
    auto v_0 = polyhedron.edges[polyhedron.faces[f][0]][0];
    faces_d[f] = faces_x0[f][0] * polyhedron.coords[v_0][0]
      + faces_x0[f][1] * polyhedron.coords[v_0][1]
      + faces_x0[f][2] * polyhedron.coords[v_0][2];
    faces_x0[f][0] *= faces_d[f];
    faces_x0[f][1] *= faces_d[f];
    faces_x0[f][2] *= faces_d[f];
    faces_edges_d[f].resize(number_of_face_edges);
    for (size_t e{}; e < number_of_face_edges; ++e)
    {
      auto e_curr = polyhedron.faces[f][e];
      faces_edges_d[f][e] =
        (edges_x0[e_curr][0] - faces_x0[f][0]) * face_edges_n[e][0]
        + (edges_x0[e_curr][1] - faces_x0[f][1]) * face_edges_n[e][1]
        + (edges_x0[e_curr][2] - faces_x0[f][2]) * face_edges_n[e][2];
    }
  }
  // Loop over monomial orders ////////////////////////////////////////////////
  for (int o{}; o < static_cast<int>(max_order) + 1; ++o)
  {
    // Loop over powers of x //////////////////////////////////////////////////
    for (int xp{o}; xp >= 0; --xp)
    {
      // Loop over powers of z ////////////////////////////////////////////////
      for (int zp{}; zp < (o - xp + 1); ++zp)
      {
        auto entry = (o - xp) * (o - xp + 1) / 2 + zp;
        auto yp = o - xp - zp;
        // Compute monomial values at vertices ////////////////////////////////
        std::vector<double> vert_vals(number_of_verts);
        for (size_t v{}; v < polyhedron.coords.size(); ++v)
        {
          vert_vals[v] =
            std::pow(polyhedron.coords[v][0], static_cast<double>(xp)) *
            std::pow(polyhedron.coords[v][1], static_cast<double>(yp)) *
            std::pow(polyhedron.coords[v][2], static_cast<double>(zp));
        }
        // Integrate monomial on edges ////////////////////////////////////////
        for (size_t e{}; e < number_of_edges; ++e)
        {
          edge_ints[o][entry][e] =
            (vert_vals[polyhedron.edges[e][0]] * edges_d1[e]
            + vert_vals[polyhedron.edges[e][1]] * edges_d2[e])
            / (static_cast<double>(o) + 1.0);
        }
        if (xp > 0)
        {
          auto prev_entry = entry;
          for (size_t e{}; e < number_of_edges; ++e)
          {
            edge_ints[o][entry][e] += static_cast<double>(xp) * edges_x0[e][0]
              * edge_ints[o - 1][prev_entry][e]
              / (static_cast<double>(o) + 1.0);
          }
        }
        if (yp > 0)
        {
          auto prev_entry = entry - o + xp;
          for (size_t e{}; e < number_of_edges; ++e)
          {
            edge_ints[o][entry][e] += static_cast<double>(yp) * edges_x0[e][1]
              * edge_ints[o - 1][prev_entry][e]
              / (static_cast<double>(o) + 1.0);
          }
        }
        if (zp > 0)
        {
          auto prev_entry = entry - o + xp - 1;
          for (size_t e{}; e < number_of_edges; ++e)
          {
            edge_ints[o][entry][e] += static_cast<double>(zp) * edges_x0[e][2]
              * edge_ints[o - 1][prev_entry][e]
              / (static_cast<double>(o) + 1.0);
          }
        }
        // Integrate monomial on faces ////////////////////////////////////////
        for (size_t f{}; f < number_of_faces; ++f)
        {
          auto number_of_face_edges = polyhedron.faces[f].size();
          for (size_t e{}; e < number_of_face_edges; ++e)
          {
            auto e_curr = polyhedron.faces[f][e];
            face_ints[o][entry][f] +=
              faces_edges_d[f][e] * edge_ints[o][entry][e_curr]
              / (static_cast<double>(o) + 2.0);
          }
        }
        if (xp > 0)
        {
          auto prev_entry = entry;
          for (size_t f{}; f < number_of_faces; ++f)
          {
            face_ints[o][entry][f] +=
              static_cast<double>(xp) * faces_x0[f][0]
              * face_ints[o - 1][prev_entry][f]
              / (static_cast<double>(o) + 2.0);
          }
        }
        if (yp > 0)
        {
          auto prev_entry = entry - o + xp;
          for (size_t f{}; f < number_of_faces; ++f)
          {
            face_ints[o][entry][f] +=
              static_cast<double>(yp) * faces_x0[f][1]
              * face_ints[o - 1][prev_entry][f]
              / (static_cast<double>(o) + 2.0);
          }
        }
        if (zp > 0)
        {
          auto prev_entry = entry - o + xp - 1;
          for (size_t f{}; f < number_of_faces; ++f)
          {
            face_ints[o][entry][f] +=
              static_cast<double>(zp) * faces_x0[f][2]
              * face_ints[o - 1][prev_entry][f]
              / (static_cast<double>(o) + 2.0);
          }
        }
        // Integrate monomial over polyhedron /////////////////////////////////
        for (size_t f{}; f < number_of_faces; ++f)
        {
          integrals[o][entry] += faces_d[f] * face_ints[o][entry][f]
            / (static_cast<double>(o) + 3.0);
        }
      }
    }
  }
  return integrals;
}

void edge_check(unsigned int a, unsigned int b, short &sign, int &idx, std::vector<std::array<unsigned int, 2>> &edges) {
    if(a<b) {
        sign = 1;
    } else {
        sign = -1;
        std::swap(a,b);
    }
    idx = -1;
    std::array<unsigned int, 2> edge_candidate = {a,b};
    for(int i = 0; i < static_cast<int>(edges.size()); i++) {
        if(edge_candidate==edges[i]) {
            idx = i;
        }
    }
    if (idx==-1) {
        edges.push_back(edge_candidate);
        idx = edges.size() - 1;
    }
}

Polyhedron OFF2polyInt(std::string file) {
    int n_coords = 0;
    int n_faces = 0;
    int n_edges = 0;

    Polyhedron polyhedron;

    int verts_per_face;

    std::ifstream infile(file);
    std::string line;

    int i = 0;
    while (std::getline(infile,line))
    {
        std::istringstream iss(line);
        if (i==1)
        {
            if(!(iss >> n_coords >> n_faces >> n_edges)) {break;}
            n_edges = n_coords + n_faces - 2; // n_edges usually not given in .off, calculate with euler polytope theorem
            break;
        }
        i++;
    }

    i = 0;
    while (std::getline(infile,line) && i<n_coords)
    {
        std::istringstream iss(line);
        double a, b, c;
        if (i<n_coords)
        {
            if(!(iss >> a >> b >> c))
            {std::cout << "error" << std::endl; break;}
            polyhedron.coords.push_back({a,b,c});
        }
        if(i == n_coords-1) { break; }
        i++;
    }

    i = 0;
    while (std::getline(infile,line))
    {
        std::istringstream iss(line);
        unsigned int a, b, c; int idx;
        short sign;
        if (i<n_faces)
        {
            std::vector<unsigned int> curr_face;
            std::vector<short> curr_face_dir;
            if(!(iss >> verts_per_face >> a >> b >> c))
            {std::cout << "error" << std::endl; break;}
            edge_check(a,b,sign,idx,polyhedron.edges); curr_face.push_back(idx); curr_face_dir.push_back(sign);
            edge_check(b,c,sign,idx,polyhedron.edges); curr_face.push_back(idx); curr_face_dir.push_back(sign);
            edge_check(c,a,sign,idx,polyhedron.edges); curr_face.push_back(idx); curr_face_dir.push_back(sign);
            polyhedron.faces.push_back(curr_face);
            polyhedron.faces_dir.push_back(curr_face_dir);
        }

        i++;
    }

    return polyhedron;
}

Eigen::VectorXd compute_pressure_constraint(std::string OFF_file, const mm::DomainDiscretization<mm::Vec3d>& d_p, int max_order)
{
    Eigen::VectorXd constraint_old = Eigen::VectorXd::Constant(d_p.size(),1.0);
    Eigen::VectorXd constraint_new = Eigen::VectorXd::Constant(d_p.size(),0.0);
    int current_order = 0; 
    bool non_negative = true; 
    Polyhedron polyhedron = OFF2polyInt(OFF_file);

    do
    {
        // set up least squares
        mm::Monomials<mm::Vec3d> mon(current_order);
        if (mon.size() > d_p.size())
            break;

        mm::WLS<mm::Monomials<mm::Vec3d>, mm::NoWeight<mm::Vec3d>, mm::NoScale, Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd>> appr(mon);

        // compute integrals of polynomials over polyhedron domain
        auto integrals = computeIntegrals(polyhedron, current_order);

        // set up rhs
        Eigen::VectorXd rhs_powers(mon.size());
        Eigen::MatrixXi Powers = mon.powers();
        int order, entry;
        for (int i = 0; i < mon.size(); i++) {
            order = Powers(0,i) + Powers(1,i) + Powers(2,i);
            entry = (order - Powers(0,i))*(order - Powers(0,i) +1) / 2 + Powers(2,i);
            rhs_powers(i) = integrals[order][entry];
        }

        // assemble Vandermonde matrix
        appr.compute({0.0, 0.0, 0.0}, d_p.positions());

        constraint_new = appr.getShapeFromRhs(rhs_powers);

        for (int i = 0; i < constraint_new.size(); i++) {
            if (constraint_new(i) < -Eigen::NumTraits<double>::epsilon()) {
                non_negative = false;
                break;
            }
        }

        if (non_negative == true)
            constraint_old = constraint_new;

        current_order++;
        
    } while (non_negative && (current_order < max_order));

    return constraint_old;
}
