#pragma once

#include "Eigen/Eigen"
#include "MeshUtilities.hpp"
#include "MeshMatricesDAO.hpp"
#include "mesh_data.hpp"


namespace Functions
{
  // ***************************************************************************
inline Gedim::GeometryUtilities::Polyhedron create_cube(const Eigen::Vector3d origin,
                                 const double edgeLength)
  {
      Gedim::GeometryUtilitiesConfig geometryConfing;
      geometryConfing.Tolerance1D = 1.0e-8;
      Gedim::GeometryUtilities geometryUtilities(geometryConfing);

      const auto cube = geometryUtilities.CreateCubeWithOrigin(origin,
                                                               edgeLength);

    return cube;
  }
  // ***************************************************************************
}
