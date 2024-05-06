#include <gtest/gtest.h>
#include "test_example.hpp"
#include "test_export_points_to_paraview.hpp"
#include "test_export_segments_to_paraview.hpp"
#include "test_export_polygons_to_paraview.hpp"
#include "test_export_polyhedra_to_paraview.hpp"
#include "test_create_triangular_mesh.hpp"
#include "test_export_triangular_mesh.hpp"
#include "test_create_tetrahedral_mesh.hpp"
#include "test_export_tetrahedral_mesh.hpp"

int main(int argc, char *argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
