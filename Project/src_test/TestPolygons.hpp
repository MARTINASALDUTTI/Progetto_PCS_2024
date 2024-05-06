#pragma once

#include "MatlabDataArray.hpp"
#include "Eigen/Eigen"
#include <iostream>
#include <vector>
#include "MatlabUtilities.hpp"
#include "mesh_data.hpp"
#include "MeshUtilities.hpp"
#include "MeshMatricesDAO.hpp"

namespace MatlabFunctions
{
  struct test_struct
  {
      Eigen::MatrixXi a;
      std::string b;
  };

  test_struct ConvertStructFROMMatlabData(const matlab::data::StructArray& MATLABStruct)
    {
      const matlab::data::TypedArray<double> matrix = MATLABStruct[0]["a"];

      return
      {
        MatlabFunctions::ConvertMatlabMatrixToEigenMatrix<int, double>(matrix),
            MatlabFunctions::ConvertMatlabStringToString(MATLABStruct[0]["b"])[0]
      };
    }

  matlab::data::StructArray ConvertToMatlabStruct(const test_struct& StructToMatlab)
    {
      matlab::data::ArrayFactory factory;
      matlab::data::StructArray structArray = factory.createStructArray({ 1 }, { "a", "b" });

      structArray[0]["a"] = MatlabFunctions::ConvertEigenMatrixToMatlabMatrix<int, double>(StructToMatlab.a);
      structArray[0]["b"] = MatlabFunctions::ConvertStringToMatlabString({ StructToMatlab.b });

      return structArray;
    }

  matlab::data::StructArray ConvertMesh2DToMatlabStruct(const Functions::Mesh2D& StructToMatlab)
  {
        matlab::data::ArrayFactory factory;
        matlab::data::StructArray structArray = factory.createStructArray({ 1 }, { "Cell0Ds", "Cell1Ds", "Cell2Ds"});

        structArray[0]["Cell0Ds"] = MatlabFunctions::ConvertEigenMatrixToMatlabMatrix<double, double>(StructToMatlab.Cell0Ds);
        structArray[0]["Cell1Ds"] = MatlabFunctions::ConvertEigenMatrixToMatlabMatrix<int, double>(StructToMatlab.Cell1Ds);
        structArray[0]["Cell2Ds"] = MatlabFunctions::ConvertVectorOfVectorToMatlabCell(StructToMatlab.Cell2Ds);

        return structArray;
  }

  matlab::data::StructArray ConvertMesh3DToMatlabStruct(const Functions::Mesh3D& StructToMatlab)
  {
      matlab::data::ArrayFactory factory;
      matlab::data::StructArray structArray = factory.createStructArray({ 1 }, { "Cell0Ds", "Cell1Ds", "Cell2Ds", "Cell3DsVertices", "Cell3DsEdges", "Cell3DsFaces"});

      structArray[0]["Cell0Ds"] = MatlabFunctions::ConvertEigenMatrixToMatlabMatrix<double, double>(StructToMatlab.Cell0Ds);
      structArray[0]["Cell1Ds"] = MatlabFunctions::ConvertEigenMatrixToMatlabMatrix<int, double>(StructToMatlab.Cell1Ds);
      structArray[0]["Cell2Ds"] = MatlabFunctions::ConvertVectorOfVectorToMatlabCell(StructToMatlab.Cell2Ds);
      structArray[0]["Cell3DsVertices"] = MatlabFunctions::ConvertVectorOfVectorToMatlabCell(StructToMatlab.Cell3DsVertices);
      structArray[0]["Cell3DsEdges"] = MatlabFunctions::ConvertVectorOfVectorToMatlabCell(StructToMatlab.Cell3DsEdges);
      structArray[0]["Cell3DsFaces"] = MatlabFunctions::ConvertVectorOfVectorToMatlabCell(StructToMatlab.Cell3DsFaces);

      return structArray;
  }

  Functions::Mesh2D ConvertMatlabStructToMesh2D(matlab::data::StructArray& MATLABStruct)
  {
      const matlab::data::TypedArray<double> point_matrix = MATLABStruct[0]["Cell0Ds"];
      const matlab::data::TypedArray<double> segment_matrix = MATLABStruct[0]["Cell1Ds"];
      const matlab::data::CellArray polygons_cell = MATLABStruct[0]["Cell2Ds"];

      Functions::Mesh2D mesh;
      mesh.Cell0Ds = MatlabFunctions::ConvertMatlabMatrixToEigenMatrix<double, double>(point_matrix);
      mesh.Cell1Ds = MatlabFunctions::ConvertMatlabMatrixToEigenMatrix<int, double>(segment_matrix);
      mesh.Cell2Ds = MatlabFunctions::ConvertCellToVectorOfVector(polygons_cell);

      return mesh;
  }

  Functions::Mesh3D ConvertMatlabStructToMesh3D(matlab::data::StructArray& MATLABStruct)
  {
      const matlab::data::TypedArray<double> point_matrix = MATLABStruct[0]["Cell0Ds"];
      const matlab::data::TypedArray<double> segment_matrix = MATLABStruct[0]["Cell1Ds"];
      const matlab::data::CellArray polygons_cell = MATLABStruct[0]["Cell2Ds"];
      const matlab::data::CellArray Vertices_cell3Ds = MATLABStruct[0]["Cell3DsVertices"];
      const matlab::data::CellArray Edges_cell3Ds = MATLABStruct[0]["Cell3DsEdges"];
      const matlab::data::CellArray Faces_cell3Ds = MATLABStruct[0]["Cell3DsFaces"];

      Functions::Mesh3D mesh;
      mesh.Cell0Ds = MatlabFunctions::ConvertMatlabMatrixToEigenMatrix<double, double>(point_matrix);
      mesh.Cell1Ds = MatlabFunctions::ConvertMatlabMatrixToEigenMatrix<int, double>(segment_matrix);
      mesh.Cell2Ds = MatlabFunctions::ConvertCellToVectorOfVector(polygons_cell);
      mesh.Cell3DsVertices = MatlabFunctions::ConvertCellToVectorOfVector(Vertices_cell3Ds);
      mesh.Cell3DsEdges = MatlabFunctions::ConvertCellToVectorOfVector(Edges_cell3Ds);
      mesh.Cell3DsFaces = MatlabFunctions::ConvertCellToVectorOfVector(Faces_cell3Ds);

      return mesh;
  }

  matlab::data::StructArray ConvertPolyhedronToMatlabStruct(const Gedim::GeometryUtilities::Polyhedron& cube)
  {
      matlab::data::ArrayFactory factory;
      matlab::data::StructArray structArray = factory.createStructArray({ 1 }, { "Vertices", "Edges", "Faces"});

      structArray[0]["Vertices"] = MatlabFunctions::ConvertEigenMatrixToMatlabMatrix<double, double>(cube.Vertices);
      structArray[0]["Edges"] = MatlabFunctions::ConvertEigenMatrixToMatlabMatrix<int, double>({ cube.Edges });
      structArray[0]["Faces"] = MatlabFunctions::ConvertVectorEigenMatrixToMatlabCell<int,double >({ cube.Faces });

      return structArray;
  }
}
