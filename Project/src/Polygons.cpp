/* ExportPolyhedraToParaview
 * Export Polyhedra to paraview
 * [] = ExportPoligonsToParaview(points, polyhedra_vertices, file_path);
*/

#include "export_polyhedra_to_paraview.hpp"
#include "mex.hpp"
#include "mexAdapter.hpp"
#include "MatlabUtilities.hpp"
#include <iostream>
#include "Eigen/Eigen"

class MexFunction : public matlab::mex::Function
{
private:
    unsigned int counter;

    // Pointer to MATLAB engine
    std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();

    // Factory (modello di progettazione che crea oggetti senza specificarne la classe) to create MATLAB data arrays
    matlab::data::ArrayFactory factory;

    void displayOnMATLAB(const std::string& stream)
    {
        matlabPtr->feval(u"fprintf", 0,
                         std::vector<matlab::data::Array>({ factory.createScalar(stream) }));
    }

    void errorOnMATLAB(const std::string& stream)
    {
        matlabPtr->feval(u"error",
                         0,
                         std::vector<matlab::data::Array>({ factory.createScalar(stream) }));
    }

public:
    MexFunction()
    {
        counter = 0;
    }

    ~MexFunction()
    {
    }

    void operator() (matlab::mex::ArgumentList outputs,
                    matlab::mex::ArgumentList inputs)
    {
        checkArguments(outputs, inputs);

        matlab::data::TypedArray<double> points = std::move(inputs[0]);
        matlab::data::CellArray polyhedra_vertices = std::move(inputs[1]);
        matlab::data::TypedArray<matlab::data::MATLABString> file_path= std::move(inputs[2]);

        std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();

        Functions::export_polyhedra_to_paraview(MatlabFunctions::ConvertMatlabMatrixToEigenMatrix<double, double>(points),
                                               MatlabFunctions::ConvertCellToVectorOfVector(polyhedra_vertices),
                                               MatlabFunctions::ConvertMatlabStringToString(file_path)[0]);
    }

    void checkArguments(matlab::mex::ArgumentList& outputs,
                        matlab::mex::ArgumentList& inputs)
    {
        //check the number of input
        if (inputs.size()!=3)
            throw std::runtime_error("the function needs three inputs");

        // Check the first input argument: points
        if (inputs[0].getType()!= matlab::data::ArrayType::DOUBLE ||
            inputs[0].getDimensions()[0]!=3      )
            throw std::runtime_error("the first input must be a double array with three row");

        // Check the second input argument: polyhedra vertices
        if (inputs[1].getType() != matlab::data::ArrayType::CELL)
        {
            throw std::runtime_error("the second input must be a cell");
        }
        else
        {
            size_t rows = inputs[1].getDimensions()[0];
            for (size_t i = 0; i < rows; ++i)
            {
                matlab::data::TypedArray<double> cellArray = inputs[1][i][0];
                size_t DimInnerRows = cellArray.getDimensions()[1];
                if (DimInnerRows < 4)
                {
                    throw std::runtime_error("the polyhedra must have at least four vertices");
                }
            }
        }

        if (inputs[2].getType()!= matlab::data::ArrayType::MATLAB_STRING)
            throw std::runtime_error("the third input must be a string");

        //check outputs
        if (outputs.size()!=0)
            throw std::runtime_error("the function doesn't return an output");
    }
};
