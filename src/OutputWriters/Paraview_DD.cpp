/*
  Copyright (c) 2014 - 2019 University of Bergen
  
  This file is part of the BROOMStyx project.

  BROOMStyx is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  BROOMStyx is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with BROOMStyx.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the AUTHORS file
  for the list of copyright holders.
*/

#include "Paraview_DD.hpp"

#include <chrono>
#include <cstdlib>
#include <stdexcept>
#include <vector>
#include <omp.h>

#include "Core/AnalysisModel.hpp"
#include "Core/ObjectFactory.hpp"
#include "Core/DomainManager.hpp"
#include "Numerics/Numerics.hpp"
#include "Util/readOperations.hpp"

using namespace broomstyx;

registerBroomstyxObject(OutputWriter, Paraview_DD)

// Constructor
Paraview_DD::Paraview_DD()
{
    _pvdFile = nullptr;
    
    // Initialize counters
    _vtuFileCount = 0;
    _writeCounter = 0;
}  

// Destructor
Paraview_DD::~Paraview_DD()
{
    if ( _pvdFile )
    {
        std::fprintf(_pvdFile, "\t</Collection>\n");
        std::fprintf(_pvdFile, "</VTKFile>\n");

        std::fclose(_pvdFile);
    }
}

// Public methods
void Paraview_DD::initialize()
{
    // Create output directory
    int dir_err = system("mkdir -p Output_Paraview");
    if ( dir_err == -1 )
        throw std::runtime_error("Error creating directory 'Output_Paraview'!\n");
    
    // Create .pvd filename
    std::string pvdFilename = "./Output_Paraview/" + _outputFilename + ".pvd";
    
    // Open .pvd file for writing
    _pvdFile = std::fopen(pvdFilename.c_str(), "w");
    
    std::fprintf(_pvdFile, "<?xml version=\"1.0\"?>\n");
    std::fprintf(_pvdFile, "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
    std::fprintf(_pvdFile, "\t<Collection>\n");
}

void Paraview_DD::readDataFrom(FILE* fp)
{
    std::string str, key, src = "Paraview (OutputWriter)";
    
    // Output filename
    verifyDeclaration(fp, key = "FILENAME", src);
    _outputFilename = getStringInputFrom(fp, "Failed to read paraview output filename from input file!", src);
    
    // Point Data
    verifyDeclaration(fp, key = "POINT_DATA", src);
    _nPointData = getIntegerInputFrom(fp, "Failed to read number of point data output from input file!", src);
    _pointData.assign(_nPointData, OutputData());
    
    for ( int i = 0; i < _nPointData; i++)
    {
        str = getStringInputFrom(fp, "Failed to read point data type from input file!", src);
        if ( str == "SCALAR" )
        {
            _pointData[i].dataType = scalar;
            _pointData[i].field.assign(1, 0);
        }
        else if ( str == "VECTOR" )
        {
            _pointData[i].dataType = vector;
            _pointData[i].field.assign(3, 0);
        }
        else if ( str == "TENSOR" )
        {
            _pointData[i].dataType = tensor;
            _pointData[i].field.assign(6, 0);
        }
        else if ( str == "CN_SCALAR" )
        {
            _pointData[i].dataType = cnScalar;
            _pointData[i].field.assign(1, 0);
        }
        else if ( str == "CN_VECTOR" )
        {
            _pointData[i].dataType = cnVector;
            _pointData[i].field.assign(3, 0);
        }
        else if ( str == "CN_TENSOR" )
        {
            _pointData[i].dataType = cnTensor;
            _pointData[i].field.assign(6, 0);
        }
        else
            throw std::runtime_error("Invalid point data type '" + str + "' encountered in input file!\nSource: " + src);
        
        // Data name to use in Paraview
        _pointData[i].name = getStringInputFrom(fp, "Failed reading point data label from input file!", src);
        
        // Nodal fields comprising components of output data
        for (int j = 0; j < (int)_pointData[i].field.size(); j++)
            _pointData[i].field[j] = getIntegerInputFrom(fp, "Failed to read nodal field number for point data from input file.", src);
    }
    
    // Cell Data
    verifyDeclaration(fp, key = "CELL_DATA", src);
    _nCellData = getIntegerInputFrom(fp, "Failed reading number of cell data output from input file!", src);
    
    _cellData.assign(_nCellData, OutputData());
    
    for ( int i = 0; i < _nCellData; i++)
    {
        str = getStringInputFrom(fp, "Failed reading cell data type from input file!", src);
        if ( str == "SCALAR" )
        {
            _cellData[i].dataType = scalar;
            _cellData[i].field.assign(1, 0);
        }
        else if ( str == "VECTOR" )
        {
            _cellData[i].dataType = vector;
            _cellData[i].field.assign(3, 0);
        }
        else if ( str == "TENSOR" )
        {
            _cellData[i].dataType = tensor;
            _cellData[i].field.assign(6, 0);
        }
        else if ( str == "PHYSTAG" )
        {
            _cellData[i].dataType = physTag;
        }
        else
            throw std::runtime_error("Invalid point data type '" + str + "' encountered in input file!\nSource: " + src);
        
        // Data name to use in Paraview
        _cellData[i].name = getStringInputFrom(fp, "Failed reading cell data label from input file!", src);
        
        // Element fields comprising components of output data
        for (int j = 0; j < (int)_cellData[i].field.size(); j++)
            _cellData[i].field[j] = getIntegerInputFrom(fp, "Failed to read cell field number for cell data from input file!", src);
    }
}

void Paraview_DD::writeOutput( double time )
{
    int nCells = analysisModel().domainManager().giveNumberOfDomainCells();
    int nCellNodes = this->giveNumberOfCellNodes();
    
    // build complete string for vtu filename
    std::string vtuFilenameInPvd, vtuFilename;
    vtuFilenameInPvd = _outputFilename + "_" + std::to_string(_vtuFileCount) + ".vtu";
    
    vtuFilename = "./Output_Paraview/" + vtuFilenameInPvd;
    
    // Write .vtu file
    std::printf("\n  %-40s", "Writing results to file ...");
    std::fflush(stdout);
    
    std::chrono::time_point<std::chrono::system_clock> tic, toc;
    std::chrono::duration<double> tictoc;
    tic = std::chrono::high_resolution_clock::now();
    
    FILE *vtuFile = fopen(vtuFilename.c_str(), "w");

    std::fprintf(vtuFile, "<?xml version=\"1.0\"?>\n");
    std::fprintf(vtuFile, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1");
    std::fprintf(vtuFile, "\" byte_order=\"LittleEndian\">\n");
    std::fprintf(vtuFile, "\t<UnstructuredGrid>\n");

    // --- nodal coordinates
    std::fprintf(vtuFile, "\t\t<Piece NumberOfPoints=\"%d\" ", nCellNodes);
    std::fprintf(vtuFile, "NumberOfCells=\"%d\">\n", nCells);
    std::fprintf(vtuFile, "\t\t\t<Points>\n");
    std::fprintf(vtuFile, "\t\t\t\t<DataArray type=\"Float32\" ");
    std::fprintf(vtuFile, "NumberOfComponents=\"3\" format=\"ascii\">\n");

    for ( int i = 0; i < nCells; i++ )
    {
        Cell* curCell = analysisModel().domainManager().giveDomainCell(i);
        std::vector<Node*> curCellNode = analysisModel().domainManager().giveNodesOf(curCell);
        for ( int j = 0; j < (int)curCellNode.size(); j++ )
        {
            RealVector coor = analysisModel().domainManager().giveCoordinatesOf(curCellNode[j]);
            if ( coor.dim() == 2 )
                std::fprintf(vtuFile, "\t\t\t\t\t%25.15e%25.15e%25.15e\n", coor(0), coor(1), 0.0);
            else if ( coor.dim() == 3 )
                std::fprintf(vtuFile, "\t\t\t\t\t%25.15e%25.15e%25.15e\n", coor(0), coor(1), coor(2));
            else {
                std::printf("\n\tERROR: Incorrect size of coordinate vector!");
                std::printf("\n\t--> In Paraview::writeOutput()\n\n");
                std::exit(1);
            }
        }    
    }
    std::fprintf(vtuFile, "\t\t\t\t</DataArray>\n");
    std::fprintf(vtuFile, "\t\t\t</Points>\n");

    // --- cell connectivities
    std::fprintf(vtuFile, "\t\t\t<Cells>\n");
    std::fprintf(vtuFile, "\t\t\t\t<DataArray type=\"Int32\" ");
    std::fprintf(vtuFile, "Name=\"connectivity\" Format=\"ascii\">\n");

    int cellNodeCount = 0;
    for ( int i = 0; i < nCells; i++)
    {
        Cell* curCell = analysisModel().domainManager().giveDomainCell(i);
        int nCellNodes = analysisModel().domainManager().giveNumberOfNodesOf(curCell);

        fprintf(vtuFile, "\t\t\t\t\t");
        for ( int j = 0; j < nCellNodes; j++) {
            std::fprintf(vtuFile, "%10d", cellNodeCount + j);
        }
        std::fprintf(vtuFile, "\n");
        cellNodeCount += nCellNodes;
    }
    std::fprintf(vtuFile, "\t\t\t\t</DataArray>\n");

    // --- cell offsets
    int offset = 0;
    std::fprintf(vtuFile, "\t\t\t\t<DataArray type=\"Int32\" ");
    std::fprintf(vtuFile, "Name=\"offsets\" format=\"ascii\">\n");
    for ( int i = 0; i < nCells; i++) {
        Cell* curCell = analysisModel().domainManager().giveDomainCell(i);
        offset += analysisModel().domainManager().giveNumberOfNodesOf(curCell);
        std::fprintf(vtuFile, "\t\t\t\t\t%d\n", offset);
    }
    std::fprintf(vtuFile, "\t\t\t\t</DataArray>\n");

    // --- cell types
    std::fprintf(vtuFile, "\t\t\t\t<DataArray type=\"UInt8\" Name=\"types\" ");
    std::fprintf(vtuFile, "format=\"ascii\">\n");
    for ( int i = 0; i < nCells; i++)
    {
        Cell* curCell = analysisModel().domainManager().giveDomainCell(i);
        Numerics* numerics = analysisModel().domainManager().giveNumericsFor(curCell);
        int dim = numerics->giveSpatialDimension();
        int nCellNodes = analysisModel().domainManager().giveNumberOfNodesOf(curCell);
        
        int cellType;
        if ( dim == 2 && nCellNodes == 3 )
            cellType = 5;
        else if ( dim == 3 && nCellNodes == 4 )
            cellType = 10;
        else if ( dim == 2 && nCellNodes == 6 )
            cellType = 22;
        else
            throw std::runtime_error("Cells of dim = " + std::to_string(dim) + " and nNodes = " + std::to_string(nCellNodes) + " not yet programmed in Paraview output writer!");
        
        std::fprintf(vtuFile, "\t\t\t\t\t%d\n", cellType);
    }
    std::fprintf(vtuFile, "\t\t\t\t</DataArray>\n");
    std::fprintf(vtuFile, "\t\t\t</Cells>\n");

    // ---- Field Data
    // A. Point Data
    if ( _nPointData > 0 )
    {
        std::fprintf(vtuFile, "\t\t\t<PointData>\n");

        // --- Individual data arrays
        for ( int i = 0; i < _nPointData; i++)
        {
            if ( _pointData[i].dataType == scalar )
            {
                // Scalar data
                std::fprintf(vtuFile, "\t\t\t\t<DataArray type=\"Float32");
                std::fprintf(vtuFile, "\" Name=\"%s\" format=\"ascii\">\n", _pointData[i].name.c_str());
                for ( int j = 0; j < nCells; j++)
                {
                    Cell* curCell = analysisModel().domainManager().giveDomainCell(j);
                    std::vector<Node*> curCellNode = analysisModel().domainManager().giveNodesOf(curCell);
                    
                    for ( int k = 0; k < (int)curCellNode.size(); k++ )
                    {
                        std::fprintf(vtuFile, "\t\t\t\t\t%25.15e\n",
                                analysisModel().domainManager().giveFieldValueAt(curCellNode[k], _pointData[i].field[0]));
                    }
                    
                }
                std::fprintf(vtuFile, "\t\t\t\t</DataArray>\n");
            }
            else if ( _pointData[i].dataType == vector )
            {
                // Vector data
                std::fprintf(vtuFile, "\t\t\t\t<DataArray type=\"Float32");
                std::fprintf(vtuFile, "\" NumberOfComponents=\"3\" Name=");
                std::fprintf(vtuFile, "\"%s\" format=\"ascii\">\n", _pointData[i].name.c_str());
                for ( int j = 0; j < nCells; j++)
                {
                    Cell* curCell = analysisModel().domainManager().giveDomainCell(j);
                    std::vector<Node*> curCellNode = analysisModel().domainManager().giveNodesOf(curCell);
                    
                    for ( int k = 0; k < (int)curCellNode.size(); k++ )
                    {
                        std::fprintf(vtuFile, "\t\t\t\t\t%25.15e%25.15e%25.15e\n",
                                analysisModel().domainManager().giveFieldValueAt(curCellNode[k], _pointData[i].field[0]),
                                analysisModel().domainManager().giveFieldValueAt(curCellNode[k], _pointData[i].field[1]),
                                analysisModel().domainManager().giveFieldValueAt(curCellNode[k], _pointData[i].field[2]));
                    }
                }
                std::fprintf(vtuFile, "\t\t\t\t</DataArray>\n");
            }
            else if ( _pointData[i].dataType == tensor )
            {
                // Tensor data
                std::fprintf(vtuFile, "\t\t\t\t<DataArray type=\"Float32\" ");
                std::fprintf(vtuFile, "NumberOfComponents=\"6\" Name=");
                std::fprintf(vtuFile, "\"%s\" format=\"ascii\">\n", _pointData[i].name.c_str());
                for ( int j = 0; j < nCells; j++)
                {
                    Cell* curCell = analysisModel().domainManager().giveDomainCell(j);
                    std::vector<Node*> curCellNode = analysisModel().domainManager().giveNodesOf(curCell);
                    
                    for ( int k = 0; k < (int)curCellNode.size(); k++ )
                    {
                        std::fprintf(vtuFile, "\t\t\t\t\t%25.15e%25.15e%25.15e%25.15e%25.15e%25.15e\n",
                                analysisModel().domainManager().giveFieldValueAt(curCellNode[k], _pointData[i].field[0]),
                                analysisModel().domainManager().giveFieldValueAt(curCellNode[k], _pointData[i].field[1]),
                                analysisModel().domainManager().giveFieldValueAt(curCellNode[k], _pointData[i].field[2]),
                                analysisModel().domainManager().giveFieldValueAt(curCellNode[k], _pointData[i].field[3]),
                                analysisModel().domainManager().giveFieldValueAt(curCellNode[k], _pointData[i].field[4]),
                                analysisModel().domainManager().giveFieldValueAt(curCellNode[k], _pointData[i].field[5]));
                    }
                }
                std::fprintf(vtuFile, "\t\t\t\t</DataArray>\n");
            }
            else if ( _pointData[i].dataType == cnScalar )
            {
                // discontinuous scalar data
                std::fprintf(vtuFile, "\t\t\t\t<DataArray type=\"Float32");
                std::fprintf(vtuFile, "\" Name=\"%s\" format=\"ascii\">\n", _pointData[i].name.c_str());
                for ( int j = 0; j < nCells; j++)
                {
                    Cell* curCell = analysisModel().domainManager().giveDomainCell(j);
                    Numerics* numerics = analysisModel().domainManager().giveNumericsFor(curCell);
                    RealVector cellNodeVals = numerics->giveCellNodeFieldValuesAt(curCell, _pointData[i].field[0]);
                    
                    for ( int k = 0; k < cellNodeVals.dim(); k++ )
                        std::fprintf(vtuFile, "\t\t\t\t\t%25.15e\n", cellNodeVals(k));
                }
                std::fprintf(vtuFile, "\t\t\t\t</DataArray>\n");
            }
            else if ( _pointData[i].dataType == cnVector )
            {
                // discontinuous scalar data
                std::fprintf(vtuFile, "\t\t\t\t<DataArray type=\"Float32");
                std::fprintf(vtuFile, "\" Name=\"%s\" format=\"ascii\">\n", _pointData[i].name.c_str());
                for ( int j = 0; j < nCells; j++)
                {
                    Cell* curCell = analysisModel().domainManager().giveDomainCell(j);
                    Numerics* numerics = analysisModel().domainManager().giveNumericsFor(curCell);
                    std::vector<RealVector> cellNodeVals;
                    cellNodeVals.assign(3, RealVector());
                    for (  int p = 0; p < 3; p++ )
                        cellNodeVals[p] = numerics->giveCellNodeFieldValuesAt(curCell, _pointData[i].field[p]);
            
                    for (  int k = 0; k < cellNodeVals[0].dim(); k++ )
                        std::fprintf(vtuFile, "\t\t\t\t\t%25.15e%25.15e%25.15e\n", 
                                cellNodeVals[0](k), cellNodeVals[1](k), cellNodeVals[2](k));
                }
                std::fprintf(vtuFile, "\t\t\t\t</DataArray>\n");
            }
            else if ( _pointData[i].dataType == cnTensor )
            {
                // Tensor data
                std::fprintf(vtuFile, "\t\t\t\t<DataArray type=\"Float32\" ");
                std::fprintf(vtuFile, "NumberOfComponents=\"6\" Name=");
                std::fprintf(vtuFile, "\"%s\" format=\"ascii\">\n", _pointData[i].name.c_str());
                for ( int j = 0; j < nCells; j++)
                {
                    Cell* curCell = analysisModel().domainManager().giveDomainCell(j);
                    Numerics* numerics = analysisModel().domainManager().giveNumericsFor(curCell);
                    std::vector<RealVector> cellNodeVals;
                    cellNodeVals.assign(6, RealVector());
                    for (  int p = 0; p < 6; p++ )
                        cellNodeVals[p] = numerics->giveCellNodeFieldValuesAt(curCell, _pointData[i].field[p]);
                    
                    for ( int k = 0; k < cellNodeVals[0].dim(); k++ )
                        std::fprintf(vtuFile, "\t\t\t\t\t%25.15e%25.15e%25.15e%25.15e%25.15e%25.15e\n",
                                cellNodeVals[0](k), cellNodeVals[1](k), cellNodeVals[2](k), 
                                cellNodeVals[3](k), cellNodeVals[4](k), cellNodeVals[5](k));
                }
                std::fprintf(vtuFile, "\t\t\t\t</DataArray>\n");
            }
        }
        std::fprintf(vtuFile, "\t\t\t</PointData>\n");
    }
    
    // B. Cell data
    if ( _nCellData > 0 )
    {
        std::fprintf(vtuFile, "\t\t\t<CellData>\n");

        // --- Individual data arrays
        for ( int i = 0; i < _nCellData; i++)
        {
                // Physical number assignment
            if ( _cellData[i].dataType == physTag )
            {
                std::fprintf(vtuFile, "\t\t\t\t<DataArray type=\"Float32\" ");
                std::fprintf(vtuFile, "Name=\"%s\" format=\"ascii\">\n", _cellData[i].name.c_str());
                for ( int j = 0; j < nCells; j++)
                {
                    Cell* curCell = analysisModel().domainManager().giveDomainCell(j);
                    int label = analysisModel().domainManager().giveLabelOf(curCell);
                    std::fprintf(vtuFile, "\t\t\t\t\t%25.15e\n", (double)label);
                }
                    
                std::fprintf(vtuFile, "\t\t\t\t</DataArray>\n");        
            }
            else if ( _cellData[i].dataType == scalar )
            {
                // Scalar data
                std::fprintf(vtuFile, "\t\t\t\t<DataArray type=\"Float32\" ");
                std::fprintf(vtuFile, "Name=\"%s\" format=\"ascii\">\n", _cellData[i].name.c_str());
                for ( int j = 0; j < nCells; j++)
                {
                    Cell* curCell = analysisModel().domainManager().giveDomainCell(j);
                    int label = analysisModel().domainManager().giveLabelOf(curCell);
                    Numerics* numerics = analysisModel().domainManager().giveNumericsForDomain(label);
                    std::fprintf(vtuFile, "\t\t\t\t\t%25.15e\n", numerics->giveCellFieldValueAt(curCell, _cellData[i].field[0]));
                }
                    
                std::fprintf(vtuFile, "\t\t\t\t</DataArray>\n");
            }
            else if ( _cellData[i].dataType == vector )
            {
                // Vector data
                std::fprintf(vtuFile, "\t\t\t\t<DataArray type=\"Float32\" ");
                std::fprintf(vtuFile, "NumberOfComponents=\"3\" Name=");
                std::fprintf(vtuFile, "\"%s\" format=\"ascii\">\n", _cellData[i].name.c_str());
                
                for ( int j = 0; j < nCells; j++) {
                    Cell* curCell = analysisModel().domainManager().giveDomainCell(j);
                    Numerics* numerics = analysisModel().domainManager().giveNumericsFor(curCell);
                    
                    std::fprintf(vtuFile, "\t\t\t\t\t%25.15e%25.15e%25.15e\n",
                            numerics->giveCellFieldValueAt(curCell, _cellData[i].field[0]),
                            numerics->giveCellFieldValueAt(curCell, _cellData[i].field[1]),
                            numerics->giveCellFieldValueAt(curCell, _cellData[i].field[2]));
                }
                std::fprintf(vtuFile, "\t\t\t\t</DataArray>\n");
            }
            else
            {
                // Tensor data
                std::fprintf(vtuFile, "\t\t\t\t<DataArray type=\"Float32\" ");
                std::fprintf(vtuFile, "NumberOfComponents=\"6\" Name=");
                std::fprintf(vtuFile, "\"%s\" format=\"ascii\">\n", _cellData[i].name.c_str());
                for ( int j = 0; j < nCells; j++)
                {
                    Cell* curCell = analysisModel().domainManager().giveDomainCell(j);
                    Numerics* numerics = analysisModel().domainManager().giveNumericsFor(curCell);
                    
                    std::fprintf(vtuFile, "\t\t\t\t\t%25.15e%25.15e%25.15e%25.15e%25.15e%25.15e\n",
                            numerics->giveCellFieldValueAt(curCell, _cellData[i].field[0]),
                            numerics->giveCellFieldValueAt(curCell, _cellData[i].field[1]),
                            numerics->giveCellFieldValueAt(curCell, _cellData[i].field[2]),
                            numerics->giveCellFieldValueAt(curCell, _cellData[i].field[3]),
                            numerics->giveCellFieldValueAt(curCell, _cellData[i].field[4]),
                            numerics->giveCellFieldValueAt(curCell, _cellData[i].field[5]));
                }
                std::fprintf(vtuFile, "\t\t\t\t</DataArray>\n");
            }
        }
        std::fprintf(vtuFile, "\t\t\t</CellData>\n");
    }
    
    std::fprintf(vtuFile, "\t\t</Piece>\n");
    std::fprintf(vtuFile, "\t</UnstructuredGrid>\n");
    std::fprintf(vtuFile, "</VTKFile>\n");

    std::fclose(vtuFile);
    
    // Increment vtu file count
    _vtuFileCount += 1;

    // Write entry for .vtu file in .pvd file
    // --------------------------------------
    std::fprintf(_pvdFile, "\t\t<DataSet timestep=\"%.15f\" file=\"%s\"/>\n", time, vtuFilenameInPvd.c_str());
    std::fflush(_pvdFile);
    
    toc = std::chrono::high_resolution_clock::now();
    tictoc = toc - tic;
    std::printf("done (time = %f sec.)\n", tictoc.count());
    std::string fullVtuFilename = "./Output_Paraview/" + vtuFilenameInPvd;
    std::printf("  --> %s\n", fullVtuFilename.c_str());
    
    _writeCounter++;
}

// Private methods
// -------------------------------------------------------------------------------
int Paraview_DD::giveNumberOfCellNodes()
{
    int nDomainCells = analysisModel().domainManager().giveNumberOfDomainCells();

    int nCellNodes = 0;
#ifdef _OPENMP
#pragma omp parallel for reduction ( + : nCellNodes )
#endif
    for ( int i = 0; i < nDomainCells; i++ )
    {
        Cell* curCell = analysisModel().domainManager().giveDomainCell(i);
        nCellNodes += analysisModel().domainManager().giveNumberOfNodesOf(curCell);
    }

    return nCellNodes;
}