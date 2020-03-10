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

#include "OutputManager.hpp"
#include <chrono>
#include <stdexcept>
#include <string>

#include "ObjectFactory.hpp"
#include "Diagnostics.hpp"
#include "OutputQuantities/OutputQuantity.hpp"
#include "OutputWriters/OutputWriter.hpp"
#include "Util/readOperations.hpp"

using namespace broomstyx;

// Constructor
OutputManager::OutputManager()
    : _outputWriter(nullptr)
    , _nCsvOutput(0)
    , _csvFile(nullptr)
{}

// Destructor
OutputManager::~OutputManager() 
{
#ifdef VERBOSE_DESTRUCTION
    std::printf("\n  Destroying OutputManager... ");
    std::fflush(stdout);
#endif

    if ( _csvFile )
        std::fclose( _csvFile );
    for ( int i = 0; i < _nCsvOutput; i++)
        if ( _csvOutput[i] )
            delete _csvOutput[i];
    
    if ( _outputWriter )
        delete _outputWriter;
    
#ifdef VERBOSE_DESTRUCTION
    std::printf("done.");
    std::fflush(stdout);
#endif

}

// Public methods
void OutputManager::initializeCSVOutput()
{
    if ( _nCsvOutput > 0 )
    {
        // Create directory
        int dir_err = system("mkdir -p Output_CSV");
        if ( dir_err == -1 )
            throw std::runtime_error("Error creating directory 'Output_CSV'!\n");
        
        // Create CSV File
        _csvFile = std::fopen(_csvFilename.c_str(), "w");
        
        // Initialize output quantities and write column labels to CSV file
        std::fprintf(_csvFile, "Time, ");
        for ( int i = 0; i < (int)_csvOutput.size(); i++)
        {
            _csvOutput[i]->initialize();
            
            std::string quantityLabel = _csvOutput[i]->giveLabel();
            std::fprintf(_csvFile, "%s", quantityLabel.c_str());
            if ( i < (int)_csvOutput.size() - 1 )
                std::fprintf(_csvFile, ", ");
            else
                std::fprintf(_csvFile, "\n");
        }
    }
}

void OutputManager::initializeOutputWriter()
{
    _outputWriter->initialize();
}

void OutputManager::readOutputWriterFromFile( FILE* fp )
{
    std::string str, src = "OutputManager";
    
    str = getStringInputFrom(fp, "Failed to read output format from input file!", src);
    _outputWriter = objectFactory().instantiateOutputWriter(str);
    _outputWriter->readDataFrom(fp);
}

void OutputManager::readDataForCSVOutputFrom( FILE* fp )
{
    std::string label, str, src = "OutputManager";
    
    _nCsvOutput = getIntegerInputFrom(fp, "Failed to read number of CSV output from input file!", src);
    
    if ( _nCsvOutput > 0 )
    {
        // Create directory
        int dir_err = system("mkdir -p Output_CSV");
        if ( dir_err == -1 )
            throw std::runtime_error("Error creating directory 'Output_CSV'!\n");
        
        _csvOutput.assign(_nCsvOutput, nullptr);
        for ( int i = 0; i < _nCsvOutput; i++)
        {
            label = getStringInputFrom(fp, "Failed to read CSV output quantity label from input file!", "OutputManager");
            str = getStringInputFrom(fp, "Failed to read CSV output type from input file!", src);
            _csvOutput[i] = objectFactory().instantiateOutputQuantity(str);
            _csvOutput[i]->setOutputLabelTo(label);
            _csvOutput[i]->readDataFrom(fp);
        }

        verifyDeclaration(fp, "CSV_FILE", src);
        str = getStringInputFrom(fp, "Failed to read CSV output filename from input file!", src);

        _csvFilename = "./Output_CSV/" + std::string(str) + ".csv";
    }
}

void OutputManager::writeOutput( double time )
{
    std::chrono::time_point<std::chrono::system_clock> tic, toc;
    std::chrono::duration<double> tictoc;
    
    tic = std::chrono::high_resolution_clock::now();
    _outputWriter->writeOutput(time);
    toc = std::chrono::high_resolution_clock::now();
    tictoc = toc - tic;
    diagnostics().addOutputWriteTime(tictoc.count());
}

void OutputManager::writeOutputQuantities( double time )
{
    std::chrono::time_point<std::chrono::system_clock> tic, toc;
    std::chrono::duration<double> tictoc;
    
    tic = std::chrono::high_resolution_clock::now();

    if ( _nCsvOutput > 0 )
    {
        std::fprintf( _csvFile, "%.15e, ", time);
        for ( int i = 0; i < (int)_csvOutput.size(); i++)
        {
            double result = _csvOutput[i]->computeOutput();
            std::fprintf(_csvFile, "%.15e", result);
            if (i < (int)_csvOutput.size() - 1)
                std::fprintf(_csvFile, ", ");
            else
                std::fprintf(_csvFile, "\n");
        }
        std::fflush(_csvFile);
        std::printf("  --> %s\n", _csvFilename.c_str());
    }

    toc = std::chrono::high_resolution_clock::now();
    tictoc = toc - tic;
    diagnostics().addOutputWriteTime(tictoc.count());
}