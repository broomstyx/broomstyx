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

#include <algorithm>
#include <chrono>
#include <cstring>
#include "Core/AnalysisModel.hpp"
#include "Core/Diagnostics.hpp"
#include "Core/ObjectFactory.hpp"

// Test source files
#include "../tests/test.hpp"

using namespace broomstyx;

int main( int argc, char **argv )
{
    if ( argc != 2 )
    {
        std::printf("\n\tError in program call: insufficient input!");
        std::printf("\n\tSample invocation: \"broomstyx <inputFile>\"\n\n");
        return 0;
    }
    
    if ( std::strcmp(argv[1],"--test") == 0 )
    	perform_tests();
    else
    {
		if ( objectFactory().hasError() )
		{
			std::printf("\n\tRegistration error detected in class factory.\n\n");
		}
		else
		{
			std::chrono::time_point<std::chrono::system_clock> tic, toc;
			std::chrono::duration<double> tictoc;
			tic = std::chrono::high_resolution_clock::now();

			try
			{
				analysisModel().initializeYourself(argv[1]);
				analysisModel().solveYourself();

				std::printf("\n\nRun successful ---> program will now terminate.");
			}
			catch (std::exception& e)
			{
				std::printf("\n%s\nException caught in main.cpp\n\n", e.what());
				std::fflush(stdout);
			}

			toc = std::chrono::high_resolution_clock::now();
			tictoc = toc - tic;
			std::printf("\nTotal runtime = %f seconds.\n\n", tictoc.count());

			diagnostics().outputDiagnostics();
		}
    }
    return 0;
}
