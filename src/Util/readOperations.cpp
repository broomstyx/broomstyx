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

#include "readOperations.hpp"
#include <cstdio>
#include <stdexcept>

namespace broomstyx
{
    std::string getDeclarationFrom( FILE* fp )
    {
        // Find first '*' in input file
        char c = 0;
        do
        {
            c = fgetc(fp);
        } while ( c != '*' );


        char buf[100];
        if ( !std::fscanf(fp, "%s", buf) )
            throw std::runtime_error("Failed to read declaration in input file!");

        std::string str(buf);
        return str;
    }

    int getIntegerInputFrom( FILE* fp, const std::string& errmsg, const std::string& src )
    {
        int val;
        if ( !std::fscanf(fp, "%d", &val) )
            throw std::runtime_error(errmsg + "\nError source: " + src);

        return val;
    }

    double getRealInputFrom( FILE* fp, const std::string& errmsg, const std::string& src )
    {
        double val;
        if ( !std::fscanf(fp, "%lf", &val) )
            throw std::runtime_error(errmsg + "\nError source: " + src);

        return val;
    }

    std::string getStringInputFrom( FILE* fp, const std::string& errmsg, const std::string& src )
    {
        char buf[100];
        if ( !std::fscanf(fp, "%s", buf) )
            throw std::runtime_error(errmsg + "\nError source: " + src);

        std::string str(buf);
        return str;
    }

    void verifyDeclaration( FILE* fp, const std::string& str, const std::string& src )
    {
        char test[100];
        if ( !std::fscanf(fp, "%s", test) )
            throw std::runtime_error("Failed reading declaration '" + str
                    + "' from input file!\nError source: " + src);

        std::string test_str = test;
        if ( test_str != str )
            throw std::runtime_error("Declaration '" + str + "' expected in input file,\n\tstring '" + test_str + "' found.\nError source: " + src);
    }

    void verifyKeyword( FILE* fp, const std::string& str, const std::string& src )
    {
        char test[100];
        if ( !std::fscanf(fp, "%s", test) )
            throw std::runtime_error("Failed reading keyword '" + str + "' from input file!\nError source: " + src);

        std::string test_str = test;
        if ( test_str != str )
            throw std::runtime_error("Keyword '" + str + "' expected in input file,\n\tstring '" + test_str + "' found.\nError source: " + src);
    }
    
    void debug( const std::string& str)
    {
        std::printf("\n%s\n", str.c_str());
        std::fflush(stdout);
    }
}