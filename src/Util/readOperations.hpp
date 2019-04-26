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

#ifndef READOPERATIONS_HPP
#define	READOPERATIONS_HPP

#include <cstdio>
#include <string>

namespace broomstyx
{
    std::string
    getDeclarationFrom( FILE* fp );
    
    int         
    getIntegerInputFrom( FILE* fp
                       , const std::string& errmsg
                       , const std::string& src );
    
    double
    getRealInputFrom( FILE* fp
                    , const std::string& errmsg
                    , const std::string& src );
    
    std::string
    getStringInputFrom( FILE* fp
                      , const std::string& errmsg
                      , const std::string& src );
    
    void
    verifyDeclaration( FILE* fp
                     , const std::string& str
                     , const std::string& src );
    
    void
    verifyKeyword( FILE* fp
                 , const std::string& str
                 , const std::string& src );
    
    void debug( const std::string& str );
}

#endif	/* READOPERATIONS_HPP */