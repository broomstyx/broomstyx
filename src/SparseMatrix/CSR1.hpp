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

#ifndef CSR1_HPP
#define	CSR1_HPP

#include "SparseMatrix.hpp"

namespace broomstyx
{
    class CSR1 : public SparseMatrix
    {
    public:
        CSR1();
        virtual ~CSR1();

        void addToComponent( int rowNum, int colNum, double val ) override;
        void atomicAddToComponent( int rowNum, int colNum, double val ) override;
        void finalizeProfile() override;
        std::tuple< int*,int* > giveProfileArrays() override;
        double* giveValArray() override;
        void initializeProfile( int dim1, int dim2 ) override;
        void initializeValues() override;
        void insertNonzeroComponentAt( int rowIdx, int colIdx) override;
        RealVector lumpRows() override;
        void       printTo( FILE* fp, int n ) override;
        RealVector times( const RealVector& x ) override;
    
    private:
        std::vector<int> _prfVec1;
        std::vector<int> _prfVec2;
        RealVector _val;
        std::vector< std::set<int> > _nz;
    };
}

#endif	/* CSR0_HPP */