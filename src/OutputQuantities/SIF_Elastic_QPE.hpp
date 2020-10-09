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

#ifndef SIF_ELASTIC_QPE_HPP
#define	SIF_ELASTIC_QPE_HPP

#include "OutputQuantity.hpp"
#include <vector>
#include "Util/RealVector.hpp"

namespace broomstyx
{
    class Dof;
    class Node;

    class SIF_Elastic_QPE final : public OutputQuantity {
    public:
        SIF_Elastic_QPE();
        virtual ~SIF_Elastic_QPE();

        void   initialize() override;
        void   readDataFrom(FILE* fp) override;
        double computeOutput() override;

    private:
        std::string _analysisMode;
        double _E;
        double _nu;

        std::string _crackTipLabel;
        Node* _crackTipNode;
        std::string _crackFaceLabel[2];
        int _dispDofNum[2];

        double _h;
        RealVector _crackTangent;
        RealVector _crackNormal;
        Node* _nodePairA[2];
        Node* _nodePairB[2];
        
        Dof* _dof[4][2];

        int _mode;
    };
}
    
#endif	/* SIF_ELASTIC_QPE_HPP */