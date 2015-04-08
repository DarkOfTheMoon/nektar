////////////////////////////////////////////////////////////////////////////////
//
//  File: BoundaryConditionBase.cpp
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  License for the specific language governing rights and limitations under
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description:
//
//
////////////////////////////////////////////////////////////////////////////////
#include <loki/Singleton.h>
#include <SpatialDomains/BoundaryConditionBase.h>

namespace Nektar
{
namespace SpatialDomains
{

    BoundaryConditionsFactory& GetBoundaryConditionsFactory()
    {
        typedef Loki::SingletonHolder<BoundaryConditionsFactory,
                                      Loki::CreateUsingNew,
                                      Loki::NoDestroy > Type;
        return Type::Instance();
    }
    
    BoundaryConditionBase::BoundaryConditionBase(
        BoundaryConditionType type,
        const std::string &userDefined):
        m_boundaryConditionType(type)
    {
        std::map<const std::string, BndUserDefinedType>  known_type;
        known_type["H"] = eHigh;
        known_type["HOutflow"] = eHighOutflow;
        known_type["I"] = eI;
        known_type["MG"] = eMG;
        known_type["Wall"] = eWall;
        known_type["WallViscous"] = eWallViscous;
        known_type["ArtificialVisc"] = eArtificialViscosity;
        known_type["Q-inflow"] = eQinflow;
        known_type["Terminal"] = eTerminal;
        known_type["R-terminal"] = eRterminal;
        known_type["CR-terminal"] = eCRterminal;
        known_type["RCR-terminal"] = eRCRterminal;
        known_type["CalcBC"] = eCalcBC;
        known_type["RinglebFlow"] = eRinglebFlow;
        known_type["Symmetry"] = eSymmetry;
        known_type["TimeDependent"] = eTimeDependent;
        known_type["Radiation"] = eRadiation;
        known_type["IsentropicVortex"] = eIsentropicVortex;
        known_type["RiemannInvariant"] = eRiemannInvariant;
        known_type["ExtrapOrder0"]     = eExtrapOrder0;
        known_type["NoUserDefined"]    = eNoUserDefined;
        
        std::map<const std::string, BndUserDefinedType>::
            const_iterator it = known_type.find(userDefined);
        if (it != known_type.end())
        {
            m_userDefined = it->second;
        }
        else
        {
            //ASSERTL0(false, std::string("Unknown boundary condition "
            //"user defined type [") + userDefined + std::string("]"));
            m_userDefined = eNoUserDefined;
        }
    }
    
}
} 
