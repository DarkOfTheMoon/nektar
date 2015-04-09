////////////////////////////////////////////////////////////////////////////////
//
//  File: BoundaryConditionBase.h
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
#ifndef NEKTAR_SPATIALDOMAINS_BOUNDARYCONDITIONSBASE_H
#define NEKTAR_SPATIALDOMAINS_BOUNDARYCONDITIONSBASE_H

#include <string>
#include <map>

#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/tuple/tuple_io.hpp>

#include <SpatialDomains/SpatialDomainsDeclspec.h>
#include <LibUtilities/BasicUtils/NekFactory.hpp>

class TiXmlElement;

namespace Nektar {

namespace LibUtilities {
class SessionReader;
typedef boost::shared_ptr<SessionReader> SessionReaderSharedPtr;
}

namespace SpatialDomains {

typedef boost::tuple<std::string, std::string> BCKey;

class BoundaryConditionBase;

/// Declaration of the boundary condition factory
typedef LibUtilities::NekFactory<
            BCKey,
            BoundaryConditionBase,
            const LibUtilities::SessionReaderSharedPtr&,
            TiXmlElement*> BoundaryConditionsFactory;

SPATIAL_DOMAINS_EXPORT BoundaryConditionsFactory& GetBoundaryConditionsFactory();


        enum BoundaryConditionType
        {
            eDirichlet,
            eNeumann,
            eRobin,
            ePeriodic,
            eNotDefined
        };

        enum BndUserDefinedType
        {
            eI,
            eMG,
            eHigh,
            eHighOutflow,
            eWall_Forces,
            eWall,
            eWallViscous,
            eArtificialViscosity,
            eSymmetry,
            eRinglebFlow,
            eMovingBody,
            eTimeDependent,
            eRadiation,
            eIsentropicVortex,
            eCalcBC,
            eQinflow,
            eTerminal,
            eRterminal,
            eCRterminal,
            eRCRterminal,
            eInflowCFS,
            eOutflowCFS,
            eRiemannInvariant,
            eExtrapOrder0,
            eNoUserDefined
        };

        const char* const BndUserDefinedTypeMap[] =
        {
            "I",
            "MG",
            "High",
            "HighOutflow",
            "Wall_Forces",
            "Wall",
            "WallViscous",
            "ArtificialVisc",
            "Symmetry",
            "RinglebFlow",
            "TimeDependent",
            "MovingBody",
            "Radiation",
            "IsentropicVortex",
            "CalcBC",
            "Qinflow",
            "Terminal",
            "Rterminal",
            "CRterminal",
            "RCRterminal",
            "InflowCFS",
            "OutflowCFS",
            "RiemannInvariant",
            "ExtrapOrder0",
            "NoUserDefined"
        };

        class BoundaryConditionBase
        {
            public:
                BoundaryConditionBase(
                    BoundaryConditionType type,
                    const std::string &userDefined = std::string("NoUserDefined"));
            
                virtual ~BoundaryConditionBase() {}

                BoundaryConditionType GetBoundaryConditionType() const
                {
                    return m_boundaryConditionType;
                }

                void SetBoundaryConditionType(BoundaryConditionType boundaryType) 
                {
                    m_boundaryConditionType = boundaryType;
                }

                void SetUserDefinedType(std::string type)
                {
                    m_userDefinedType = type;
                }

                std::string GetUserDefinedType() const
                {
                    return m_userDefinedType;
                }

                const std::string GetBndTypeAsString(BndUserDefinedType type)
                {
                    return BndUserDefinedTypeMap[type];
                }

                void SetIsTimeDependent(bool val)
                {
                    m_isTimeDependent = val;
                }

                bool IsTimeDependent(void)
                {
                    return m_isTimeDependent;
                }

            protected:
                BoundaryConditionType m_boundaryConditionType;
                std::string           m_userDefinedType;
                bool                  m_isTimeDependent; 
        };

        typedef boost::shared_ptr<BoundaryConditionBase> BoundaryConditionBaseSharedPtr;
}
}

#endif
