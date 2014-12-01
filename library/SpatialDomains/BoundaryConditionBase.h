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

                void SetUserDefined(BndUserDefinedType type)
                {
                    m_userDefined = type;
                }

                BndUserDefinedType GetUserDefined() const
                {
                    return m_userDefined;
                }

                const std::string GetBndTypeAsString(BndUserDefinedType type)
                {
                    return BndUserDefinedTypeMap[type];
                }

            protected:
                BoundaryConditionType m_boundaryConditionType;
                BndUserDefinedType    m_userDefined;
        };

        typedef boost::shared_ptr<BoundaryConditionBase> BoundaryConditionBaseSharedPtr;
}
}

#endif
