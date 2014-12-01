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
