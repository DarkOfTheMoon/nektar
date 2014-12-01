#ifndef NEKTAR_SPATIALDOMAINS_BOUNDARYCONDITIONTYPES_H
#define NEKTAR_SPATIALDOMAINS_BOUNDARYCONDITIONTYPES_H

#include <LibUtilities/BasicUtils/Equation.h>
#include <SpatialDomains/BoundaryConditionBase.h>
#include <tinyxml.h>

namespace Nektar {
namespace SpatialDomains {

class DirichletBoundaryCondition : public BoundaryConditionBase
{
    public:
        static BoundaryConditionBaseSharedPtr create(
                const LibUtilities::SessionReaderSharedPtr &pSession,
                TiXmlElement* pConditions) {
            return MemoryManager<DirichletBoundaryCondition>::AllocateSharedPtr(pSession, pConditions);
        }

            DirichletBoundaryCondition(const LibUtilities::SessionReaderSharedPtr &pSession,
                                        TiXmlElement* pBoundaryConditions);
            
            DirichletBoundaryCondition(
                const LibUtilities::SessionReaderSharedPtr &pSession,
                const std::string& eqn,
                const std::string& userDefined = std::string("NoUserDefined"),
                const std::string& filename=std::string("")):
                    BoundaryConditionBase(eDirichlet, userDefined),
                    m_dirichletCondition(pSession, eqn),
                    m_filename(filename)
            {
            }

        
            LibUtilities::EquationSharedPtr m_dirichletConditionPtr;
            LibUtilities::Equation m_dirichletCondition;
            std::string m_filename;
            static BCKey m_type;
        };


}
}

#endif
