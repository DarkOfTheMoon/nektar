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

    private:
        LibUtilities::EquationSharedPtr m_dirichletConditionPtr;
        LibUtilities::Equation m_dirichletCondition;
        std::string m_filename;
        static BCKey m_type;
};

class NeumannBoundaryCondition : public BoundaryConditionBase
{
    public:
        static BoundaryConditionBaseSharedPtr create(
                const LibUtilities::SessionReaderSharedPtr &pSession,
                TiXmlElement* pConditions) {
            return MemoryManager<NeumannBoundaryCondition>::AllocateSharedPtr(pSession, pConditions);
        }

        NeumannBoundaryCondition(const LibUtilities::SessionReaderSharedPtr &pSession,
                                       TiXmlElement* pBoundaryConditions);

        NeumannBoundaryCondition(
            const LibUtilities::SessionReaderSharedPtr &pSession,
            const std::string& eqn,
            const std::string& userDefined = std::string("NoUserDefined"),
            const std::string& filename=std::string("")):
                BoundaryConditionBase(eNeumann, userDefined),
                m_neumannCondition(pSession, eqn),
                m_filename(filename)
        {

        }

    private:
        LibUtilities::Equation m_neumannCondition;
        std::string m_filename;
        static BCKey m_type;
};

class RobinBoundaryCondition : public BoundaryConditionBase
{
    public:
        static BoundaryConditionBaseSharedPtr create(
                const LibUtilities::SessionReaderSharedPtr &pSession,
                TiXmlElement* pConditions) {
            return MemoryManager<NeumannBoundaryCondition>::AllocateSharedPtr(pSession, pConditions);
        }


        RobinBoundaryCondition(const LibUtilities::SessionReaderSharedPtr &pSession,
                               TiXmlElement* pBoundaryConditions);

        RobinBoundaryCondition(
            const LibUtilities::SessionReaderSharedPtr &pSession,
            const std::string &a,
            const std::string &b,
            const std::string &userDefined = std::string("NoUserDefined"),
            const std::string& filename=std::string("")):
                BoundaryConditionBase(eRobin, userDefined),
                m_robinFunction(pSession, a),
                m_robinPrimitiveCoeff(pSession, b),
                m_filename(filename)
        {
        }

    private:
        // \frac{\partial {u}}{\partial{n}} +
        // m_robinPrimativeCoeff(x,y,z)*u = m_robinFunction(x,y,z)
        LibUtilities::Equation m_robinFunction;
        LibUtilities::Equation m_robinPrimitiveCoeff;
        std::string m_filename;
        static BCKey m_type;
};


//Periodic need to set up extra parameters...from xml file

class PeriodicBoundaryCondition : public BoundaryConditionBase
{
    public:
        static BoundaryConditionBaseSharedPtr create(
                const LibUtilities::SessionReaderSharedPtr &pSession,
                TiXmlElement* pConditions) {
            return MemoryManager<PeriodicBoundaryCondition>::AllocateSharedPtr(pSession, pConditions);
        }


        PeriodicBoundaryCondition(const LibUtilities::SessionReaderSharedPtr &pSession,
                                      TiXmlElement* pBoundaryConditions);

        PeriodicBoundaryCondition(const unsigned int n):
            BoundaryConditionBase(ePeriodic),
            m_connectedBoundaryRegion(n)
        {
        }

    private:
        unsigned int m_connectedBoundaryRegion;
        static BCKey m_type;
};

class NotDefinedBoundaryCondition : public BoundaryConditionBase
{
    public:
        static BoundaryConditionBaseSharedPtr create(
                const LibUtilities::SessionReaderSharedPtr &pSession,
                TiXmlElement* pConditions) {
            return MemoryManager<NotDefinedBoundaryCondition>::AllocateSharedPtr(pSession, pConditions);
        }

        NotDefinedBoundaryCondition(const LibUtilities::SessionReaderSharedPtr &pSession,
                                      TiXmlElement* pBoundaryConditions);

        NotDefinedBoundaryCondition(
            const LibUtilities::SessionReaderSharedPtr &pSession,
            const std::string& eqn,
            const std::string& userDefined = std::string("NoUserDefined"),
            const std::string& filename=std::string("")):
                BoundaryConditionBase(eNotDefined, userDefined),
        m_notDefinedCondition(pSession, eqn),
        m_filename(filename)
        {
        }

    private:
        LibUtilities::Equation m_notDefinedCondition;
        std::string m_filename;
        static BCKey m_type;
};

}
}

#endif
