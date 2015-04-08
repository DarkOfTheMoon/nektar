////////////////////////////////////////////////////////////////////////////////
//
//  File: BoundaryConditionTypes.h
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
            m_dirichletCondition(MemoryManager<LibUtilities::Equation>::AllocateSharedPtr(pSession, eqn)),
            m_filename(filename)
        {
        }

        const LibUtilities::EquationSharedPtr &GetDirichletCondition(void) const
        {
            return m_dirichletCondition;
        }

        std::string  GetFilename(void)
        {
            return m_filename; 
        }

    private:
        LibUtilities::EquationSharedPtr m_dirichletCondition;
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

        NeumannBoundaryCondition(const LibUtilities::SessionReaderSharedPtr &pSession, TiXmlElement* pBoundaryConditions);

       NeumannBoundaryCondition(
              const LibUtilities::SessionReaderSharedPtr &pSession,
              const std::string& eqn,
              const std::string& userDefined = std::string("NoUserDefined"),
              const std::string& filename=std::string("")):
        BoundaryConditionBase(eNeumann, userDefined),
            m_neumannCondition(MemoryManager<LibUtilities::Equation>::AllocateSharedPtr(pSession, eqn)),
            m_filename(filename)
            {
            }
        
        const LibUtilities::EquationSharedPtr &GetNeumannCondition(void) const
        {
            return m_neumannCondition;
        }
        
        std::string  GetFilename(void)
        {
            return m_filename; 
        }
    private:
        LibUtilities::EquationSharedPtr m_neumannCondition;
        std::string m_filename;
        static BCKey m_type;
};

class RobinBoundaryCondition : public BoundaryConditionBase
{
    public:
        static BoundaryConditionBaseSharedPtr create(
                const LibUtilities::SessionReaderSharedPtr &pSession,
                TiXmlElement* pConditions) {
            return MemoryManager<RobinBoundaryCondition>::AllocateSharedPtr(pSession, pConditions);
        }

        
        RobinBoundaryCondition(const LibUtilities::SessionReaderSharedPtr &pSession, TiXmlElement* pBoundaryConditions);

        RobinBoundaryCondition(
            const LibUtilities::SessionReaderSharedPtr &pSession,
            const std::string &a,
            const std::string &b,
            const std::string &userDefined = std::string("NoUserDefined"),
            const std::string& filename=std::string("")):
        BoundaryConditionBase(eRobin, userDefined),
            m_robinFunction(MemoryManager<LibUtilities::Equation>::AllocateSharedPtr(pSession, a)),
            m_robinPrimitiveCoeff(MemoryManager<LibUtilities::Equation>::AllocateSharedPtr(pSession, b)),
            m_filename(filename)
            {
            }
        
        const LibUtilities::EquationSharedPtr GetRobinFunction(void) const
        {
            return m_robinFunction;
        }


        const LibUtilities::EquationSharedPtr &GetRobinPrimitiveCoeff(void) const
        {
            return m_robinPrimitiveCoeff;
        }
        
        std::string  GetFilename(void)
        {
            return m_filename; 
        }
        
        private:
        // \frac{\partial {u}}{\partial{n}} +
        // m_robinPrimativeCoeff(x,y,z)*u = m_robinFunction(x,y,z)
        LibUtilities::EquationSharedPtr m_robinFunction;
        LibUtilities::EquationSharedPtr m_robinPrimitiveCoeff;
        std::string m_filename;
        static BCKey m_type;
};


//Periodic need to set up extra parameters...from xml file

class PeriodicBoundaryCondition : public BoundaryConditionBase
{
public:
    static BoundaryConditionBaseSharedPtr create(
                const LibUtilities::SessionReaderSharedPtr &pSession,
                TiXmlElement* pConditions) 
    {
        return MemoryManager<PeriodicBoundaryCondition>::AllocateSharedPtr(pSession, pConditions);
    }
    
    
    PeriodicBoundaryCondition(const LibUtilities::SessionReaderSharedPtr &pSession, TiXmlElement* pBoundaryConditions);
    
    PeriodicBoundaryCondition(const unsigned int n):
    BoundaryConditionBase(ePeriodic),
        m_connectedBoundaryRegion(n)
        {
        }
    
    int GetConnectedBoundaryRegion(void)
    {
        return m_connectedBoundaryRegion; 
    }
    
    private:
        unsigned int m_connectedBoundaryRegion;
        static BCKey m_type;
};

// not sure this is required
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
            m_notDefinedCondition(MemoryManager<LibUtilities::Equation>::AllocateSharedPtr(pSession, eqn)),
            m_filename(filename)
            {
            }
                
       private:
                LibUtilities::EquationSharedPtr m_notDefinedCondition;
        std::string m_filename;
        static BCKey m_type;
};

}
}

#endif
