////////////////////////////////////////////////////////////////////////////////
//
//  File: BoundaryConditionsTypes.cpp
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
#include <vector>
#include <LibUtilities/BasicUtils/ParseUtils.hpp>
#include <SpatialDomains/BoundaryConditionTypes.h>

namespace Nektar {
namespace SpatialDomains {

    // Register boundary conditions
    BCKey DirichletBoundaryCondition::m_type = GetBoundaryConditionsFactory().RegisterCreatorFunction(
        BCKey("D","NoUserDefined"),
        DirichletBoundaryCondition::create,
        "Dirichlet");


    BCKey NeumannBoundaryCondition::m_type = GetBoundaryConditionsFactory().RegisterCreatorFunction(
        BCKey("N","NoUserDefined"),
        NeumannBoundaryCondition::create,
        "Neumann");


    BCKey RobinBoundaryCondition::m_type = GetBoundaryConditionsFactory().RegisterCreatorFunction(
        BCKey("R","NoUserDefined"),
        RobinBoundaryCondition::create,
        "Robin");

    BCKey PeriodicBoundaryCondition::m_type = GetBoundaryConditionsFactory().RegisterCreatorFunction(
        BCKey("P","NoUserDefined"),
        PeriodicBoundaryCondition::create,
        "Periodic");
    
    DirichletBoundaryCondition::DirichletBoundaryCondition(const LibUtilities::SessionReaderSharedPtr &pSession, TiXmlElement* pBoundaryConditions):
        BoundaryConditionBase(eDirichlet, std::string("")),
        m_dirichletCondition(MemoryManager<LibUtilities::Equation>::AllocateSharedPtr(pSession, std::string("")))
    {

    TiXmlAttribute *attr = pBoundaryConditions->FirstAttribute();

     std::vector<std::string>::iterator iter;
     std::string attrName;
     std::vector<std::string> vars = pSession->GetVariables();
     std::string attrData = pBoundaryConditions->Attribute("VAR");

    if (!attrData.empty())
    {
        iter = std::find(vars.begin(), vars.end(), attrData);
        ASSERTL0(iter != vars.end(), (std::string("Cannot find variable: ") + attrData).c_str());
    }

    attr = attr->Next();
    while(attr)
    {
        attrName = attr->Name();
        if(attrName=="VALUE")
        {
            ASSERTL0(attrName == "VALUE", (std::string("Unknown attribute: ") + attrName).c_str());

            attrData = attr->Value();
            ASSERTL0(!attrData.empty(), "VALUE attribute must have associated value.");
            m_dirichletCondition = MemoryManager<LibUtilities::Equation>::AllocateSharedPtr(pSession, attrData);
            
            pSession->SubstituteExpressions(attrData);
        }
        else if(attrName=="FILE")
        {

            ASSERTL0(attrName == "FILE", (std::string("Unknown attribute: ") + attrName).c_str());

            attrData = attr->Value();
            ASSERTL0(!attrData.empty(), "FILE attribute must be specified.");

            pSession->SubstituteExpressions(attrData);

            m_filename = attrData;
        }

    attr = attr->Next();
    }

}

NeumannBoundaryCondition::NeumannBoundaryCondition(
        const LibUtilities::SessionReaderSharedPtr &pSession,
              TiXmlElement* pBoundaryConditions):
    BoundaryConditionBase(eNeumann, std::string("")),
                             m_neumannCondition(MemoryManager<LibUtilities::Equation>::AllocateSharedPtr(pSession, std::string("")))
{

    TiXmlAttribute *attr = pBoundaryConditions->FirstAttribute();

    std::vector<std::string>::iterator iter;
    std::string attrName;
    std::vector<std::string> vars = pSession->GetVariables();
    std::string attrData = pBoundaryConditions->Attribute("VAR");

    if (!attrData.empty())
    {
        iter = std::find(vars.begin(), vars.end(), attrData);
        ASSERTL0(iter != vars.end(), (std::string("Cannot find variable: ") + attrData).c_str());
    }

    attr = attr->Next();
    while(attr)
    {
        attrName = attr->Name();
        if(attrName=="VALUE")
        {
            ASSERTL0(attrName == "VALUE", (std::string("Unknown attribute: ") + attrName).c_str());

            attrData = attr->Value();
            ASSERTL0(!attrData.empty(), "VALUE attribute must have associated value.");
            m_neumannCondition = MemoryManager<LibUtilities::Equation>::AllocateSharedPtr(pSession, attrData);

            pSession->SubstituteExpressions(attrData);
        }
        else if(attrName=="FILE")
        {

            ASSERTL0(attrName == "FILE", (std::string("Unknown attribute: ") + attrName).c_str());

            attrData = attr->Value();
            ASSERTL0(!attrData.empty(), "FILE attribute must be specified.");

            pSession->SubstituteExpressions(attrData);

            m_filename = attrData;
        }

        attr = attr->Next();
    }
}

RobinBoundaryCondition::RobinBoundaryCondition(
        const LibUtilities::SessionReaderSharedPtr &pSession,
              TiXmlElement* pBoundaryConditions):
    BoundaryConditionBase(eRobin, std::string("")),
                             m_robinFunction(MemoryManager<LibUtilities::Equation>::AllocateSharedPtr(pSession, std::string(""))),
                             m_robinPrimitiveCoeff(MemoryManager<LibUtilities::Equation>::AllocateSharedPtr(pSession, std::string("")))
{
    TiXmlAttribute *attr = pBoundaryConditions->FirstAttribute();


    std::vector<std::string>::iterator iter;
    std::string attrName;
    std::vector<std::string> vars = pSession->GetVariables();
    std::string attrData = pBoundaryConditions->Attribute("VAR");

    attr = attr->Next();

    if (attr)
    {
        std::string attrName1;
        std::string attrData1;

        while(attr)
        {
            attrName1 = attr->Name();

            if(attrName1 == "VALUE")
            {
                attrData1 = attr->Value();
                ASSERTL0(!attrData1.empty(), "VALUE attributes must have associated values.");

                pSession->SubstituteExpressions(attrData1);

                m_robinFunction = MemoryManager<LibUtilities::Equation>::AllocateSharedPtr(pSession, attrData1);

                attr = attr->Next();
                ASSERTL0(attr, "Unable to read PRIMCOEFF attribute.");

                attrName1= attr->Name();
                ASSERTL0(attrName1 == "PRIMCOEFF", (std::string("Unknown attribute: ") + attrName1).c_str());

                attrData1 = attr->Value();
                ASSERTL0(!attrData1.empty(), "PRIMCOEFF attributes must have associated values.");

                pSession->SubstituteExpressions(attrData1);

                m_robinPrimitiveCoeff = MemoryManager<LibUtilities::Equation>::AllocateSharedPtr(pSession, attrData1);
            }
            else if(attrName1=="FILE")
            {

                ASSERTL0(attrName1 == "FILE", (std::string("Unknown attribute: ") + attrName1).c_str());

                attrData1 = attr->Value();
                ASSERTL0(!attrData1.empty(), "FILE attribute must be specified.");

                pSession->SubstituteExpressions(attrData1);
                m_filename = attrData1;
            }
            attr = attr->Next();

        }
    }
}


PeriodicBoundaryCondition::PeriodicBoundaryCondition(
        const LibUtilities::SessionReaderSharedPtr &pSession,
              TiXmlElement* pBoundaryConditions):
    BoundaryConditionBase(ePeriodic, std::string(""))
{

    TiXmlAttribute *attr = pBoundaryConditions->FirstAttribute();


    std::vector<std::string>::iterator iter;
    std::string attrName;
    std::vector<std::string> vars = pSession->GetVariables();
    std::string attrData = pBoundaryConditions->Attribute("VAR");

    attr = attr->Next();
    if(attr)
    {
        attrName = attr->Name();
        
        ASSERTL0(attrName == "VALUE", (std::string("Unknown attribute: ") + attrName).c_str());
        
        attrData = attr->Value();
        ASSERTL0(!attrData.empty(), "VALUE attribute must have associated value.");
        
        int beg = attrData.find_first_of("[");
        int end = attrData.find_first_of("]");
        std::string periodicBndRegionIndexStr = attrData.substr(beg+1,end-beg-1);
        ASSERTL0(beg < end, (std::string("Error reading periodic boundary region definition using") + 
                             attrData).c_str());

        std::vector<unsigned int> periodicBndRegionIndex;
        bool parseGood = ParseUtils::GenerateSeqVector(periodicBndRegionIndexStr.c_str(), periodicBndRegionIndex);
        
        ASSERTL0(parseGood && (periodicBndRegionIndex.size()==1), (std::string("Unable to identify periodic boundary condition in string: ")  + attrData).c_str());
        
        m_connectedBoundaryRegion=periodicBndRegionIndex[0];
    }
    else
    {
        ASSERTL0(false, "Periodic boundary conditions should be explicitely defined");
    }

}

NotDefinedBoundaryCondition::NotDefinedBoundaryCondition(
        const LibUtilities::SessionReaderSharedPtr &pSession,
              TiXmlElement* pBoundaryConditions):
    BoundaryConditionBase(eNotDefined, std::string("")),
    m_notDefinedCondition(MemoryManager<LibUtilities::Equation>::AllocateSharedPtr(pSession, std::string("")))
{
}

}
}
