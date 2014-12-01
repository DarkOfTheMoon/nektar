#include <vector>

#include <SpatialDomains/BoundaryConditionTypes.h>

namespace Nektar {
namespace SpatialDomains {

BCKey DirichletBoundaryCondition::m_type = GetBoundaryConditionsFactory().RegisterCreatorFunction(
        BCKey("D",""),
        DirichletBoundaryCondition::create,
        "Dirichlet");

DirichletBoundaryCondition::DirichletBoundaryCondition(const LibUtilities::SessionReaderSharedPtr &pSession,
                            TiXmlElement* pBoundaryConditions):
        BoundaryConditionBase(eDirichlet, std::string("")),
        m_dirichletCondition(pSession, std::string(""))
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
            LibUtilities::Equation m_dirichletCondition(pSession,attrData);
            m_dirichletConditionPtr= MemoryManager<LibUtilities::Equation>::AllocateSharedPtr(pSession,attrData);
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
    m_neumannCondition(pSession, std::string(""))
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
            LibUtilities::Equation m_neumannCondition(pSession,attrData);
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
    m_robinFunction(pSession, std::string("")),
    m_robinPrimitiveCoeff(pSession, std::string(""))
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

                // here I need to instantiate a with attrData1;
                LibUtilities::Equation m_robinFunction(pSession,attrData1);


                attr = attr->Next();
                ASSERTL0(attr, "Unable to read PRIMCOEFF attribute.");

                attrName1= attr->Name();
                ASSERTL0(attrName1 == "PRIMCOEFF", (std::string("Unknown attribute: ") + attrName1).c_str());

                attrData1 = attr->Value();
                ASSERTL0(!attrData1.empty(), "PRIMCOEFF attributes must have associated values.");

                pSession->SubstituteExpressions(attrData1);

                LibUtilities::Equation m_robinPrimitiveCoeff(pSession,attrData1);
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
    //ASSERTL0(beg < end, (std::string("Error reading periodic boundary region definition for boundary region: ")
    //                     + boundaryRegionIDStrm.str()).c_str());

    std::vector<unsigned int> periodicBndRegionIndex;
    //bool parseGood = ParseUtils::GenerateSeqVector(periodicBndRegionIndexStr.c_str(), periodicBndRegionIndex);

    //ASSERTL0(parseGood && (periodicBndRegionIndex.size()==1), (std::string("Unable to read periodic boundary condition for boundary region: ")
    //                                                           + boundaryRegionIDStrm.str()).c_str());

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
    m_notDefinedCondition(pSession, std::string(""))
{
}

}
}
