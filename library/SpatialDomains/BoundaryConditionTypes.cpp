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

}
}
