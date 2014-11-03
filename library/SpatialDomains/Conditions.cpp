////////////////////////////////////////////////////////////////////////////////
//
//  File: BoundaryConditions.cpp
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

#include <LibUtilities/BasicUtils/ParseUtils.hpp>
#include <SpatialDomains/Conditions.h>
#include <tinyxml.h>

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
		
		
        BoundaryConditions::BoundaryConditions(const LibUtilities::SessionReaderSharedPtr &pSession, const MeshGraphSharedPtr &meshGraph)
            : m_meshGraph(meshGraph), 
              m_session  (pSession)
              
        {
            Read(m_session->GetElement("Nektar/Conditions"));
        }


        BoundaryConditions::BoundaryConditions(void)
        {
        }

        BoundaryConditions::~BoundaryConditions(void)
        {
        }


        /**
         *
         */
        void BoundaryConditions::Read(TiXmlElement *conditions)
        {
            ASSERTL0(conditions, "Unable to find CONDITIONS tag in file.");

            TiXmlElement *boundaryRegions = conditions->FirstChildElement("BOUNDARYREGIONS");

            if(boundaryRegions)
            {
                ReadBoundaryRegions(conditions);

                ReadBoundaryConditions(conditions);
            }
        }
		
		std::ostream &operator<<(std::ostream &os, boost::tuple<std::string,std::string> &p)
        {
            os << boost::get<0>(p)  << " "
			   << boost::get<1>(p);
            return os;
        }


        /**
         *
         */
        void BoundaryConditions::ReadBoundaryRegions(TiXmlElement *conditions)
        {
            // ensure boundary regions only read once per class definition
            if(m_boundaryRegions.size() != 0)
            {
                return;
            }

            TiXmlElement *boundaryRegions = conditions->FirstChildElement("BOUNDARYREGIONS");
            ASSERTL0(boundaryRegions, "Unable to find BOUNDARYREGIONS block.");

            // See if we have boundary regions defined.
            TiXmlElement *boundaryRegionsElement = boundaryRegions->FirstChildElement("B");

            while (boundaryRegionsElement)
            {
                /// All elements are of the form: "<B ID="#"> ... </B>", with
                /// ? being the element type.
                int indx;
                int err = boundaryRegionsElement->QueryIntAttribute("ID", &indx);
                ASSERTL0(err == TIXML_SUCCESS, "Unable to read attribute ID.");

                TiXmlNode* boundaryRegionChild = boundaryRegionsElement->FirstChild();
                // This is primarily to skip comments that may be present.
                // Comments appear as nodes just like elements.
                // We are specifically looking for text in the body
                // of the definition.
                while(boundaryRegionChild && boundaryRegionChild->Type() != TiXmlNode::TEXT)
                {
                    boundaryRegionChild = boundaryRegionChild->NextSibling();
                }

                ASSERTL0(boundaryRegionChild, "Unable to read variable definition body.");
                std::string boundaryRegionStr = boundaryRegionChild->ToText()->ValueStr();

                std::string::size_type indxBeg = boundaryRegionStr.find_first_of('[') + 1;
                std::string::size_type indxEnd = boundaryRegionStr.find_last_of(']') - 1;

                ASSERTL0(indxBeg <= indxEnd, (std::string("Error reading boundary region definition:") + boundaryRegionStr).c_str());

                std::string indxStr = boundaryRegionStr.substr(indxBeg, indxEnd - indxBeg + 1);

                if (!indxStr.empty())
                {
                    // Extract the composites from the string and return them in a list.
                    BoundaryRegionShPtr boundaryRegion(MemoryManager<BoundaryRegion>::AllocateSharedPtr());

                    ASSERTL0(m_boundaryRegions.count(indx) == 0,
                             "Boundary region "+indxStr+ " defined more than "
                             "once!");
                    
                    m_meshGraph->GetCompositeList(indxStr, *boundaryRegion);
                    m_boundaryRegions[indx] = boundaryRegion;
                }

                boundaryRegionsElement = boundaryRegionsElement->NextSiblingElement("B");
            }
        }


        /**
         *
         */
        void BoundaryConditions::ReadBoundaryConditions(TiXmlElement *conditions)
        {
            // Protect against multiple reads.
            if(m_boundaryConditions.size() != 0)
            {
                return;
            }

            // Read REGION tags
            TiXmlElement *boundaryConditionsElement = conditions->FirstChildElement("BOUNDARYCONDITIONS");
            ASSERTL0(boundaryConditionsElement, "Boundary conditions must be specified.");

            TiXmlElement *regionElement = boundaryConditionsElement->FirstChildElement("REGION");

            // Read R (Robin), D (Dirichlet), N (Neumann), P (Periodic) C(Cauchy) tags
            while (regionElement)
            {
                BoundaryConditionMapShPtr boundaryConditions = MemoryManager<BoundaryConditionMap>::AllocateSharedPtr();

                int boundaryRegionID;
                int err = regionElement->QueryIntAttribute("REF", &boundaryRegionID);
                ASSERTL0(err == TIXML_SUCCESS, "Error reading boundary region reference.");

                // Find the boundary region corresponding to this ID.
                std::string boundaryRegionIDStr;
                std::ostringstream boundaryRegionIDStrm(boundaryRegionIDStr);
                boundaryRegionIDStrm << boundaryRegionID;

                ASSERTL0(m_boundaryRegions.count(boundaryRegionID) == 1,
                         "Boundary region " + boost::lexical_cast<
                         string>(boundaryRegionID)+ " not found");

                TiXmlElement *conditionElement = regionElement->FirstChildElement();
                std::vector<std::string> vars = m_session->GetVariables();
				
                while (conditionElement)
                {
                    // Check type.
                    std::string conditionType = conditionElement->Value();
                    std::string attrData;
					//for the inner map
					std::vector<std::string>::iterator iter;
                    std::string attrName;
					
                    attrData = conditionElement->Attribute("VAR");
					
                    if (!attrData.empty())
                    {
                        iter = std::find(vars.begin(), vars.end(), attrData);
                        ASSERTL0(iter != vars.end(), (std::string("Cannot find variable: ") + attrData).c_str());
                    }
					//end for the inner map
					
                    // All have var specified, or else all variables are zero.
                    TiXmlAttribute *attr = conditionElement->FirstAttribute();
					std::string userdefined = conditionElement->Attribute("USERDEFINEDTYPE");
					
					boost::tuple<std::string,std::string> bnd_pair=boost::make_tuple(conditionType,userdefined);
					BoundaryConditionShPtr bnd=GetBoundaryConditionsFactory().CreateInstance(bnd_pair,m_session, conditionElement);
					//I need to put bnd it inside a map 
					(*boundaryConditions)[*iter]=bnd;
				}						
					
				//outer mapping to regions ID
                m_boundaryConditions[boundaryRegionID] = boundaryConditions;
                regionElement = regionElement->NextSiblingElement("REGION");
            }
        }
    }
}
