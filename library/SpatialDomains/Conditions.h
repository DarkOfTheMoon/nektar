////////////////////////////////////////////////////////////////////////////////
//
//  File:  Conditions.h
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
#ifndef NEKTAR_SPATIALDOMAINS_BOUNDARYCONDITIONS_H
#define NEKTAR_SPATIALDOMAINS_BOUNDARYCONDITIONS_H

#include <string>
#include <map>

#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/Equation.h>
#include <SpatialDomains/MeshGraph.h>


namespace Nektar
{
    struct OneD;

    namespace SpatialDomains
    {
		class BoundaryConditionBase;
		
		/// Declaration of the boundary condition factory
		typedef LibUtilities::NekFactory< 
		std::pair<std::string,std::string>,
		BoundaryConditionBase,
		const LibUtilities::SessionReaderSharedPtr&,
		const TiXmlElement*> BoundaryConditionsFactory;
	
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

        struct BoundaryConditionBase
        {
            BoundaryConditionBase(
                BoundaryConditionType type,
                const std::string &userDefined = std::string("NoUserDefined")):
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

            virtual ~BoundaryConditionBase()
            {};

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
            BoundaryConditionType  m_boundaryConditionType;
            BndUserDefinedType     m_userDefined;
        };


        struct DirichletBoundaryCondition : public BoundaryConditionBase
        {
			
			DirichletBoundaryCondition(const LibUtilities::SessionReaderSharedPtr &pSession,
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
        };

        struct NeumannBoundaryCondition : public BoundaryConditionBase
        {
			
			
			
			NeumannBoundaryCondition(const LibUtilities::SessionReaderSharedPtr &pSession,
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

            LibUtilities::Equation m_neumannCondition;
            std::string m_filename;
        };

        struct RobinBoundaryCondition : public BoundaryConditionBase
        {
            
			RobinBoundaryCondition(const LibUtilities::SessionReaderSharedPtr &pSession,
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
            // \frac{\partial {u}}{\partial{n}} +
            // m_robinPrimativeCoeff(x,y,z)*u = m_robinFunction(x,y,z)
            LibUtilities::Equation m_robinFunction;
            LibUtilities::Equation m_robinPrimitiveCoeff;
            std::string m_filename;
        };


        struct PeriodicBoundaryCondition : public BoundaryConditionBase
        {
			
			PeriodicBoundaryCondition(const LibUtilities::SessionReaderSharedPtr &pSession,
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
				//					 + boundaryRegionIDStrm.str()).c_str());
				
				vector<unsigned int> periodicBndRegionIndex;
				//bool parseGood = ParseUtils::GenerateSeqVector(periodicBndRegionIndexStr.c_str(), periodicBndRegionIndex);
				
				//ASSERTL0(parseGood && (periodicBndRegionIndex.size()==1), (std::string("Unable to read periodic boundary condition for boundary region: ")
				//														   + boundaryRegionIDStrm.str()).c_str());

			    m_connectedBoundaryRegion=periodicBndRegionIndex[0];
				}
				else
				{
					ASSERTL0(false, "Periodic boundary conditions should be explicitely defined");
				}
                   				
			}	

			
			PeriodicBoundaryCondition(const unsigned int n):
                BoundaryConditionBase(ePeriodic),
                m_connectedBoundaryRegion(n)
            {
            }
			
            unsigned int m_connectedBoundaryRegion;
        };

        struct NotDefinedBoundaryCondition : public BoundaryConditionBase
        {
            
		
			
			NotDefinedBoundaryCondition(const LibUtilities::SessionReaderSharedPtr &pSession,
									  const TiXmlElement* pBoundaryConditions):
			BoundaryConditionBase(eNotDefined, std::string("")),
			m_notDefinedCondition(pSession, std::string(""))			
			{
				
			}
			
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

            LibUtilities::Equation m_notDefinedCondition;
            std::string m_filename;
			 
        };


        typedef std::map<int, Composite> BoundaryRegion;
        typedef boost::shared_ptr<BoundaryRegion> BoundaryRegionShPtr;
        typedef boost::shared_ptr<const BoundaryRegion> ConstBoundaryRegionShPtr;
        typedef std::map<int, BoundaryRegionShPtr> BoundaryRegionCollection;
		typedef boost::shared_ptr<LibUtilities::Equation> EquationSharedPtr;


        typedef boost::shared_ptr<BoundaryConditionBase> BoundaryConditionShPtr;
        typedef boost::shared_ptr<DirichletBoundaryCondition> DirichletBCShPtr;
        typedef boost::shared_ptr<NeumannBoundaryCondition>   NeumannBCShPtr;
        typedef boost::shared_ptr<RobinBoundaryCondition>     RobinBCShPtr;
		
        typedef std::map<std::string,BoundaryConditionShPtr>  BoundaryConditionMap;
        typedef boost::shared_ptr<BoundaryConditionMap>  BoundaryConditionMapShPtr;
        typedef std::map<int, BoundaryConditionMapShPtr> BoundaryConditionCollection;

        const static Array<OneD, BoundaryConditionShPtr> NullBoundaryConditionShPtrArray;
		

        class BoundaryConditions
        {
        public:
            SPATIAL_DOMAINS_EXPORT BoundaryConditions(const LibUtilities::SessionReaderSharedPtr &pSession, const MeshGraphSharedPtr &meshGraph);

            SPATIAL_DOMAINS_EXPORT BoundaryConditions(void);
            SPATIAL_DOMAINS_EXPORT ~BoundaryConditions(void);

            const BoundaryRegionCollection &GetBoundaryRegions(void) const
            {
                return m_boundaryRegions;
            }

            void AddBoundaryRegions(const int regionID, BoundaryRegionShPtr &bRegion)
            {
                m_boundaryRegions[regionID] = bRegion;
            }

            const BoundaryConditionCollection &GetBoundaryConditions(void) const
            {
                return m_boundaryConditions;
            }


            void AddBoundaryConditions(const int regionID, BoundaryConditionMapShPtr &bCond)
            {
                m_boundaryConditions[regionID] = bCond; 
            }

            const std::string GetVariable(unsigned int indx)
            {
                return m_session->GetVariable(indx);
            }

        protected:
            /// The mesh graph to use for referencing geometry info.
            MeshGraphSharedPtr                      m_meshGraph;
            LibUtilities::SessionReaderSharedPtr    m_session;

            BoundaryRegionCollection                m_boundaryRegions;
            BoundaryConditionCollection             m_boundaryConditions;

        private:

            /// Read segments (and general MeshGraph) given TiXmlDocument.
            void Read(TiXmlElement *conditions);

            void ReadBoundaryRegions(TiXmlElement *regions);
            void ReadBoundaryConditions(TiXmlElement *conditions);
        };

        typedef boost::shared_ptr<BoundaryConditions> 
            BoundaryConditionsSharedPtr;
    }
}

#endif //NEKTAR_SPATIALDOMAINS_BOUNDARYCONDITIONS_H

