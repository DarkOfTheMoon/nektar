///////////////////////////////////////////////////////////////////////////////
//
// File: CollectionOptimisation.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: Collection top class definition
//
///////////////////////////////////////////////////////////////////////////////

#include <Collections/CollectionOptimisation.h>
#include <LibUtilities/BasicUtils/ParseUtils.hpp>

namespace Nektar 
{
    namespace Collections 
    {        

        // static manager for Operator ImplementationMap
        map<OpImpTimingKey,OperatorImpMap> CollectionOptimisation::m_opImpMap;

        CollectionOptimisation::CollectionOptimisation(
                               LibUtilities::SessionReaderSharedPtr pSession,
                               ImplementationType defaultType)
        {
            map<ElmtOrder, ImplementationType> defaults;
            map<ElmtOrder, ImplementationType>::iterator it;
            
            m_setByXml = false;
            
            // Default all elements to eStdMat
            defaults[ElmtOrder(LibUtilities::eSegment,       -1)] = defaultType;
            defaults[ElmtOrder(LibUtilities::eTriangle,      -1)] = defaultType;
            defaults[ElmtOrder(LibUtilities::eQuadrilateral, -1)] = defaultType;
            defaults[ElmtOrder(LibUtilities::eTetrahedron,   -1)] = defaultType;
            defaults[ElmtOrder(LibUtilities::ePyramid,       -1)] = defaultType;
            defaults[ElmtOrder(LibUtilities::ePrism,         -1)] = defaultType;
            defaults[ElmtOrder(LibUtilities::eHexahedron,    -1)] = defaultType;
            
            
            map<string, LibUtilities::ShapeType> elTypes;
            map<string, LibUtilities::ShapeType>::iterator it2;
            elTypes["S"] = LibUtilities::eSegment;
            elTypes["T"] = LibUtilities::eTriangle;
            elTypes["Q"] = LibUtilities::eQuadrilateral;
            elTypes["A"] = LibUtilities::eTetrahedron;
            elTypes["P"] = LibUtilities::ePyramid;
            elTypes["R"] = LibUtilities::ePrism;
            elTypes["H"] = LibUtilities::eHexahedron;
            
            map<string, OperatorType> opTypes;
            for (int i = 0; i < SIZE_OperatorType; ++i)
            {
                opTypes[OperatorTypeMap[i]] = (OperatorType)i;
                m_global[(OperatorType)i] = defaults;
            }
            
            map<string, ImplementationType> impTypes;
            for (int i = 0; i < SIZE_ImplementationType; ++i)
            {
                impTypes[ImplementationTypeMap[i]] = (ImplementationType)i;
            }
            
            if(pSession.get()) // turn off file reader if dummy pointer is given
            {
                TiXmlDocument &doc = pSession->GetDocument();
                TiXmlHandle docHandle(&doc);
                TiXmlElement *master
                    = docHandle.FirstChildElement("NEKTAR").Element();
                ASSERTL0(master, "Unable to find NEKTAR tag in file.");
                
                TiXmlElement *xmlCol = master->FirstChildElement("COLLECTIONS");
                
                if (xmlCol)
                {
                    TiXmlElement *elmt = xmlCol->FirstChildElement();
                    m_setByXml = true;
                    
                    while (elmt)
                    {
                        string tagname = elmt->ValueStr();
                        
                        ASSERTL0(boost::iequals(tagname, "OPERATOR"),
                                 "Only OPERATOR tags are supported inside the "
                                 "COLLECTIONS tag.");
                        
                        const char *attr = elmt->Attribute("TYPE");
                        ASSERTL0(attr, "Missing TYPE in OPERATOR tag.");
                        string opType(attr);
                        
                        ASSERTL0(opTypes.count(opType) > 0,
                                 "Unknown OPERATOR type " + opType + ".");
                        
                        OperatorType ot = opTypes[opType];
                        
                        TiXmlElement *elmt2 = elmt->FirstChildElement();
                        
                        while (elmt2)
                        {
                            string tagname = elmt2->ValueStr();
                            ASSERTL0(boost::iequals(tagname, "ELEMENT"),
                                     "Only ELEMENT tags are supported inside the "
                                     "OPERATOR tag.");
                            
                            const char *attr = elmt2->Attribute("TYPE");
                            ASSERTL0(attr, "Missing TYPE in ELEMENT tag.");
                            
                            string elType(attr);
                            it2 = elTypes.find(elType);
                            ASSERTL0(it2 != elTypes.end(),
                                     "Unknown element type "+elType+" in ELEMENT "
                                     "tag");
                            
                            const char *attr2 = elmt2->Attribute("IMPTYPE");
                            ASSERTL0(attr2, "Missing IMPTYPE in ELEMENT tag.");
                            string impType(attr2);
                            ASSERTL0(impTypes.count(impType) > 0,
                                     "Unknown IMPTYPE type " + impType + ".");
                            
                            const char *attr3 = elmt2->Attribute("ORDER");
                            ASSERTL0(attr3, "Missing ORDER in ELEMENT tag.");
                            string order(attr3);
                            
                            if (order == "*")
                            {
                                m_global[ot][ElmtOrder(it2->second, -1)] 
                                    = impTypes[impType];
                            }
                            else
                            {
                                vector<unsigned int> orders;
                                ParseUtils::GenerateSeqVector(order.c_str(), orders);
                                
                                for (int i = 0; i < orders.size(); ++i)
                                {
                                    m_global[ot][ElmtOrder(it2->second, orders[i])] = impTypes[impType];
                                }
                            }
                            
                            elmt2 = elmt2->NextSiblingElement();
                        }
                        
                        elmt = elmt->NextSiblingElement();
                    }

		    // Print out operator map
		    if(pSession->DefinesCmdLineArgument("verbose"))
		    {
		        if(pSession->GetComm()->GetRank() == 0)
			{
			  map<OperatorType, map<ElmtOrder, ImplementationType> >::iterator mIt;
			  map<ElmtOrder, ImplementationType>::iterator eIt;
			  for (mIt = m_global.begin(); mIt != m_global.end(); mIt++)
			    {
			      cout << "Operator " << OperatorTypeMap[mIt->first] << ":"
				   << endl;
			      
			      for (eIt = mIt->second.begin(); eIt != mIt->second.end(); eIt++)
				{
				  cout << "- " << LibUtilities::ShapeTypeMap[eIt->first.first]
				       << " order " << eIt->first.second
				       << " -> " << ImplementationTypeMap[eIt->second]
				       << endl;
				}
			    }
			}
		    }
		}
		
            }
        }


        OperatorImpMap  CollectionOptimisation::GetOperatorImpMap(StdRegions::StdExpansionSharedPtr pExp)
        {
            map<OperatorType, map<ElmtOrder, ImplementationType> >::iterator it;
            map<ElmtOrder, ImplementationType>::iterator it2;

            OperatorImpMap ret;
            ElmtOrder searchKey(pExp->DetShapeType(),
                                pExp->GetBasisNumModes(0));
            ElmtOrder defSearch(pExp->DetShapeType(), -1);

            for (it = m_global.begin(); it != m_global.end(); ++it)
            {
                ImplementationType impType;

                it2 = it->second.find(searchKey);

                if (it2 == it->second.end())
                {
                    it2 = it->second.find(defSearch);
                    if (it2 == it->second.end())
                    {
                        cout << "shouldn't be here..." << endl;
                        impType = eStdMat;
                    }
                    else
                    {
                        impType = it2->second;
                    }
                }
                else
                    {
                        impType = it2->second;
                    }
                    
                    ret[it->first] = impType;
                }

                return ret;
            }
            
        OperatorImpMap CollectionOptimisation::SetWithTimings(
                            StdRegions::StdExpansionSharedPtr pExp,
                            vector<SpatialDomains::GeometrySharedPtr> pGeom,
                            OperatorImpMap &impTypes,
                            bool verbose )
            {
                OperatorImpMap ret;

                // check to see if already defined for this expansion
                OpImpTimingKey OpKey(pExp,pGeom.size(),pExp->GetNumBases());
                if(m_opImpMap.count(OpKey) != 0)
                {
                    ret = m_opImpMap[OpKey];
                    return ret;
                }
#if 0 
		map<OpImpTimingKey,OperatorImpMap>::iterator it;
	        cout << "OpKey: " << OpKey.m_nbasis << " " << OpKey.m_ngeoms << endl;
	  	int cnt;
		cout << m_opImpMap.size() << endl;
	        for(cnt = 1, it = m_opImpMap.begin(); it != m_opImpMap.end(); ++it, ++cnt)
		{
			cout <<  cnt << ": " << it->first.m_nbasis << " " << it->first.m_ngeoms << endl;
		}
#endif
	
                int maxsize = pGeom.size()*max(pExp->GetNcoeffs(),pExp->GetTotPoints());
                Array<OneD, NekDouble> inarray(maxsize,1.0);
                Array<OneD, NekDouble> outarray1(maxsize);
                Array<OneD, NekDouble> outarray2(maxsize);
                Array<OneD, NekDouble> outarray3(maxsize);

                Timer t;

                if(verbose)
                {
                    cout << "Collection Implemenation for " << LibUtilities::ShapeTypeMap[pExp->DetShapeType()];
                    cout << " ( ";
                    for(int i = 0; i < pExp->GetNumBases(); ++i)
                    {
                        cout << pExp->GetBasis(i)->GetNumModes() <<" ";
                    }
                    cout << ")" <<  " for ngeoms = " << pGeom.size() << endl;
                }
                // set  up an array of collections
                CollectionVector coll; 
                for(int imp = 1; imp < SIZE_ImplementationType; ++imp)
                {
                    OperatorImpMap impTypes = SetFixedImpType((ImplementationType) imp);
                                                                       
                    Collection collloc(pExp,pGeom,impTypes);
                    coll.push_back(collloc);
                }
                    
                // Determine the number of tests to do in one second
                Array<OneD, int> Ntest(SIZE_OperatorType);
                for(int i = 0; i < SIZE_OperatorType; ++i)
                {
                    OperatorType OpType = (OperatorType)i;

                    t.Start();
                    coll[0].ApplyOperator(OpType,
                                       inarray,
                                       outarray1,
                                       outarray2,
                                       outarray3);
                    t.Stop();
                    
                    NekDouble oneTest = t.TimePerTest(1);
                    
                    Ntest[i] = max((int)(0.25/oneTest),1);
                }

                Array<OneD, NekDouble> timing(SIZE_ImplementationType);
                // loop over all operators and determine fastest implementation
                for(int i = 0; i < SIZE_OperatorType; ++i)
                {
                    OperatorType OpType = (OperatorType)i;

                    // call collection implementation in thorugh ExpList. 
                    for (int imp = 0; imp < coll.size(); ++imp)
                    {
                        t.Start();
                        for(int n = 0; n < Ntest[i]; ++n)
                        {
                            coll[imp].ApplyOperator(OpType,
                                                  inarray,
                                                  outarray1,
                                                  outarray2,
                                                  outarray3);
                        }
                        t.Stop();
                        timing[imp] = t.TimePerTest(Ntest[i]);
                    }
                    // determine optimal implementation. Note +1 to
                    // remove NoImplementationType flag
                    int minImp = Vmath::Imin(coll.size(),timing,1)+1;

                    if(verbose)
                    {
                        cout << "\t " << OperatorTypeMap[i] << ": " << ImplementationTypeMap[minImp]; 
                        cout << "\t (";
                        for(int j = 0; j < coll.size(); ++j)
                        {
                            cout << timing[j] ;
                            if(j != coll.size()-1)
                            {
                                cout <<", ";
                            }
                        }
                        cout << ")" <<endl;
                    }
                    // could reset global map if reusing  method? 
                    //m_global[OpType][pExp->DetShapeType()] = (ImplementationType)minImp;
                    // set up new map
                    ret[OpType] = (ImplementationType)minImp;
                }

                // store map for use by another expansion. 
                m_opImpMap[OpKey] = ret;
                return ret;
            }
    }
}
