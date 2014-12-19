///////////////////////////////////////////////////////////////////////////////
//
// File: Collection.h 
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

#ifndef NEKTAR_LIBRARY_COLLECTIONS_COLLECTIONOPTIMISATION_H
#define NEKTAR_LIBRARY_COLLECTIONS_COLLECTIONOPTIMISATION_H

#include <Collections/Collection.h>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <tinyxml.h>

namespace Nektar 
{
    namespace Collections 
    {        
        class OpImpTimingKey
        {
        public:
            /// Constructor
            OpImpTimingKey(StdRegions::StdExpansionSharedPtr pExp, int ngeoms, int nbases):
                    m_exp(pExp),
                    m_ngeoms(ngeoms),
                    m_nbasis(nbases)
            {
            }
            
            /// Destructor
                ~OpImpTimingKey(void)
            {
            }
            

            bool operator<(const OpImpTimingKey &rhs) const
            {

                if(m_nbasis <   rhs.m_nbasis) 
                {
                    return true; 
                }

                if(m_nbasis > rhs.m_nbasis)
                {
                    return false; 
                }
                
                for(int i = 0; i < m_nbasis; ++i)
                {
                    if((m_exp->GetBasis(i)->GetBasisKey() 
                        != rhs.m_exp->GetBasis(i)->GetBasisKey()))
                    {
		      return (m_exp->GetBasis(i)->GetBasisKey() < rhs.m_exp->GetBasis(i)->GetBasisKey());
                    }
                }
                
                if((m_ngeoms < 100)&&(rhs.m_ngeoms < 100))
                {

                    if(m_ngeoms < rhs.m_ngeoms)
                    {
                        return true;
                    }
                    else
                    {
                        return false;
                    }
                }
                
                return false;
            }

            StdRegions::StdExpansionSharedPtr m_exp; 
            int m_ngeoms; 
            int m_nbasis; 
        private:
        };
        
        
        class CollectionOptimisation
        {
            typedef pair<LibUtilities::ShapeType, int> ElmtOrder;

        public:
            // Constuctor 
            CollectionOptimisation( LibUtilities::SessionReaderSharedPtr pSession,
                                    ImplementationType defaultType = eStdMat);

            ~CollectionOptimisation(){};

            
            /// Get Operator Implementation Map from XMl or using default; 
            OperatorImpMap  GetOperatorImpMap(StdRegions::StdExpansionSharedPtr pExp);

            // Get Map by doing autotuning testing. 
            OperatorImpMap SetWithTimings(vector<StdRegions::StdExpansionSharedPtr> pGeom,
                                          OperatorImpMap &impTypes,
                                          bool verbose = true);


                
            bool SetByXml(void)
            {
                return m_setByXml;
            }

        private:
            static map<OpImpTimingKey,OperatorImpMap> m_opImpMap;
            map<OperatorType, map<ElmtOrder, ImplementationType> > m_global;
            bool m_setByXml;
        };
    }
}
#endif
