///////////////////////////////////////////////////////////////////////////////
//
// File: Operator.h
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
// Description: Operator top class definition
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBRARY_COLLECTIONS_OPERATOR_H
#define NEKTAR_LIBRARY_COLLECTIONS_OPERATOR_H

#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <StdRegions/StdExpansion.h>
#include <SpatialDomains/Geometry.h> 

#define OPERATOR_CREATE(cname)                                  \
    static OperatorKey m_type;                                  \
    static OperatorKey m_typeArr[];                             \
    static OperatorSharedPtr create(                            \
        StdRegions::StdExpansionSharedPtr pExp,                 \
        vector<SpatialDomains::GeometrySharedPtr> pGeom,        \
        boost::shared_ptr<CoalescedGeomData> GeomData)          \
    {                                                           \
        return MemoryManager<cname>                             \
            ::AllocateSharedPtr(pExp, pGeom, GeomData);         \
    }
 
namespace Nektar {
    namespace Collections {

        class CoalescedGeomData;

        enum OperatorType
        {
            eBwdTrans,
            eIProductWRTBase,
            eIProductWRTDerivBase,
            ePhysDeriv,
            SIZE_OperatorType
        };
        
        const char* const OperatorTypeMap[] =
        {
            "BwdTrans",
            "IProductWRTBase",
            "IProductWRTDerivBase",
            "PhysDeriv"
        };

        enum ImplementationType
        {
            eNoImpType,
            eIterPerExp,
            eStdMat,
            eSumFac,
            SIZE_ImplementationType
        };
        
        const char* const ImplementationTypeMap[] =
        {
            "NoImplementationType",
            "IterPerExp",
            "StdMat",
            "SumFac"
        };

        typedef map<OperatorType, ImplementationType> OperatorImpMap;
        
        /// simple Operator Implementation Map generator
        OperatorImpMap SetFixedImpType(ImplementationType defaultType);
        
        class Operator;
        typedef boost::shared_ptr<Operator> OperatorSharedPtr;
        
        typedef boost::tuple<
            LibUtilities::ShapeType, OperatorType, ImplementationType, bool> OperatorKey;
        bool operator< (OperatorKey const &p1, OperatorKey const &p2);
        std::ostream &operator<<(std::ostream &os, OperatorKey const &p);

        typedef Nektar::LibUtilities::NekFactory<
            OperatorKey,
            Operator,
            StdRegions::StdExpansionSharedPtr,
            vector<SpatialDomains::GeometrySharedPtr>,
            boost::shared_ptr<CoalescedGeomData> > OperatorFactory;
        OperatorFactory& GetOperatorFactory();
        
        class Operator
        {
        public:
        Operator(StdRegions::StdExpansionSharedPtr pExp,
                 vector<SpatialDomains::GeometrySharedPtr> pGeom,
                 boost::shared_ptr<CoalescedGeomData> GeomData)
            : m_stdExp (pExp),
                m_numElmt(pGeom.size()),
                m_wspSize(0)
                {
                }
            
            virtual void operator()(const Array<OneD, const NekDouble> &input,
                                    Array<OneD,       NekDouble> &output0,
                                    Array<OneD,       NekDouble> &output1,
                                    Array<OneD,       NekDouble> &output2,
                                    Array<OneD,       NekDouble> &wsp = NullNekDouble1DArray) = 0;

            virtual ~Operator(void);
            
            int GetWspSize()
            {
                return m_wspSize;
            }
            
        protected:
            StdRegions::StdExpansionSharedPtr m_stdExp;
            unsigned int m_numElmt;
            unsigned int m_wspSize;
        };
    }
}
#endif
