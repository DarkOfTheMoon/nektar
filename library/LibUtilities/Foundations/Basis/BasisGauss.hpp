///////////////////////////////////////////////////////////////////////////////
//
// File: BasisGauss.cpp
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
// Description: Specialisations for Gauss-type bases
//
///////////////////////////////////////////////////////////////////////////////


#ifndef LIBRARY_LIBUTILITIES_FOUNDATIONS_BASIS_BASISGAUSS_HPP_
#define LIBRARY_LIBUTILITIES_FOUNDATIONS_BASIS_BASISGAUSS_HPP_

#include <LibUtilities/Foundations/Basis/BasisTemplate.hpp>
#include <LibUtilities/Polylib/Polylib.h>
#include <LibUtilities/LinearAlgebra/Blas.hpp>

namespace Nektar
{
namespace LibUtilities
{
namespace Foundations
{

/**
 * @brief Specialisation for Modified Legendre basis.
 */
template<typename TData, typename TShape, typename TPts>
class Basis<TData, TShape, TPts, ModifiedLegendre> : public BasisBase<TData>
{
        BASIS_DEFINES_VALIDATION(ModifiedLegendre)
        BASIS_CORE_FUNCTIONS(ModifiedLegendre)
        BASIS_CORE_DATA

    public:
        Basis() : BaseType() {}

        Basis(const BasisKey& pKey) : BaseType(pKey), m_points(pKey.m_ptsKey)
        {
            Populate();
        }

        // Populates the data members of the base class for this basis.
        void Populate() {
            BaseType::m_bdata = Array<OneD, TData>(
                    BaseType::m_key.m_ptsKey.m_numpoints[0]*
                    BaseType::m_key.m_nummodes[0]);

            // Get the point distribution and integration weights
            Array<OneD, const NekDouble> z;
            Array<OneD, const NekDouble> w;
            m_points.GetZW(z,w);

            int i         = 0;
            int p         = 0;
            int numPoints = m_points.GetNumPoints();
            int numModes  = BaseType::m_key.m_nummodes[0];

            for(i = 0; i < numPoints; ++i)
            {
                BaseType::m_bdata[i] = 0.5*(1-z[i]);
                BaseType::m_bdata[numPoints + i] = 0.5*(1+z[i]);
            }

            TData* mode = BaseType::m_bdata.data() + 2*numPoints;

            for(p = 2; p < numModes; ++p, mode += numPoints)
            {
                Polylib::jacobfd(numPoints, z.data(), mode, NULL, p-2,1.0,1.0);

                for(i = 0; i < numPoints; ++i)
                {
                    mode[i] *= BaseType::m_bdata[i] *
                               BaseType::m_bdata[numPoints+i];
                }
            }

            // define derivative basis
//            const NekDouble* D = &(m_points.GetD()->GetRawPtr())[0];
//            Blas::Dgemm('n', 'n', numPoints, numModes,  numPoints,
//                        1.0, D,   numPoints,
//                        BaseType::m_bdata.data(),       numPoints,
//                        0.0, BaseType::m_dbdata.data(), numPoints);
        }

};

}
}
}


#endif /* LIBRARY_LIBUTILITIES_FOUNDATIONS_BASIS_BASISGAUSS_HPP_ */
