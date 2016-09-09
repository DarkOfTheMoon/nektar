///////////////////////////////////////////////////////////////////////////////
//
// File: PointsFekete.cpp
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
// Description: Specialisations for Fekete point distribution
//
///////////////////////////////////////////////////////////////////////////////

#ifndef LIBUTILITIES_FOUNDATIONS_POINTS_POINTSFEKETE
#define LIBUTILITIES_FOUNDATIONS_POINTS_POINTSFEKETE

#include <iostream>
#include <type_traits>

#include <LibUtilities/Foundations/Points/PointsTemplate.hpp>

namespace Nektar
{
namespace LibUtilities
{
namespace Foundations
{

/**
 * @brief Specialisation for Fekete
 */
template<typename TData, typename TShape>
class Points<TData, TShape, Fekete> : public PointsBase<TData>
{
        typedef Points<TData, TShape, Fekete> ThisType;
        typedef PointsBase<TData> BaseType;
        typedef Fekete PointsType;

        static_assert(traits::points_traits<PointsType>::dimension == traits::shape_traits<TShape>::dimension,
                "Points dimension and shape dimension do not agree.");
        static_assert(traits::distribution_traits<TShape, PointsType>::is_valid,
                "Not a valid combination of shape and points type.");

    public:
        typedef traits::points_traits<PointsType> get_traits;

        static PointsSharedPtr<TData> create(const PointsKey& pKey)
        {
            PointsSharedPtr<TData> p = MemoryManager<ThisType>::AllocateSharedPtr(pKey);
            return p;
        }

        Points() : PointsBase<TData>() {}

        Points(const PointsKey& pKey) : PointsBase<TData>(pKey)
        {
            Populate(pKey);
        }

        void Populate(const PointsKey& p) {
            BaseType::m_key = p;
            const int n = p.m_numpoints[0];
            BaseType::template AllocateArrays<PointsType>();
            std::cout << "Need to generate Fekete points" << std::endl;
        }

    private:
        virtual Array<OneD, TData> v_GetInterpMatrixData(const int& npts, const Array<OneD, const TData>& pts)
        {
            Array<OneD, TData> interp(npts * BaseType::m_key.m_numpoints[0]);
            std::cout << "Need to generate Fekete interp" << std::endl;
            return interp;
        }

};

}
}
}

#endif
