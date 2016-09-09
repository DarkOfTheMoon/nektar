///////////////////////////////////////////////////////////////////////////////
//
// File: PointsGauss.cpp
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
// Description: Specialisations for Gauss-type point distributions
//
///////////////////////////////////////////////////////////////////////////////

#ifndef LIBUTILITIES_FOUNDATIONS_POINTS_POINTSGAUSS
#define LIBUTILITIES_FOUNDATIONS_POINTS_POINTSGAUSS

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
 * @brief Specialisation for GaussGaussLegendre
 */
template<typename TData, typename TShape>
class Points<TData, TShape, GaussGaussLegendre> : public PointsBase<TData>
{
        typedef Points<TData, TShape, GaussGaussLegendre> ThisType;
        typedef PointsBase<TData> BaseType;
        typedef GaussGaussLegendre PointsType;

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

        virtual ~Points() {}

        void Populate(const PointsKey& p) {
            BaseType::m_key = p;
            const int n = p.m_numpoints[0];
            BaseType::template AllocateArrays<PointsType>();
            Polylib::zwgj(BaseType::m_points[0].data(),
                          BaseType::m_weights.data(),
                          n,0.0,0.0);
        }

    private:
        virtual Array<OneD, TData> v_GetInterpMatrixData(const int& npts, const Array<OneD, const TData>& pts)
        {
            Array<OneD, TData> interp(npts * BaseType::m_key.m_numpoints[0]);
            Polylib::Imgj(interp.data(), BaseType::m_points[0].data(),
                          pts.data(), BaseType::m_key.m_numpoints[0], npts, 0.0, 0.0);
            return interp;
        }
};


/**
 * Specialisation for GaussRadauMLegendre
 */
template<typename TData, typename TShape>
class Points<TData, TShape, GaussRadauMLegendre> : public PointsBase<TData>
{
        typedef Points<TData, TShape, GaussRadauMLegendre> ThisType;
        typedef PointsBase<TData> BaseType;
        typedef GaussGaussLegendre PointsType;

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
            Polylib::zwgrjm(BaseType::m_points[0].data(),
                          BaseType::m_weights.data(),
                          n,0.0,0.0);
        }

    private:
        virtual Array<OneD, TData> v_GetInterpMatrixData(const int& npts, const Array<OneD, const TData>& pts)
        {
            Array<OneD, TData> interp(npts * BaseType::m_key.m_numpoints[0]);
            Polylib::Imgrjm(interp.data(), BaseType::m_points[0].data(),
                            pts.data(), BaseType::m_key.m_numpoints[0], npts, 0.0, 0.0);
            return interp;
        }

};

}
}
}

#endif
