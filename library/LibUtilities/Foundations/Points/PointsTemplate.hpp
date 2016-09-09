///////////////////////////////////////////////////////////////////////////////
//
// File: PointsTemplate.cpp
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
// Description: Primary template for Points class
//
///////////////////////////////////////////////////////////////////////////////

#ifndef LIBUTILITIES_FOUNDATIONS_POINTS_POINTSTEMPLATE
#define LIBUTILITIES_FOUNDATIONS_POINTS_POINTSTEMPLATE

#include <iostream>
#include <type_traits>

#include <LibUtilities/Foundations/Points.hpp>


namespace Nektar
{
namespace LibUtilities
{
namespace Foundations
{

/**
 * @brief Primary template for Points classes.
 * @tparam TData Data type of the point coordinates and weight values.
 * @tparam TShape The native shape of this point distribution
 * @tparam TPtsTuple A std::tuple of data types containing one or more point
 *                   distributions which form the composite point distribution.
 */
template<typename TData, typename TShape, typename TPtsTuple>
class Points;

/**
 * @brief Composite Points classes.
 * @tparam TData Data type of the point coordinates and weight values.
 * @tparam TShape The native shape of this point distribution
 * @tparam TPts A list of data types containing one or more points distributions
 *                   which form the composite points distribution.
 */
template<typename TData, typename TShape, typename... TPts>
class Points<TData, TShape, std::tuple<TPts...>> : public PointsBase<TData>
{
        /// The type of this class
        typedef Points<TData, TShape, std::tuple<TPts...>> ThisType;
        /// The type of the base class
        typedef PointsBase<TData> BaseType;
        /// The type of the tuple for holding constituent point types
        typedef std::tuple<Points<TData,
                           typename traits::points_traits<TPts>::native_shape,
                           TPts>...
                          > TupleType;


        static_assert(sizeof...(TPts) > 0, "No point type given.");
        static_assert(sizeof...(TPts) < 4, "Too many point types given.");
        static_assert(traits::points_traits<TPts...>::dimension == traits::shape_traits<TShape>::dimension,
                "Points dimension and shape dimension do not agree.");

    public:
        typedef traits::points_traits<TPts...> get_traits;

        static PointsSharedPtr<TData> create(const PointsKey& pKey)
        {
            PointsSharedPtr<TData> p = MemoryManager<ThisType>::AllocateSharedPtr(pKey);
            return p;
        }

        Points() : PointsBase<TData>()  {}

        Points(const PointsKey& pKey) : BaseType(pKey) {
            BaseType::template AllocateArrays<TPts...>();
            InitialiseSubType<sizeof...(TPts)-1>(pKey, x);
        }

        virtual ~Points() {}

    private:
        TupleType x;

        template<size_t i>
        typename std::enable_if<i!=0, void>::type
        InitialiseSubType(const PointsKey& p, TupleType& pTuple) {
            PointsKey tmpKey;
            tmpKey.m_numpoints[0] = p.m_numpoints[i];
            tmpKey.m_params = p.m_params;
            std::get<i>(pTuple).Populate(tmpKey);
            InitialiseSubType<i-1>(p, pTuple);
        }

        template<size_t i>
        typename std::enable_if<i==0, void>::type
        InitialiseSubType(const PointsKey& p, TupleType& pTuple) {
            PointsKey tmpKey;
            tmpKey.m_numpoints[0] = p.m_numpoints[i];
            tmpKey.m_params = p.m_params;
            std::get<i>(pTuple).Populate(tmpKey);
        }

        template<size_t i>
        inline typename std::enable_if<i < sizeof...(TPts), PointsBase<TData>*>::type
        GetTupleEntryPtr()
        {
            return &std::get<i>(x);
        }

        template<size_t i>
        inline typename std::enable_if<i >= sizeof...(TPts), PointsBase<TData>*>::type
        GetTupleEntryPtr()
        {
            throw std::logic_error("Constituent points index out of range.");
        }

        virtual const Array<OneD, const TData> v_GetZ(const int& i) {
                switch (i)
                {
                    case 0: return GetTupleEntryPtr<0>()->GetZ();
                    case 1: return GetTupleEntryPtr<1>()->GetZ();
                    case 2: return GetTupleEntryPtr<2>()->GetZ();
                    default: throw i;
                }
        }

        virtual Array<OneD, TData> v_GetInterpMatrixData(const int& npts, const Array<OneD, const TData>& pts)
        {
            Array<OneD, TData> interp(0);
            return interp;
        }
};

}
}
}

#endif
