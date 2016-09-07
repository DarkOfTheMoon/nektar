/*
 * BasisGauss.hpp
 *
 *  Created on: 2 Sep 2016
 *      Author: cc
 */

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
 * @class Basis
 * @brief Specialisation for Modified Legendre basis.
 */
template<typename TData, typename TShape, typename TPts>
class Basis<TData, TShape, TPts, ModifiedLegendre> : public BasisBase<TData>
{
        typedef Basis<TData, TShape, TPts, ModifiedLegendre> ThisType;
        typedef BasisBase<TData> BaseType;

        static_assert(is_not_tuple<TPts>::value, "Is a tuple.");
        static_assert(traits::points_traits<TPts>::dimension == traits::shape_traits<TShape>::dimension,
                "Points dimension does not match shape dimension,");
        static_assert(traits::basis_traits<ModifiedLegendre>::dimension == traits::shape_traits<TShape>::dimension,
                "Basis dimension and shape dimension do not agree.");
        static_assert(traits::expansion_traits<TShape, ModifiedLegendre>::is_valid,
                "Not a valid combination of shape and basis type.");

    public:
        /// Creates an instance of this basis
        static BasisSharedPtr<TData> create(const BasisKey& pKey)
        {
            BasisSharedPtr<TData> p = MemoryManager<ThisType>::AllocateSharedPtr(pKey);
            return p;
        }

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

    protected:
        virtual int v_GetNumConstituentBases() const
        {
            return traits::basis_traits<ModifiedLegendre>::num_constituent_bases;
        }

        virtual const PointsBase<TData>& v_GetPoints() const
        {
            return m_points;
        }

        virtual LibUtilities::Foundations::ShapeType v_GetShapeType() const
        {
            return traits::shape_traits<TShape>::type;
        }

        virtual std::string v_GetShapeName() const
        {
            return std::string(traits::shape_traits<TShape>::name);
        }

    private:
        Points<TData, typename traits::points_traits<TPts>::native_shape, TPts> m_points;

};

}
}
}


#endif /* LIBRARY_LIBUTILITIES_FOUNDATIONS_BASIS_BASISGAUSS_HPP_ */
