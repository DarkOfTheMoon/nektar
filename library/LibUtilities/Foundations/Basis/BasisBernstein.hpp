#ifndef LIBRARY_LIBUTILITIES_FOUNDATIONS_BASIS_BASISBERNSTEIN_HPP_
#define LIBRARY_LIBUTILITIES_FOUNDATIONS_BASIS_BASISBERNSTEIN_HPP_

#include <LibUtilities/Foundations/Basis/BasisTemplate.hpp>

namespace Nektar
{
namespace LibUtilities
{
namespace Foundations
{


/**
 *
 */
template<typename TData, typename TShape, typename TPts>
class Basis<TData, TShape, TPts, BernsteinTriangle> : public BasisBase<TData>
{
        typedef Basis<TData, TShape, TPts, BernsteinTriangle> ThisType;
        typedef BasisBase<TData> BaseType;

        static_assert(traits::points_traits<TPts>::dimension == traits::shape_traits<TShape>::dimension,
                "Points dimension does not match shape dimension,");
        static_assert(traits::basis_traits<BernsteinTriangle>::dimension == traits::shape_traits<TShape>::dimension,
                "Basis dimension and shape dimension do not agree.");
        static_assert(traits::expansion_traits<TShape, BernsteinTriangle>::is_valid,
                "Not a valid combination of shape and basis type.");

    public:
        static BasisSharedPtr<TData> create(const BasisKey& pKey)
        {
            BasisSharedPtr<TData> p = MemoryManager<ThisType>::AllocateSharedPtr(pKey);
            return p;
        }

        Basis() : BaseType() {}

        Basis(const BasisKey& pKey) : BaseType(pKey)
        {
            Populate(pKey);
        }
        void Populate(const BasisKey& p) {
            BaseType::m_key = p;
            //const int n = p.m_nummodes[0];
            //BasisBase<TData>::template AllocateArrays<ModifiedLegendre>();
        }

    protected:
        virtual int v_GetNumConstituentBases() const
        {
            return traits::basis_traits<BernsteinTriangle>::num_constituent_bases;
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
#endif
