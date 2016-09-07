#ifndef LIBRARY_LIBUTILITIES_FOUNDATIONS_BASIS_BASISTEMPLATE_HPP_
#define LIBRARY_LIBUTILITIES_FOUNDATIONS_BASIS_BASISTEMPLATE_HPP_

#include <LibUtilities/Foundations/Basis.hpp>

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/LibUtilitiesDeclspec.h>
#include <LibUtilities/Foundations/Foundations.hpp>
#include <LibUtilities/Foundations/ShapeTypes.hpp>
#include <LibUtilities/Foundations/Points/PointsTypes.hpp>
#include <LibUtilities/Foundations/Basis/BasisTypes.hpp>
#include <LibUtilities/Foundations/Points.hpp>

namespace Nektar
{
namespace LibUtilities
{
namespace Foundations
{

/**
 * @class Basis
 * @brief Primary template for composite Points classes.
 *
 * No implementation is provided, since this non-specialised template should
 * never be instantiated.
 */
template<typename TData, typename TShape, typename... TTuple>
class Basis;






/**
 * @class Basis
 * @brief Specialisation for composite basis.
 *
 * A composite basis is a multi-dimensional basis made up of two or more
 * basis types. For each constituent basis type, the corresponding points type
 * of matching dimension must also be provided.
 */
template<typename TData, typename TShape, typename... TPts, typename... TBasis>
class Basis<TData, TShape, std::tuple<TPts...>, std::tuple<TBasis...>> : public BasisBase<TData>
{
        typedef Basis<TData, TShape, std::tuple<TPts...>, std::tuple<TBasis...>> ThisType;
        typedef BasisBase<TData> BaseType;
        typedef std::tuple<Basis<TData, typename traits::basis_traits<TBasis>::native_shape, TPts, TBasis>...> TupleType;

        static_assert(sizeof...(TBasis) > 0, "No basis type given.");
        static_assert(sizeof...(TBasis) < 4, "Too many basis types given. Must be 3 or less.");
        static_assert(sizeof...(TPts) == sizeof...(TBasis), "Basis and Points types do not match.");
        static_assert(traits::basis_traits<TBasis...>::dimension == traits::shape_traits<TShape>::dimension,
                "Basis dimension and shape dimension do not agree.");
        static_assert(traits::points_traits<TPts...>::dimension == traits::shape_traits<TShape>::dimension,
                "Points dimension does not match shape dimension,");

    public:
        static BasisSharedPtr<TData> create(const BasisKey& pKey)
        {
            BasisSharedPtr<TData> b = MemoryManager<ThisType>::AllocateSharedPtr(pKey);
            return b;
        }

        Basis(const BasisKey& pKey) : BasisBase<TData>(pKey) {
            //BasisBase<TData>::template AllocateArrays<TBasis...>();
        }

        virtual int v_GetNumConstituentBases() const
        {
            return traits::basis_traits<TBasis...>::num_constituent_bases;
        }

        virtual const BasisBase<TData>& v_GetConstituentBasis(int i) const
        {
            switch (i)
            {
                case 0: return GetTupleEntry<0>();
                case 1: return GetTupleEntry<1>();
                case 2: return GetTupleEntry<2>();
                default: throw i;
            }
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

    protected:
        Points<TData, TShape, std::tuple<TPts...>> m_points;
        TupleType x;

    private:
        template<size_t i>
        inline typename std::enable_if<i < sizeof...(TBasis), const BasisBase<TData>&>::type
        GetTupleEntry() const
        {
            return std::get<i>(x);
        }

        template<size_t i>
        inline typename std::enable_if<i >= sizeof...(TBasis), const BasisBase<TData>&>::type
        GetTupleEntry() const
        {
            throw std::logic_error("Out of range.");
        }
};

}
}
}
#endif
