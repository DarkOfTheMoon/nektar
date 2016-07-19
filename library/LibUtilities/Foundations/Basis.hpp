#ifndef LIBUTILITIES_FOUNDATIONS_BASIS
#define LIBUTILITIES_FOUNDATIONS_BASIS

#include <iostream>
#include <type_traits>

#include <loki/Singleton.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/LibUtilitiesDeclspec.h>
#include <LibUtilities/Polylib/Polylib.h>
#include <LibUtilities/Foundations/ShapeTypes.hpp>
#include <LibUtilities/Foundations/PointsTypes.hpp>
#include <LibUtilities/Foundations/BasisTypes.hpp>
#include <LibUtilities/Foundations/Points.hpp>

template <class T>
struct is_not_tuple { static const bool value = true; };

template <typename... Args>
struct is_not_tuple<std::tuple<Args...>> { static const bool value = false; };

template <typename... Args>
struct is_not_tuple<const std::tuple<Args...>> {
        static const bool value = false;
};

namespace Nektar
{
namespace LibUtilities
{
namespace Foundations
{

typedef double NekDouble;
typedef std::string BasisParamKey;
typedef NekDouble BasisParamValue;
typedef std::map<BasisParamKey, BasisParamValue> BasisParamList;

class BasisKey
{
    public:
      unsigned int m_dim;
      unsigned int m_nummodes[3];
      BasisParamList m_params;
      PointsKey m_ptsKey;
};


/**
 * @class BasisBase
 * @brief Basis base class defining interface and data members
 */
template<typename TData>
class BasisBase
{
    public:
        /// Destructor.
        virtual ~BasisBase()
        {
        };

        /// Return order of basis from the basis specification.
        inline int GetNumModes() const
        {
            return m_key.m_nummodes[0];
        }

        /// Return total number of modes from the basis specification.
        virtual int GetTotNumModes() const
        {
            return m_key.m_nummodes[0];
        }

    protected:
        BasisKey                m_key;

        BasisBase(const BasisKey& pKey) {m_key = pKey;}
        BasisBase() {}

        template<typename TBasis1, typename TBasis2, typename... TBasisOther>
        unsigned int GetNumberOfModes(int i = 0)
        {
            return traits::basis_traits<TBasis1>::get_total_modes(m_key.m_nummodes[i])
                    * GetNumberOfModes<TBasis2, TBasisOther...>(i + 1);
        }

        template<typename TBasis>
        unsigned int GetNumberOfModes(int i = 0)
        {
            return traits::basis_traits<TBasis>::get_total_modes(m_key.m_nummodes[i]);
        }

};


template<typename TData>
using BasisSharedPtr = boost::shared_ptr<BasisBase<TData>>;

template<typename TData>
using BasisFactory = LibUtilities::NekFactory<
        std::string, BasisBase<TData>, const BasisKey&>;

template<typename TData>
LIB_UTILITIES_EXPORT BasisFactory<TData>& GetBasisFactory()
{
    typedef Loki::SingletonHolder<BasisFactory<TData>,
                                  Loki::CreateUsingNew,
                                  Loki::NoDestroy> Type;
    // Putting Loki::ClassLevelLockable causes an assertion!
    return Type::Instance();
}

/**
 * @class Basis
 * @brief Primary template for composite Points classes.
 */
template<typename TData, typename TShape, typename... TTuple>
class Basis : public BasisBase<TData> {};

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

    protected:
        Points<TData, TShape, std::tuple<TPts...>> m_points;
        TupleType x;
};


/**
 * Bi-Basis specialisation
 */
//template<typename TData, typename TShape, template<typename, typename, typename...> class TPts, typename TBasis1, typename TBasis2>
//class Basis<TData, TShape, TPts, TBasis1, TBasis2> : public BasisBase<TData>
//{
//        typedef Points<TData, TShape, TPts, TBasis1, TBasis2> ThisType;
//        typedef PointsBase<TData> BaseType;
//
//        static_assert(TPts::get_traits::dimension == traits::shape_traits<TShape>::dimension,
//                "Points dimension does not match shape dimension,");
//        static_assert(traits::basis_traits<TBasis1, TBasis2>::dimension == traits::shape_traits<TShape>::dimension,
//                "Basis dimension and shape dimension do not agree.");
//
//    public:
//        static BasisSharedPtr<TData> create(const BasisKey& pKey)
//        {
//            BasisSharedPtr<TData> p = MemoryManager<ThisType>::AllocateSharedPtr(pKey);
//            return p;
//        }
//
//        Basis(const BasisKey& pKey) : BaseType(pKey) {
//            x1.Populate(pKey);
//            BasisKey tmpKey;
//            tmpKey.m_nummodes[0] = pKey.m_nummodes[1];
//            tmpKey.m_params = pKey.m_params;
//            x2.Populate(tmpKey);
//        }
//
//        /// Return total number of modes from the basis specification.
//        virtual int GetTotNumModes() const
//        {
//            return BaseType::template GetNumberOfModes<TBasis1, TBasis2>();
//        }
//
//
//    private:
//        TPts m_points;
//        Basis<TData, typename traits::basis_traits<TBasis1>::native_shape, TBasis1> x1;
//        Basis<TData, typename traits::basis_traits<TBasis2>::native_shape, TBasis2> x2;
//};


/**
 * Tri-Basis specialisation
 */
//template<typename TData, typename TShape, template<typename, typename, typename...> class TPts, typename TBasis1, typename TBasis2, typename TBasis3>
//class Basis<TData, TShape, TPts, TBasis1, TBasis2, TBasis3> : public BasisBase<TData>
//{
//        typedef Basis<TData, TShape, TPts, TBasis1, TBasis2, TBasis3> ThisType;
//        typedef BasisBase<TData> BaseType;
//
//        static_assert(TPts::get_traits::dimension == traits::shape_traits<TShape>::dimension,
//                "Points dimension does not match shape dimension,");
//        static_assert(traits::basis_traits<TBasis1, TBasis2, TBasis3>::dimension == traits::shape_traits<TShape>::dimension,
//                "Basis dimension and shape dimension do not agree.");
//
//    public:
//        static BasisSharedPtr<TData> create(const BasisKey& pKey)
//        {
//            BasisSharedPtr<TData> p = MemoryManager<ThisType>::AllocateSharedPtr(pKey);
//            return p;
//        }
//
//        Basis(const BasisKey& pKey) : BaseType(pKey) {
//            //BasisBase<TData>::template AllocateArrays<TPts1, TPts2, TPts3>();
//            BasisKey tmpKey;
//            x1.Populate(pKey);
//            tmpKey.m_nummodes[0] = pKey.m_nummodes[1];
//            tmpKey.m_params = pKey.m_params;
//            x2.Populate(tmpKey);
//            tmpKey.m_nummodes[0] = pKey.m_nummodes[2];
//            tmpKey.m_params = pKey.m_params;
//            x3.Populate(tmpKey);
//        }
//
//    private:
//        TPts m_points;
//        Basis<TData, typename traits::basis_traits<TBasis1>::native_shape, TBasis1> x1;
//        Basis<TData, typename traits::basis_traits<TBasis2>::native_shape, TBasis2> x2;
//        Basis<TData, typename traits::basis_traits<TBasis3>::native_shape, TBasis3> x3;
//};


/**
 *
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

    private:
        Points<TData, typename traits::points_traits<TPts>::native_shape, TPts> m_points;

};


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

    private:
        Points<TData, typename traits::points_traits<TPts>::native_shape, TPts> m_points;

};
}
}
}

#endif
