#include <iostream>
#include <type_traits>

#include <loki/Singleton.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/LibUtilitiesDeclspec.h>
#include <LibUtilities/Polylib/Polylib.h>
#include <LibUtilities/FoundationsNew/ShapeTypes.hpp>
#include <LibUtilities/FoundationsNew/PointsTypes.hpp>

namespace Nektar
{
namespace LibUtilities
{
namespace Foundations
{

typedef double NekDouble;
typedef std::string PointsParamKey;
typedef NekDouble PointsParamValue;
typedef std::map<PointsParamKey, PointsParamValue> PointsParamList;

class PointsKey
{
    public:
      unsigned int m_numpoints[3];
      PointsParamList m_params;
};


/**
 * @class PointsBase
 * @brief Points base class defining interface and data members
 */
template<typename TData>
class PointsBase
{
    public:
        virtual ~PointsBase()
        {

        }

        inline unsigned int GetNumPoints() const
        {
            return m_key.m_numpoints[0];
        }

        virtual const Array<OneD, const TData> GetZ(const int& i = 0)
        {
            return m_points[0];
        }

        inline const Array<OneD, const TData>& GetW()
        {
            return m_weights;
        }

        inline void GetZW(Array<OneD, const TData> &z,
            Array<OneD, const TData> &w) const
        {
            z = m_points[0];
            w = m_weights;
        }

        inline void GetPoints(Array<OneD, const TData> &x) const
        {
            x = m_points[0];
        }

        inline void GetPoints(Array<OneD, const TData> &x,
                              Array<OneD, const TData> &y) const
        {
            x = m_points[0];
            y = m_points[1];
        }

        inline void GetPoints(Array<OneD, const TData> &x,
                              Array<OneD, const TData> &y,
                              Array<OneD, const TData> &z) const
        {
            x = m_points[0];
            y = m_points[1];
            z = m_points[2];
        }

    protected:
        Array<OneD, TData> m_points[3]; // x, y, z coordinates
        Array<OneD, TData> m_weights;
        PointsKey m_key;

        PointsBase(const PointsKey& pKey) {m_key = pKey;}
        PointsBase() {}

        template<typename TPts1, typename TPts2, typename... TPtsOther>
        unsigned int GetNumberOfPoints(int i = 0)
        {
            return traits::points_traits<TPts1>::get_total_points(m_key.m_numpoints[i])
                    * GetNumberOfPoints<TPts2, TPtsOther...>(i + 1);
        }

        template<typename TPts>
        unsigned int GetNumberOfPoints(int i = 0)
        {
            return traits::points_traits<TPts>::get_total_points(m_key.m_numpoints[i]);
        }

        template<typename... TPts>
        void AllocateArrays() {
            const unsigned int npts = GetNumberOfPoints<TPts...>();
            std::cout << "Allocating storage of size: " << npts << std::endl;
            for (unsigned int i = 0; i < traits::points_traits<TPts...>::dimension; ++i)
            {
                m_points[i] = Array<OneD, TData>(npts);
            }
            m_weights = Array<OneD, TData>(npts);
        }

};

template<typename TData>
using PointsSharedPtr = boost::shared_ptr<PointsBase<TData>>;

template<typename TData>
using PointsFactory = LibUtilities::NekFactory<
        std::string, PointsBase<TData>, const PointsKey&>;

template<typename TData>
LIB_UTILITIES_EXPORT PointsFactory<TData>& GetPointsFactory()
{
    typedef Loki::SingletonHolder<PointsFactory<TData>,
                                  Loki::CreateUsingNew,
                                  Loki::NoDestroy> Type;
    // Putting Loki::ClassLevelLockable causes an assertion!
    return Type::Instance();
}

//typedef LibUtilities::NekFactory<std::string, PointsBase<double>, const PointsKey&> PointsFactory;
//LIB_UTILITIES_EXPORT PointsFactory& GetPointsFactory()
//{
//    typedef Loki::SingletonHolder<PointsFactory,
//                                  Loki::CreateUsingNew,
//                                  Loki::NoDestroy> Type;
//    return Type::Instance();
//}




/**
 * @class Points
 * @brief Primary template for composite Points classes.
 */
template<typename TData, typename TShape, typename... TPts>
class Points : public PointsBase<TData>
{
        typedef Points<TData, TShape, TPts...> ThisType;
        typedef PointsBase<TData> BaseType;

        static_assert(sizeof...(TPts) > 0, "No point type given.");
        static_assert(sizeof...(TPts) < 4, "Too many point types given.");
        static_assert(traits::points_traits<TPts...>::dimension == traits::shape_traits<TShape>::dimension,
                "Points dimension and shape dimension do not agree.");

    public:
        typedef traits::points_traits<TPts...> traits_info;

        Points(const PointsKey& pKey) : BaseType(pKey) {
            BaseType::template AllocateArrays<TPts...>();
//            InitialiseSubType<sizeof...(TPts)-1>(pKey, x);
        }
//
//
//        inline const Array<OneD, const TData> GetZ(const int& i) {
//            switch (i)
//            {
//                case 0: return std::get<0>(x).GetZ();
//                case 1: return std::get<1>(x).GetZ();
//                default: throw;
//            }
//        }
//
//
//    private:
//        std::tuple<Points<TData, TShape, TPts>...> x;
//
//        template<size_t i>
//        typename std::enable_if<i!=0, void>::type
//        InitialiseSubType(const PointsKey& p, std::tuple<Points<TData, TShape, TPts>...>& pTuple) {
//            PointsKey tmpKey;
//            tmpKey.m_numpoints[0] = p.m_numpoints[i];
//            tmpKey.m_params = p.m_params;
//            cout << "Init subtype " << i << " with numpoints=" << tmpKey.m_numpoints[0] << endl;
//            std::get<i>(pTuple).Populate(p);
//            InitialiseSubType<i-1>(p, pTuple);
//        }
//
//        template<size_t i>
//        typename std::enable_if<i==0, void>::type
//        InitialiseSubType(const PointsKey& p, std::tuple<Points<TData, TShape, TPts>...>& pTuple) {
//            PointsKey tmpKey;
//            tmpKey.m_numpoints[0] = p.m_numpoints[i];
//            tmpKey.m_params = p.m_params;
//            cout << "Init subtype " << i << " with numpoints=" << tmpKey.m_numpoints[0] << endl;
//            std::get<i>(pTuple).Populate(tmpKey);
//        }
};


/**
 * Bi-Points specialisation
 */
template<typename TData, typename TShape, typename TPts1, typename TPts2>
class Points<TData, TShape, TPts1, TPts2> : public PointsBase<TData>
{
        typedef Points<TData, TShape, TPts1, TPts2> ThisType;
        typedef PointsBase<TData> BaseType;

        static_assert(traits::points_traits<TPts1, TPts2>::dimension == traits::shape_traits<TShape>::dimension,
                "Points dimension and shape dimension do not agree.");

    public:
        typedef traits::points_traits<TPts1, TPts2> get_traits;

        Points(const PointsKey& pKey) : BaseType(pKey) {
            //PointsBase<TData>::template AllocateArrays<TPts1, TPts2>();
            x1.Populate(pKey);
            PointsKey tmpKey;
            tmpKey.m_numpoints[0] = pKey.m_numpoints[1];
            tmpKey.m_params = pKey.m_params;
            x2.Populate(tmpKey);
        }

        virtual const Array<OneD, const TData> GetZ(const int& i) {
            switch (i)
            {
                case 0: return x1.GetZ();
                case 1: return x2.GetZ();
                default: throw;
            }
        }

    private:
        Points<TData, typename traits::points_traits<TPts1>::native_shape, TPts1> x1;
        Points<TData, typename traits::points_traits<TPts2>::native_shape, TPts2> x2;
};


/**
 * Tri-Points specialisation
 */
template<typename TData, typename TShape, typename TPts1, typename TPts2, typename TPts3>
class Points<TData, TShape, TPts1, TPts2, TPts3> : public PointsBase<TData>
{
        typedef Points<TData, TShape, TPts1, TPts2, TPts3> ThisType;
        typedef PointsBase<TData> BaseType;

        static_assert(traits::points_traits<TPts1, TPts2, TPts3>::dimension == traits::shape_traits<TShape>::dimension,
                "Points dimension and shape dimension do not agree.");

    public:
        typedef traits::points_traits<TPts1, TPts2, TPts3> get_traits;

        Points(const PointsKey& pKey) : BaseType(pKey) {
            //PointsBase<TData>::template AllocateArrays<TPts1, TPts2, TPts3>();
            PointsKey tmpKey;
            x1.Populate(pKey);
            tmpKey.m_numpoints[0] = pKey.m_numpoints[1];
            tmpKey.m_params = pKey.m_params;
            x2.Populate(tmpKey);
            tmpKey.m_numpoints[0] = pKey.m_numpoints[2];
            tmpKey.m_params = pKey.m_params;
            x3.Populate(tmpKey);
        }

        virtual const Array<OneD, const TData> GetZ(const int& i) {
            switch (i)
            {
                case 0: return x1.GetZ();
                case 1: return x2.GetZ();
                case 2: return x3.GetZ();
                default: throw;
            }
        }

    private:
        Points<TData, typename traits::points_traits<TPts1>::native_shape, TPts1> x1;
        Points<TData, typename traits::points_traits<TPts2>::native_shape, TPts2> x2;
        Points<TData, typename traits::points_traits<TPts3>::native_shape, TPts3> x3;
};


namespace detail {


}

/**
 * Specialisation for GaussGaussLegendre
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
        void Populate(const PointsKey& p) {
            BaseType::m_key = p;
            const int n = p.m_numpoints[0];
            BaseType::template AllocateArrays<PointsType>();
            Polylib::zwgj(BaseType::m_points[0].data(),
                          BaseType::m_weights.data(),
                          n,0.0,0.0);
        }

};


/**
 * Specialisation for GaussGaussLegendre
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

};
}
}
}
