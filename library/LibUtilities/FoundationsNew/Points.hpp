#include <iostream>
#include <type_traits>

#include <loki/Singleton.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/LibUtilitiesDeclspec.h>
#include <LibUtilities/Polylib/Polylib.h>
#include "PointsTypes.hpp"

namespace Nektar
{
namespace LibUtilities
{

//// Dummy Array class for the purpose of isolated development
//struct OneD {};
//template<typename dim, typename datatype>
//class Array {};

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

        inline const Array<OneD, const TData>& GetZ()
        {
            return m_points;
        }

        inline const Array<OneD, const TData>& GetW()
        {
            return m_weights;
        }

    protected:
        Array<OneD, TData> m_points[3];
        Array<OneD, TData> m_weights;
        PointsKey m_key;

        PointsBase(const PointsKey& pKey) : m_points(), m_weights(), m_key(pKey) {}

        template<typename TPts>
        void AllocateArrays() {
            const unsigned int npts = traits::points_traits<TPts>::get_total_points(m_key.m_numpoints[0]);
            for (unsigned int i = 0; i < traits::points_traits<TPts>::dimension; ++i)
            {
                m_points[i] = Array<OneD, TData>(npts);
            }
        }

};

template<typename TData>
using PointsFactory = LibUtilities::NekFactory<
        std::string, PointsBase<TData>, const PointsKey&>;

template<typename TData>
LIB_UTILITIES_EXPORT PointsFactory<TData>& GetPointsFactory()
{
    typedef Loki::SingletonHolder<PointsFactory<TData>,
                                  Loki::CreateUsingNew,
                                  Loki::NoDestroy,
                                  Loki::ClassLevelLockable> Type;
    return Type::Instance();
}


/**
 * @class Points
 * @brief Primary template for Points classes.
 */
template<typename TShape, typename TPts, typename TData>
class Points : public PointsBase<TData>
{
    //static_assert(false, "No implementation provided for this points type.");

    public:
        Points(const PointsKey& pKey) : PointsBase<TData>(pKey) {}
};


/**
 *
 */
template<typename TShape, typename TData>
class Points<TShape, GaussGaussLegendre, TData> : public PointsBase<TData>
{
    static_assert(traits::distribution_traits<TShape, GaussGaussLegendre>::is_valid,
            "Not a valid combination of shape and points type.");
    public:
        Points(const PointsKey& pKey) : PointsBase<TData>(pKey)
        {
            PointsBase<TData>::template AllocateArrays<GaussGaussLegendre>();
            std::cout << PointsBase<TData>::m_points[0].num_elements();
//            Polylib::zwgj(PointsBase<TData>::m_points[0].data(),
//                          PointsBase<TData>::m_weights.data(),
//                          pKey.m_numpoints,0.0,0.0);
        }
};

}
}
}
