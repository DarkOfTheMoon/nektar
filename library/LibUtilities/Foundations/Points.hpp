#ifndef LIBUTILITIES_FOUNDATIONS_POINTS
#define LIBUTILITIES_FOUNDATIONS_POINTS

#include <iostream>
#include <type_traits>

#include <loki/Singleton.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/LibUtilitiesDeclspec.h>
#include <LibUtilities/Polylib/Polylib.h>
#include <LibUtilities/LinearAlgebra/NekMatrixFwd.hpp>
#include <LibUtilities/Foundations/ShapeTypes.hpp>
#include <LibUtilities/Foundations/PointsTraits.hpp>

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

template<typename TData>
class PointsBase;

template<typename TData>
using PointsSharedPtr = boost::shared_ptr<PointsBase<TData>>;

template<typename TData>
using PointsFactory = LibUtilities::NekFactory<
        std::string, PointsBase<TData>, const PointsKey&>;

// For some reason get undefined symbols when using a templated form...
//template<typename TData>
//LIB_UTILITIES_EXPORT PointsFactory<TData>& GetPointsFactory();

LIB_UTILITIES_EXPORT PointsFactory<NekDouble>& GetPointsFactory();

/**
 * @class PointsBase
 * @brief Points base class defining interface and data members.
 *
 * A points distribution may be multi-dimensional and may be a composite of two
 * or three constituent point distributions, combined in a tensor-product
 * fashion, but conditional that the overall dimension must be less than 3.
 */
template<typename TData>
class PointsBase
{
        typedef boost::shared_ptr<NekMatrix<TData> > MatrixType;

    public:
        /**
         * @brief Destructor
         */
        virtual ~PointsBase() {}


        /**
         * @brief Get the number of points in this distribution
         */
        inline unsigned int GetNumPoints() const
        {
            return m_key.m_numpoints[0];
        }

        /**
         * @brief Get the coordinates of one of the component point distributions
         */
        inline const Array<OneD, const TData> GetZ(const int& i = 0)
        {
            return v_GetZ(i);
        }

        /**
         * @brief Get the quadrature weights
         */
        inline const Array<OneD, const TData>& GetW()
        {
            return m_weights;
        }

        /**
         *
         */
        inline void GetZW(Array<OneD, const TData> &z,
            Array<OneD, const TData> &w) const
        {
            z = m_points[0];
            w = m_weights;
        }

        /**
         *
         */
        inline void GetPoints(Array<OneD, const TData> &x) const
        {
            x = m_points[0];
        }

        /**
         *
         */
        inline void GetPoints(Array<OneD, const TData> &x,
                              Array<OneD, const TData> &y) const
        {
            x = m_points[0];
            y = m_points[1];
        }

        /**
         *
         */
        inline void GetPoints(Array<OneD, const TData> &x,
                              Array<OneD, const TData> &y,
                              Array<OneD, const TData> &z) const
        {
            x = m_points[0];
            y = m_points[1];
            z = m_points[2];
        }

        /**
         *
         */
        inline MatrixType GetInterpMatrix(
                const Array<OneD, const TData>& pts, const int npts = 1 )
        {
            NekDouble* t = v_GetInterpMatrix(npts, pts).data();
            MatrixType returnval(MemoryManager<NekMatrix<NekDouble> >::AllocateSharedPtr(npts,GetNumPoints(),t));

            return returnval;
        }

    protected:
        Array<OneD, TData> m_points[3]; /// x, y, z coordinates of points
        Array<OneD, TData> m_weights;   /// Integration weights
        PointsKey m_key;                /// Defines number of points

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

    private:
        virtual const Array<OneD, const TData> v_GetZ(const int& i)
        {
            ASSERTL1(i == 0, "Only 0 is valid for fundamental types.");
            return m_points[0];
        }

        virtual Array<OneD, TData> v_GetInterpMatrixData(const int& npts, const Array<OneD, const TData>& pts) = 0;
};



/**
 * @class Points
 * @brief Primary template for composite Points classes.
 */
template<typename TData, typename TShape, typename TPtsTuple>
class Points {};

template<typename TData, typename TShape, typename... TPts>
class Points<TData, TShape, std::tuple<TPts...>> : public PointsBase<TData>
{
        typedef Points<TData, TShape, std::tuple<TPts...>> ThisType;
        typedef PointsBase<TData> BaseType;
        typedef std::tuple<Points<TData, typename traits::points_traits<TPts>::native_shape, TPts>...> TupleType;

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

    private:
        TupleType x;

        template<size_t i>
        typename std::enable_if<i!=0, void>::type
        InitialiseSubType(const PointsKey& p, TupleType& pTuple) {
            PointsKey tmpKey;
            tmpKey.m_numpoints[0] = p.m_numpoints[i];
            tmpKey.m_params = p.m_params;
            std::cout << "Init subtype " << i << " with numpoints=" << tmpKey.m_numpoints[0] << std::endl;
            std::get<i>(pTuple).Populate(p);
            InitialiseSubType<i-1>(p, pTuple);
        }

        template<size_t i>
        typename std::enable_if<i==0, void>::type
        InitialiseSubType(const PointsKey& p, TupleType& pTuple) {
            PointsKey tmpKey;
            tmpKey.m_numpoints[0] = p.m_numpoints[i];
            tmpKey.m_params = p.m_params;
            std::cout << "Init subtype " << i << " with numpoints=" << tmpKey.m_numpoints[0] << std::endl;
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
            throw std::logic_error("Out of range.");
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

    private:
        virtual Array<OneD, TData> v_GetInterpMatrixData(const int& npts, const Array<OneD, const TData>& pts)
        {
            Array<OneD, TData> interp(npts * BaseType::m_key.m_numpoints[0]);
            Polylib::Imgrjm(interp.data(), BaseType::m_points[0].data(),
                            pts.data(), BaseType::m_key.m_numpoints[0], npts, 0.0, 0.0);
            return interp;
        }

};

/**
 * Specialisation for Fekete
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
