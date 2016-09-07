#ifndef LIBUTILITIES_FOUNDATIONS_POINTS
#define LIBUTILITIES_FOUNDATIONS_POINTS

#include <iostream>
#include <type_traits>

#include <loki/Singleton.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/LibUtilitiesDeclspec.h>
#include <LibUtilities/Polylib/Polylib.h>
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>
#include <LibUtilities/Foundations/ShapeTypes.hpp>
#include <LibUtilities/Foundations/Points/PointsTraits.hpp>

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
         * @brief Get the derivative matrix
         */
        inline const MatrixType& GetD(int dir = 0) const
        {
            return m_ddata[dir];
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
        MatrixType         m_ddata[3];  /// derivative

        PointsKey m_key;                /// Defines number of points

        PointsBase(const PointsKey& pKey) {m_key = pKey;}
        PointsBase() {}

        template<typename TPts1, typename TPts2, typename... TPtsOther>
        unsigned int GetNumberOfPoints(int i = 0)
        {
            return traits::points_traits<TPts1>::get_total_points(m_key.m_numpoints[i])
                    * GetNumberOfPoints<TPts2, TPtsOther...>(i+1);
        }

        template<typename TPts>
        unsigned int GetNumberOfPoints(int i = 0)
        {
            return traits::points_traits<TPts>::get_total_points(m_key.m_numpoints[i]);
        }

        template<typename... TPts>
        void AllocateArrays() {
            const unsigned int npts = GetNumberOfPoints<TPts...>();
            unsigned int i = 0;
            for (i = 0; i < traits::points_traits<TPts...>::dimension; ++i)
            {
                m_points[i] = Array<OneD, TData>(npts);
            }
            m_weights = Array<OneD, TData>(npts);
            for (i = 0; i < traits::points_traits<TPts...>::dimension; ++i)
            {
                m_ddata[i] = MemoryManager<NekMatrix<TData> >::AllocateSharedPtr(npts,npts);
            }
        }

    private:
        virtual const Array<OneD, const TData> v_GetZ(const int& i)
        {
            if (i != 0)
            {
                throw std::logic_error("Constituent points index out of range.");
            }
            return m_points[0];
        }

        virtual Array<OneD, TData> v_GetInterpMatrixData(const int& npts, const Array<OneD, const TData>& pts) = 0;
};







}
}
}

#endif
