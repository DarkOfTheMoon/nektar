///////////////////////////////////////////////////////////////////////////////
//
// File: Points.hpp
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
// Description: Base class for Points distributions
//
///////////////////////////////////////////////////////////////////////////////

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

/// Data type for parameter key name in PointsKey
typedef std::string PointsParamKey;
/// Data type for parameter value in PointsKey
typedef NekDouble PointsParamValue;
/// List of parameters key/value pairs
typedef std::map<PointsParamKey, PointsParamValue> PointsParamList;

/**
 * @brief Information about a Points object
 */
class PointsKey
{
    public:
        /// Number of points for each of the constituent point distributions
        unsigned int m_numpoints[3];
        /// List of zero or more parameters
        PointsParamList m_params;
};


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
         * @brief Get the points key
         */
        inline const PointsKey& GetKey() const
        {
            return m_key;
        }

        /**
         * @brief Get the coordinates of one of the component point distributions
         */
        inline const Array<OneD, const TData>& GetZ(const int& i = 0) const
        {
            return v_GetZ(i);
        }

        /**
         * @brief Get the quadrature weights
         */
        inline const Array<OneD, const TData>& GetW() const
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
         * @brief Get both the coordinates and quadrature weights.
         * @param z  Output array in which to put the coordinates.
         * @param w  Output array in which to put the quadrature weights.
         */
        inline void GetZW(Array<OneD, const TData> &z,
                          Array<OneD, const TData> &w) const
        {
            z = m_points[0];
            w = m_weights;
        }

        /**
         * @brief Get the coordinates of 1D point distributions
         */
        inline void GetPoints(Array<OneD, const TData> &x) const
        {
            x = m_points[0];
        }

        /**
         * @brief Get the coordinates of 2D point distributions
         */
        inline void GetPoints(Array<OneD, const TData> &x,
                              Array<OneD, const TData> &y) const
        {
            x = m_points[0];
            y = m_points[1];
        }

        /**
         * @brief Get the coordinates of 3D point distributions
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
         * @brief Get the interpolation matrix from this points distribution
         *        to the points distribution provided by \pts.
         */
        inline MatrixType GetInterpMatrix(
                const Array<OneD, const TData>& pts, const int npts = 1 ) const
        {
            NekDouble* t = v_GetInterpMatrixData(npts, pts).data();
            MatrixType returnval(MemoryManager<NekMatrix<NekDouble> >::AllocateSharedPtr(npts,GetNumPoints(),t));

            return returnval;
        }

    protected:
        Array<OneD, TData> m_points[3]; ///< x, y, z coordinates of points
        Array<OneD, TData> m_weights;   ///< Integration weights
        MatrixType         m_ddata[3];  ///< derivative

        PointsKey m_key;                ///< Defines number of points

        PointsBase(const PointsKey& pKey) {m_key = pKey;}
        PointsBase() {}

        /**
         * @brief Allocate necessary storage for a composite point distribution
         */
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
        /**
         * @brief Calculate the total number of points in a composite points
         *        distribution using induction.
         */
        template<typename TPts1, typename TPts2, typename... TPtsOther>
        unsigned int GetNumberOfPoints(int i = 0)
        {
            return traits::points_traits<TPts1>::get_total_points(m_key.m_numpoints[i])
                    * GetNumberOfPoints<TPts2, TPtsOther...>(i+1);
        }

        /**
         * @brief Calculate the total number of points in a composite points
         *        distribution using induction. Terminating case.
         */
        template<typename TPts>
        unsigned int GetNumberOfPoints(int i = 0)
        {
            return traits::points_traits<TPts>::get_total_points(m_key.m_numpoints[i]);
        }

        /**
         * @copydoc PointsBase::GetZ(const int&)
         */
        virtual const Array<OneD, const TData>& v_GetZ(const int& i) const
        {
            if (i != 0)
            {
                throw std::logic_error("Constituent points index out of range.");
            }
            return m_points[0];
        }

        virtual Array<OneD, TData> v_GetInterpMatrixData(
                const int& npts,
                const Array<OneD, const TData>& pts) const = 0;

};

/// A shared pointer to a PointsBase object
template<typename TData>
using PointsSharedPtr = boost::shared_ptr<PointsBase<TData>>;

/// A factory for generating Points objects
template<typename TData>
using PointsFactory = LibUtilities::NekFactory<
        std::string, PointsBase<TData>, const PointsKey&>;

LIB_UTILITIES_EXPORT PointsFactory<NekDouble>& GetPointsFactory();


}
}
}

#endif
