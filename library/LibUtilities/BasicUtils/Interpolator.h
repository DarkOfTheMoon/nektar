///////////////////////////////////////////////////////////////////////////////
//
// File Interpolator.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2015 Kilian Lackhove
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
// Description: Interpolator
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_BASIC_UTILS_INTERPOLATOR_H
#define NEKTAR_LIB_UTILITIES_BASIC_UTILS_INTERPOLATOR_H

#include <vector>

#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/VmathArray.hpp>
#include <LibUtilities/BasicUtils/PtsField.h>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

using namespace std;

namespace Nektar
{
namespace LibUtilities
{

enum InterpMethod{
    eNoMethod,
    eNearestNeighbour,
    eQuadratic,
    eShepard,
    eGauss,
};

class Interpolator
{
    public:

        LIB_UTILITIES_EXPORT Interpolator(
            InterpMethod method = eNoMethod,
            short int coordId = -1,
            NekDouble filtWidth = 0.0):
            m_method(method),
            m_filtWidth(filtWidth),
            m_coordId(coordId)
        {
        };

        LIB_UTILITIES_EXPORT void CalcWeights(
            const PtsFieldSharedPtr inField,
            PtsFieldSharedPtr &outField);

        LIB_UTILITIES_EXPORT void Interpolate(
            const PtsFieldSharedPtr inField,
            PtsFieldSharedPtr &outField);

        LIB_UTILITIES_EXPORT int GetDim() const;

        LIB_UTILITIES_EXPORT NekDouble GetFiltWidth() const;

        LIB_UTILITIES_EXPORT int GetCoordId() const;

        LIB_UTILITIES_EXPORT InterpMethod GetInterpMethod() const;

        LIB_UTILITIES_EXPORT PtsFieldSharedPtr GetInField() const;

        LIB_UTILITIES_EXPORT PtsFieldSharedPtr GetOutField() const;


        template<typename FuncPointerT, typename ObjectPointerT>
        void SetProgressCallback(FuncPointerT func,
                ObjectPointerT obj)
        {
            m_progressCallback = boost::bind(func, obj, _1, _2);
        }

    private:

        class PtsPoint
        {
            public:

                int                                         idx;
                Array<OneD, NekDouble>                      coords;
                NekDouble                                   dist;

                LIB_UTILITIES_EXPORT PtsPoint():
                    idx(-1),
                    coords(Array<OneD, NekDouble>(3)),
                    dist(1E30)
                {
                };

                LIB_UTILITIES_EXPORT PtsPoint(int idx, Array<OneD, NekDouble> coords,
                                            NekDouble dist):
                    idx(idx),
                    coords(coords),
                    dist(dist)
                {
                };

                LIB_UTILITIES_EXPORT bool operator < (const PtsPoint &comp) const
                {
                    return (dist < comp.dist);
                };
        };

        static const int                            m_dim = 3;
        typedef bg::model::point<NekDouble, m_dim, bg::cs::cartesian> BPoint;
        typedef std::pair<BPoint, unsigned int>                       PtsPointPair;
        typedef bgi::rtree< PtsPointPair, bgi::rstar<16> >            PtsRtree;

        PtsFieldSharedPtr                           m_inField;
        PtsFieldSharedPtr                           m_outField;

        /// Interpolation Method
        InterpMethod                                m_method;
        /// A tree structure to speed up the neighbour search.
        /// Note that we fill it with an iterator, so instead of rstar, the
        /// packing algorithm is used.
        boost::shared_ptr<PtsRtree>                 m_rtree;
        /// Interpolation weights for each neighbour.
        /// Structure: m_weights[physPtIdx][neighbourIdx]
        Array<OneD, Array<OneD, float> >            m_weights;
        /// Indices of the relevant neighbours for each physical point.
        /// Structure: m_neighInds[ptIdx][neighbourIdx]
        Array<OneD, Array<OneD, unsigned int> >     m_neighInds;
        NekDouble                                   m_filtWidth;
        short int                                   m_coordId;

        boost::function<void (const int position, const int goal)> m_progressCallback;

        LIB_UTILITIES_EXPORT void CalcW_Gauss(
            const PtsPoint &searchPt,
            const NekDouble sigma);

        LIB_UTILITIES_EXPORT void CalcW_Linear(const PtsPoint &searchPt, int coordId);

        LIB_UTILITIES_EXPORT void CalcW_NNeighbour(const PtsPoint &searchPt);

        LIB_UTILITIES_EXPORT void CalcW_Shepard(const PtsPoint &searchPt);

        LIB_UTILITIES_EXPORT void CalcW_Quadratic(const PtsPoint &searchPt, int coordId);

        LIB_UTILITIES_EXPORT void FindNeighbours(const PtsPoint &searchPt,
                vector<PtsPoint > &neighbourPts,
                const NekDouble dist);

        LIB_UTILITIES_EXPORT void FindNNeighbours(const PtsPoint &searchPt,
                vector<PtsPoint > &neighbourPts,
                const unsigned int numPts = 1);
};

typedef boost::shared_ptr<Interpolator> InterpolatorSharedPtr;
static InterpolatorSharedPtr NullInterpolator;


}
}



#endif

