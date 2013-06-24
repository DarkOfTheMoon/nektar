////////////////////////////////////////////////////////////////////////////////
//
//  File:  Spline.h
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  License for the specific language governing rights and limitations under
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description:  Builds and evaluates cubic splines through given set of coordinates in 2D
//
//
////////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_LIBUTILITIES_SPLINE_H
#define NEKTAR_LIBUTILITIES_SPLINE_H

#include <string>
#include <LibUtilities/LinearAlgebra/Lapack.hpp>
#include <LibUtilities/BasicUtils/VmathArray.hpp>
#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
#include <LibUtilities/LibUtilitiesDeclspec.h>
#include <MultiRegions/ContField2D.h>
#include <iostream>
#include <fstream>

using namespace std;

namespace Nektar
{
  namespace LibUtilities
  {
	  enum SplineBoundaryType
	  {
		eSplineBoundary,
		eNatural,
		eClamped,
		eClampedQuadratic,
		eClampedQuadraticLagrange,
		eNotaKnot,
		eSplineBoundaryTypeSize
	  };

  	  class CubicSpline
  	  {
  	  	  public:

		  CubicSpline()
		  {

		  };

  		  CubicSpline(Array<OneD,NekDouble> &x, Array<OneD,NekDouble> &y, SplineBoundaryType bcleft, SplineBoundaryType bcright,
  				NekDouble bcvalueleft, NekDouble bcvalueright);

  		  CubicSpline(const CubicSpline &In);

  		  /// Build cubic spline through coordinates
  		  void UpdateCubicSpline(Array<OneD,NekDouble> &x, Array<OneD,NekDouble> &y);

  		  void SolveTriangularMatrixForSlope();

  		  void SolveTriMatrix (int n, Array<OneD,NekDouble> &a, Array<OneD,NekDouble> &b, Array<OneD,NekDouble> &c, Array<OneD,NekDouble> &v, Array<OneD,NekDouble> &x);

  		  void ComputeCoefficients();

  		  void ComputeSplineNormals();

  		  void ComputeSplineCurvature();

  		  NekDouble GetSplineValue(NekDouble x);

  		  NekDouble GetSplineDerivativeValue(NekDouble x);

  		  NekDouble GetSplineDiff(int id);

  		  void WriteMatlabFiles(string fname, int cnt);

  		  void WriteMatlabFilesInterpolation(string fname, int cnt, int nintpoints);

  		  boost::shared_ptr<CubicSpline> SetUpSpline(MultiRegions::ExpListSharedPtr &pBoundaryExpansion,
  				  SplineBoundaryType bcleft, SplineBoundaryType bcright,
    				NekDouble bcvalueleft, NekDouble bcvalueright);

  		  void UpdateSpline(MultiRegions::ExpListSharedPtr &pBoundaryExpansion);

  		  int GetSplineCoordOffset(NekDouble xcoord);

  		  void GetSplineLocalNormals(NekDouble xcoord, int n_quad, Array<OneD, Array<OneD, NekDouble> > &locnormals);

  		  int GetSplineNormal(NekDouble xcoord, NekDouble &nx, NekDouble &ny);
  		  int GetSplineCurvature(NekDouble xcoord, NekDouble &curv);

  		  void UpdateSplineCornerPoints(MultiRegions::ExpListSharedPtr &pBoundaryExpansion);
  		  boost::shared_ptr<CubicSpline>  SetUpSplineCornerPoints(MultiRegions::ExpListSharedPtr &pBoundaryExpansion);

  		  void ComputeGLLPointsAlongEdge(Array<OneD,NekDouble> &xedge, Array<OneD,NekDouble> &yedge, int &npoints);

  		  NekDouble ComputeArcLengthIntegralWithSimpsonsRule(NekDouble &a, NekDouble &b);
  		  NekDouble ComputeArcLengthIntegralWithGaussianQuadrature(NekDouble &a, NekDouble &b);
  		  NekDouble ComputeSplineTotalArcLength();
  		  NekDouble ComputeSplineArcLength(NekDouble &a, NekDouble &b);

  		  inline int GetTotPoints()
		  {
			return m_totpoints;
		  }

  		  inline Array<OneD, Array<OneD,NekDouble> > GetCoords()
		  {
			return m_coords;
		  }

  		  inline Array<OneD, Array<OneD,NekDouble> > GetNormals()
		  {
			return m_normals;
		  }

  		  inline Array<OneD,NekDouble> GetCurv()
		  {
			return m_curv;
		  }

  		  inline NekDouble GetArcLengthFunction(NekDouble &x)
		  {
			return sqrt(1+pow(GetSplineDerivativeValue(x),2.0));
		  }


          void SetCoords(Array<OneD,NekDouble> &x, Array<OneD,NekDouble> &y);

  	  	  protected:

          	  /// Slope of the spline, i.e. second derivative of spline s''
              Array<OneD,NekDouble> m_M;

          	  //total number of coordinate points
          	  int m_totpoints;

          	  // Boundary conditions left and right
          	  SplineBoundaryType m_bcleft, m_bcright;
          	  NekDouble m_bcleftvalue, m_bcrightvalue;

          	  /// coefficients of spline
          	  Array<OneD,NekDouble> m_a;
          	  Array<OneD,NekDouble> m_b;
          	  Array<OneD,NekDouble> m_c;
          	  Array<OneD,NekDouble> m_d;

  		  	  /// coordinates that are used to build the spline
  		  	  Array<OneD, Array<OneD,NekDouble> > m_coords;

  		  	  /// Normals along spline
  		  	  Array<OneD, Array<OneD,NekDouble> > m_normals;

  		  	  /// Curvature along spline
  		  	  Array<OneD,NekDouble>  m_curv;

  	  };

  	typedef boost::shared_ptr<CubicSpline> CubicSplineSharedPtr;

  }
}

#endif /* SPLINE_H_ */
