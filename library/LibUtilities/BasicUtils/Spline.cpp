////////////////////////////////////////////////////////////////////////////////
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
#ifndef NEKTAR_LIBUTILITIES_SPLINE_CPP
#define NEKTAR_LIBUTILITIES_SPLINE_CPP

#include <LibUtilities/BasicUtils/Spline.h>

//using namespace std;

namespace Nektar
{
    namespace LibUtilities
    {
    	CubicSpline::CubicSpline()
    	{}

    	CubicSpline::CubicSpline(Array<OneD,NekDouble> &x, Array<OneD,NekDouble> &y, SplineBoundaryType bcleft, SplineBoundaryType bcright, NekDouble bcvalueleft, NekDouble bcvalueright)
		{
    		int spacedim = 2;
    		int i,j;
    		int cnt=0;
    		Array<OneD,NekDouble> RHS;
    		NekDouble h;
    		Array<OneD,NekDouble> xtmp(x.num_elements());
    		Array<OneD,NekDouble> ytmp(y.num_elements());

    		Vmath::Vcopy(x.num_elements(),x,1,xtmp,1);
    		Vmath::Vcopy(y.num_elements(),y,1,ytmp,1);

    		//bubble sort coordinates after increasing x
			for(i = 0; i < xtmp.num_elements(); i++)
			{
				for(j = 1; j < xtmp.num_elements(); j++)
				{
					if(xtmp[j] < xtmp[j-1])
					{
						std::swap(xtmp[j], xtmp[j-1]);
						std::swap(ytmp[j], ytmp[j-1]);
					}
				}
			}

    		// Allocate Memory
    		// Find out number of distinct points in terms of x because for h=0 slope computation fails
    		for(i=0;i<x.num_elements();i++)
    		{
    			if(i==x.num_elements()-1)
				{
					h = abs(xtmp[i] - xtmp[i-1]);
				}
				else
				{
					h = abs(xtmp[i+1] - xtmp[i]);
				}
    			if(h>=1e-10)
    			{
    				cnt++;
    			}
    		}

    		m_totpoints = cnt;

    	/*	std::cout << "Original number of x values: " << x.num_elements() << std::endl;

    		std::cout << "Total number of points to be interpolated: " << m_totpoints << std::endl; */

    		m_M = Array<OneD,NekDouble> (m_totpoints);

    		m_a = Array<OneD,NekDouble> (m_totpoints);
    		m_b = Array<OneD,NekDouble> (m_totpoints);
    		m_c = Array<OneD,NekDouble> (m_totpoints);
    		m_d = Array<OneD,NekDouble> (m_totpoints);

    		m_curv = Array<OneD,NekDouble> (m_totpoints);

    		m_coords = Array<OneD, Array<OneD, NekDouble> > (spacedim);
    		m_normals = Array<OneD, Array<OneD, NekDouble> > (spacedim);

    		m_coords[0] = Array<OneD, NekDouble>(m_totpoints*spacedim);
    		m_normals[0] = Array<OneD, NekDouble>(m_totpoints*spacedim);
    		for(i = 1; i < spacedim; ++i)
			{
				m_coords[i] = m_coords[i-1] + m_totpoints;
				m_normals[i] = m_normals[i-1] + m_totpoints;
			}

    		// Copy coords x,y into m_coords
    		SetCoords(xtmp,ytmp);

    		/*std::cout << "Spline from (x,y)=(" << m_coords[0][0]  << "," << m_coords[1][0]  << ") to (x,y)=("
    				<< m_coords[0][m_totpoints-1]  << "," << m_coords[1][m_totpoints-1]  << ")" << std::endl; */

    		m_bcleft = bcleft;
    		m_bcright = bcright;

    		m_bcleftvalue = bcvalueleft;
    		m_bcrightvalue = bcvalueright;

    		if(m_bcleft==eClampedQuadratic)
    		{
    			m_bcleftvalue = 2*((m_coords[1][1]-m_coords[1][0])/(m_coords[0][1]-m_coords[0][0]))-(m_coords[1][2]-m_coords[1][1])/(m_coords[0][2]-m_coords[0][1]);
    		}

    		if(m_bcright==eClampedQuadratic)
			{
				m_bcrightvalue = 2*((m_coords[1][m_totpoints-1]-m_coords[1][m_totpoints-2])/(m_coords[0][m_totpoints-1]-m_coords[0][m_totpoints-2]))-(m_coords[1][m_totpoints-2]-m_coords[1][m_totpoints-3])/(m_coords[0][m_totpoints-2]-m_coords[0][m_totpoints-3]);
			}

    		if(m_bcleft==eClampedQuadraticLagrange)
			{
				m_bcleftvalue = (((2*m_coords[0][0])-(m_coords[0][1]+m_coords[0][2]))/((m_coords[0][0]-m_coords[0][1])*(m_coords[0][0]-m_coords[0][2])))*m_coords[1][0]
				              + ((m_coords[0][0]-m_coords[0][2])/((m_coords[0][1]-m_coords[0][0])*(m_coords[0][1]-m_coords[0][2])))*m_coords[1][1]
							  + ((m_coords[0][0]-m_coords[0][1])/((m_coords[0][2]-m_coords[0][0])*(m_coords[0][2]-m_coords[0][1])))*m_coords[1][2];
			}

    		if(m_bcright==eClampedQuadraticLagrange)
			{
				m_bcrightvalue = (((2*m_coords[0][m_totpoints-1])-(m_coords[0][m_totpoints-2]+m_coords[0][m_totpoints-3]))/((m_coords[0][m_totpoints-1]-m_coords[0][m_totpoints-2])*(m_coords[0][m_totpoints-1]-m_coords[0][m_totpoints-3])))*m_coords[1][m_totpoints-1]
							  + ((m_coords[0][m_totpoints-1]-m_coords[0][m_totpoints-3])/((m_coords[0][m_totpoints-2]-m_coords[0][m_totpoints-1])*(m_coords[0][m_totpoints-2]-m_coords[0][m_totpoints-3])))*m_coords[1][m_totpoints-2]
							  + ((m_coords[0][m_totpoints-1]-m_coords[0][m_totpoints-2])/((m_coords[0][m_totpoints-3]-m_coords[0][m_totpoints-1])*(m_coords[0][m_totpoints-3]-m_coords[0][m_totpoints-2])))*m_coords[1][m_totpoints-3];
			}

    		// Fill the rest of the Slope Array m_M
    		SolveTriangularMatrixForSlope();

    		// Compute coefficients a,b,c,d
    		ComputeCoefficients();

    		/*cout << "Slope values = " << endl;
    		for(i=0;i<m_totpoints;i++)
			{
    			 cout << "M[" << i << "]=" << m_M[i] << endl;
    			 cout << "(x,y)=(" << m_coords[0][i] << "," << m_coords[1][i] << ")" << endl;
    			 cout << "a[" << i << "]=" << m_a[i] << endl;
    			 cout << "b[" << i << "]=" << m_b[i] << endl;
    			 cout << "c[" << i << "]=" << m_c[i] << endl;
    			 cout << "d[" << i << "]=" << m_d[i] << endl;
			}*/

    		// Compute normals at given coordinates x,y from computed spline
    		ComputeSplineNormals();

    		// Compute curvature at given coordinates x,y from computed spline
    		ComputeSplineCurvature();

		}

    	CubicSpline::CubicSpline(const CubicSpline &In):
				m_M(In.m_M),
				m_totpoints(In.m_totpoints),
				m_bcleft(In.m_bcleft),
				m_bcright(In.m_bcright),
				m_a(In.m_a),
				m_b(In.m_b),
				m_c(In.m_c),
				m_d(In.m_d),
				m_coords(In. m_coords),
				m_normals(In.m_normals),
				m_curv(In.m_curv)
			{
			}

    	void CubicSpline::SolveTriangularMatrixForSlope()
    	{
    		// Prepare Right Hand Side and Slope Values to call tridiagonal matrix solver
			Array<OneD,NekDouble> h;
			Array<OneD,NekDouble> dl,d,du,b,outarray;
			// Number of matrix entries
			int N = m_totpoints;
			int istart,iend;

    		// Set boundary conditions in Slope Array
			if(m_bcleft == eNatural)
			{
				 m_M[0] = 0.0;
			     N = N-1;
			}
			if(m_bcright == eNatural)
			{
				m_M[m_totpoints-1] = 0.0;
				N = N-1;
			}


			int i;
			int info;

    		//cout << "Left BC value" << m_bcleftvalue;
    		// cout << "Right BC value" << m_bcrightvalue << endl;

			// all have dimension N-2
    		h = Array<OneD,NekDouble> (m_totpoints-1);
    		dl = Array<OneD,NekDouble> (N-1);
    		d = Array<OneD,NekDouble> (N);
    		du = Array<OneD,NekDouble> (N-1);
    		b = Array<OneD,NekDouble> (N);
    		outarray = Array<OneD,NekDouble> (N);

    		for(i=0;i<m_totpoints-1;i++)
    		{
    			h[i] = m_coords[0][i+1]-m_coords[0][i];
    			if(!h[i])
    			{
    				cout << "Warning distance from one point to the other is zero in " << i << endl;
    			}
    			//cout << "h[" << i << "]=" << h[i] << endl;
    		}

    		switch(m_bcleft)
    		{
    			case eNatural: istart = 0;
    						   break;
    			case eClampedQuadratic:
    			case eClampedQuadraticLagrange:
    			case eClamped: istart = 1;
    						   d[0] = 2*h[0];
    						   du[0] = h[0];
    						   dl[0] = h[0];
    						   b[0] = 6.0*(((m_coords[1][1]-m_coords[1][0])/h[0])-m_bcleftvalue);
    						   break;
    			case eNotaKnot: istart =1;
    							dl[0] = h[0]-(h[1]*h[1])/h[0];
							    d[0] = 2*h[0]+(h[1]*h[1])/h[0]+3*h[1];
							    du[0] = h[0];
							    b[0] = 6*(((m_coords[1][2]-m_coords[1][1])/h[1])-((m_coords[1][1]-m_coords[1][0])/h[0]));
    							break;
    			default:
    				ASSERTL0(false,"SplineBCTYPE not found");
    				break;
    		}

    		switch(m_bcright)
			{
				case eNatural: iend = N;
							   break;
				case eClampedQuadratic:
				case eClampedQuadraticLagrange:
				case eClamped: iend = N-1;
							   d[N-1] = 2*h[m_totpoints-2];
							   du[N-2] =  h[m_totpoints-2];
							   dl[N-2] =  h[m_totpoints-2];
							   b[N-1] = 6.0*(m_bcrightvalue-((m_coords[1][m_totpoints-1]-m_coords[1][m_totpoints-2])/h[m_totpoints-2]));
							   break;
    			case eNotaKnot: iend =N-1;
    							dl[N-2] = h[m_totpoints-2]-(h[m_totpoints-3]*h[m_totpoints-3])/h[m_totpoints-2];
							    d[N-1] = 2*h[m_totpoints-2]+(h[m_totpoints-3]*h[m_totpoints-3])/h[m_totpoints-2]+3*h[m_totpoints-3];
							    du[N-2] = h[m_totpoints-2];
							    b[N-1] = 6*(((m_coords[1][m_totpoints-1]-m_coords[1][m_totpoints-2])/h[m_totpoints-2])-((m_coords[1][m_totpoints-2]-m_coords[1][m_totpoints-3])/h[m_totpoints-3]));
							    break;
				default:
					ASSERTL0(false,"SplineBCTYPE not found");
					break;
			}


    		for(i=istart;i<iend;i++)
		    {
    			if(istart==0)
				{
    				d[i] = 2*(h[i]+h[i+1]);
				}
    			else
    			{
    				d[i] = 2*(h[i-1]+h[i]);
    			}
		    }

    		for(i=istart;i<iend-1;i++)
			{
    			if(istart==0)
				{
    				du[i] = h[i+1];
    				dl[i] = h[i+1];
				}
    			else
    			{
    				du[i] = h[i];
    				dl[i] = h[i];
    			}
			}
    		for(i=istart;i<iend;i++)
    		{
    			if(istart==0)
    			{
    				b[i] = 6*(((m_coords[1][i+2]-m_coords[1][i+1])/h[i+1])-((m_coords[1][i+1]-m_coords[1][i])/h[i]));
    			}
    			else
    			{
    				b[i] = 6*(((m_coords[1][i+1]-m_coords[1][i])/h[i])-((m_coords[1][i]-m_coords[1][i-1])/h[i-1]));
    			}
    		}
    	/*	for(i=0;i<N;i++)
    		{
    			cout << "d[" << i << "]=" << d[i] << endl;
    			cout << "b[" << i << "]=" << b[i] << endl;

    		}
    		for(i=0;i<N-1;i++)
			{
				cout << "du[" << i << "]=" << du[i] << endl;

			} */

    		// Include boundary conditions in b
    		//b[0] -= h[0]*m_M[0];
    		//b[N-1] -= h[m_totpoints-2]*m_M[m_totpoints-1];

    		//Solve triangular system
    		//Lapack::dgtsv_(N,1,&dl[0],&d[0],&du[0],&b[0],N,info);

    		SolveTriMatrix (N, dl, d, du, b, outarray);

    		if(info==0)
    		{
    			//cout << "Tridiagonal Matrix successfully inverted" << endl;
    		}
    		else if(info<0)
    		{
    			cout << "Illegal value in "<< info << endl;
    		}
    		else
			{
				cout << "Value in "<< info <<" is zero and has not been computed" << endl;
			}

    		switch(m_bcleft)
    		{
    			case eNatural: for(i=1;i<m_totpoints-1;i++)
							   {
									m_M[i] = b[i-1];
							   }
    						   break;
    			case eNotaKnot:
    			case eClampedQuadratic:
    			case eClampedQuadraticLagrange:
    			case eClamped: for(i=0;i<m_totpoints-1;i++)
							   {
									m_M[i] = b[i];
							   }
    						   break;
    			default:
    				ASSERTL0(false,"SplineBCTYPE not found");
    				break;
    		}

    		switch(m_bcright)
			{
				case eNatural: break;
				case eNotaKnot:
				case eClampedQuadratic:
				case eClampedQuadraticLagrange:
				case eClamped: m_M[m_totpoints-1] = b[m_totpoints-1];
							   break;
				default:
					ASSERTL0(false,"SplineBCTYPE not found");
					break;
			}

    	}

    	void CubicSpline::SolveTriMatrix (int n, Array<OneD,NekDouble> &a, Array<OneD,NekDouble> &b, Array<OneD,NekDouble> &c, Array<OneD,NekDouble> &v, Array<OneD,NekDouble> &x)
    	{
    	        /**
    	         * n - number of equations
    	         * a - sub-diagonal (means it is the diagonal below the main diagonal) -- indexed from 1..n-1
    	         * b - the main diagonal
    	         * c - sup-diagonal (means it is the diagonal above the main diagonal) -- indexed from 0..n-2
    	         * v - right part
    	         * x - the answer
    	         */
    	        for (int i = 1; i < n; i++)
    	        {
    	                double m = a[i]/b[i-1];
    	                b[i] = b[i] - m*c[i-1];
    	                v[i] = v[i] - m*v[i-1];
    	        }

    	        x[n-1] = v[n-1]/b[n-1];

    	        for (int i = n - 2; i >= 0; i--)
    	                x[i]=(v[i]-c[i]*x[i+1])/b[i];
    	}

    	void CubicSpline::ComputeCoefficients()
    	{
    		NekDouble h;
    		int i;
    		for(i=0;i<m_totpoints-1;i++)
    		{
    			h = m_coords[0][i+1] - m_coords[0][i];
    			m_a[i] = (m_M[i+1] - m_M[i]) / (6.0*h);
    			m_b[i] = (m_M[i]) / (2.0);
    			m_c[i] = (m_coords[1][i+1]-m_coords[1][i])/h - ((m_M[i+1] + 2.0*m_M[i]) /6.0) * h;
    			m_d[i] = m_coords[1][i];
    		}
    		i=m_totpoints-1;
    		h = m_coords[0][i] - m_coords[0][i-1];
    		m_a[i] = (m_M[i] - m_M[i-1]) / (6.0*h);
			m_b[i] = (m_M[i]) / (2.0);
			m_c[i] = (m_coords[1][i]-m_coords[1][i-1])/h + ((m_M[i-1] + 2.0*m_M[i]) /6.0) * h;
			m_d[i] = m_coords[1][i];
    	}

    	void CubicSpline::ComputeSplineNormals()
    	{
    		NekDouble norm;
    		NekDouble splinediff;
    		for(int i=0;i<m_totpoints;i++)
    		{
    			splinediff = m_c[i];
    			norm = sqrt(splinediff*splinediff+1);
    			m_normals[0][i]=-splinediff/norm;
    			m_normals[1][i]=1/norm;
    		}
    	}

    	void CubicSpline::ComputeSplineCurvature()
		{
    		NekDouble splinediff,splinediff2;
    		for(int i=0;i<m_totpoints;i++)
			{
    			splinediff = m_c[i];
    			splinediff2 = 2*m_b[i];
    			m_curv[i]=splinediff2/sqrt(pow(1+splinediff*splinediff,3));
			}
		}

    	NekDouble CubicSpline::ComputeArcLengthIntegralWithSimpsonsRule(NekDouble &a, NekDouble &b)
		{
    		NekDouble fa,f2ab,fa2b,fb;
    		NekDouble result = 0;
    		fa = GetArcLengthFunction(a);
    		result = (2*a+b)/3.0;
    		f2ab = GetArcLengthFunction(result);
    		result = (a+2*b)/3.0;
    		fa2b = GetArcLengthFunction(result);
    		fb = GetArcLengthFunction(b);

    		//cout << "(fa,f2ab,fa2b,fb)=(" << fa << "," << f2ab << "," << fa2b << "," << fb << ")" << endl;

    		result = ((b-a)/8.0)*(fa+3.0*f2ab+3.0*fa2b+fb);
    		return result;
		}

    	NekDouble CubicSpline::ComputeArcLengthIntegralWithGaussianQuadrature(NekDouble &a, NekDouble &b)
    	{
    		LibUtilities::PointsType quadPointsType = LibUtilities::eGaussLobattoLegendre;
    		int nQuadPoints = 6;
    		const LibUtilities::PointsKey quadPointsKey(nQuadPoints, quadPointsType);

			// Declare a BasisKey which uniquely defines the (discrete) one-dimensional
			// spectral/hp expansion basis
			// As we will use the basis only to evaluate an integral, the type and order of the basis
			// can be chosen arbitrarily
			const LibUtilities::BasisKey basisKey(LibUtilities::eModified_A,2,quadPointsKey);

			// Using these keys, now define an spectral/hp expansion of corresponding type on a
			// standard region
			StdRegions::StdSegExpSharedPtr segExpansion =
				MemoryManager<StdRegions::StdSegExp>::AllocateSharedPtr(basisKey);

			// Calculate the coordinates xi of all the quadrature points of this (discrete)
			// expansion
			int nTotQuadPoints = segExpansion->GetTotPoints();
			Array<OneD,NekDouble> xi(nTotQuadPoints);
			segExpansion->GetCoords(xi);

			// Analytically evaluate the Jacobian.
			NekDouble jacobian;
			jacobian = 0.5*(b-a);

			// Calculate the values of the function to be integrated
			int i;
			Array<OneD,NekDouble> integrand(nTotQuadPoints);
			for(i = 0; i < nTotQuadPoints; i++)
			{
				// Calculate the local coordinate of quadrature zeros
				double x = a * (1-xi[i])/2 + b * (1+xi[i])/2;

				integrand[i] = GetArcLengthFunction(x);

			}

			// Do the numerical integration by calling the ecorresponding Nektar++ routine of
			// the StdExpansion class
			NekDouble result = jacobian * segExpansion->Integral(integrand);

			return result;
    	}

    	NekDouble CubicSpline::ComputeSplineTotalArcLength()
		{
    		NekDouble result = 0;
    		NekDouble a = m_coords[0][0];
    		NekDouble b = m_coords[0][m_totpoints-1];
    		//cout << "(a,b)=( " << a << ", " << b << ")" << endl;
    		//result = ComputeArcLengthIntegralWithSimpsonsRule(a, b);
    		result = ComputeArcLengthIntegralWithGaussianQuadrature(a, b);
    		return result;
		}

    	NekDouble CubicSpline::ComputeSplineArcLength(NekDouble &a, NekDouble &b)
		{
			NekDouble result = 0;
			//result = ComputeArcLengthIntegralWithSimpsonsRule(a, b);
			result = ComputeArcLengthIntegralWithGaussianQuadrature(a, b);
			return result;
		}

    	// For a given set of points xedge, compute the GLL point locations along an edge and return them in xedge,yedge
    	void CubicSpline::ComputeGLLPointsAlongEdge(Array<OneD,NekDouble> &xedge, Array<OneD,NekDouble> &yedge, int &npoints)
		{
    		int nQuadPoints=npoints;
    		int i;
    		NekDouble g,Dg,xn;
    		LibUtilities::PointsType quadPointsType = LibUtilities::eGaussLobattoLegendre;
			//Decalre variables to hold the quadrature zeros and weights
			Array<OneD, NekDouble> xi;
			const PointsKey quadPointsKey(nQuadPoints, quadPointsType);
			/* obtain quadrature points on standard element*/
			xi = PointsManager()[quadPointsKey]->GetZ();

			/* Compute Total Arc-Length of Edge */
			NekDouble L = ComputeSplineArcLength(xedge[0], xedge[npoints-1]);
			NekDouble diff = 1.0;

			/* Solve root finding problem with Newton method to find point x_i at GLL location */
			/* xedge[0] and xedge[npoints-1] remains unchanged (they are mapped onto -1 and 1) */
			for(i=1;i<npoints-1;i++)
			{
				//cout << "Before: " <<  i << ": " << xedge[i] << endl;
				diff = 1.0;
				/*initial guess should be x-value given */
				while(diff>1e-5)
				{
					xn = xedge[i];
					g = (1.0/L)*ComputeSplineArcLength(xedge[0], xn)-0.5*(xi[i]+1);
					Dg = 1.0/L*GetArcLengthFunction(xn);
					xedge[i] = xn - g/Dg;
					diff = abs(xedge[i]-xn);
					//cout << diff << endl;
				}
				//cout  << "After: " << i << ": " << xedge[i] << endl;
			}

			for(i=0;i<npoints;i++)
			{
				yedge[i] = GetSplineValue(xedge[i]);
			}
		}

    	void CubicSpline::UpdateCubicSpline(Array<OneD,NekDouble> &x, Array<OneD,NekDouble> &y)
    	{
       		//bubble sort coordinates after increasing x
    			for(int i = 0; i < x.num_elements(); i++)
    			{
    				for(int j = 1; j < x.num_elements(); j++)
    				{
    					if(x[j] < x[j-1])
    					{
    						std::swap(x[j], x[j-1]);
    						std::swap(y[j], y[j-1]);
    					}
    				}
    			}

           		// Copy coords x,y into m_coords
    			SetCoords(x,y);

        		if(m_bcleft==eClampedQuadratic)
        		{
        			m_bcleftvalue = 2*((m_coords[1][1]-m_coords[1][0])/(m_coords[0][1]-m_coords[0][0]))-(m_coords[1][2]-m_coords[1][1])/(m_coords[0][2]-m_coords[0][1]);
        		}

        		if(m_bcright==eClampedQuadratic)
    			{
    				m_bcrightvalue = 2*((m_coords[1][m_totpoints-1]-m_coords[1][m_totpoints-2])/(m_coords[0][m_totpoints-1]-m_coords[0][m_totpoints-2]))-(m_coords[1][m_totpoints-2]-m_coords[1][m_totpoints-3])/(m_coords[0][m_totpoints-2]-m_coords[0][m_totpoints-3]);
    			}

        		if(m_bcleft==eClampedQuadraticLagrange)
    			{
    				m_bcleftvalue = (((2*m_coords[0][0])-(m_coords[0][1]+m_coords[0][2]))/((m_coords[0][0]-m_coords[0][1])*(m_coords[0][0]-m_coords[0][2])))*m_coords[1][0]
    				              + ((m_coords[0][0]-m_coords[0][2])/((m_coords[0][1]-m_coords[0][0])*(m_coords[0][1]-m_coords[0][2])))*m_coords[1][1]
    							  + ((m_coords[0][0]-m_coords[0][1])/((m_coords[0][2]-m_coords[0][0])*(m_coords[0][2]-m_coords[0][1])))*m_coords[1][2];
    			}

        		if(m_bcright==eClampedQuadraticLagrange)
    			{
    				m_bcrightvalue = (((2*m_coords[0][m_totpoints-1])-(m_coords[0][m_totpoints-2]+m_coords[0][m_totpoints-3]))/((m_coords[0][m_totpoints-1]-m_coords[0][m_totpoints-2])*(m_coords[0][m_totpoints-1]-m_coords[0][m_totpoints-3])))*m_coords[1][m_totpoints-1]
    							  + ((m_coords[0][m_totpoints-1]-m_coords[0][m_totpoints-3])/((m_coords[0][m_totpoints-2]-m_coords[0][m_totpoints-1])*(m_coords[0][m_totpoints-2]-m_coords[0][m_totpoints-3])))*m_coords[1][m_totpoints-2]
    							  + ((m_coords[0][m_totpoints-1]-m_coords[0][m_totpoints-2])/((m_coords[0][m_totpoints-3]-m_coords[0][m_totpoints-1])*(m_coords[0][m_totpoints-3]-m_coords[0][m_totpoints-2])))*m_coords[1][m_totpoints-3];
    			}

    			// Fill the rest of the Slope Array m_M
    			SolveTriangularMatrixForSlope();

    			// Compute coefficients a,b,c,d
    			ComputeCoefficients();

    			// Compute normals at given coordinates x,y from computed spline
    			ComputeSplineNormals();

    			// Compute curvature at given coordinates x,y from computed spline
    			ComputeSplineCurvature();

    			//int nintpoints = 10;
        		//CubicSpline->WriteMatlabFiles("Splinefreesurface",nintpoints);
    	}

        void CubicSpline::SetCoords(Array<OneD,NekDouble> &x, Array<OneD,NekDouble> &y)
        {
        	NekDouble h;
        	int i,j;
        	int cnt=0;

        	// Copy nonzero distance x into coords
    		for(i=0;i<x.num_elements();i++)
    		{
    			if(i==x.num_elements()-1)
    			{
    				h = abs(x[i] - x[i-1]);
    			}
    			else
    			{
    				h = abs(x[i+1] - x[i]);
    			}
    			if(h>=1e-10)
    			{
    				m_coords[0][cnt] = x[i];
    				m_coords[1][cnt] = y[i];
    				//cout << "(i,cnt)="<< i << ","<< cnt << ":(xc,yc)=(" << m_coords[0][cnt] << "," << m_coords[1][cnt] << ")" << endl;
    				cnt++;
    			}
    		}
        }

        NekDouble CubicSpline::GetSplineValue(NekDouble x)
        {
        	NekDouble value;
        	int i;
        	NekDouble h;
        	// Find which part of spline x coordinate belongs to
        	for(i=0;i<m_totpoints-1;i++)
        	{
        		if(x >= m_coords[0][i] && x<=m_coords[0][i+1])
        		{
        			h = x-m_coords[0][i];
        			value = m_a[i]*h*h*h + m_b[i]*h*h + m_c[i]*h + m_d[i];
        			//cout << "(x,y) = (" << x << "," << value << "), (a,b,c,d)=(" << m_a[i] << "," << m_b[i] << "," << m_c[i] << "," << m_d[i] << ")" << endl;
        		}
        	}

        	return value;

        }

        NekDouble CubicSpline::GetSplineDerivativeValue(NekDouble x)
		{
			NekDouble value;
			int i;
			NekDouble h;
			// Find which part of spline x coordinate belongs to
			for(i=0;i<m_totpoints-1;i++)
			{
				if(x >= m_coords[0][i] && x<=m_coords[0][i+1])
				{
					h = x-m_coords[0][i];
					value = 3*m_a[i]*h*h + 2*m_b[i]*h + m_c[i];
					//cout << "(x,y) = (" << x << "," << value << "), (a,b,c,d)=(" << m_a[i] << "," << m_b[i] << "," << m_c[i] << "," << m_d[i] << ")" << endl;
				}
			}

			return value;

		}

        NekDouble CubicSpline::GetSplineDiff(int id)
		{
			return m_c[id];
		}

        // if index=0 means index could not be found
        int CubicSpline::GetSplineCoordOffset(NekDouble xcoord)
	    {
        	int index=0;
        	int i;
			// Find index in m_coords/m_normals to which xstart belongs
			for(i=0;i<m_totpoints;i++)
			{
				if(xcoord == m_coords[0][i])
				{
					index = i;
				}
			}

			return index;
	    }

        //xcoord is xmin of edge, n_quad number of points along edge, locnormals
        void CubicSpline::GetSplineLocalNormals(NekDouble xcoord, int n_quad, Array<OneD, Array<OneD, NekDouble> > &locnormals)
	    {
        	int offset = GetSplineCoordOffset(xcoord);
        	int i;
        	cout << "offset=" << offset << endl;
			// Find index in m_coords/m_normals to which xstart belongs
			for(i=0;i<n_quad;i++)
			{
				locnormals[0][i] = m_normals[0][offset+i];
				locnormals[1][i] = m_normals[1][offset+i];
				cout << "(x,y)=" << "(" << m_coords[0][offset+i] << "," << m_coords[1][offset+i] << "),\t(nx,ny)=(" << m_normals[0][offset+i] << "," << m_normals[1][offset+i] << ")" << endl;
			}
	    }

        // if index=0 means index could not be found
        int CubicSpline::GetSplineNormal(NekDouble xcoord, NekDouble &nx, NekDouble &ny)
 	    {
         	int index=0;
         	int i;
 			// Find index in m_coords/m_normals to which xstart belongs
 			for(i=0;i<m_totpoints;i++)
 			{
 				if(xcoord == m_coords[0][i])
 				{
 					nx=m_normals[0][i];
 					ny=m_normals[1][i];
 				}
 			}

 			return index;
 	    }

        // if index=0 means index could not be found
		int CubicSpline::GetSplineCurvature(NekDouble xcoord, NekDouble &curv)
		{
			int index=0;
			int i;
			// Find index in m_coords/m_normals to which xstart belongs
			for(i=0;i<m_totpoints;i++)
			{
				if(xcoord == m_coords[0][i])
				{
					curv=m_curv[i];
				}
			}

			return index;
		}

    	void CubicSpline::WriteMatlabFiles(string fname, int cnt)
    	{
    		stringstream convert;
    		convert <<  fname << "_" << cnt << ".txt";
    		string filename = convert.str();

			cout << "Writing... " << filename << endl;

    		ofstream outfile(filename.c_str());
    		int nvariables = 2;


			outfile << "x" << "\t";
			outfile.width(10);
			outfile << "y" << "\t";
			outfile.width(10);
			outfile << "nx" << "\t";
			outfile.width(10);
			outfile << "ny" << "\t";
			outfile.width(10);
			outfile << "curv" << "\n";
			for(int i=0;i<m_totpoints;i++)
			{
				outfile.width(10);
				outfile << m_coords[0][i] << "\t";
				outfile.width(10);
				outfile << m_coords[1][i] << "\t";
				outfile.width(10);
				outfile << m_normals[0][i] << "\t";
				outfile.width(10);
				outfile << m_normals[1][i] << "\t";
				outfile.width(10);
				outfile << m_curv[i] << "\n";
			}

			outfile.close();
    	}

    	void CubicSpline::WriteMatlabFilesInterpolation(string fname, int cnt, int nintpoints)
		{
    		int points = nintpoints*(m_totpoints-1)+1;

			string filename;
			stringstream convert;
			convert <<  fname << "_" << cnt << "interpolated.txt";
			filename = convert.str();
			cout << "Writing... " << filename << endl;

			ofstream outfile(filename.c_str());

			Array<OneD,NekDouble> x(points);
			Array<OneD,NekDouble> y(points);
			int k = 0;
			int j;
			NekDouble h;

			for(int i=0;i<m_totpoints-1;i++)
			{
				x[k]=m_coords[0][i];
				h = m_coords[0][i+1] - m_coords[0][i];

				for(j=0;j<nintpoints;j++)
				{
					x[k]=m_coords[0][i]+h*j/nintpoints;
					y[k] = GetSplineValue(x[k]);
					k++;
				}
				//k = k+nintpoints;
				//k++;
			}
			//k++;
			x[k]=m_coords[0][m_totpoints-1];
			y[k]= GetSplineValue(x[k]);

			outfile.width(10);
			outfile << "xint" << "\t";
			outfile.width(10);
			outfile << "yint" << "\n";
			for(int i=0;i<points;i++)
			{
				outfile.width(10);
				outfile << x[i] << "\t";
				outfile.width(10);
				outfile << y[i] << "\n";
			}

			outfile.close();
		}

    	CubicSplineSharedPtr CubicSpline::SetUpSpline(MultiRegions::ExpListSharedPtr &pBoundaryExpansion, SplineBoundaryType bcleft, SplineBoundaryType bcright, NekDouble bcvalueleft, NekDouble bcvalueright)
    	{
    		CubicSplineSharedPtr CubicSpline;
    		int nq= pBoundaryExpansion->GetTotPoints();

    		Array<OneD,NekDouble> x(nq,0.0);
    		Array<OneD,NekDouble> y(nq,0.0);
    		Array<OneD,NekDouble> z(nq,0.0);

    		pBoundaryExpansion->GetCoords(x,y,z);

    		CubicSpline = MemoryManager<LibUtilities::CubicSpline>::AllocateSharedPtr(x,y,bcleft,bcright,bcvalueleft,bcvalueright);

    		return CubicSpline;
    	}

    	void CubicSpline::UpdateSpline(MultiRegions::ExpListSharedPtr &pBoundaryExpansion)
    	{
    		int i,j;
    		int nq= pBoundaryExpansion->GetTotPoints();

    		Array<OneD,NekDouble> x(nq,0.0);
    		Array<OneD,NekDouble> y(nq,0.0);
    		Array<OneD,NekDouble> z(nq,0.0);

    		pBoundaryExpansion->GetCoords(x,y,z);

    		UpdateCubicSpline(x,y);

    	}

    	void CubicSpline::UpdateSplineCornerPoints(MultiRegions::ExpListSharedPtr &pBoundaryExpansion)
    	{
    		int i,j;
    		int el;
    		int nel = pBoundaryExpansion->GetExpSize();
    		StdRegions::StdExpansion1DSharedPtr EdgeExp;

    		Array<OneD,NekDouble> x(nel,0.0);
			Array<OneD,NekDouble> y(nel,0.0);
			Array<OneD,NekDouble> z(nel,0.0);

    		for(el = 0; el < pBoundaryExpansion->GetExpSize(); ++el)
    		{
    			EdgeExp =  boost::dynamic_pointer_cast<StdRegions::StdExpansion1D> (pBoundaryExpansion->GetExp(el));
        		int nquad_e = EdgeExp->GetNumPoints(0);
				Array<OneD,NekDouble> xedge(nquad_e,0.0);
				Array<OneD,NekDouble> yedge(nquad_e,0.0);
				Array<OneD,NekDouble> zedge(nquad_e,0.0);

				EdgeExp->GetCoords(xedge,yedge,zedge);
				int middle = nquad_e/2;
				x[el] = xedge[middle];
				y[el] = yedge[middle];
    		}

    		//bubble sort coordinates after increasing x
			for(i = 0; i < x.num_elements(); i++)
			{
				for(j = 1; j < x.num_elements(); j++)
				{
					if(x[j] < x[j-1])
					{
						std::swap(x[j], x[j-1]);
						std::swap(y[j], y[j-1]);
					}
				}
			}

       		// Copy coords x,y into m_coords
			SetCoords(x,y);

			// Fill the rest of the Slope Array m_M
			SolveTriangularMatrixForSlope();

			// Compute coefficients a,b,c,d
			ComputeCoefficients();

			// Compute normals at given coordinates x,y from computed spline
			ComputeSplineNormals();

			// Compute curvature at given coordinates x,y from computed spline
			ComputeSplineCurvature();

			int nintpoints = 10;
    		//CubicSpline->WriteMatlabFiles("Splinefreesurface",nintpoints);
    	}

    	LibUtilities::CubicSplineSharedPtr CubicSpline::SetUpSplineCornerPoints(MultiRegions::ExpListSharedPtr &pBoundaryExpansion)
    	{
    		SplineBoundaryType bcleft,bcright;
    		CubicSplineSharedPtr CubicSpline;

    		bcleft = LibUtilities::eNatural;
    		bcright = LibUtilities::eNatural;

       		int nel = pBoundaryExpansion->GetExpSize();
			StdRegions::StdExpansion1DSharedPtr EdgeExp;

			Array<OneD,NekDouble> x(nel,0.0);
			Array<OneD,NekDouble> y(nel,0.0);
			Array<OneD,NekDouble> z(nel,0.0);

			for(int el = 0; el < pBoundaryExpansion->GetExpSize(); ++el)
			{
				EdgeExp =  boost::dynamic_pointer_cast<StdRegions::StdExpansion1D> (pBoundaryExpansion->GetExp(el));
				int nquad_e = EdgeExp->GetNumPoints(0);
				Array<OneD,NekDouble> xedge(nquad_e,0.0);
				Array<OneD,NekDouble> yedge(nquad_e,0.0);
				Array<OneD,NekDouble> zedge(nquad_e,0.0);

				EdgeExp->GetCoords(xedge,yedge,zedge);
				int middle = nquad_e/2;
				x[el] = xedge[middle];
				y[el] = yedge[middle];
			}

			NekDouble bcr=0,bcl=0;

    		CubicSpline = MemoryManager<LibUtilities::CubicSpline>::AllocateSharedPtr(x,y,bcleft,bcright,bcl,bcr);

    		//CubicSpline->WriteMatlabFiles("Splinefreesurface",nintpoints);

    		return CubicSpline;
    	}

    }
}


#endif //NEKTAR_LIBUTILITIES_SPLINE_CPP
