////////////////////////////////////////////////////////////////////////////////
//
// File LogConformation.cpp
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
// Description: Viscoelastic constitutive equations in log conformation formulation
//
///////////////////////////////////////////////////////////////////////////////

#include <ViscoelasticSolver/EquationSystems/LogConformation.h>


namespace Nektar
{
    /**
     * Constructor. Creates ...
     *
     * \param
     * \param
     */

LogConformation::LogConformation(const LibUtilities::SessionReaderSharedPtr &pSession)
    :Viscoelastic(pSession)
    {
    }

void LogConformation::v_InitObject()
 {

 			Viscoelastic::v_InitObject();
    	        int i,j;
    	        int numfields = m_fields.num_elements();
    	        if(m_equationType==eViscoelasticLogc)
    	        {
    			std::string velids[] = {"sxx","sxy","syy"};

    	        // Set up stress components to point to the last three m_fields;
    			// pressure has to be the last field
    			// Be aware implementation is for 2D only so far!!!
    			// Set up int value to identify which m_fields contain the stress tensor components
    			int StressDim = m_stress.num_elements();
    	        m_logc = Array<OneD,int>(StressDim);

    	        for(i = 0; i < StressDim; ++i)
    	        {
    	            for(j = 0; j < numfields; ++j)
    	            {
    	                std::string var = m_boundaryConditions->GetVariable(j);
    	                if(NoCaseStringCompare(velids[i],var) == 0)
    	                {
    	                    m_logc[i] = j;
    	                    break;
    	                }

    	                if(j == numfields)
    	                {
    	                    std::string error = "Failed to find field: " + var;
    	                    ASSERTL0(false,error.c_str());
    	                }
    	            }
    	        }

    	        int phystot = m_fields[0]->GetTotPoints();

    	        /* Change number of Convective fields because elastic extra stress tensor components t_xx
    	         * t_xy, and t_yy are not convected. Instead the logconformation tensor is convected.
    	         */
    	       // m_nConvectiveFields = numfields-(1+StressDim);

//TODO: Set boundary conditions here from Inflow values of the stresses
    	        SetBoundaryConditionsLogc(0.0);

    	        }
    }

/*!
 Computes the analytical solution of the
 spectral decomposition of an symmetric positive definite tensor T in 2D.
 @param[in] phystot number of physical points for which the spectral decomposition is computed
 @param[in] T in terms of the tensor components T[0]=Txx, T[1]=Txy, T[2]=Tyy
 @param[out]  EigVal  give back the Eigenvalues EigValue[0]=lambda1 and EigValue[1]=lambda2
 @param[out]  EigVec   contains the information about the eigenvector tensor R in terms of
					   EigVec[0]=e1^2, EigVec[1]=e1*e2, EigVec[2]=e2^2
 */

  void LogConformation::CompSpectralDecomp2D(int phystot, const Array<OneD, Array<OneD, NekDouble> > &T, Array<OneD, Array<OneD, NekDouble> > &EigVal,
                       Array<OneD, Array<OneD, NekDouble> > &EigVec)
  {
	  Array<OneD, NekDouble> tmp1(phystot,0.0);
	  Array<OneD, NekDouble> tmp2(phystot,0.0);
	  int i;
	  NekDouble Diff;
	  //TODO: change this back!!!
	  NekDouble E1,E2;

	  //calc eigenvalues
	  // tmp1=0.5*sqrt((Txx-Tyy)^2+4*Txy^2)
	  Vmath::Vsub(phystot,T[0],1,T[2],1,tmp1,1);
	  Vmath::Smul(phystot,2.0,T[1],1,tmp2,1);
	  Vmath::Vvtvvtp (phystot,tmp1,1,tmp1,1,tmp2,1, tmp2,1,tmp1, 1);
	  Vmath::Vsqrt(phystot, tmp1, 1, tmp1, 1);
	  Vmath::Smul(phystot,0.5,tmp1,1,tmp1,1);

	  // tmp2=0.5*(Txx+Tyy)
	  Vmath::Vadd(phystot,T[0],1,T[2],1,tmp2,1);
	  Vmath::Smul(phystot,0.5,tmp2,1,tmp2,1);

	  // lambda1=tmp1 + tmp2, lambda2=tmp1 - tmp2
	  Vmath::Vadd(phystot,tmp2,1,tmp1,1,EigVal[0],1);
	  Vmath::Vsub(phystot,tmp2,1,tmp1,1,EigVal[1],1);

	  //calc Eigenvector in form of e1^2, e1e2, e2^2
	  // lambda1 - Txx
	  Vmath::Vsub(phystot,EigVal[0],1,T[0],1,tmp1,1);

	  for(i=0;i<phystot;i++)
	  {
		 // E1 = fabs(EigVal[0][i]);
		 // E2 = fabs(EigVal[1][i]);
		  Diff = EigVal[0][i]-EigVal[1][i];
		 // Diff = E1 -E2;
		  Diff = fabs(Diff);
		  if(Diff<1e-10)
		  {
			  EigVec[0][i]=1.0;
			  EigVec[1][i]=0.0;
			  EigVec[2][i]=0.0;
			  // also force eigenvalues to be the same in that case
			  // set second equals first because second has the danger of being less than zero
			  EigVal[0][i]=T[0][i];
			  EigVal[1][i]=EigVal[0][i];
		  }
		  else
		  {
			  tmp2[i] = T[1][i]*T[1][i]+tmp1[i]*tmp1[i];
			  EigVec[0][i]=(T[1][i]*T[1][i])/tmp2[i];
			  EigVec[1][i]=T[1][i]*tmp1[i]/tmp2[i];
			  EigVec[2][i]=tmp1[i]*tmp1[i]/tmp2[i];
		  }
	  }
/*	  // e1^2 = Txy*Txy
	  Vmath::Vmul(phystot, T[1], 1, T[1], 1,  EigVec[0], 1);
	  // e1e2 = tmp1*Txy
	  Vmath::Vmul(phystot, tmp1, 1, T[1], 1,  EigVec[1], 1);
	  // e2^2 = tmp1*tmp1
	  Vmath::Vmul(phystot, tmp1, 1, tmp1, 1,  EigVec[2], 1);

	  // tmp1 = Txy*Txy + tmp1*tmp1
	  Vmath::Vadd(phystot,EigVec[0],1,EigVec[2],1,tmp1,1);

	  for(i=0;i<3;i++)
	  {
		  Vmath::Vdiv(phystot, EigVec[i], 1 , tmp1 , 1 , EigVec[i], 1);
	  }*/

  }

	/*!
	 Computes the transformation of a 2D symmetric positive definite tensor from its principal frame
	 into the global frame:
	 T = R*Tp*R^t (R: matrix of eigenvectors, Tp: tensor T in prin)
	 @param[in] phystot number of physical points
	 @param[in] Tp tensor in its principal frame in terms of the tensor components Tp[0]=Tpxx, Tp[1]=Tpxy, Tp[2]=Tpyy
	 @param[in] EigVec Eigenvector tensor R in terms of EigVec[0]=e1^2, EigVec[1]=e1*e2, EigVec[2]=e2^2
	 @param[out]  T  tensor components T[0]=Txx, T[1]=Txy, T[2]=Tyy in the global frame
	 */

	  void LogConformation::PrincipalToGlobalTrans(int phystot, const Array<OneD, Array<OneD, NekDouble> > &Tp, const Array<OneD, Array<OneD, NekDouble> > &EigVec,
                           Array<OneD, Array<OneD, NekDouble> > &T)
	  {
		  // Txx

		  // Tyy = 2.0*e1e2*Tpxy
		  Vmath::Vmul(phystot, Tp[1], 1, EigVec[1], 1,  T[2], 1);
		  Vmath::Smul(phystot,2.0,T[2],1,T[2],1);

		  // Txx = Tpxx*e1^2-(2.0*e1e2*Tpxy) + Tpyy*e2^2
		  Vmath::Vvtvm(phystot, Tp[0], 1, EigVec[0],1, T[2], 1,T[0],1);
		  Vmath::Vvtvp(phystot, Tp[2], 1, EigVec[2],1, T[0], 1,T[0],1);

		  // Tyy= Tpxx*e2^2+(2.0*e1e2*Tpxy) + Tpyy*e1^2
		  Vmath::Vvtvp(phystot, Tp[0], 1, EigVec[2],1, T[2], 1,T[2],1);
		  Vmath::Vvtvp(phystot, Tp[2], 1, EigVec[0],1, T[2], 1,T[2],1);

		  // Txy
		  Vmath::Vsub(phystot,EigVec[0],1,EigVec[2],1,T[1],1);
		  Vmath::Vmul(phystot, Tp[1], 1, T[1], 1,  T[1], 1);
		  Vmath::Vvtvp(phystot, Tp[0], 1, EigVec[1],1, T[1], 1,T[1],1);
		  Vmath::Svvtvp(phystot,-1.0,Tp[2],1,EigVec[1],1,T[1], 1, T[1], 1);

	  }

	  void LogConformation::GlobalToPrincipalTrans(int phystot, const Array<OneD, Array<OneD, NekDouble> > &T, const Array<OneD, Array<OneD, NekDouble> > &EigVec,
	  							   Array<OneD, Array<OneD, NekDouble> > &Tp)
	  {
		  int components=T.num_elements();

		  /* symmetric tensor with 3 relevant components T[0]=Txx, T[1]=Txy, T[2]=Tyy */
		 if(components==3)
		 {
		  // Tpyy = 2.0*e1e2*Txy
		  Vmath::Vmul(phystot, T[1], 1, EigVec[1], 1,  Tp[2], 1);
		  Vmath::Smul(phystot,2.0,Tp[2],1,Tp[2],1);

		  // Tpxx = Txx*e1^2+(2.0*e1e2*Txy) + Tyy*e2^2
		  Vmath::Vvtvp(phystot, T[0], 1, EigVec[0],1, Tp[2], 1,Tp[0],1);
		  Vmath::Vvtvp(phystot, T[2], 1, EigVec[2],1, Tp[0], 1,Tp[0],1);

		  // Tpyy= Txx*e2^2-(2.0*e1e2*Txy) + Tyy*e1^2
		  Vmath::Vvtvm(phystot, T[0], 1, EigVec[2],1, Tp[2], 1,Tp[2],1);
		  Vmath::Vvtvp(phystot, T[2], 1, EigVec[0],1, Tp[2], 1,Tp[2],1);

		  // Tpxy
		  Vmath::Vsub(phystot,EigVec[0],1,EigVec[2],1,Tp[1],1);
		  Vmath::Vmul(phystot, T[1], 1, Tp[1], 1,  Tp[1], 1);
		  Vmath::Vvtvp(phystot, T[2], 1, EigVec[1],1, Tp[1], 1,Tp[1],1);
		  Vmath::Svvtvp(phystot,-1.0,T[0],1,EigVec[1],1,Tp[1], 1, Tp[1], 1);
		 }
		 /* asymmetric tensor with 4 relevant components T[0]=Txx, T[1]=Txy, T[2]=Tyx, T[3]=Tyy */
		 else
		 {
			 // Tpyy = e1e2*Txy + Tyx*e1e2
		  Vmath::Vmul(phystot, T[1], 1, EigVec[1], 1,  Tp[3], 1);
		  Vmath::Vvtvp(phystot, T[2], 1, EigVec[1],1, Tp[3], 1,Tp[3],1);

		  // Tpxx = Txx*e1^2+e1e2*Txy + Tyx*e1e2 + Tyy*e2^2
		  Vmath::Vvtvp(phystot, T[0], 1, EigVec[0],1, Tp[3], 1,Tp[0],1);
		  Vmath::Vvtvp(phystot, T[3], 1, EigVec[2],1, Tp[0], 1,Tp[0],1);

		  // Tpyy= Txx*e2^2-(e1e2*Txy + Tyx*e1e2) + Tpyy*e1^2
		  Vmath::Vvtvm(phystot, T[0], 1, EigVec[2],1, Tp[3], 1,Tp[3],1);
		  Vmath::Vvtvp(phystot, T[3], 1, EigVec[0],1, Tp[3], 1,Tp[3],1);

		  // Tpxy=-Txx*e1e2+ e1^2*Txy - Tyx*e2^2 + Tyy*e1e2
		  Vmath::Vvtvvtp(phystot,T[0],1,EigVec[1],1,T[2],1,EigVec[2],1,Tp[1], 1);
		  Vmath::Vvtvm(phystot, T[1], 1, EigVec[0],1, Tp[1], 1,Tp[1],1);
		  Vmath::Vvtvp(phystot, T[3], 1, EigVec[1],1, Tp[1], 1,Tp[1],1);

		  // Tpyx=-Txx*e1e2 - e2^2*Txy + Tyx*e1^2 + Tyy*e1e2
		  Vmath::Vvtvvtp(phystot,T[0],1,EigVec[1],1,T[1],1,EigVec[2],1,Tp[2], 1);
		  Vmath::Vvtvm(phystot, T[2], 1, EigVec[0],1, Tp[2], 1,Tp[2],1);
		  Vmath::Vvtvp(phystot, T[3], 1, EigVec[1],1, Tp[2], 1,Tp[2],1);

		 }
	  }

	  /*!
	  	 Computes conformation tensor from the elastic stress tensor components:
	  	 c = Wi/(1-beta)*tau + I
	  	 @param[in] phystot number of physical points
	  	 @param[in] Elastic stress tensor in tensor components Elastic[0]=tauxx, Elastic[1]=tauxy, Elastic[2]=tauyy
	  	 @param[out]  C conformation tensor components C[0]=Cxx, C[1]=Txy, C[2]=Tyy
	  	 */
	  void LogConformation::ElasticToConformationTensor(int phystot, const Array<OneD, Array<OneD, NekDouble> > &Elastic, Array<OneD, Array<OneD, NekDouble> > &C)
	  {
		  NekDouble G = m_weissenberg/(1-m_beta);
		  int StressDim = m_stress.num_elements();

		  for(int i=0;i<StressDim;i++)
		  {
			  Vmath::Smul(phystot, G, Elastic[i],1, C[i], 1);
		  }

		  Vmath::Sadd(phystot,1.0,C[0],1, C[0], 1);
		  Vmath::Sadd(phystot,1.0,C[2],1, C[2], 1);
	  }

	  /*!
	  	 Computes log-conformation tensor from the elastic stress tensor components:
	  	 1. c = Wi/(1-beta)*tau + I, 2. c=R*Lambda*R^t, s=R*log(Lambda)*R^t
	  	 @param[in] phystot number of physical points
	  	 @param[in] Elastic stress tensor in tensor components Elastic[0]=tauxx, Elastic[1]=tauxy, Elastic[2]=tauyy
	  	 @param[out]  C conformation tensor components C[0]=Cxx, C[1]=Txy, C[2]=Tyy
	  	 */
	  void LogConformation::ElasticToLogConformationTensor(int phystot, const Array<OneD, Array<OneD, NekDouble> > &Elastic, Array<OneD, Array<OneD, NekDouble> > &Logc)
	  {
		  /* Allocate memory for:
		   * Conformation Tensor 3 Component
		   * EigenVector Matrix 3 components
		   * EigenValues 2
		   * all of them need phystot entries*/
		  int StressDim = m_stress.num_elements();
		  int i;
		  Array<OneD, Array<OneD, NekDouble> > C(StressDim);
		  Array<OneD, Array<OneD, NekDouble> > EigVec(StressDim);
		  Array<OneD, Array<OneD, NekDouble> > EigVal(m_spacedim);

		  C[0] = Array<OneD, NekDouble>(phystot*StressDim);
		  EigVec[0] = Array<OneD, NekDouble>(phystot*StressDim);
		  for(i = 1; i < StressDim; ++i)
		  {
		      C[i] = C[i-1] + phystot;
		      EigVec[i] = EigVec[i-1] + phystot;
		  }

		  EigVal[0] = Array<OneD, NekDouble>(phystot*m_spacedim);
		  for(i = 1; i < m_spacedim; ++i)
		  {
			  EigVal[i] = EigVal[i-1] + phystot;
		  }

		  ElasticToConformationTensor(phystot,Elastic,C);

		  CompSpectralDecomp2D(phystot,C, EigVal,EigVec);

		  //take logarithm of EigenValues to get eigenvalues of s
		  for(i = 0; i < m_spacedim; ++i)
		  {
			  Vmath::Vlog(phystot,EigVal[i],1,EigVal[i],1);
		  }

		  //transform s into global frame use C to create matrix components of principal frame
		  // in terms of eigenvalues because C is no longer needed
		  C[0] = EigVal[0];
		  Vmath::Zero(phystot,C[1],1);
		  C[2] = EigVal[1];

		  PrincipalToGlobalTrans(phystot,C,EigVec,Logc);

	/*	  cout << "logcxxintrans=" ;
		  				 		  for(i=0;i<phystot;i++)
		  				 		  {

		  				 				cout << Logc[0][i] << endl;
		  				 		  }
	 cout << "logcxyintrans=" ;
									  for(i=0;i<phystot;i++)
									  {

											cout << Logc[1][i] << endl;
									  }
		 cout << "logcyyintrans=" ;
										  for(i=0;i<phystot;i++)
										  {

												cout << Logc[2][i] << endl;
										  }
*/

	  }

	  void LogConformation::LogConformationTensorToElastic(int phystot, const Array<OneD, Array<OneD, NekDouble> > &Logc, Array<OneD, Array<OneD, NekDouble> > &Elastic)
	  {
		  /* Allocate memory for:
		   * Conformation Tensor 3 Component
		   * EigenVector Matrix 3 components
		   * EigenValues 2
		   * all of them need phystot entries*/
		  int i;
		  NekDouble G = (1-m_beta)/m_weissenberg;
		  int StressDim = m_logc.num_elements();

		  Array<OneD, Array<OneD, NekDouble> > C(StressDim);

		  C[0] = Array<OneD, NekDouble>(phystot*StressDim);

		  for(i = 1; i < StressDim; ++i)
		  {
			  C[i] = C[i-1] + phystot;
		  }

		  LogConformationToConformationTensor(phystot,Logc,C);

		  for(int i=0;i<StressDim;i++)
		  {
			  Vmath::Smul(phystot, G, C[i],1, Elastic[i], 1);
		  }

		  Vmath::Sadd(phystot,-G,Elastic[0],1, Elastic[0], 1);
		  Vmath::Sadd(phystot,-G,Elastic[2],1, Elastic[2], 1);

	/*	  cout << "logcxxintrans=" ;
		  				 		  for(i=0;i<phystot;i++)
		  				 		  {

		  				 				cout << Logc[0][i] << endl;
		  				 		  }
	 cout << "logcxyintrans=" ;
									  for(i=0;i<phystot;i++)
									  {

											cout << Logc[1][i] << endl;
									  }
		 cout << "logcyyintrans=" ;
										  for(i=0;i<phystot;i++)
										  {

												cout << Logc[2][i] << endl;
										  }
*/

	  }

	  void LogConformation::LogConformationToConformationTensor(int phystot, const Array<OneD, Array<OneD, NekDouble> > &Logc, Array<OneD, Array<OneD, NekDouble> > &C)
	  {
		  int StressDim = m_stress.num_elements();
		  int i;
		  Array<OneD, Array<OneD, NekDouble> > Cp(StressDim);
		  Array<OneD, Array<OneD, NekDouble> > EigVec(StressDim);
		  Array<OneD, Array<OneD, NekDouble> > EigVal(m_spacedim);

		  Cp[0] = Array<OneD, NekDouble>(phystot*StressDim);
		  EigVec[0] = Array<OneD, NekDouble>(phystot*StressDim);
		  for(i = 1; i < StressDim; ++i)
		  {
		      Cp[i] = Cp[i-1] + phystot;
		      EigVec[i] = EigVec[i-1] + phystot;
		  }

		  EigVal[0] = Array<OneD, NekDouble>(phystot*m_spacedim);
		  for(i = 1; i < m_spacedim; ++i)
		  {
			  EigVal[i] = EigVal[i-1] + phystot;
		  }

		  CompSpectralDecomp2D(phystot,Logc, EigVal,EigVec);

		  //take logarithm of EigenValues to get eigenvalues of s
		  for(i = 0; i < m_spacedim; ++i)
		  {
			  Vmath::Vexp(phystot,EigVal[i],1,EigVal[i],1);
		  }

		  //transform s into global frame use C to create matrix components of principal frame
		  // in terms of eigenvalues because C is no longer needed
		  Cp[0] = EigVal[0];
		  Vmath::Zero(phystot,Cp[1],1);
		  Cp[2] = EigVal[1];

		  PrincipalToGlobalTrans(phystot,Cp,EigVec,C);

	  }

		void LogConformation::SetBoundaryConditionsLogc(NekDouble time)
		{
			int nvariables = m_fields.num_elements();
			int StressDim = m_logc.num_elements();
			int cnt = 0;

			//cout << "time" << time << endl;

			// loop over Boundary Regions for logc component sxx
			for(int n = 0; n < m_fields[m_logc[0]]->GetBndConditions().num_elements(); ++n)
			{
				// Set Dirichlet Inflow Boundary Values for the logarithm of the conformation tensor
				// computed from the elastic stress inflow values
				if (m_fields[m_logc[0]]->GetBndConditions()[n]->GetUserDefined() == SpatialDomains::eStressDirichlet)
				{
					if(m_projectionType == MultiRegions::eGalerkin)
					{
						SetBoundaryInflowLogc(n,cnt,time);
					}
					else
					{
						SetBoundaryInflowLogcDG(n,cnt,time);

						/*int npoints = m_fields[m_logc[0]]->GetBndCondExpansions()[n]->GetNpoints();
				        Array<OneD,NekDouble> x0(npoints,0.0);
				        Array<OneD,NekDouble> x1(npoints,0.0);
				        Array<OneD,NekDouble> x2(npoints,0.0);

				        m_fields[m_logc[0]]->GetBndCondExpansions()[n]->GetCoords(x0,x1,x2);


						for(int j=0;j<npoints;j++)
						{
							cout << "(x,y)=" << x0[j] << "," << x1[j] << "\t";
							for(int i=0; i<StressDim; i++)
							{
								cout<<m_fields[m_logc[i]]->GetBndCondExpansions()[n]->GetPhys()[j] << ",";
							}
							cout << endl;
						} */
					}

					// forward transform to fill the modal coeffs
					for(int i=0; i<StressDim; i++)
					{
						m_fields[m_logc[i]]->GetBndCondExpansions()[n]->FwdTrans_BndConstrained(m_fields[m_logc[i]]->GetBndCondExpansions()[n]->GetPhys(),m_fields[m_logc[i]]->GetBndCondExpansions()[n]->UpdateCoeffs());

					}

				}

				cnt +=m_fields[m_logc[0]]->GetBndCondExpansions()[n]->GetExpSize();
			}
		}
/*!
 * Set Inflow values of the logc tensor to fitting dirichlet values derived from the elastic stress
 * tensor.
 */
		void LogConformation::SetBoundaryInflowLogcDG(int bcRegion, int cnt, NekDouble time)
		{
			int nvariables      = m_stress.num_elements(); //only stress components tauxx and tauxy have to be computed
		    int nTraceNumPoints = m_fields[m_stress[0]]->GetTrace()->GetNpoints();
		    int i;
		    int static count=0;
	/*	    Array<OneD, Array<OneD, NekDouble> > Fwd(nvariables);
	        Fwd[0] = Array<OneD, NekDouble>(nTraceNumPoints*nvariables);
	        Fwd[1] = Fwd[0] + nTraceNumPoints;
	        Fwd[2] = Fwd[1] + nTraceNumPoints; */

		    // get physical values of the forward trace (from exp to phys)
	        // to still have other bcs
	        // Fill Fwd with values of the elastic stress tensor
		/*    for (int i = 0; i < nvariables; ++i)
		    {
		    	m_fields[m_stress[i]]->ExtractTracePhys(Fwd[i]);
		    } */

		    /* id1: id of the start point of the boundary expansion
		     * id2: id of the start point for the trace values
		     */
		    for(int e = 0; e < m_fields[m_stress[0]]->GetBndCondExpansions()[bcRegion]->GetExpSize(); ++e)
		    {
		    	int npoints = m_fields[m_stress[0]]->GetBndCondExpansions()[bcRegion]->GetExp(e)->GetNumPoints(0);
		    	int id1  = m_fields[m_stress[0]]->GetBndCondExpansions()[bcRegion]->GetPhys_Offset(e) ;

		    //	cout << "id1" << id1;
		    //	int id2  = m_fields[m_stress[0]]->GetTrace()->GetPhys_Offset(m_fields[m_stress[0]]->GetTraceMap()->GetBndCondCoeffsToGlobalCoeffsMap(cnt+e));

		   /* 	Array<OneD,NekDouble> x0(npoints,0.0);
		    	Array<OneD,NekDouble> x1(npoints,0.0);
		    	Array<OneD,NekDouble> x2(npoints,0.0);

		    	m_fields[m_stress[0]]->GetBndCondExpansions()[bcRegion]->GetExp(e)->GetCoords(x0,x1,x2);

		    	NekDouble y; */
			    Array<OneD, Array<OneD, NekDouble> > Elastic(nvariables);
		        Elastic[0] = Array<OneD, NekDouble>(npoints*nvariables);
		        for(i=1;i<nvariables;i++)
		        {
		        	Elastic[i] = Elastic[i-1] + npoints;
		        }

		    	/* Copy elastic stress values from boundary expansion to Fwd trace
		    	 * for edge e */
		        for(i=0;i<nvariables;i++)
		        {
		        	Vmath::Vcopy(npoints,&(m_fields[m_stress[i]]->GetBndCondExpansions()[bcRegion]->GetPhys())[id1], 1,&Elastic[i][0],1);
		        }

		       // for(i=0;i<nvariables;i++)
		     	/*	      		       { cout << "elasticbefore"; //<< i ;
		     		      			  for(int j=0;j<npoints;j++)
		     		      			  {

		     		      					cout<< Elastic[0][j] << endl;
		     		      			  }
		     		      		       } */

		        ElasticToLogConformationTensor(npoints,Elastic,Elastic);

		        for(i=0;i<nvariables;i++)
		        {
		        	Vmath::Vcopy(npoints,&Elastic[i][0], 1,&(m_fields[m_logc[i]]->GetBndCondExpansions()[bcRegion]->UpdatePhys())[id1],1);
		        }
	        Array<OneD,NekDouble> x0(npoints,0.0);
	        Array<OneD,NekDouble> x1(npoints,0.0);
	        Array<OneD,NekDouble> x2(npoints,0.0);

	        m_fields[m_stress[0]]->GetBndCondExpansions()[bcRegion]->GetExp(e)->GetCoords(x0,x1,x2);
	        string outname = "Logcanalytical.txt";
	        ofstream outfileA(outname.c_str(),ios::app);

	        if(count==0)
	        {
	        // Writes header of Matlab file
	        outfileA.width(10);
	        outfileA << "xana" << "\t";
	        outfileA.width(10);
	        outfileA << "yana" << "\t";
	        outfileA.width(10);
	        outfileA << "sxxana" << "\t";
	        outfileA.width(10);
	        outfileA << "sxyana" << "\t";
	       // outfile.close();
	        outfileA.width(10);
	        outfileA << "syyana" << "\t"  << "\n";
	       // outfile.close();
	        cout << "Writing outfile: Logcanalytical.txt"  << endl;
	        count++;
	        }

	        for(i=0;i<npoints;i++)
	        {
	        	outfileA.width(10);
	        	outfileA.precision(8);
	        	outfileA << x0[i]<< "\t";
	        	outfileA.width(10);
	        	outfileA.precision(8);
	        	outfileA << x1[i]<< "\t";
	        	for(int j=0;j<nvariables;j++)
	        	{
	        		outfileA.width(10);
	        		outfileA.precision(8);
	        		outfileA<<(m_fields[m_logc[j]]->GetBndCondExpansions()[bcRegion]->GetPhys())[id1+i] << "\t";
	        	}
	        	outfileA << "\n";
	        }

	      /*  for(i=0;i<npoints;i++)
	        {
	        	cout<<"(x,y)=" << x0[i] << "," << x1[i] << "\t";
	        	cout<<"(sxx,sxy,syy)=";
	        	for(int j=0;j<nvariables;j++)
	        	{
	        		cout<<(m_fields[m_logc[j]]->GetBndCondExpansions()[bcRegion]->GetPhys())[id1+i] << ",";
	        	}
	        	cout << endl;
	        } */

	/*	       LogConformationTensorToElastic(npoints,Elastic,Elastic);

		       for(int j=0;j<nvariables;j++)
		      		       { cout << "elasticafter" << j << endl;
		      			  for(int i=0;i<npoints;i++)
		      			  {

		      					cout<< Elastic[j][i] << endl;
		      			  }
		      		       } */
/*

 Compute conformation tensor from these edge values

				// Loop on all the points of that edge
				for(int j = 0; j < npoints; j++)
				{
	// Fwd[1][kk] and Fwd[2][kk] and copy them separately into into boundary expansion
					// y=x1[j];
					int idcnt = id1+j;
					int kk = id2+j;
					m_fields[m_stress[0]]->GetBndCondExpansions()[bcRegion]->GetPhys())[id1]
					// CalcDirichletWatersTAU(y, Fwd[0][kk], Fwd[1][kk],time);
					// Test it first with steady state solution!! :)
					// SetDirichletWatersSteadyState(y, Fwd[0][kk], Fwd[1][kk]);

					//cout << "(x,y)=(" << x0[j] << "," << y << ")" << endl;
					//cout << "tauxx=" << Fwd[0][kk] << " ,tauxy=" << Fwd[1][kk] << endl;
			    }

			    Vmath::Vcopy(npoints,&Fwd[0][id2], 1,&(m_fields[m_stress[0]]->GetBndCondExpansions()[bcRegion]->UpdatePhys())[id1],1);
			    Vmath::Vcopy(npoints,&Fwd[1][id2], 1,&(m_fields[m_stress[1]]->GetBndCondExpansions()[bcRegion]->UpdatePhys())[id1],1);*/
	/*			   m_fields[0]->GetBndCondExpansions()[bcRegion]->FwdTrans_BndConstrained(m_fields[0]->GetBndCondExpansions()[bcRegion]->GetPhys(),
				    		    		m_fields[0]->GetBndCondExpansions()[bcRegion]->UpdateCoeffs());
				    m_fields[m_stress[0]]->GetBndCondExpansions()[bcRegion]->FwdTrans_BndConstrained(m_fields[m_stress[0]]->GetBndCondExpansions()[bcRegion]->GetPhys(),
				    		    		m_fields[m_stress[0]]->GetBndCondExpansions()[bcRegion]->UpdateCoeffs());
				    m_fields[m_stress[1]]->GetBndCondExpansions()[bcRegion]->FwdTrans_BndConstrained(m_fields[m_stress[1]]->GetBndCondExpansions()[bcRegion]->GetPhys(),
				    		    		m_fields[m_stress[1]]->GetBndCondExpansions()[bcRegion]->UpdateCoeffs()); */

		    }

		  }

		void LogConformation::SetBoundaryInflowLogc(int bcRegion, int cnt, NekDouble time)
		{
			int npoints = m_fields[m_stress[0]]->GetBndCondExpansions()[bcRegion]->GetNpoints();
			int nvariables      = m_stress.num_elements();
			int i;

		    Array<OneD, Array<OneD, NekDouble> > Elastic(nvariables);
	        Elastic[0] = Array<OneD, NekDouble>(npoints*nvariables);
	        for(i=1;i<nvariables;i++)
	        {
	        	Elastic[i] = Elastic[i-1] + npoints;
	        }

	    	/* Copy elastic stress values from boundary expansion to Fwd trace
	    	 * for edge e */
	        for(i=0;i<nvariables;i++)
	        {
	        	Vmath::Vcopy(npoints,&(m_fields[m_stress[i]]->GetBndCondExpansions()[bcRegion]->GetPhys())[0], 1,&Elastic[i][0],1);
	        }

	       // for(i=0;i<nvariables;i++)
	     	/*	      		       { cout << "elasticbefore"; //<< i ;
	     		      			  for(int j=0;j<npoints;j++)
	     		      			  {

	     		      					cout<< Elastic[0][j] << endl;
	     		      			  }
	     		      		       } */
	        for(int j=0;j<nvariables;j++)
	        		       { cout << "logccontinuous" << j << endl;
	        			  for(i=0;i<npoints;i++)
	        			  {
	        					cout<<(m_fields[m_logc[j]]->GetBndCondExpansions()[bcRegion]->GetPhys())[i] << endl;
	        			  }
	        		       }
	        ElasticToLogConformationTensor(npoints,Elastic,Elastic);

	        for(i=0;i<nvariables;i++)
	        {
	        	Vmath::Vcopy(npoints,&Elastic[i][0], 1,&(m_fields[m_logc[i]]->GetBndCondExpansions()[bcRegion]->UpdatePhys())[0],1);
	        }
		}

		void LogConformation::EvaluateRhsViscoelasticEquationLogc(const Array<OneD, const Array<OneD, NekDouble> > &inarray,
			    	                Array<OneD, Array<OneD, NekDouble> > &outarray)
		{
			if(m_projectionType == MultiRegions::eGalerkin)
			{
					int vstart=m_logc[0];
					int phystot = m_fields[0]->GetTotPoints();
					int ncoeffs = m_fields[0]->GetNcoeffs();
					//int nvariables = m_fields.num_elements();
					int StressDim  = m_logc.num_elements();
					int i;
					int k=15;

					Array<OneD, Array<OneD, NekDouble> > stressoutarray(StressDim);
					// Vmath::Zero(nCoeffs, outarray, 1);

					stressoutarray[0] = Array<OneD, NekDouble>(phystot*StressDim);
					for(i = 1; i < StressDim; ++i)
					{
						stressoutarray[i] = stressoutarray[i-1] + phystot;
					}

					cout << "Warning: Advection has to be reimplemented" << endl;
				    //StressAdvection(inarray,outarray,vstart);

				/*    cout << "\tLogcxxconvection: " << outarray[m_logc[0]][k] << endl;
				    cout << "\tLogcxxin: " << inarray[m_logc[0]][k] << endl;
				    cout << "\tuin: " << inarray[m_velocity[0]][k] << endl; */
					if(m_stressTreatment==eExplicit)
					{
					   // Adding source terms to outarray
					   CalcSourceTermLogc(inarray,stressoutarray);
					   for(i = 0; i < StressDim; ++i)
					   {
						   Vmath::Vadd(phystot,stressoutarray[i],1,outarray[m_logc[i]],1,outarray[m_logc[i]],1);
					   }
					}
					else
					{
						// source term is calculated implicitely
					}
					// cout << "\tLogcxxtotal: " << outarray[m_logc[0]][k] << endl;
			}
			else
			{
				/*
				 * Create an Array to store the local expansion modes coeffsout
				 * Evaluate the weak advection (maybe inverse with mass matrix not sure yet)
				 * Evaluate the source term and cast it into weak form have a physfield for that purpose
				 * or use outarray as tmp save
				 * maybe then compute inverse of mass matrix and apply BwdTrans to calculate outarray
				 * (similar to UnsteadyAdvection solver)
				 */
				int ncoeffs = m_fields[m_logc[0]]->GetNcoeffs();
				int StressDim  = m_logc.num_elements();
				int nvariables = m_fields.num_elements();
				int i;
				int k=15;
				int vstart = m_logc[0];

				/* Store results of the weak form (dimension is of N_eof number of all local expansions)
				 * in Weakoutarray
				 * and use stressinarray to calculate new stress values
				 */
		        Array<OneD, Array<OneD, NekDouble> > Weakoutarray(StressDim);
		        Weakoutarray[0] = Array<OneD, NekDouble>(ncoeffs*StressDim);
		        for(i = 1; i < StressDim; ++i)
		        {
		           Weakoutarray[i] = Weakoutarray[i-1] + ncoeffs;
		        }

		        StressWeakDGAdvection(inarray,Weakoutarray,vstart);

		    	/*for(i = 0; i < StressDim; ++i)
				{
					 m_fields[m_logc[i]]->MultiplyByElmtInvMass(Weakoutarray[i],
														Weakoutarray[i]);
					 m_fields[m_logc[i]]->BwdTrans(Weakoutarray[i],outarray[m_logc[i]]);
				}*/
#ifdef OUTPUT
		    	cout << "\tLogcxx in convection: " << outarray[m_logc[0]][k] << endl;
#endif

		    /*    for(i = 0; i < StressDim; ++i)
		        Vmath::Zero(ncoeffs,Weakoutarray[i],1); */

		        if(m_stressTreatment==eExplicit)
		        {
		          CalcWeakDGSourceTermLogc(inarray,Weakoutarray,outarray);
		        }
		        else
		        {
		        	// source term is calculated implicitely
		        }
#ifdef OUTPUT
				cout << "\tLogcxx: " << outarray[m_logc[0]][k] << endl;
#endif

		        /* This gives logc in terms of physical values */
				for(i = 0; i < StressDim; ++i)
				{
					/* Calculate d\hat{txx}/dt
					 * Coefficients of the RHS */

				     m_fields[m_logc[i]]->MultiplyByElmtInvMass(Weakoutarray[i],
				                                        Weakoutarray[i]);
				     /* Calculate dtxx/dt physical values at quadrature points
				      * physical values of RHS */

				     m_fields[m_logc[i]]->BwdTrans(Weakoutarray[i],outarray[m_logc[i]]);
				}

	/*			int nq = m_fields[0]->GetNpoints();

						                 Array<OneD,NekDouble> x0(nq);
						                 Array<OneD,NekDouble> x1(nq);
						                 Array<OneD,NekDouble> x2(nq);

						                 // get the coordinates (assuming all fields have the same
						                 // discretisation)
						                 m_fields[0]->GetCoords(x0,x1,x2);

				for(int i=0;i<nq;i++)
				{
					if(x0[i]==0.5)
					{
						cout << "i:" << i;
						cout << " y=" << x1[i] << "\tLogcxxin: " << outarray[m_logc[0]][i];
						cout << "\tLogcxyin: " << outarray[m_logc[1]][i];
						cout << "\tLogcyyin: " << outarray[m_logc[2]][i];
						cout << endl;
					}
				}
				*/





			}
		}

		void LogConformation::CalcWeakDGSourceTermLogc(const Array<OneD, const Array<OneD, NekDouble> > &inarray,
	            Array<OneD, Array<OneD, NekDouble> > &Weakoutarray, Array<OneD, Array<OneD, NekDouble> > &outarray)
		{
			int phystot = m_fields[m_logc[0]]->GetTotPoints();
			int ncoeffs = m_fields[m_logc[0]]->GetNcoeffs();
			//int nvariables = m_fields.num_elements();
			int StressDim  = m_logc.num_elements();
			int i;
			int k=15;


			Array<OneD, Array<OneD, NekDouble> > stressoutarray(StressDim);
	        Array<OneD, NekDouble> iprod(ncoeffs);
	        Array<OneD, Array<OneD, NekDouble> > C(StressDim);
	        // Vmath::Zero(nCoeffs, outarray, 1);

	        stressoutarray[0] = Array<OneD, NekDouble>(phystot*StressDim);
	        C[0] = Array<OneD, NekDouble>(phystot*StressDim);
	        for(i = 1; i < StressDim; ++i)
	        {
	        	stressoutarray[i] = stressoutarray[i-1] + phystot;
				  C[i] = C[i-1] + phystot;
	        }

	        CalcSourceTermLogc(inarray,stressoutarray);


	    /*    int nq = m_fields[m_logc[0]]->GetNpoints();

	        						                 Array<OneD,NekDouble> x0(nq);
	        						                 Array<OneD,NekDouble> x1(nq);
	        						                 Array<OneD,NekDouble> x2(nq);

	        						                 // get the coordinates (assuming all fields have the same
	        						                 // discretisation)
	        						                 m_fields[m_logc[0]]->GetCoords(x0,x1,x2);

	       				for(int i=0;i<nq;i++)
	        				{
	        					if(x1[i]==0.5)
	        					{
	        						for(int j=0;j<StressDim;j++)
	        						{
	        						stressoutarray[j][i]=0.0;
	        						}
	        						/*cout << "i:" << i;
	        						cout << " y=" << x1[i] << "\tLogcxxin: " << outarray[m_logc[0]][i];
	        						cout << "\tLogcxyin: " << outarray[m_logc[1]][i];
	        						cout << "\tLogcyyin: " << outarray[m_logc[2]][i];
	        						cout << endl;
	        					}
	        				} */
#ifdef OUTPUT
	        cout << "\toutarraybeforemultiply by inverse: " << stressoutarray[0][k] << endl;
#endif
//TODO:CHANGE THIS BACK
	        for(i = 0; i < StressDim; ++i)
	        {
	        	// m_fields[m_logc[i]]->IProductWRTBase(stressoutarray[i],Weakoutarray[i]);
	           m_fields[m_logc[i]]->IProductWRTBase(stressoutarray[i],iprod);
	           Vmath::Vadd(ncoeffs,iprod,1,Weakoutarray[i],1,Weakoutarray[i],1);
	        }
#ifdef OUTPUT
	        cout << "\tWeakoutarray of RARt: " << Weakoutarray[0][k] << endl;
#endif
	     /*   for(i = 0; i < StressDim; ++i)
			{
				 m_fields[m_logc[i]]->MultiplyByElmtInvMass(Weakoutarray[i],
													Weakoutarray[i]);
				 m_fields[m_logc[i]]->BwdTrans(Weakoutarray[i],stressoutarray[i]);
				 Vmath::Vadd(phystot,stressoutarray[i],1,outarray[m_logc[i]],1,outarray[m_logc[i]],1);
				 //outarray[m_logc[i]] = stressoutarray[i]; //erases the influence of convection
			} */
#ifdef OUTPUT
	        cout << "\tstressoutarray of RARt: " << stressoutarray[0][k] << endl;
#endif
		}

		void LogConformation::CalcSourceTermLogc(const Array<OneD, const Array<OneD, NekDouble> > &inarray,
		                Array<OneD, Array<OneD, NekDouble> > &outarray)
		{
			int phystot = m_fields[m_logc[0]]->GetTotPoints();
			int nvariables = m_fields.num_elements();
			int StressDim  = m_logc.num_elements();
			int VelDim  = m_velocity.num_elements();
			NekDouble Diff=0.0;
			NekDouble E1,E2;
			int i;
			int k=15;

			NekDouble factor = 0.0;
			int nvar = VelDim*VelDim;

			 Array<OneD,NekDouble> x0(phystot);
			 Array<OneD,NekDouble> x1(phystot);
			 Array<OneD,NekDouble> x2(phystot);

			 // get the coordinates (assuming all fields have the same
			 // discretisation)
			 m_fields[m_logc[0]]->GetCoords(x0,x1,x2);

			 Array<OneD, Array<OneD, NekDouble> > S(StressDim);
			 Array<OneD, Array<OneD, NekDouble> > Ap(StressDim);
			 Array<OneD, Array<OneD, NekDouble> > EigVec(StressDim);
			 Array<OneD, Array<OneD, NekDouble> > EigValS(m_spacedim);
			 Array<OneD, Array<OneD, NekDouble> > EigValC(m_spacedim);
			 Array<OneD, Array<OneD, NekDouble> > F(m_spacedim);
			 Array<OneD, Array<OneD, NekDouble> > L(nvar);
			 Array<OneD, Array<OneD, NekDouble> > Lp(nvar);
			 Array<OneD, NekDouble> ftr(phystot,1.0);
			 Array<OneD, NekDouble> tmp(phystot,0.0);

			 Ap[0] = Array<OneD, NekDouble>(phystot*StressDim);
			 S[0] = Array<OneD, NekDouble>(phystot*StressDim);
			 EigVec[0] = Array<OneD, NekDouble>(phystot*StressDim);
			  for(i = 1; i < StressDim; ++i)
			  {
				  Ap[i] = Ap[i-1] + phystot;
				  S[i] = S[i-1] + phystot;

				  EigVec[i] = EigVec[i-1] + phystot;
			  }

			  EigValC[0] = Array<OneD, NekDouble>(phystot*m_spacedim);
			  EigValS[0] = Array<OneD, NekDouble>(phystot*m_spacedim);
			  F[0] = Array<OneD, NekDouble>(phystot*m_spacedim);
			  for(i = 1; i < m_spacedim; ++i)
			  {
				  EigValC[i] = EigValC[i-1] + phystot;
				  EigValS[i] = EigValS[i-1] + phystot;
				  F[i] = F[i-1] + phystot;
			  }

			 Lp[0] = Array<OneD, NekDouble>(phystot*nvar);
			 L[0] = Array<OneD, NekDouble>(phystot*nvar);
			  for(i = 1; i < nvar; ++i)
			  {
				  Lp[i] = Lp[i-1] + phystot;
				  L[i] = L[i-1] + phystot;
			  }

			// Calculate velocity gradient tensor components
			m_fields[0]->PhysDeriv(0,inarray[m_velocity[0]], L[0]); //dUdx
			m_fields[0]->PhysDeriv(1,inarray[m_velocity[0]], L[1] ); //dUdy
			m_fields[0]->PhysDeriv(0,inarray[m_velocity[1]], L[2]); //dVdx
			m_fields[0]->PhysDeriv(1,inarray[m_velocity[1]], L[3]); //dVdy
#ifdef OUTPUT
			cout << "\tDudx: " << L[0][k] << endl;
			cout << "\tDudy: " << L[1][k] << endl;
			cout << "\tDvdx: " << L[2][k] << endl;
			cout << "\tDvdy: " << L[3][k] << endl;
#endif
			//1. Compute spectral decomposition of S:
			for(i=0;i<StressDim;i++)
			{
				S[i]=inarray[m_logc[i]];
			}
#ifdef OUTPUT
			cout << "\tSxx: " << S[0][k] << endl;
			cout << "\tSxy: " << S[1][k] << endl;
			cout << "\tSyy: " << S[2][k] << endl;
#endif
			CompSpectralDecomp2D(phystot,S, EigValS,EigVec);
#ifdef OUTPUT

			cout << "\teigcvals1: " << EigValS[0][k] << endl;
			cout << "\teigcvals2: " << EigValS[1][k] << endl;

			cout << "\te1: " << EigVec[0][k] << endl;
			cout << "\te1e2: " << EigVec[1][k] << endl;
			cout << "\te2: " << EigVec[1][k] << endl;
#endif
			//2. Compute eigenvalues of C
		  for(i = 0; i < m_spacedim; ++i)
		  {
			  Vmath::Vexp(phystot,EigValS[i],1,EigValC[i],1);
		  }
#ifdef OUTPUT
		  cout << "\teigcvalc1: " << EigValC[0][k] << endl;
		  cout << "\teigcvalc2: " << EigValC[1][k] << endl;
#endif
		  //3. Transform L into principal frame
		  GlobalToPrincipalTrans(phystot, L, EigVec,Lp);

#ifdef OUTPUT
			cout << "\tDudxp: " << Lp[0][k] << endl;
			cout << "\tDudyp: " << Lp[1][k] << endl;
			cout << "\tDvdxp: " << Lp[2][k] << endl;
			cout << "\tDvdyp: " << Lp[3][k] << endl;
#endif
		  //4. Compute F
		  //compute f(trace)
		  if(m_eps)
		  {
			  Vmath::Vadd(phystot, EigValC[0] , 1, EigValC[1] ,1,  ftr, 1);
			  Vmath::Sadd(phystot, -3.0, ftr,1, ftr, 1);
			  Vmath::Smul(phystot, m_eps,ftr, 1,ftr, 1);
			  if(m_lptt)
			  {
				  Vmath::Sadd(phystot, 1.0, ftr,1, ftr, 1);
			  }
			  else
			  {
				  Vmath::Vexp(phystot, ftr, 1, ftr, 1);
			  }
		  }
		  //4. compute F
		  // Use L as temporary variable because its not used anymore
		 for(i=0;i<m_spacedim;i++)
		 {
		  Vmath::Sadd(phystot, -1.0, EigValC[i],1, L[i], 1);
		  Vmath::Vmul(phystot,ftr, 1, L[i],1,  F[i], 1);
		  if(m_alpha)
		  {
		  Vmath::Vmul(phystot,L[i], 1, L[i],1,  L[i], 1);
		  Vmath::Smul(phystot, m_alpha,L[i],1,L[i],1);
		  Vmath::Vadd(phystot, L[i], 1, F[i],1, F[i], 1);
		  }
		  Vmath::Smul(phystot,(-1.0/m_weissenberg),F[i],1,F[i],1);
		 }

#ifdef OUTPUT
		  cout << "\tf1: " << F[0][k] << endl;
		  cout << "\tf2: " << F[1][k] << endl;
#endif
		 //5. Compute Ap
		 Vmath::Vdiv(phystot, F[0], 1,EigValC[0],1,Ap[0], 1);
		 Vmath::Svtvp(phystot, 2.0, Lp[0],1, Ap[0], 1,Ap[0], 1);

		 Vmath::Vdiv(phystot, F[1], 1,EigValC[1],1,Ap[2], 1);
		 Vmath::Svtvp(phystot, 2.0, Lp[3],1, Ap[2], 1,Ap[2], 1);

		  for(i=0;i<phystot;i++)
		  {
			  Diff = EigValC[0][i]-EigValC[1][i];
			  //E1 = fabs(EigValC[0][i]);
			  //E2 = fabs(EigValC[1][i]);
			 // Diff = EigVal[0][i]-EigVal[1][i];
			 // Diff = E1 -E2;
			  Diff = fabs(Diff);
			 // cout << "Diff at (x,y)=\t" << x0[i] << "\t," << x1[i] << "\t=" << Diff << endl;
			  if(Diff<1e-10)
			  {
				  Ap[1][i]= Lp[1][i]+Lp[2][i];
#ifdef OUTPUT
				  cout << "EV degenerated at (x,y)=\t" << x0[i] << "\t," << x1[i] << "\tid=" << i << endl;
#endif
				  //cout << "Ap" << Ap[1][i]<<endl;
			  }
			  else
			  {
				  Ap[1][i]=((EigValS[0][i]-EigValS[1][i])/(EigValC[0][i]-EigValC[1][i]))*(EigValC[1][i]*Lp[1][i]+EigValC[0][i]*Lp[2][i]);
			  }
		  }

#ifdef OUTPUT
		  cout << "\tApxx: " << Ap[0][k] << endl;
		  cout << "\tApxy: " << Ap[1][k] << endl;
		  cout << "\tApyy: " << Ap[2][k] << endl;
#endif
		  //6. Transform Ap into the global frame and save the result in outarray
		  PrincipalToGlobalTrans(phystot, Ap, EigVec,outarray);

		  //TODO: Experiment calc logc from tau for critical values
	/*	  for(i=0;i<phystot;i++)
		  {
			  Diff = EigValC[0][i]-EigValC[1][i];
			  //E1 = fabs(EigValC[0][i]);
			  //E2 = fabs(EigValC[1][i]);
			 // Diff = EigVal[0][i]-EigVal[1][i];
			 // Diff = E1 -E2;
			  Diff = fabs(Diff);
			 // cout << "Diff at (x,y)=\t" << x0[i] << "\t," << x1[i] << "\t=" << Diff << endl;
			  if(Diff<1e-5)
			  {
				  CalcElasticToLogConformationTensorforCriticalPoints(i, S);
				  for(int j = 0 ; j < m_stress.num_elements(); j++)
				{
				  outarray[j][i] = S[j][i];
				  cout << "outarray from tau" << outarray[j][i] << endl;
				}
				  cout << "EV degenerated at (x,y)=\t" << x0[i] << "\t," << x1[i] << "\tid=" << i << endl;

				  //cout << "Ap" << Ap[1][i]<<endl;
			  }
			  else
			  {

			  }
		  } */

/*		  int nq = m_fields[0]->GetNpoints();

		                 Array<OneD,NekDouble> x0(nq);
		                 Array<OneD,NekDouble> x1(nq);
		                 Array<OneD,NekDouble> x2(nq);

		                 // get the coordinates (assuming all fields have the same
		                 // discretisation)
		                 m_fields[0]->GetCoords(x0,x1,x2); */

		  //7. Use S as temporary variable to transform C into global frame
	/*	  S[0] = EigValC[0];
		  Vmath::Zero(phystot,S[1],1);
		  S[2] = EigValC[1];

		  PrincipalToGlobalTrans(phystot, S, EigVec,C); */

	/* For weak discontinuous version outarray is just for storing the stress values to cast it into the weak form
	 * and not the outarray directly
	 */
		}

	  void LogConformation::LogcToElasticL2Projection(const Array<OneD, const Array<OneD, NekDouble> > &inarray,
			  Array<OneD, Array<OneD, NekDouble> > &outarray)
	  {
		  int StressDim=m_stress.num_elements();
		  int phystot = m_fields[m_logc[0]]->GetTotPoints();
		  int ncoeffs = m_fields[m_logc[0]]->GetNcoeffs();
		  int i;

		  Array<OneD, Array<OneD, NekDouble> > Weakoutarray(StressDim);
		  Array<OneD, Array<OneD, NekDouble> > Logc(StressDim);
		  Array<OneD, Array<OneD, NekDouble> > Elastic(StressDim);

		  Weakoutarray[0] = Array<OneD, NekDouble>(ncoeffs*StressDim);
		  for(i = 1; i < StressDim; ++i)
		  {
			  Weakoutarray[i] = Weakoutarray[i-1] + ncoeffs;
		  }

		  Logc[0] = Array<OneD, NekDouble>(phystot*StressDim);
		  Elastic[0] = Array<OneD, NekDouble>(phystot*StressDim);
		  for(i = 1; i < StressDim; ++i)
		  {
			  Logc[i] = Logc[i-1] + phystot;
			  Elastic[i] = Elastic[i-1] + phystot;
		  }

		  for(i = 0; i < StressDim; ++i)
		  {
			  Logc[i]=outarray[m_logc[i]];
		  }

		  LogConformationTensorToElastic(phystot,Logc, Elastic);


		  /* TODO: Do I need L2 projection yes or no?! */
		 /* for(i = 0; i < StressDim; ++i)
		  {
			  m_fields[m_stress[i]]->UpdatePhys() = Elastic[i];
		  } */

		  /* This gives logc in terms of physical values */
		  for(i = 0; i < StressDim; ++i)
		  {
			  // m_fields[m_stress[i]]->IProductWRTBase(Elastic[i],Weakoutarray[i]);
				/* Calculate d\hat{txx}/dt
				 * Coefficients of the RHS
				 */
			 // m_fields[m_stress[i]]->MultiplyByElmtInvMass(Weakoutarray[i],Weakoutarray[i]);
			  /* Calculate dtxx/dt physical values at quadrature points
			   * physical values of RHS
			   */
			  //m_fields[m_stress[i]]->BwdTrans(Weakoutarray[i],outarray[m_stress[i]]);
      		m_fields[m_stress[i]]->FwdTrans(Elastic[i],m_fields[m_stress[i]]->UpdateCoeffs());
      		m_fields[m_stress[i]]->BwdTrans(m_fields[m_stress[i]]->GetCoeffs(),m_fields[m_stress[i]]->UpdatePhys());
		  }
	  }

      void LogConformation::AddStressDivergenceLogc(const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                                                     Array<OneD, Array<OneD, NekDouble> > &outarray)
        {
            int i,j;
            int nvariables = inarray.num_elements();
            int StressDim  = m_stress.num_elements();
            int VelDim     = m_velocity.num_elements();
            int physTot = m_fields[0]->GetTotPoints();
            Array<OneD, NekDouble> wk(physTot);

            // divergence for first velocity component
            Array<OneD, NekDouble> div0(2*physTot, 0.0);

            // divergence for second velocity component v
            Array<OneD, NekDouble> div1;

            div1 = div0 + physTot;

           // Array<OneD, Array<OneD, NekDouble> > stress(StressDim);

//            for(i = 0; i < StressDim; ++i)
  //          {
    //            stress[i] = inarray[m_stress[i]];
     //       }

            // Calculate Divergence of stress tensor
            if(m_spacedim==2)
            {
                // PhysDeriv calculates the derivative of the stress tensor, values of the input array
                // must be in physical space
//inarray[m_stress[i]] or m_fields[m_stress[0]]->GetPhys()
                   // wk= dtxx/dx
                       m_fields[m_stress[0]]->PhysDeriv(0,m_fields[m_stress[0]]->GetPhys() , wk);
                       // div0 = wk + div0
                       Vmath::Vadd(physTot,wk,1,div0,1,div0,1);
                   // wk = dtxy/dy
                       m_fields[m_stress[1]]->PhysDeriv(1,m_fields[m_stress[1]]->GetPhys() , wk);
                   // div0 = dtxx/dx + dtxy/dy
                       Vmath::Vadd(physTot,wk,1,div0,1,div0,1);


                   // div1 = dtxy/dx
                       m_fields[m_stress[1]]->PhysDeriv(0,m_fields[m_stress[1]]->GetPhys() , wk);
                       Vmath::Vadd(physTot,wk,1,div1,1,div1,1);
                   // div1 = dtxy/dx + dtyy/dy
                       m_fields[m_stress[2]]->PhysDeriv(1,m_fields[m_stress[2]]->GetPhys() , wk);
                       Vmath::Vadd(physTot,wk,1,div1,1,div1,1);


                   // multiply m_kinvis with divergence of stress
                     //  Vmath::Smul(physTot,m_kinvis,div0,1,div0,1);
                     //  Vmath::Smul(physTot,m_kinvis,div1,1,div1,1);


                   // Add result (kinvis*div(stress)) to the velocity components of the outarray
                   Vmath::Vadd(physTot,div0,1,outarray[m_velocity[0]],1,outarray[m_velocity[0]],1);
                   Vmath::Vadd(physTot,div1,1,outarray[m_velocity[1]],1,outarray[m_velocity[1]],1);

            }
            else
            {
            	ASSERTL0(false,"Calculation of Divergence of Stress is not yet implemented in 3D");
            }

        }

      void LogConformation::AddStressDivergenceWeakFormLogc(const Array<OneD, const Array<OneD, NekDouble> > &inarray,
													  Array<OneD, Array<OneD, NekDouble> > &forcing)
        {
        	int  nCoeffs = m_fields[0]->GetNcoeffs();
        	 int physTot = m_fields[m_logc[0]]->GetTotPoints();
        	 int StressDim = m_logc.num_elements();
        	 int VelDim = m_velocity.num_elements();
        	Array<OneD, Array<OneD, NekDouble> > Weakoutarray(VelDim);
        	Weakoutarray[0] = Array<OneD, NekDouble>(nCoeffs*VelDim);
        	for(int i = 1; i < VelDim; ++i)
        	{
        		Weakoutarray[i] = Weakoutarray[i-1] + nCoeffs;
        	}

        	Array<OneD, Array<OneD, NekDouble> > stressin(StressDim);
        	stressin[0] = Array<OneD, NekDouble>(physTot*StressDim);
			for(int i = 1; i < StressDim; ++i)
			{
				stressin[i] = stressin[i-1] + physTot;
			}

			for(int i = 0; i < StressDim; ++i)
			{
				stressin[i] = m_fields[m_stress[i]]->GetPhys();
			}

        	// Compute Weakoutarray=(tau^n,grad(phi_u))
        	CalcStressDivergenceWeakForm(stressin,Weakoutarray);

        	// -(tau^n,grad(phi_u))
            for(int i = 0; i < VelDim; ++i)
            {
            	Vmath::Neg(nCoeffs,Weakoutarray[i],1);
            }

            // Weakoutarray = Weakoutarray + <tau^n*n,phi_u>
            AddStressTimesNormalToVelocityNeumannBC(stressin,Weakoutarray);

        	Vmath::Vadd(nCoeffs,forcing[0],1,Weakoutarray[0],1,forcing[0],1);
        	Vmath::Vadd(nCoeffs,forcing[1],1,Weakoutarray[1],1,forcing[1],1);

        }

      void LogConformation::SetInitialConditionsLogc(NekDouble initialtime)
       {
           std::string restartstr = "RESTART";

           cout << "Initial Conditions:" << endl;

           // Check for restart file.
        	   int npoints = m_fields[0]->GetNpoints();
        	   int nvariables = m_stress.num_elements();

               Array<OneD, Array<OneD, NekDouble> > Elastic(nvariables);
               Elastic[0] = Array<OneD, NekDouble>(npoints*nvariables);
               for(int i=1;i<m_stress.num_elements();i++)
               {
            	   Elastic[i] = Elastic[i-1] + npoints;
               }

               for(int i = 0 ; i < m_stress.num_elements(); i++)
			 {
				 for(int j = 0; j < npoints; j++)
				 {
					 Elastic[i][j] = (m_fields[m_stress[i]]->GetPhys())[j];
				 }
			 }

               ElasticToLogConformationTensor(npoints,Elastic,Elastic);

               for(int i = 0 ; i < m_logc.num_elements(); i++)
               {

                   for(int j = 0; j < npoints; j++)
                   {
                       (m_fields[m_logc[i]]->UpdatePhys())[j]
                               = Elastic[i][j];
                   }
                   m_fields[m_logc[i]]->SetPhysState(true);
                   m_fields[m_logc[i]]->FwdTrans(m_fields[m_logc[i]]->GetPhys(),
                                         m_fields[m_logc[i]]->UpdateCoeffs());
                   cout << "\tField "<< m_boundaryConditions->GetVariable(m_logc[i])
                        <<": " << m_boundaryConditions->GetVariable(m_stress[i]) << endl;
               }

       }

	  void LogConformation::CalcElasticToLogConformationTensorforCriticalPoints(int id, Array<OneD, Array<OneD, NekDouble> > &Logc)
	  {
		  int nvariables = m_stress.num_elements();
          Array<OneD, Array<OneD, NekDouble> > Elastic(nvariables);
          Elastic[0] = Array<OneD, NekDouble>(nvariables);
          for(int i=1;i<m_stress.num_elements();i++)
          {
       	   Elastic[i] = Elastic[i-1] + 1;
          }

          for(int i = 0 ; i < m_stress.num_elements(); i++)
          {
        	  Elastic[i][0] = (m_fields[m_stress[i]]->GetPhys())[id];
          }

          ElasticToLogConformationTensor(1,Elastic,Elastic);

          for(int i = 0 ; i < m_stress.num_elements(); i++)
          {
        	  Logc[i][id] = Elastic[i][0];
          }
	  }

} //end of namespace


