////////////////////////////////////////////////////////////////////////////////
//
// File Viscoelastic.cpp
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
// Description: Viscoelastic class definition built on
// Navier-Stokes class
//
///////////////////////////////////////////////////////////////////////////////

#include <ViscoelasticSolver/EquationSystems/Viscoelastic.h>

namespace Nektar
{
    Viscoelastic::Viscoelastic(const LibUtilities::SessionReaderSharedPtr &pSession):
            IncNavierStokes(pSession),
    m_alpha(0.0),
    m_xi(0.0),
    m_eps(0.0),
    m_lptt(false),
    m_gamma(1.0)
    {
    }

    void Viscoelastic::v_InitObject()
    {

		IncNavierStokes::v_InitObject();
		int i,j;
		int numfields = m_fields.num_elements();

		// standard values
		m_alpha = 0.0;
		m_xi = 0.0;
		m_eps = 0.0;
		m_lptt = false;
		m_gamma = 1.0;

		if(m_equationType==eUnsteadyViscoelastic || m_equationType==eViscoelasticDEVSSG || m_equationType==eViscoelasticDEVSSGALE || m_equationType==eViscoelasticLogc)
		{
		std::string velids[] = {"txx","txy","tyy"};

		// Set up stress components to point to the last three m_fields;
		// pressure has to be the last field
		// Be aware implementation is for 2D only so far!!!
		// Set up int value to identify which m_fields contain the stress tensor components
		int StressDim = 3;
		m_stress = Array<OneD,int>(StressDim);

		for(i = 0; i < StressDim; ++i)
		{
			for(j = 0; j < numfields; ++j)
			{
				std::string var = m_boundaryConditions->GetVariable(j);
				if(NoCaseStringCompare(velids[i],var) == 0)
				{
					m_stress[i] = j;
					break;
				}

				if(j == numfields)
				{
					std::string error = "Failed to find field: " + var;
					ASSERTL0(false,error.c_str());
				}
			}
		}

		if(m_equationType==eViscoelasticDEVSSG || m_equationType==eViscoelasticDEVSSGALE )
		{
			int GDim = m_spacedim*m_spacedim;
			/*std::string velids[] = {"G11","G12","G21","G22" };

			// Set up stress components to point to the last three m_fields;
			// pressure has to be the last field
			// Be aware implementation is for 2D only so far!!!
			// Set up int value to identify which m_fields contain the stress tensor components


			m_Gindex = Array<OneD,int>(GDim);

			for(i = 0; i < GDim; ++i)
			{
				for(j = 0; j < numfields; ++j)
				{
					std::string var = m_boundaryConditions->GetVariable(j);
					if(NoCaseStringCompare(velids[i],var) == 0)
					{
						m_Gindex[i] = j;
						break;
					}

					if(j == numfields)
					{
						std::string error = "Failed to find field: " + var;
						ASSERTL0(false,error.c_str());
					}
				}
			} */

			 m_G  = Array<OneD, MultiRegions::ExpListSharedPtr>(GDim);

				for(i = 0; i < GDim; ++i)
				{
					m_G[i] = MemoryManager<MultiRegions::ExpList2D>::AllocateSharedPtr(m_session,m_graph->GetExpansions("txx"));
				}

		  /*  m_G  = Array<OneD, MultiRegions::ExpListSharedPtr>(GDim);

			for(i = 0; i < GDim; ++i)
			{
				m_G[i] = m_fields[m_Gindex[i]];
			} */

				  // Read in shear thinning parameter xi for the PTT model

				m_session->LoadParameter("theta", m_theta, (1-m_beta));

		}

		int phystot = m_fields[m_stress[0]]->GetTotPoints();
		m_stressold   = Array<OneD, Array<OneD, NekDouble> >(StressDim);
		m_stressold[0] = Array<OneD, NekDouble>(phystot*StressDim,0.0);

		for(i = 1; i < StressDim; ++i)
		{
			m_stressold[i] = m_stressold[i-1] + phystot;
		}

		const std::string stressStr = m_session->GetSolverInfo("STRESSTREATMENT");
		 for(i = 0; i < (int) eStressTreatmentSize; ++i)
		 {
			 if(nocase_cmp(kStressTreatmentStr[i],stressStr) == 0 )
			 {
				 m_stressTreatment = (StressTreatment)i;
				 break;
			 }
		 }
		 ASSERTL0(i != eStressTreatmentSize,"StressTreatment not found in SOLVERINFO section");

		 m_session->LoadParameter("Weissenberg", m_weissenberg,1.0);
		 m_relaxationtime = m_weissenberg;

		// Read in anisotropy parameter alpha for the Giesekus model
		 m_session->LoadParameter("alpha", m_alpha,0.0);

		// Read in shear thinning parameter xi for the PTT model
		 m_session->LoadParameter("xi", m_xi,0.0);

		// Read in elongational parameter epsilon for the PTT model
		 m_session->LoadParameter("eps", m_eps,0.0);

		// Read in if linear PTT model is required otherwise exponential PTT is calculated
		 m_session->LoadParameter("LPTT", m_lptt,1.0);

		// Read in upwind parameter gamma=1.0 fully upwind
		 m_session->LoadParameter("gamma", m_gamma,0.0);
		}
    }

    void Viscoelastic::SetInitialConditionsStress(NekDouble initialtime)
    {

    	/* initialize tau^n=tau(t=0.0) */
    	for(int i=0;i<m_stress.num_elements();i++)
		{
			m_fields[m_stress[i]]->BwdTrans(m_fields[m_stress[i]]->GetCoeffs(),m_stressold[i]);
		}

/* Compute initial velocity gradient tensor */
// TODO: undo this
 if(m_equationType==eViscoelasticDEVSSG || m_equationType==eViscoelasticDEVSSGALE )
    	{
    		m_fields[m_velocity[0]]->PhysDeriv(0,m_fields[m_velocity[0]]->GetPhys(), m_G[0]->UpdatePhys());
    		m_fields[m_velocity[0]]->PhysDeriv(1,m_fields[m_velocity[0]]->GetPhys(), m_G[1]->UpdatePhys());
    		m_fields[m_velocity[1]]->PhysDeriv(0,m_fields[m_velocity[1]]->GetPhys(), m_G[2]->UpdatePhys());
    		m_fields[m_velocity[1]]->PhysDeriv(1,m_fields[m_velocity[1]]->GetPhys(), m_G[3]->UpdatePhys());
    	}

    	for(int i = 0; i < m_G.num_elements(); i++)
    	{
    		m_G[i]->FwdTrans_IterPerExp(m_G[i]->GetPhys(),
    										m_G[i]->UpdateCoeffs());
    	}

    }

    // Adds the divergence of the stress tensor to outarray momentum equation
        void Viscoelastic::AddStressDivergence(const Array<OneD, const Array<OneD, NekDouble> > &inarray,
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

                   // wk= dtxx/dx
                       m_fields[m_stress[0]]->PhysDeriv(0,inarray[m_stress[0]], wk);
                       // div0 = wk + div0
                       Vmath::Vadd(physTot,wk,1,div0,1,div0,1);
                   // wk = dtxy/dy
                       m_fields[m_stress[1]]->PhysDeriv(1,inarray[m_stress[1]], wk);
                   // div0 = dtxx/dx + dtxy/dy
                       Vmath::Vadd(physTot,wk,1,div0,1,div0,1);


                   // div1 = dtxy/dx
                       m_fields[m_stress[1]]->PhysDeriv(0,inarray[m_stress[1]], wk);
                       Vmath::Vadd(physTot,wk,1,div1,1,div1,1);
                   // div1 = dtxy/dx + dtyy/dy
                       m_fields[m_stress[2]]->PhysDeriv(1,inarray[m_stress[2]], wk);
                       Vmath::Vadd(physTot,wk,1,div1,1,div1,1);


                   // TODO: change this back! multiply m_kinvis with divergence of stress
                       Vmath::Smul(physTot,m_kinvis,div0,1,div0,1);
                       Vmath::Smul(physTot,m_kinvis,div1,1,div1,1);


                   // Add result (kinvis*div(stress)) to the velocity components of the outarray
                   Vmath::Vadd(physTot,div0,1,outarray[m_velocity[0]],1,outarray[m_velocity[0]],1);
                   Vmath::Vadd(physTot,div1,1,outarray[m_velocity[1]],1,outarray[m_velocity[1]],1);

            }
            else
            {
            	ASSERTL0(false,"Calculation of Divergence of Stress is not yet implemented in 3D");
            }

        }

        // Adds the divergence of the stress tensor to outarray momentum equation
 void Viscoelastic::SubGGTDivergence(const Array<OneD, const Array<OneD, NekDouble> > &inarray,
											  Array<OneD, Array<OneD, NekDouble> > &outarray)
 {
	 int i,j;
	 int nvariables = inarray.num_elements();
	 int StressDim  = m_stress.num_elements();
	 int VelDim     = m_velocity.num_elements();
	 int physTot = m_fields[0]->GetTotPoints();
	 Array<OneD, NekDouble> wk(physTot);
	 Array<OneD, NekDouble> tmp(physTot);

	 // divergence for first velocity component
	 Array<OneD, NekDouble> div0(2*physTot, 0.0);

	 // divergence for second velocity component v
	 Array<OneD, NekDouble> div1;

	 div1 = div0 + physTot;

	 // Calculate Divergence of stress tensor
	 if(m_spacedim==2)
	 {
		 // PhysDeriv calculates the derivative of the stress tensor, values of the input array
		 // must be in physical space

				Vmath::Smul(physTot,2.0,m_G[0]->GetPhys(),1,tmp,1);
			// wk= dtxx/dx
				m_fields[m_stress[0]]->PhysDeriv(0,tmp, wk);
				// div0 = wk + div0
				Vmath::Vadd(physTot,wk,1,div0,1,div0,1);

			    Vmath::Vadd(physTot,m_G[1]->GetPhys(),1,m_G[2]->GetPhys(),1,tmp,1);
			// wk = dtxy/dy
				m_fields[m_stress[1]]->PhysDeriv(1,tmp, wk);
			// div0 = dtxx/dx + dtxy/dy
				Vmath::Vadd(physTot,wk,1,div0,1,div0,1);

			// div1 = dtxy/dx
				m_fields[m_stress[1]]->PhysDeriv(0,tmp, wk);
				Vmath::Vadd(physTot,wk,1,div1,1,div1,1);

				Vmath::Smul(physTot,2.0,m_G[3]->GetPhys(),1,tmp,1);
			// div1 = dtxy/dx + dtyy/dy
				m_fields[m_stress[2]]->PhysDeriv(1,tmp, wk);
				Vmath::Vadd(physTot,wk,1,div1,1,div1,1);


			// TODO: change this back! multiply m_kinvis with divergence of stress
				Vmath::Smul(physTot,m_kinvis*(-m_theta),div0,1,div0,1);
				Vmath::Smul(physTot,m_kinvis*(-m_theta),div1,1,div1,1);

		  /*  	cout <<  "div(G+GT)L2-projection: " << endl;


					  cout << "Div0:" << ": Max :" << Vmath::Vmax(physTot,div0,1);
					  cout << " Div0: " << Vmath::Vmin(physTot,div0,1)<< "\t";
					  cout << "Div1:" << ": Max :" << Vmath::Vmax(physTot,div1,1);
					  cout << " Div1: " << Vmath::Vmin(physTot,div1,1)<< "\t";

					   cout << endl; */

			// Add result (kinvis*div(stress)) to the velocity components of the outarray
			Vmath::Vadd(physTot,outarray[m_velocity[0]],1,div0,1,outarray[m_velocity[0]],1);
			Vmath::Vadd(physTot,outarray[m_velocity[1]],1,div1,1,outarray[m_velocity[1]],1);

	 }
	 else
	 {
		ASSERTL0(false,"Calculation of Divergence of Stress is not yet implemented in 3D");
	 }

 }

 void Viscoelastic::SubGDivergence(const Array<OneD, const Array<OneD, NekDouble> > &inarray,
											  Array<OneD, Array<OneD, NekDouble> > &outarray)
 {
	 int i,j;
	 int nvariables = inarray.num_elements();
	 int StressDim  = m_stress.num_elements();
	 int VelDim     = m_velocity.num_elements();
	 int physTot = m_fields[0]->GetTotPoints();
	 Array<OneD, NekDouble> wk(physTot);
	 Array<OneD, NekDouble> tmp(physTot);

	 // divergence for first velocity component
	 Array<OneD, NekDouble> div0(2*physTot, 0.0);

	 // divergence for second velocity component v
	 Array<OneD, NekDouble> div1;

	 div1 = div0 + physTot;

	 // Calculate Divergence of stress tensor
	 if(m_spacedim==2)
	 {
		 // PhysDeriv calculates the derivative of the stress tensor, values of the input array
		 // must be in physical space

				//Vmath::Smul(physTot,2.0,m_G[0]->GetPhys(),1,tmp,1);
				tmp = m_G[0]->GetPhys();
			// wk= dtxx/dx
				m_fields[m_stress[0]]->PhysDeriv(0,tmp, wk);
				// div0 = wk + div0
				Vmath::Vadd(physTot,wk,1,div0,1,div0,1);

			  //  Vmath::Vadd(physTot,m_G[1]->GetPhys(),1,m_G[2]->GetPhys(),1,tmp,1);
				tmp = m_G[1]->GetPhys();
			// wk = dtxy/dy
				m_fields[m_stress[1]]->PhysDeriv(1,tmp, wk);
			// div0 = dtxx/dx + dtxy/dy
				Vmath::Vadd(physTot,wk,1,div0,1,div0,1);


				tmp = m_G[2]->GetPhys();
			// div1 = dtxy/dx
				m_fields[m_stress[1]]->PhysDeriv(0,tmp, wk);
				Vmath::Vadd(physTot,wk,1,div1,1,div1,1);

			 //   Vmath::Smul(physTot,2.0,m_G[3]->GetPhys(),1,tmp,1);
				tmp = m_G[3]->GetPhys();
			// div1 = dtxy/dx + dtyy/dy
				m_fields[m_stress[2]]->PhysDeriv(1,tmp, wk);
				Vmath::Vadd(physTot,wk,1,div1,1,div1,1);


			// TODO: change this back! multiply m_kinvis with divergence of stress
				Vmath::Smul(physTot,m_kinvis*(-m_theta),div0,1,div0,1);
				Vmath::Smul(physTot,m_kinvis*(-m_theta),div1,1,div1,1);

		  /*  	cout <<  "div(G+GT)L2-projection: " << endl;


					  cout << "Div0:" << ": Max :" << Vmath::Vmax(physTot,div0,1);
					  cout << " Div0: " << Vmath::Vmin(physTot,div0,1)<< "\t";
					  cout << "Div1:" << ": Max :" << Vmath::Vmax(physTot,div1,1);
					  cout << " Div1: " << Vmath::Vmin(physTot,div1,1)<< "\t";

					   cout << endl; */

			// Add result (kinvis*div(stress)) to the velocity components of the outarray
			Vmath::Vadd(physTot,outarray[m_velocity[0]],1,div0,1,outarray[m_velocity[0]],1);
			Vmath::Vadd(physTot,outarray[m_velocity[1]],1,div1,1,outarray[m_velocity[1]],1);

	 }
	 else
	 {
		ASSERTL0(false,"Calculation of Divergence of Stress is not yet implemented in 3D");
	 }

 }
/* Computes the viscoelastic constitutive equation for the elastic stress
 *
 * inarray: physical values at the quadrature points of the values from timestep n tau_xx^n
 * outarray: values at quadrature points
 */

void Viscoelastic::AddStressDivergenceWeakFormExplicit(const Array<OneD, const Array<OneD, NekDouble> > &inarray,
											  Array<OneD, Array<OneD, NekDouble> > &forcing)
 {
	static int ncalls = 1; // Number of time this function has been called.
	Array<OneD, NekDouble> tmp;
	int  n=0;
	int VelDim = m_velocity.num_elements();
	int  nint    = min(ncalls++,m_intSteps);
	int  nlevels = m_stressdivergenceweaku.num_elements();
	int  nCoeffs = m_fields[0]->GetNcoeffs();
	Array<OneD, Array<OneD, NekDouble> > Weakoutarray(VelDim);
	Weakoutarray[0] = Array<OneD, NekDouble>(nCoeffs*VelDim);
	for(int i = 1; i < VelDim; ++i)
	{
		Weakoutarray[i] = Weakoutarray[i-1] + nCoeffs;
	}

	// Reshuffle Storage vector
	tmp = m_stressdivergenceweaku[nlevels-1];
	for(n = nlevels-1; n > 0; --n)
	{
		m_stressdivergenceweaku[n] = m_stressdivergenceweaku[n-1];
	}
	m_stressdivergenceweaku[0] = tmp;

	tmp = m_stressdivergenceweakv[nlevels-1];
	for(n = nlevels-1; n > 0; --n)
	{
		m_stressdivergenceweakv[n] = m_stressdivergenceweakv[n-1];
	}
	m_stressdivergenceweakv[0] = tmp;

	//Calculate (tau,grad(phi))
    // tau^n is given in m_fields(m_stress[i])




//Maybe here if m_fields[m_stress[i]]->GetPhysState == false ??? -> FwdTrans to update physical space

	// Compute Weakoutarray=(tau^n,grad(phi_u))
	CalcStressDivergenceWeakForm(m_stressold,Weakoutarray);

	// -(tau^n,grad(phi_u))
    for(int i = 0; i < VelDim; ++i)
    {
    	Vmath::Neg(nCoeffs,Weakoutarray[i],1);
    }

    // Weakoutarray = Weakoutarray + <tau^n*n,phi_u>
    AddStressTimesNormalToVelocityNeumannBC(m_stressold,Weakoutarray);

    if(m_equationType==eViscoelasticDEVSSG || m_equationType==eViscoelasticDEVSSGALE)
    {
    	AddGDivergenceWeakForm(Weakoutarray);
    	AddGTimesNormalToVelocityNeumannBC(Weakoutarray);
    }

    m_stressdivergenceweaku[0] = Weakoutarray[0];
    m_stressdivergenceweakv[0] = Weakoutarray[1];

	// Extrapolate to n+1
	Vmath::Smul(nCoeffs,kStressDivergenceExtrapolation[nint-1][nint-1],
			m_stressdivergenceweaku[nint-1],1,m_stressdivergenceweaku[nlevels-1],1);
	Vmath::Smul(nCoeffs,kStressDivergenceExtrapolation[nint-1][nint-1],
			m_stressdivergenceweakv[nint-1],1,m_stressdivergenceweakv[nlevels-1],1);

	for(n = 0; n < nint-1; ++n)
	{
		Vmath::Svtvp(nCoeffs,kStressDivergenceExtrapolation[nint-1][n],
				m_stressdivergenceweaku[n],1,m_stressdivergenceweaku[nlevels-1],1,
				m_stressdivergenceweaku[nlevels-1],1);
		Vmath::Svtvp(nCoeffs,kStressDivergenceExtrapolation[nint-1][n],
						m_stressdivergenceweakv[n],1,m_stressdivergenceweakv[nlevels-1],1,
						m_stressdivergenceweakv[nlevels-1],1);
	}

	Vmath::Vadd(nCoeffs,forcing[0],1,m_stressdivergenceweaku[nlevels-1],1,forcing[0],1);
	Vmath::Vadd(nCoeffs,forcing[1],1,m_stressdivergenceweakv[nlevels-1],1,forcing[1],1);
 }

void Viscoelastic::CalcStressDivergenceWeakForm(const Array<OneD, const Array<OneD, NekDouble> > &stressin,
												  Array<OneD, Array<OneD, NekDouble> > &Weakoutarray)
 {
	int  nCoeffs = m_fields[0]->GetNcoeffs();
	Array<OneD, NekDouble> iprod(nCoeffs);

    m_fields[m_velocity[0]]->IProductWRTDerivBase(0, stressin[0], iprod);
    m_fields[m_velocity[0]]->IProductWRTDerivBase(1, stressin[1], Weakoutarray[0]);
    Vmath::Vadd(nCoeffs, iprod, 1, Weakoutarray[0], 1, Weakoutarray[0], 1);

    m_fields[m_velocity[1]]->IProductWRTDerivBase(0, stressin[1], iprod);
    m_fields[m_velocity[1]]->IProductWRTDerivBase(1, stressin[2], Weakoutarray[1]);
    Vmath::Vadd(nCoeffs, iprod, 1, Weakoutarray[1], 1, Weakoutarray[1], 1);
 }

void Viscoelastic::AddStressTimesNormalToVelocityNeumannBC(const Array<OneD, const Array<OneD, NekDouble> > &stressin,
											  Array<OneD, Array<OneD, NekDouble> > &Weakoutarray)
{

    int i,el,n,cnt, offset, phys_offset,nq;
    int edge, elmtid;
    Array<OneD, NekDouble> e_outarray;
    Array<OneD, const NekDouble> txx,txy,tyy,P;
   // Array<OneD, Array<OneD, StdRegions::StdExpansion1DSharedPtr> >
    //    elmtToTrace = m_traceMap->GetElmtToTrace();
    StdRegions::StdExpansionSharedPtr EdgeExp;

    Array<OneD, int> ElmtID,EdgeID;

    //Fills ElmtID and EdgeID with global ids of elements id number and edges id number of the boundary expansion
    m_fields[m_velocity[0]]->GetBoundaryToElmtMap(ElmtID,EdgeID);

    int maxpts = 0;
    int ncoeffs = m_fields[0]->GetNcoeffs();
    int phystot = m_fields[m_velocity[0]]->GetTotPoints();
    Array<OneD, NekDouble> pressurecoeff(ncoeffs,0.0);
    Array<OneD, NekDouble> pressurephys(phystot,0.0);

    // Project pressure onto velocity space
    m_pressure->BwdTrans(m_pressure->GetCoeffs(), m_pressure->UpdatePhys());

    // project onto velocity space without enforcing BCs
    m_fields[0]->FwdTrans_IterPerExp(m_pressure->GetPhys(),pressurecoeff);
    m_fields[0]->BwdTrans(pressurecoeff,pressurephys);


     // find the maximum values of points
     for(cnt = n = 0; n < m_fields[m_velocity[0]]->GetBndConditions().num_elements(); ++n)
     {
 //Loop over all elements in boundary region
     	//cnt starts with element 0 then going through elements
         for(i = 0; i < m_fields[m_velocity[0]]->GetBndCondExpansions()[n]->GetExpSize(); ++i)
         {
             maxpts = max(maxpts, m_fields[m_velocity[0]]->GetExp(ElmtID[cnt++])->GetTotPoints());
         }
     }

     Array<OneD, NekDouble> txxedge(3*maxpts);
     Array<OneD, NekDouble> txyedge = txxedge + maxpts;
     Array<OneD, NekDouble> tyyedge = txyedge + maxpts;

    // Neumann BC for velocity component u
	for(cnt = n = 0; n < m_fields[m_velocity[0]]->GetBndConditions().num_elements(); ++n)
	{
		// Waters solution
		if((m_fields[m_velocity[0]]->GetBndConditions()[n])->GetBoundaryConditionType() == SpatialDomains::eNeumann)
		{
			// loop over elements along boundary
			for(el = 0; el < m_fields[m_velocity[0]]->GetBndCondExpansions()[n]->GetExpSize(); ++el,cnt++)
			{
				EdgeExp =  boost::dynamic_pointer_cast<StdRegions::StdExpansion1D> (m_fields[m_velocity[0]]->GetBndCondExpansions()[n]->GetExp(el));

                nq   = m_fields[m_velocity[0]]->GetExp(elmtid)->GetTotPoints();

                //cout << EdgeID.num_elements() << endl;
				// grab edge on the boundary
				edge = EdgeID[cnt];
				elmtid = ElmtID[cnt];

				offset = m_fields[m_velocity[0]]->GetPhys_Offset(elmtid);

				txx = stressin[0] + offset;
				txy = stressin[1] + offset;
				//P = (m_pressure->GetPhys())+ m_pressure->GetPhys_Offset(elmtid);
				P = pressurephys + offset;

				Array<OneD, NekDouble> cauchyxx(nq,0.0);

				Vmath::Vsub(nq,txx,1,P,1,cauchyxx,1);

				/*int nedge = EdgeExp->GetNumPoints(0);

				m_fields[m_velocity[0]]->GetExp(elmtid)->GetEdgePhysVals(edge,EdgeExp,P,tyyedge);

				for(int j=0;j<nedge;j++)
				{
					cout << "P=" << tyyedge[j] << endl;
				} */

				m_fields[m_velocity[0]]->GetExp(elmtid)->GetEdgePhysVals(edge,EdgeExp,cauchyxx,txxedge);
				m_fields[m_velocity[0]]->GetExp(elmtid)->GetEdgePhysVals(edge,EdgeExp,txy,txyedge);

				int nquad_e = EdgeExp->GetNumPoints(0);

				/* for(int j=0;j<nquad_e;j++)
				{
					cout << txxedge[j] << ", ";
				}
				cout << endl; */

				e_outarray = Weakoutarray[0] + m_fields[m_velocity[0]]->GetCoeff_Offset(elmtid);

				// in: physical values of txx and txy, out:weak values of forcing function
				m_fields[m_velocity[0]]->GetExp(elmtid)->AddEdgeNormBoundaryInt(edge,EdgeExp,txxedge,
			                                                txyedge,
			                                                e_outarray);
			}
		}

		cnt +=m_fields[m_velocity[0]]->GetBndCondExpansions()[n]->GetExpSize();
	}

	   // Neumann BC for velocity component u
		for(cnt = n = 0; n < m_fields[m_velocity[1]]->GetBndConditions().num_elements(); ++n)
		{
			// Waters solution
			if((m_fields[m_velocity[1]]->GetBndConditions()[n])->GetBoundaryConditionType() == SpatialDomains::eNeumann)
			{
				// loop over elements along boundary
				for(el = 0; el < m_fields[m_velocity[1]]->GetBndCondExpansions()[n]->GetExpSize(); ++el,cnt++)
				{
					EdgeExp =  boost::dynamic_pointer_cast<StdRegions::StdExpansion1D> (m_fields[m_velocity[1]]->GetBndCondExpansions()[n]->GetExp(el));

					nq   = m_fields[m_velocity[1]]->GetExp(elmtid)->GetTotPoints();
					// grab edge on the boundary
					edge = EdgeID[cnt];
					elmtid = ElmtID[cnt];

					offset = m_fields[m_velocity[1]]->GetPhys_Offset(elmtid);

					txy = stressin[1] + offset;
					tyy = stressin[2] + offset;
					//P = (m_pressure->GetPhys())+ m_pressure->GetPhys_Offset(elmtid);
					P = pressurephys + offset;

					Array<OneD, NekDouble> cauchyyy(nq,0.0);
					Vmath::Vsub(nq,txx,1,P,1,cauchyyy,1);

					m_fields[m_velocity[1]]->GetExp(elmtid)->GetEdgePhysVals(edge,EdgeExp,txy,txyedge);
					m_fields[m_velocity[1]]->GetExp(elmtid)->GetEdgePhysVals(edge,EdgeExp,cauchyyy,tyyedge);

					/*int nquad_e = EdgeExp->GetNumPoints(0);

					for(int j=0;j<nquad_e;j++)
					{
						cout << txyedge[j] << ", ";
					}
					cout << endl; */

					e_outarray = Weakoutarray[1] + m_fields[m_velocity[1]]->GetCoeff_Offset(elmtid);

					// in: physical values of txx and txy, out:weak values of forcing function
					m_fields[m_velocity[1]]->GetExp(elmtid)->AddEdgeNormBoundaryInt(edge,EdgeExp,txyedge,
				                                                tyyedge,
				                                                e_outarray);
				}
			}

			cnt +=m_fields[m_velocity[1]]->GetBndCondExpansions()[n]->GetExpSize();
		}

   /* for(n = 0; n < GetExpSize(); ++n)
    {
        offset = GetCoeff_Offset(n);
        for(e = 0; e < (*m_exp)[n]->GetNedges(); ++e)
        {
            t_offset = GetTrace()->GetPhys_Offset(elmtToTrace[n][e]->GetElmtId());

            (*m_exp)[n]->AddEdgeNormBoundaryInt(e,elmtToTrace[n][e],
                                                Fx + t_offset,
                                                Fy + t_offset,
                                                e_outarray = outarray+offset);
        }
    } */
}

void Viscoelastic::AddStressDivergenceWeakFormExplicit2(const Array<OneD, const Array<OneD, NekDouble> > &inarray,
											  Array<OneD, Array<OneD, NekDouble> > &forcing)
 {

	//Calculate (tau,grad(phi))
    // tau^n is given in m_fields(m_stress[i])
    int nCoeffs = m_fields[0]->GetNcoeffs();
    int StressDim = m_stress.num_elements();
    int phystot = m_fields[0]->GetTotPoints();

    Array<OneD, NekDouble> wk(nCoeffs);

    // divergence for first velocity component
    Array<OneD, NekDouble> div0(2*nCoeffs, 0.0);

    // divergence for second velocity component v
    Array<OneD, NekDouble> div1;

    div1 = div0 + nCoeffs;

    Array<OneD, Array<OneD, NekDouble> > stress(StressDim);
    stress[0] = Array<OneD, NekDouble>(phystot*StressDim,0.0);

	for(int i = 1; i < StressDim; ++i)
	{
		stress[i] = stress[i-1] + phystot;
	}

  /*  for(int i=0;i<m_stress.num_elements();i++)
    {
    	m_fields[m_stress[i]]->BwdTrans(m_fields[m_stress[i]]->GetCoeffs(),m_stressold[i]);
    } */

    // wk= dtxx/dx
		m_fields[0]->IProductWRTDerivBase(0, m_stressold[0], div0);
    // wk = dtxy/dy
		m_fields[0]->IProductWRTDerivBase(1, m_stressold[1], wk);
    // div0 = dtxx/dx + dtxy/dy
        Vmath::Vadd(nCoeffs,wk,1,div0,1,div0,1);

    // div1 = dtxy/dx
        m_fields[0]->IProductWRTDerivBase(0, m_stressold[1], div1);
    // div1 = dtxy/dx + dtyy/dy
        m_fields[0]->IProductWRTDerivBase(1, m_stressold[2], wk);
        Vmath::Vadd(nCoeffs,wk,1,div1,1,div1,1);

        //m_fields[m_stress[0]]->BwdTrans(div0,stress[0]);
        //m_fields[m_stress[0]]->BwdTrans(div1,stress[1]);
        //m_fields[m_stress[0]]->FwdTrans(stress[0],div0);
        //m_fields[m_stress[0]]->FwdTrans(stress[1],div1);

    // TODO: change this back! multiply m_kinvis with divergence of stress
      Vmath::Smul(nCoeffs,-1.0,div0,1,div0,1);
        Vmath::Smul(nCoeffs,-1.0,div1,1,div1,1);
	//Vmath::Vsub(nCoeffs,&(forcing[0])[0],1,&div0[0],1,&(forcing[0])[0],1);
	//Vmath::Vsub(nCoeffs,&(forcing[1])[0],1,&div1[0],1,&(forcing[1])[0],1);


	forcing[0] = div0;
	forcing[1] = div1;
 }

void Viscoelastic::AddStressDivergenceWeakFormExplicit3(const Array<OneD, const Array<OneD, NekDouble> > &inarray,
											  Array<OneD, Array<OneD, NekDouble> > &forcing)
 {

    int i,j;
    int StressDim  = m_stress.num_elements();
    int physTot = m_fields[0]->GetTotPoints();
    int nCoeffs = m_fields[0]->GetNcoeffs();
	static int ncalls = 1; // Number of time this function has been called.
	Array<OneD, NekDouble> tmp;
	int  n=0;
	int  nint    = min(ncalls++,m_intSteps);
	int  nlevels = m_stressdivergenceweaku.num_elements();

    Array<OneD, NekDouble> wk(physTot);

    // divergence for first velocity component
    Array<OneD, NekDouble> div0(2*physTot, 0.0);

    // divergence for second velocity component v
    Array<OneD, NekDouble> div1;

    div1 = div0 + physTot;

	// Reshuffle Storage vector
	tmp = m_stressdivergenceweaku[nlevels-1];
	for(n = nlevels-1; n > 0; --n)
	{
		m_stressdivergenceweaku[n] = m_stressdivergenceweaku[n-1];
	}
	m_stressdivergenceweaku[0] = tmp;

	tmp = m_stressdivergenceweakv[nlevels-1];
	for(n = nlevels-1; n > 0; --n)
	{
		m_stressdivergenceweakv[n] = m_stressdivergenceweakv[n-1];
	}
	m_stressdivergenceweakv[0] = tmp;

    if(m_spacedim==2)
    {
        // PhysDeriv calculates the derivative of the stress tensor, values of the input array
        // must be in physical space

        for(int i=0;i<m_stress.num_elements();i++)
        {
        	m_fields[m_stress[i]]->BwdTrans(m_fields[m_stress[i]]->GetCoeffs(),m_stressold[i]);
        }

           // wk= dtxx/dx
               m_fields[m_stress[0]]->PhysDeriv(0,m_stressold[0], wk);
               // div0 = wk + div0
               Vmath::Vadd(physTot,wk,1,div0,1,div0,1);
           // wk = dtxy/dy
               m_fields[m_stress[1]]->PhysDeriv(1,m_stressold[1], wk);
           // div0 = dtxx/dx + dtxy/dy
               Vmath::Vadd(physTot,wk,1,div0,1,div0,1);


           // div1 = dtxy/dx
               m_fields[m_stress[1]]->PhysDeriv(0,m_stressold[1], wk);
               Vmath::Vadd(physTot,wk,1,div1,1,div1,1);
           // div1 = dtxy/dx + dtyy/dy
               m_fields[m_stress[2]]->PhysDeriv(1,m_stressold[2], wk);
               Vmath::Vadd(physTot,wk,1,div1,1,div1,1);

               m_stressdivergenceweaku[0] = div0;
               m_stressdivergenceweakv[0] = div1;

           	// Extrapolate to n+1
           	Vmath::Smul(physTot,kStressDivergenceExtrapolation[nint-1][nint-1],
           			m_stressdivergenceweaku[nint-1],1,m_stressdivergenceweaku[nlevels-1],1);
           	Vmath::Smul(physTot,kStressDivergenceExtrapolation[nint-1][nint-1],
           			m_stressdivergenceweakv[nint-1],1,m_stressdivergenceweakv[nlevels-1],1);

           	for(n = 0; n < nint-1; ++n)
           	{
           		Vmath::Svtvp(physTot,kStressDivergenceExtrapolation[nint-1][n],
           				m_stressdivergenceweaku[n],1,m_stressdivergenceweaku[nlevels-1],1,
           				m_stressdivergenceweaku[nlevels-1],1);
           		Vmath::Svtvp(physTot,kStressDivergenceExtrapolation[nint-1][n],
           						m_stressdivergenceweakv[n],1,m_stressdivergenceweakv[nlevels-1],1,
           						m_stressdivergenceweakv[nlevels-1],1);
           	}

           	//Create weak form and save results in force function
               m_fields[m_stress[0]]->IProductWRTBase(m_stressdivergenceweaku[nlevels-1],forcing[0]);
               m_fields[m_stress[0]]->IProductWRTBase(m_stressdivergenceweakv[nlevels-1],forcing[1]);
    }

 }

        void Viscoelastic::EvaluateRhsViscoelasticEquation(const Array<OneD, const Array<OneD, NekDouble> > &inarray,
	    	                Array<OneD, Array<OneD, NekDouble> > &outarray)
	{
		if(m_projectionType == MultiRegions::eGalerkin)
		{
			int vstart=m_stress[0];
		    StressAdvection(inarray,outarray,vstart);
	        if(m_stressTreatment==eExplicit)
	        {
	 		   // Adding source terms to outarray
	 		   CalcSourceTerm(inarray,outarray);
	        }
	        else
	        {
	        	// source term is calculated implicitely
	        }

		}
		else
		{
			/* Follow SWE solver
			 * Create an Array to store the local expansion modes coeffsout
			 * Evaluate the weak advection (maybe inverse with mass matrix not sure yet)
			 * Evaluate the source term and cast it into weak form have a physfield for that purpose
			 * or use outarray as tmp save
			 * maybe then compute inverse of mass matrix and apply BwdTrans to calculate outarray
			 * (similar to UnsteadyAdvection solver)
			 */
			int vstart = m_stress[0];
			int ncoeffs = m_fields[m_stress[0]]->GetNcoeffs();
			int StressDim  = m_stress.num_elements();
			int nvariables = m_fields.num_elements();
			int i;


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

	        if(m_stressTreatment==eExplicit)
	        {
	          CalcWeakDGSourceTerm(inarray,Weakoutarray);
	        }
	        else
	        {
	        	// source term is calculated implicitely
	        }

			for(i = 0; i < StressDim; ++i)
			{
				/* Calculate d\hat{txx}/dt
					 * Coefficients of the RHS
					 */
				 m_fields[m_stress[i]]->MultiplyByElmtInvMass(Weakoutarray[i],
			                                        Weakoutarray[i]);
			     /* Calculate dtxx/dt physical values at quadrature points
			      * physical values of RHS
			      */
			     m_fields[m_stress[i]]->BwdTrans(Weakoutarray[i],outarray[m_stress[i]]);
			}

		}

	}

	void Viscoelastic::EvaluateViscoelasticEquationImplicitely(const Array<OneD, const Array<OneD, NekDouble> > &inarray,
		    	                Array<OneD, Array<OneD, NekDouble> > &outarray,const NekDouble aii_Dt)
	{
				/* Follow SWE solver
				 * Create an Array to store the local expansion modes coeffsout
				 * Evaluate the weak advection (maybe inverse with mass matrix not sure yet)
				 * Evaluate the source term and cast it into weak form have a physfield for that purpose
				 * or use outarray as tmp save
				 * maybe then compute inverse of mass matrix and apply BwdTrans to calculate outarray
				 * (similar to UnsteadyAdvection solver)
				 */
				int phystot = m_fields[m_stress[0]]->GetTotPoints();
				int ncoeffs = m_fields[m_stress[0]]->GetNcoeffs();
				int nvariables = m_fields.num_elements();
				int StressDim  = m_stress.num_elements();
				int VelDim  = m_velocity.num_elements();
				int InarrayDim;
			/*	if(m_equationType == eSteadyViscoelastic)
				{
					InarrayDim=m_nConvectiveFields + m_spacedim;
				}
				else
				{ */
					InarrayDim=m_nConvectiveFields;
			//	}

				int i;
				int it=1;
				static int outputcalls = 1;
				NekDouble maximum, max1;
				maximum = 1.0;
				Array<OneD, NekDouble> max(StressDim);
				for(i = 0; i<StressDim; i++)
							{
								max[i]=1.0;
							}

				Array<OneD, Array<OneD, NekDouble> > stressinarray(InarrayDim);
				Array<OneD, Array<OneD, NekDouble> > stressoutarray(StressDim);
		        Array<OneD, Array<OneD, NekDouble> > Weakoutarray(StressDim);
		        Array<OneD, NekDouble> Diff(phystot);

		        Weakoutarray[0] = Array<OneD, NekDouble>(ncoeffs*StressDim);
		        for(i = 1; i < StressDim; ++i)
		        {
		           Weakoutarray[i] = Weakoutarray[i-1] + ncoeffs;
		        }

		        stressoutarray[0] = Array<OneD, NekDouble>(phystot*StressDim,0.0);
		        for(i = 1; i < StressDim; ++i)
		        {
		        	stressoutarray[i] = stressoutarray[i-1] + phystot;
		        }

		        stressinarray[0] = Array<OneD, NekDouble>(phystot*InarrayDim);
				for(i = 1; i < InarrayDim; ++i)
				{
					stressinarray[i] = stressinarray[i-1] + phystot;
				}

				/* fill stressinarray with u^{n+1} and intermediate stress
				 * components
				 */
		/*		if(m_equationType==eSteadyViscoelastic)
				{
					for(i = 0; i<m_spacedim; i++)
					{
						m_fields[m_velocity[i]]->BwdTrans(m_fields[m_velocity[i]]->GetCoeffs(),stressinarray[m_velocity[i]]);
					}
				}
				else
				{ */
					for(i = 0; i<VelDim; i++)
					{
						stressinarray[m_velocity[i]]=outarray[m_velocity[i]];
					}
			//	}

				//TODO: after testing include this again. Less iterations are necessary with this approach
				for(i = 0; i<StressDim; i++)
				{
					stressinarray[m_stress[i]]=m_stressold[i];
					//stressinarray[m_stress[i]]=inarray[m_stress[i]];
				}

/*				cout << "Values stressinarray";
				for(i = 0; i < m_nConvectiveFields; ++i)
				{
				cout << GetVariable(i) << ": Max :" << Vmath::Vmax(phystot,stressinarray[i],1);
				cout << " Min: " << Vmath::Vmin(phystot,stressinarray[i],1) << "\t";
				}
				cout << endl; */
				Array<OneD, NekDouble> ftr(phystot,1.0);
				Array<OneD, NekDouble> trace(phystot,0.0);
				Array<OneD, NekDouble> tmp(phystot,0.0);


		        /*
		         * iterate over this function until approximation gives new time step values
		         * while(maxerror<0.01) do (it = 0 , it++)
		         * write error and iteration to screen
		         */

				//while(maximum>1e-8 || it<=5)
				while(maximum>1e-9)
				{
					CalcSourceTermImplicit(stressinarray,stressoutarray);

					/* Store results of the weak form (dimension is of N_eof number of all local expansions)
					 * in Weakoutarray
					 * and use stressinarray to calculate new stress values
					 */
				    // Calculate trace if PTT models are use in 2D: trace = tauxx+ tauyy
				    // Calculate f(tr) for the PTT models
			        if(m_eps)
			        {
			          Vmath::Vadd(phystot,stressinarray[m_stress[0]],1,stressinarray[m_stress[2]],1,trace,1);
			          // tmp = (m_eps*m_relaxationtime/(1.0-m_beta))*trace
			          Vmath::Smul(phystot,(m_eps*m_relaxationtime/(1.0-m_beta)),trace,1,tmp,1);
			          if(m_lptt)
			          {
			        	  // ftr = 1 + tmp
			        	  Vmath::Sadd(phystot,1.0,tmp,1,ftr,1);
			          }
			          else
			          {
			        	  // ftr = exp(tmp)
			        	  Vmath::Vexp(phystot, tmp, 1, ftr, 1);
			          }

			        }

					for(i = 0; i < StressDim; ++i)
					{
						/*
						* add intermediate stress field to stressoutarray
						*/
						// outarray += Wi*gamma_0/delta_t*tauxx


						/*if(m_eps)
						{
							// tyy -= 1/m_relaxtiontime*ftr*tauyy
							Vmath::Svvtvp(phystot,(-1.0/m_relaxationtime),ftr,1,inarray[m_stress[2]],1,tyy,1,tyy,1);
						}
						else
						{
							// tyy -= 1/m_relaxtiontime*tauyy
							Vmath::Svtvp(phystot,(-1.0/m_relaxationtime),inarray[m_stress[2]],1,tyy,1,tyy,1);
						} */

						/* outarray = source + gamma_0/deltat tauintermediate */
						Vmath::Svtvp(phystot,(1.0/aii_Dt),inarray[m_stress[i]],1,stressoutarray[i],1,stressoutarray[i],1);

						// tmp = ftr + gamma_0/deltat*Wi
						Vmath::Sadd(phystot,(m_weissenberg/aii_Dt),ftr,1,tmp,1);
						// tmp = Wi/tmp;
						Vmath::Sdiv(phystot,m_weissenberg, tmp,1, tmp, 1);
						Vmath::Vmul(phystot, tmp, 1, stressoutarray[i],1,  stressoutarray[i], 1);
					/* OLD:	Vmath::Svtvp(phystot,(m_weissenberg/aii_Dt),inarray[m_stress[i]],1,stressoutarray[i],1,stressoutarray[i],1);

						Vmath::Smul(phystot,1/(1+(m_weissenberg/aii_Dt)),stressoutarray[i],1,stressoutarray[i],1); */

/*						cout << "stressoutarray";

					    cout << GetVariable(m_stress[i]) << ": Max :" << Vmath::Vmax(phystot,stressoutarray[i],1);
						cout << " Min: " << Vmath::Vmin(phystot,stressoutarray[i],1) << "\t";
						cout << endl;*/
						/* Cast into weak form
						 *
						 */
					// m_fields[m_stress[i]]->IProductWRTBase(stressoutarray[i],Weakoutarray[i]);
						// Calculate d\hat{txx}/dt
						// Coefficients of the RHS

						// m_fields[m_stress[i]]->MultiplyByElmtInvMass(Weakoutarray[i],
						//									Weakoutarray[i]);

						//  Calculate dtxx/dt physical values at quadrature points
						//  physical values of RHS

						 //m_fields[m_stress[i]]->BwdTrans(Weakoutarray[i],stressoutarray[i]);

						 //m_fields[m_stress[i]]->SetPhys(stressoutarray[i]);
						 //m_fields[m_stress[i]]->SetPhysState(true);

						//Product WRT Base and MultiplyByElmtInvMass can be replaced by FwdTrans
						// Updates Coeffs field as well
						m_fields[m_stress[i]]->FwdTrans(stressoutarray[i],m_fields[m_stress[i]]->UpdateCoeffs());
						m_fields[m_stress[i]]->BwdTrans(m_fields[m_stress[i]]->GetCoeffs(),stressoutarray[i]);


						 /* Calculate maximum error */
						 Vmath::Vsub(phystot,stressoutarray[i],1,stressinarray[m_stress[i]],1,Diff,1);
						 max[i] = Vmath::Vmax(phystot,Diff,1);

						 //max[i] = m_fields[m_stress[i]]->L2(stressinarray[m_stress[i]]);

						 /* form inarray to contain physical values from old time step */
						 /* use the updated stress for input at next iteration */

						 stressinarray[m_stress[i]]=stressoutarray[i];
					}


					max1 = (max[0]<max[1]) ? max[1] : max[0];
					maximum = (max1 <max[2]) ? max[2] : max1;

					if(!((outputcalls+1)%m_infosteps))
				    {
						cout << "it:" << it << " max:" << maximum << endl;
				    }

					it++;
				}

				for(i = 0; i<StressDim; i++)
				{
					outarray[m_stress[i]]=stressoutarray[i];
				}

				outputcalls++;
	}

	// TODO: Formulate one function for each component and fill inarray that way otherwise function is too big
	void Viscoelastic::CalcSourceTerm(const Array<OneD, const Array<OneD, NekDouble> > &inarray,
	                Array<OneD, Array<OneD, NekDouble> > &outarray)
	{
		int phystot = m_fields[0]->GetTotPoints();
		int nvariables = m_fields.num_elements();
		int StressDim  = m_stress.num_elements();
		int VelDim  = m_velocity.num_elements();

		NekDouble factor = 0.0;
		int nvar = VelDim*VelDim;

		// cout << "Set stress boundary conditions with previous solution y" << endl;
		Array<OneD, NekDouble> txx(phystot,0.0);
		Array<OneD, NekDouble> txy(phystot,0.0);
		Array<OneD, NekDouble> tyy(phystot,0.0);
		Array<OneD, NekDouble> ftr(phystot,0.0);

		// number of arrays needed for calculations
		// in case of Giesekus or PTT models trace and tmp array is needed
		// if(m_alpha || m_eps || m_xi)
		// {
		  nvar +=  2;
		// }

		// VelDim*VelDim*phystot = Number of components of rate of deformation tensor
		Array<OneD, NekDouble> dUdx(nvar*phystot,0.0);
		Array<OneD, NekDouble> dUdy  = dUdx + phystot;
		Array<OneD, NekDouble> dVdx = dUdy  + phystot;
		Array<OneD, NekDouble> dVdy = dVdx + phystot;

	//	if(m_alpha || m_eps || m_xi)
	//	{
			Array<OneD, NekDouble> trace = dVdy + phystot;
			Array<OneD, NekDouble> tmp = trace + phystot;
//		}

		// Calculate velocity gradient tensor components
		m_fields[0]->PhysDeriv(0,inarray[m_velocity[0]], dUdx);
		m_fields[0]->PhysDeriv(1,inarray[m_velocity[0]], dUdy);
		m_fields[0]->PhysDeriv(0,inarray[m_velocity[1]], dVdx);
		m_fields[0]->PhysDeriv(1,inarray[m_velocity[1]], dVdy);

	    // Calculate trace if Giesekus or PTT models are use in 2D: trace = tauxx+ tauyy
	    if(m_alpha || m_eps)
	    {
	      Vmath::Vadd(phystot,inarray[m_stress[0]],1,inarray[m_stress[2]],1,trace,1);
	    }
	    // Calculate f(tr) for the PTT models
        if(m_eps)
        {
          // tmp = (m_eps*m_relaxationtime/(1.0-m_beta))*trace
          Vmath::Smul(phystot,(m_eps*m_relaxationtime/(1.0-m_beta)),trace,1,tmp,1);
          if(m_lptt)
          {
        	  // ftr = 1 + tmp
        	  Vmath::Sadd(phystot,1.0,tmp,1,ftr,1);
          }
          else
          {
        	  // ftr = exp(tmp)
        	  Vmath::Vexp(phystot, tmp, 1, ftr, 1);
          }

        }

		// Calculate txx
		// txx = 2*(1.0-m_beta)/m_relaxationtime * dUdx
		Vmath::Smul(phystot,(2*(1.0-m_beta)/m_relaxationtime),dUdx,1,txx,1);
		// txx += 2*(1-m_xi)*dUdx*tauxx
		Vmath::Svvtvp(phystot,2*(1-m_xi),dUdx,1,inarray[m_stress[0]],1,txx,1,txx,1);
		// txx += (2-m_xi)*dUdy*tauxy
		Vmath::Svvtvp(phystot,(2-m_xi),dUdy,1,inarray[m_stress[1]],1,txx,1,txx,1);

		// txx += -m_xi*dVdx*tauxy
		if(m_xi)
		{
			Vmath::Svvtvp(phystot,(-m_xi),dVdx,1,inarray[m_stress[1]],1,txx,1,txx,1);
		}
		// Vmath::Smul(phystot,2.0,txx,1,txx,1);

		if(m_eps)
		{
			// txx -= 1/m_relaxtiontime*ftr*tauxx
			Vmath::Svvtvp(phystot,(-1.0/m_relaxationtime),ftr,1,inarray[m_stress[0]],1,txx,1,txx,1);
		}
		else
		{
			// txx -= 1/m_relaxtiontime*tauxx
			Vmath::Svtvp(phystot,(-1.0/m_relaxationtime),inarray[m_stress[0]],1,txx,1,txx,1);
		}

		//Terms arising from the Giesekus model
		if(m_alpha)
		{
			// tmp = tauxx*tauxx+tauxy*tauxy
		    // z = v*w + x*y
		    Vmath::Vvtvvtp (phystot,inarray[m_stress[0]],1,inarray[m_stress[0]],1,inarray[m_stress[1]],1, inarray[m_stress[1]],1,tmp, 1);
		    // txx += -m_alpha/(1-m_beta)*tmp
		    Vmath::Svtvp(phystot,(-m_alpha/(1-m_beta)),tmp,1,txx,1,txx,1);

		}



		// Calculate txy Dtxy/Dt = RHS what we calculate here
		Vmath::Vadd(phystot,dUdy,1,dVdx,1,txy,1);
		Vmath::Smul(phystot,((1.0-m_beta)/m_relaxationtime),txy,1,txy,1);

		// those terms are zero because of div(u)=0
		// Vmath::Vvtvp(phystot,dUdx,1,inarray[m_stress[1]],1,txy,1,txy,1);
		// Vmath::Vvtvp(phystot,dVdy,1,inarray[m_stress[1]],1,txy,1,txy,1);

		Vmath::Svvtvp(phystot,(0.5*(2.0-m_xi)),dUdy,1,inarray[m_stress[2]],1,txy,1,txy,1);
		Vmath::Svvtvp(phystot,(0.5*(2.0-m_xi)),dVdx,1,inarray[m_stress[0]],1,txy,1,txy,1);

		// txx += -m_xi*dVdx*tauxy
			if(m_xi)
			{
				Vmath::Svvtvp(phystot,(-0.5*m_xi),dVdx,1,inarray[m_stress[2]],1,txy,1,txy,1);
				Vmath::Svvtvp(phystot,(-0.5*m_xi),dUdy,1,inarray[m_stress[0]],1,txy,1,txy,1);
			}
			// Vmath::Smul(phystot,2.0,txx,1,txx,1);

			if(m_eps)
			{
				// txx -= 1/m_relaxtiontime*ftr*tauxx
				Vmath::Svvtvp(phystot,(-1.0/m_relaxationtime),ftr,1,inarray[m_stress[1]],1,txy,1,txy,1);
			}
			else
			{
				// txx -= 1/m_relaxtiontime*tauxx
				Vmath::Svtvp(phystot,(-1.0/m_relaxationtime),inarray[m_stress[1]],1,txy,1,txy,1);
			}

			//Terms arising from the Giesekus model
			if(m_alpha)
			{
				// tmp = tauxx*tauxx+tauxy*tauxy
			    // z = v*w + x*y
			    Vmath::Vvtvvtp (phystot,inarray[m_stress[0]],1,inarray[m_stress[1]],1,inarray[m_stress[2]],1, inarray[m_stress[1]],1,tmp, 1);
			    // txx += -m_alpha/(1-m_beta)*tmp
			    Vmath::Svtvp(phystot,(-m_alpha/(1-m_beta)),tmp,1,txy,1,txy,1);

			}

		// Calculate tyy
		Vmath::Smul(phystot,(2*(1.0-m_beta)/m_relaxationtime),dVdy,1,tyy,1);

		Vmath::Svvtvp(phystot,2*(1-m_xi),dVdy,1,inarray[m_stress[2]],1,tyy,1,tyy,1);
		Vmath::Svvtvp(phystot,(2-m_xi),dVdx,1,inarray[m_stress[1]],1,tyy,1,tyy,1);


		// txx += -m_xi*dVdx*tauxy
		if(m_xi)
		{
			Vmath::Svvtvp(phystot,(-m_xi),dUdy,1,inarray[m_stress[1]],1,tyy,1,tyy,1);
		}
		// Vmath::Smul(phystot,2.0,txx,1,txx,1);

		if(m_eps)
		{
			// tyy -= 1/m_relaxtiontime*ftr*tauyy
			Vmath::Svvtvp(phystot,(-1.0/m_relaxationtime),ftr,1,inarray[m_stress[2]],1,tyy,1,tyy,1);
		}
		else
		{
			// tyy -= 1/m_relaxtiontime*tauyy
			Vmath::Svtvp(phystot,(-1.0/m_relaxationtime),inarray[m_stress[2]],1,tyy,1,tyy,1);
		}

		//Terms arising from the Giesekus model
		if(m_alpha)
		{
			// tmp = tauyy*tauyy+tauxy*tauxy
		    // z = v*w + x*y
		    Vmath::Vvtvvtp (phystot,inarray[m_stress[2]],1,inarray[m_stress[2]],1,inarray[m_stress[1]],1, inarray[m_stress[1]],1,tmp, 1);
		    // tyy += -m_alpha/(1-m_beta)*tmp
		    Vmath::Svtvp(phystot,(-m_alpha/(1-m_beta)),tmp,1,tyy,1,tyy,1);

		}

/* For weak discontinuous version outarray is just for storing the stress values to cast it into the weak form
 * and not the outarray directly
 */
		if(m_projectionType == MultiRegions::eGalerkin)
		{
			Vmath::Vadd(phystot,txx,1,outarray[m_stress[0]],1,outarray[m_stress[0]],1);
			Vmath::Vadd(phystot,txy,1,outarray[m_stress[1]],1,outarray[m_stress[1]],1);
			Vmath::Vadd(phystot,tyy,1,outarray[m_stress[2]],1,outarray[m_stress[2]],1);
		}
		else
		{
			outarray[0] = txx;
			outarray[1] = txy;
			outarray[2] = tyy;
		}
	}

	/* Call it 2 for test purposes,remove old version CalcImplicit later */
	void Viscoelastic::CalcSourceTermImplicit(const Array<OneD, const Array<OneD, NekDouble> > &inarray,
		                Array<OneD, Array<OneD, NekDouble> > &outarray)
		{
			int phystot = m_fields[0]->GetTotPoints();
			int nvariables = m_fields.num_elements();
			int StressDim  = m_stress.num_elements();
			int VelDim  = m_velocity.num_elements();

			NekDouble factor = 0.0;
			int nvar = VelDim*VelDim;

			// cout << "Set stress boundary conditions with previous solution y" << endl;
			Array<OneD, NekDouble> txx(phystot,0.0);
			Array<OneD, NekDouble> txy(phystot,0.0);
			Array<OneD, NekDouble> tyy(phystot,0.0);


			// number of arrays needed for calculations
			// in case of Giesekus or PTT models trace and tmp array is needed
			// if(m_alpha || m_eps || m_xi)
			// {
			//  nvar +=  2;
			// }
			nvar += 1;

			// VelDim*VelDim*phystot = Number of components of rate of deformation tensor
			Array<OneD, NekDouble> dUdx(nvar*phystot,0.0);
			Array<OneD, NekDouble> dUdy  = dUdx + phystot;
			Array<OneD, NekDouble> dVdx = dUdy  + phystot;
			Array<OneD, NekDouble> dVdy = dVdx + phystot;

			Array<OneD, NekDouble> tmp = dVdy + phystot;

		//	if(m_alpha || m_eps || m_xi)
		//	{

	//		}

			// Calculate velocity gradient tensor components
			m_fields[0]->PhysDeriv(0,inarray[m_velocity[0]], dUdx);
			m_fields[0]->PhysDeriv(1,inarray[m_velocity[0]], dUdy);
			m_fields[0]->PhysDeriv(0,inarray[m_velocity[1]], dVdx);
			m_fields[0]->PhysDeriv(1,inarray[m_velocity[1]], dVdy);


			// Calculate txx
			// txx = 2*(1.0-m_beta)/m_relaxationtime * dUdx
			Vmath::Smul(phystot,(2*(1.0-m_beta)/m_relaxationtime),dUdx,1,txx,1);
			// txx += 2*(1-m_xi)*dUdx*tauxx
			Vmath::Svvtvp(phystot,2*(1-m_xi),dUdx,1,inarray[m_stress[0]],1,txx,1,txx,1);
			// txx += (2-m_xi)*dUdy*tauxy
			Vmath::Svvtvp(phystot,(2-m_xi),dUdy,1,inarray[m_stress[1]],1,txx,1,txx,1);

			// txx += -m_xi*dVdx*tauxy
			if(m_xi)
			{
				Vmath::Svvtvp(phystot,(-m_xi),dVdx,1,inarray[m_stress[1]],1,txx,1,txx,1);
			}

			//Terms arising from the Giesekus model
			if(m_alpha)
			{
				// tmp = tauxx*tauxx+tauxy*tauxy
			    // z = v*w + x*y
			    Vmath::Vvtvvtp (phystot,inarray[m_stress[0]],1,inarray[m_stress[0]],1,inarray[m_stress[1]],1, inarray[m_stress[1]],1,tmp, 1);
			    // txx += -m_alpha/(1-m_beta)*tmp
			    Vmath::Svtvp(phystot,(-m_alpha/(1-m_beta)),tmp,1,txx,1,txx,1);

			}



			// Calculate txy Dtxy/Dt = RHS what we calculate here
			Vmath::Vadd(phystot,dUdy,1,dVdx,1,txy,1);
			Vmath::Smul(phystot,((1.0-m_beta)/m_relaxationtime),txy,1,txy,1);

			// those terms are zero because of div(u)=0
			// Vmath::Vvtvp(phystot,dUdx,1,inarray[m_stress[1]],1,txy,1,txy,1);
			// Vmath::Vvtvp(phystot,dVdy,1,inarray[m_stress[1]],1,txy,1,txy,1);

			Vmath::Svvtvp(phystot,(0.5*(2.0-m_xi)),dUdy,1,inarray[m_stress[2]],1,txy,1,txy,1);
			Vmath::Svvtvp(phystot,(0.5*(2.0-m_xi)),dVdx,1,inarray[m_stress[0]],1,txy,1,txy,1);

			// txx += -m_xi*dVdx*tauxy
				if(m_xi)
				{
					Vmath::Svvtvp(phystot,(-0.5*m_xi),dVdx,1,inarray[m_stress[2]],1,txy,1,txy,1);
					Vmath::Svvtvp(phystot,(-0.5*m_xi),dUdy,1,inarray[m_stress[0]],1,txy,1,txy,1);
				}

				//Terms arising from the Giesekus model
				if(m_alpha)
				{
					// tmp = tauxx*tauxx+tauxy*tauxy
				    // z = v*w + x*y
				    Vmath::Vvtvvtp (phystot,inarray[m_stress[0]],1,inarray[m_stress[1]],1,inarray[m_stress[2]],1, inarray[m_stress[1]],1,tmp, 1);
				    // txx += -m_alpha/(1-m_beta)*tmp
				    Vmath::Svtvp(phystot,(-m_alpha/(1-m_beta)),tmp,1,txy,1,txy,1);

				}

			// Calculate tyy
			Vmath::Smul(phystot,(2*(1.0-m_beta)/m_relaxationtime),dVdy,1,tyy,1);

			Vmath::Svvtvp(phystot,2*(1-m_xi),dVdy,1,inarray[m_stress[2]],1,tyy,1,tyy,1);
			Vmath::Svvtvp(phystot,(2-m_xi),dVdx,1,inarray[m_stress[1]],1,tyy,1,tyy,1);


			// txx += -m_xi*dVdx*tauxy
			if(m_xi)
			{
				Vmath::Svvtvp(phystot,(-m_xi),dUdy,1,inarray[m_stress[1]],1,tyy,1,tyy,1);
			}
			// Vmath::Smul(phystot,2.0,txx,1,txx,1);



			//Terms arising from the Giesekus model
			if(m_alpha)
			{
				// tmp = tauyy*tauyy+tauxy*tauxy
			    // z = v*w + x*y
			    Vmath::Vvtvvtp (phystot,inarray[m_stress[2]],1,inarray[m_stress[2]],1,inarray[m_stress[1]],1, inarray[m_stress[1]],1,tmp, 1);
			    // tyy += -m_alpha/(1-m_beta)*tmp
			    Vmath::Svtvp(phystot,(-m_alpha/(1-m_beta)),tmp,1,tyy,1,tyy,1);

			}

	/* For weak discontinuous version outarray is just for storing the stress values to cast it into the weak form
	 * and not the outarray directly
	 */
				outarray[0] = txx;
				outarray[1] = txy;
				outarray[2] = tyy;
		}


	/*void Viscoelastic::StressAdvection(const Array<OneD, const Array<OneD, NekDouble> > &inarray,
	                Array<OneD, Array<OneD, NekDouble> > &outarray, int vstart)
	{
       int i,j;
	   int nvariables = inarray.num_elements();
	   int nqtot      = m_fields[0]->GetTotPoints();
	   int VelDim     = m_velocity.num_elements();
	   int StressDim  = m_stress.num_elements();
	   Array<OneD, Array<OneD, NekDouble> > velocity(VelDim);
	   Array<OneD, NekDouble > Deriv;

	   for(i = 0; i < VelDim; ++i)
	   {
		   velocity[i] = inarray[m_velocity[i]];
	   }

	   //m_ALE->ComputeVelocityMinusMeshVelocity(m_fields,velocity);
	  // ALE::ALE Ale;
	   //ComputeVelocityMinusMeshVelocity(velocity);

	   // Set up Derivative work space;
	   if(wk.num_elements())
	   {
		   ASSERTL0(wk.num_elements() >= nqtot*VelDim,"Workspace is not sufficient");
		   Deriv = wk;
	   }
	   else
	   {
		   Deriv = Array<OneD, NekDouble> (nqtot*VelDim);
	   }

	   m_advObject->DoAdvection(m_fields,m_nConvectiveFields,
								m_velocity,inarray,outarray,Deriv);

	   for(i = 0; i < StressDim; ++i)
					 {
						AdvectionNonConservativeForm(velocity,inarray[vstart+i],outarray[vstart+i],Deriv);
						Vmath::Neg(nqtot,outarray[vstart+i],1);
					 }
	} */

	void Viscoelastic::StressAdvection(const Array<OneD, const Array<OneD, NekDouble> > &inarray,
	                Array<OneD, Array<OneD, NekDouble> > &outarray, int vstart)
	{
        int i,j;
        int nvariables = inarray.num_elements();
        int nqtot      = m_fields[0]->GetTotPoints();
        int VelDim     = m_velocity.num_elements();
        int StressDim  = m_stress.num_elements();

        Array<OneD, Array<OneD, NekDouble> > velocity(VelDim);
        Array<OneD, NekDouble > Deriv;
	/*	if(m_equationType==eSteadyViscoelastic)
		{
			for(i = 0; i<m_spacedim; i++)
			{
				m_fields[m_velocity[i]]->BwdTrans(m_fields[m_velocity[i]]->GetCoeffs(),velocity[i]);
			}
		}
		else
		{ */
        for(i = 0; i < VelDim; ++i)
        {
            velocity[i] = inarray[m_velocity[i]];
        }
	//	}

        // Set up Derivative work space;
        Deriv = Array<OneD, NekDouble> (nqtot*VelDim);

//                for(i = 0; i < m_nConvectiveFields; ++i)
//                {
//                    AdvectionNonConservativeForm(velocity,inarray[i],outarray[i],Deriv);
//                    Vmath::Neg(nqtot,outarray[i],1);
//                }

                for(i = 0; i < StressDim; ++i)
                {
                	AdvectionNonConservativeForm(velocity,inarray[vstart+i],outarray[vstart+i],Deriv);
                	Vmath::Neg(nqtot,outarray[vstart+i],1);
                }

	}

	void Viscoelastic::StressWeakDGAdvection(const Array<OneD, const Array<OneD, NekDouble> > &inarray,
		                Array<OneD, Array<OneD, NekDouble> > &Weakoutarray,int vstart)
	{
		int ncoeffs = m_fields[vstart]->GetNcoeffs();
		int StressDim  = m_stress.num_elements();

        int i;
        int VelDim     = m_velocity.num_elements();

        Array<OneD, Array<OneD, NekDouble> > velocity(VelDim);
        Array<OneD, Array<OneD, NekDouble> > stressinarray(StressDim);


	/*	if(m_equationType==eSteadyViscoelastic)
		{
			for(i = 0; i<m_spacedim; i++)
			{
				m_fields[m_velocity[i]]->BwdTrans(m_fields[m_velocity[i]]->GetCoeffs(),velocity[i]);
			}
		}
		else
		{ */
	        for(i = 0; i < VelDim; ++i)
	        {
	            velocity[i] = inarray[m_velocity[i]];
	        }

	 /*       for(int i=0;i<m_velocity.num_elements();i++)
	        {
	        	m_fields[m_velocity[i]]->BwdTrans(m_fields[m_velocity[i]]->GetCoeffs(),velocity[i]);
	        } */
	//	}

        for(i = 0; i < StressDim; ++i)
        {
           stressinarray[i] = inarray[vstart+i];
        }


		WeakDGAdvection(velocity,stressinarray, Weakoutarray,true,true,StressDim,vstart);

		// Compute -N(tau)
		for(i = 0; i < StressDim; ++i)
		{
			Vmath::Neg(ncoeffs,Weakoutarray[i],1);
		}
	}
	void Viscoelastic::WeakDGAdvection( Array<OneD, Array<OneD, NekDouble> > &velocity,
	                      const Array<OneD, Array<OneD, NekDouble> >& InField,
	                      Array<OneD, Array<OneD, NekDouble> >& OutField,
	                      bool NumericalFluxIncludesNormal,
	                      bool InFieldIsInPhysSpace,
	                      int nvariables,int vstart)
	          {
	              int i;
	              int nVelDim         = velocity.num_elements();
	              int nPointsTot      = m_fields[vstart]->GetNpoints();
	              int ncoeffs         = m_fields[vstart]->GetNcoeffs();
	              int nTracePointsTot = m_fields[vstart]->GetTrace()->GetNpoints();

	              if (!nvariables)
	              {
	                  nvariables      = m_fields.num_elements();
	              }

	              Array<OneD, Array<OneD, NekDouble> > fluxvector(nVelDim);
	              Array<OneD, Array<OneD, NekDouble> > physfield (nvariables);
	              Array<OneD, Array<OneD, NekDouble> > physoutfield (nvariables);
	              Array<OneD, NekDouble> iprod(ncoeffs);

	              for(i = 0; i < nVelDim; ++i)
	              {
	                  fluxvector[i]    = Array<OneD, NekDouble>(nPointsTot);
	              }

	              // Get the variables in physical space
	              // already in physical space
	              if(InFieldIsInPhysSpace == true)
	              {
	                  for(i = 0; i < nvariables; ++i)
	                  {
	                      physfield[i] = InField[i];
	                  }
	              }
	              // otherwise do a backward transformation
	              else
	              {
	                  for(i = 0; i < nvariables; ++i)
	                  {
	                      // Could make this point to m_fields[i]->UpdatePhys();
	                      physfield[i] = Array<OneD, NekDouble>(nPointsTot);
	                      m_fields[vstart+i]->BwdTrans(InField[i],physfield[i]);
	                  }
	              }

	              // Get the advection part (without numerical flux)
	              for(i = 0; i < nvariables; ++i)
	              {
	                  // Get the ith component of the  flux vector in (physical space)
	                  GetFluxVector(velocity, i, physfield, fluxvector, vstart);

	                  // Calculate the i^th value of (\grad_i \phi, F)
	                  WeakAdvectionGreensDivergenceForm(fluxvector,OutField[i],vstart);
	              }

	  	/*		for(i = 0; i < nvariables; ++i)
	  			{
	  			    m_fields[vstart+i]->MultiplyByElmtInvMass(OutField[i],
	  			                                        iprod);

	  			  physoutfield[i] = Array<OneD, NekDouble>(nPointsTot);
	  			     m_fields[vstart+i]->BwdTrans(iprod,physoutfield[i]);
	  			}


	              string fname = "cylinderphysfield.txt";
	              PrintValueAtCylinder(physoutfield,fname); */

	              // Get the numerical flux and add to the modal coeffs
	              // if the NumericalFluxs function already includes the
	              // normal in the output
	              if (NumericalFluxIncludesNormal == true)
	              {
	                  Array<OneD, Array<OneD, NekDouble> > numflux   (nvariables);

	                  for(i = 0; i < nvariables; ++i)
	                  {
	                      numflux[i]   = Array<OneD, NekDouble>(nTracePointsTot);
	                  }

	                  // Evaluate numerical flux in physical space which may in
	                  // general couple all component of vectors
	                  NumericalFlux(velocity, physfield, numflux,vstart);

	  	           /*   fname = "cylindernumflux0.txt";
	  	              PrintValueAtCylinderTrace(numflux[0],fname); */
	                  // Evaulate  <\phi, \hat{F}\cdot n> - OutField[i]
	                  for(i = 0; i < nvariables; ++i)
	                  {
	                      Vmath::Neg(ncoeffs,OutField[i],1);
	                      m_fields[vstart+i]->AddTraceIntegral(numflux[i],OutField[i]);
	                      m_fields[vstart+i]->SetPhysState(false);
	                  }

	  	/*  			for(i = 0; i < nvariables; ++i)
	  	  			{
	  	  			    m_fields[vstart+i]->MultiplyByElmtInvMass(OutField[i],
	  	  			                                        iprod);
	  	  			     m_fields[vstart+i]->BwdTrans(iprod,physoutfield[i]);
	  	  			}


	  	              fname = "cylinderphysfieldplustrace.txt";
	  	              PrintValueAtCylinder(physoutfield,fname); */
	              }

	              // if the NumericalFlux function does not include the
	              // normal in the output
	              else
	              {
	                  Array<OneD, Array<OneD, NekDouble> > numfluxX   (nvariables);
	                  Array<OneD, Array<OneD, NekDouble> > numfluxY   (nvariables);

	                  for(i = 0; i < nvariables; ++i)
	                  {
	                      numfluxX[i]   = Array<OneD, NekDouble>(nTracePointsTot);
	                      numfluxY[i]   = Array<OneD, NekDouble>(nTracePointsTot);
	                  }

	                  // Evaluate numerical flux in physical space which may in
	                  // general couple all component of vectors
	                  // UNDO: NumericalFlux(velocity, physfield, numfluxX, numfluxY);

	                  // Evaulate  <\phi, \hat{F}\cdot n> - OutField[i]
	                  for(i = 0; i < nvariables; ++i)
	                  {
	                      Vmath::Neg(ncoeffs,OutField[i],1);
	                      m_fields[vstart+i]->AddTraceIntegral(numfluxX[i], numfluxY[i],
	                                                    OutField[i]);
	                      m_fields[vstart+i]->SetPhysState(false);
	                  }
	              }
	          }

	void Viscoelastic::CalcWeakDGSourceTerm(const Array<OneD, const Array<OneD, NekDouble> > &inarray,
            Array<OneD, Array<OneD, NekDouble> > &Weakoutarray)
	{
		int phystot = m_fields[0]->GetTotPoints();
		int ncoeffs = m_fields[0]->GetNcoeffs();
		int nvariables = m_fields.num_elements();
		int StressDim  = m_stress.num_elements();
		int i;

		Array<OneD, Array<OneD, NekDouble> > stressoutarray(StressDim);
        Array<OneD, NekDouble> iprod(ncoeffs);
        // Vmath::Zero(nCoeffs, outarray, 1);

        stressoutarray[0] = Array<OneD, NekDouble>(phystot*StressDim);
        for(i = 1; i < StressDim; ++i)
        {
        	stressoutarray[i] = stressoutarray[i-1] + phystot;
        }

        CalcSourceTerm(inarray,stressoutarray);

        for(i = 0; i < StressDim; ++i)
        {
            m_fields[m_stress[i]]->IProductWRTBase(stressoutarray[i],iprod);
            Vmath::Vadd(ncoeffs,iprod,1,Weakoutarray[i],1,Weakoutarray[i],1);
        }
	}

    void Viscoelastic::GetFluxVector(Array<OneD, Array<OneD, NekDouble> > &velocity, const int i, Array<OneD, Array<OneD, NekDouble> > &physfield,
                           Array<OneD, Array<OneD, NekDouble> > &flux, int vstart)
    {
        ASSERTL1(flux.num_elements() == velocity.num_elements(),"Dimension of flux array and velocity array do not match");

        for(int j = 0; j < flux.num_elements(); ++j)
        {
            Vmath::Vmul(m_fields[vstart]->GetNpoints(),physfield[i],1,velocity[j],1,flux[j],1);
        }
    }

    void Viscoelastic::NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &velocity, Array<OneD, Array<OneD, NekDouble> > &physfield, Array<OneD, Array<OneD, NekDouble> > &numflux, int vstart)
    {
         int i;
         int nTraceNumPoints = m_fields[vstart]->GetTrace()->GetNpoints();
         int nvel = m_spacedim; //m_velocity.num_elements();
         Array<OneD, Array<OneD, NekDouble> > traceNormals(m_spacedim);
         //traceNormals = Array<OneD, Array<OneD, NekDouble> >(m_spacedim);

         for(i = 0; i < m_spacedim; ++i)
         {
            traceNormals[i] = Array<OneD, NekDouble> (m_fields[vstart]->GetTrace()->GetNpoints());
         }

         m_fields[vstart]->GetTrace()->GetNormals(traceNormals);

         Array<OneD, NekDouble > Fwd(nTraceNumPoints);
         Array<OneD, NekDouble > Bwd(nTraceNumPoints);
         Array<OneD, NekDouble > Vn (nTraceNumPoints,0.0);
/*int nq = m_fields[vstart]->GetNpoints();

Array<OneD,NekDouble> x0(nq);
Array<OneD,NekDouble> x1(nq);
Array<OneD,NekDouble> x2(nq);
Array<OneD,NekDouble> tx0(nTraceNumPoints);
Array<OneD,NekDouble> tx1(nTraceNumPoints);
Array<OneD,NekDouble> tx2(nTraceNumPoints);

// get the coordinates (assuming all fields have the same
// discretisation)
m_fields[vstart]->GetCoords(x0,x1,x2);

m_fields[vstart]->GetTrace()->GetCoords(tx0,tx1,tx2); */

         // Get Edge Velocity - Could be stored if time independent
         for(i = 0; i < nvel; ++i)
         {
        	 m_fields[vstart]->ExtractTracePhys(velocity[i], Fwd);

             Vmath::Vvtvp(nTraceNumPoints,traceNormals[i],1,Fwd,1,Vn,1,Vn,1);


         }



         for(i = 0; i < numflux.num_elements(); ++i)
         {
             m_fields[vstart+i]->GetFwdBwdTracePhys(physfield[i],Fwd,Bwd);
             //evaluate upwinded m_fields[i]
             // if Vn >= 0, flux = uFwd, i.e.,
             //  edge::eForward, if V*n>=0 <=> V*n_F>=0, pick uflux = uFwd
             //  edge::eBackward, if V*n>=0 <=> V*n_B<0, pick uflux = uFwd

             // else if Vn < 0, flux = uBwd, i.e.,
             //  edge::eForward, if V*n<0 <=> V*n_F<0, pick uflux = uBwd
             //  edge::eBackward, if V*n<0 <=> V*n_B>=0, pick uflux = uBwd

              //m_fields[vstart+i]->GetTrace()->Upwind(Vn,Fwd,Bwd,numflux[i]);

              for (int j = 0; j < nTraceNumPoints; ++j)
             {
            	if (Vn[j] >= 0.0)
            	{
            		numflux[i][j]  = m_gamma*Fwd[j] + (1.-m_gamma)*Bwd[j];
            	}
            	else
            	{
            		numflux[i][j]  = m_gamma*Bwd[j] + (1.-m_gamma)*Fwd[j];
            	}
             }
            /* for other values then gamma 1 */
              // Imposing weak boundary condition with flux
              // if Vn >= 0, uflux = uBwd at Neumann, i.e.,
              //  edge::eForward, if V*n>=0 <=> V*n_F>=0, pick uflux = uBwd
              //  edge::eBackward, if V*n>=0 <=> V*n_B<0, pick uflux = uBwd

              // if Vn >= 0, uflux = uFwd at Neumann, i.e.,
              //  edge::eForward, if V*n<0 <=> V*n_F<0, pick uflux = uFwd
              //  edge::eBackward, if V*n<0 <=> V*n_B>=0, pick uflux = uFwd

              if(m_fields[vstart+i]->GetBndCondExpansions().num_elements())
              {
                  WeakPenaltyforScalar(vstart+i,physfield[i],numflux[i],0.0);
              }

             // calculate m_fields[i]*Vn
              Vmath::Vmul(nTraceNumPoints,numflux[i],1,Vn,1,numflux[i],1);

             /*				for(int j=0;j<nTraceNumPoints;j++)
              				{ cout << GetVariable(vstart+i) << endl;
             					cout << "(x,y)= (" << tx0[j] << ", " << tx1[j] << ")\t";
              cout << "numflux" << j << ":" << numflux[i][j] <<endl;
              				} */
         }
     }
    void Viscoelastic::WeakAdvectionGreensDivergenceForm(
                const Array<OneD, Array<OneD, NekDouble> > &F,
                Array<OneD, NekDouble> &outarray, int vstart)
    {
        // use dimension of Velocity vector to dictate dimension of operation
        int ndim    = F.num_elements();
        int nCoeffs = m_fields[vstart]->GetNcoeffs();

        Array<OneD, NekDouble> iprod(nCoeffs);
        Vmath::Zero(nCoeffs, outarray, 1);

        for (int i = 0; i < ndim; ++i)
        {
            m_fields[vstart+i]->IProductWRTDerivBase(i, F[i], iprod);
            Vmath::Vadd(nCoeffs, iprod, 1, outarray, 1, outarray, 1);
        }
    }
    void Viscoelastic::WeakPenaltyforScalar(int var,
                 Array<OneD, NekDouble> &physfield,
                      Array<OneD, NekDouble> &penaltyflux,
                      NekDouble time)
    {
        int i, j, e, npoints, id1, id2;
        // Number of boundary regions
        int nbnd = m_fields[var]->GetBndCondExpansions().num_elements();
        int Nfps, numBDEdge;
        int nTraceNumPoints = m_fields[var]->GetTrace()->GetNpoints();
        int cnt = 0;

        Array<OneD, NekDouble > uplus(nTraceNumPoints);

        m_fields[var]->ExtractTracePhys(physfield,uplus);
        for(i = 0; i < nbnd; ++i)
        {
            // Number of boundary expansion related to that region
            numBDEdge = m_fields[var]->GetBndCondExpansions()[i]->GetExpSize();
            npoints = m_fields[var]->GetBndCondExpansions()[i]->GetNpoints();
            Array<OneD,NekDouble> x0(npoints,0.0);
            Array<OneD,NekDouble> x1(npoints,0.0);
            Array<OneD,NekDouble> x2(npoints,0.0);

            m_fields[var]->GetBndCondExpansions()[i]->GetCoords(x0,x1,x2);

            // Weakly impose boundary conditions by modifying flux values
            for (e = 0; e < numBDEdge ; ++e)
            {
                // Number of points on the expansion
                Nfps = m_fields[var]->GetBndCondExpansions()[i]->GetExp(e)->GetNumPoints(0) ;
                id1 = m_fields[var]->GetBndCondExpansions()[i]->GetPhys_Offset(e);
                id2 = m_fields[var]->GetTrace()->GetPhys_Offset(m_fields[var]->GetTraceMap()->GetBndCondTraceToGlobalTraceMap(cnt++));

                // For Dirichlet boundary condition: uflux = g_D given in the boundary expansion
                if(m_fields[var]->GetBndConditions()[i]->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                {
                	/*if(m_fields[var]->GetBndConditions()[i]->GetUserDefined().GetEquation() == "Wall")
                	{
                		Vmath::Vcopy(Nfps,&uplus[id2],1,&penaltyflux[id2],1);
                	}
                	else */
                	//{
                		Vmath::Vcopy(Nfps,&(m_fields[var]->GetBndCondExpansions()[i]->GetPhys())[id1],1,&penaltyflux[id2],1);
                	//}
                }
                // For Neumann boundary condition: uflux = u+
                else if((m_fields[var]->GetBndConditions()[i])->GetBoundaryConditionType() == SpatialDomains::eNeumann)
                {
                    Vmath::Vcopy(Nfps,&uplus[id2],1,&penaltyflux[id2],1);
                }
            }
        }
    }
	void Viscoelastic::SetBoundaryConditionsStress(const Array<OneD, const Array<OneD, NekDouble> > &physarray, NekDouble time)
	{
		int nvariables = m_fields.num_elements();
		int cnt = 0;

		// Set Boundary Condition for velocities
		for(int n = 0; n < m_fields[0]->GetBndConditions().num_elements(); ++n)
		{
			// Waters solution
			if (m_fields[0]->GetBndConditions()[n]->GetUserDefined() == SpatialDomains::eWaters)
			{

				SetBoundaryWatersU2(n,cnt,time);
				// forward transform to fill the modal coeffs
				m_fields[0]->GetBndCondExpansions()[n]->FwdTrans_BndConstrained(m_fields[0]->GetBndCondExpansions()[n]->GetPhys(),m_fields[0]->GetBndCondExpansions()[n]->UpdateCoeffs());
			}

			cnt +=m_fields[0]->GetBndCondExpansions()[n]->GetExpSize();
		}

		// reset counter to set bc for stress components
		cnt = 0;

		// loop over Boundary Regions for stress component tauxx
		for(int n = 0; n < m_fields[m_stress[0]]->GetBndConditions().num_elements(); ++n)
		{
			// cout << "n" << n << endl;
			// Time Dependent Boundary Condition (specified in meshfile)
			if (m_fields[m_stress[0]]->GetBndConditions()[n]->GetUserDefined() == SpatialDomains::eTimeDependent)
			{
				for (int i = 0; i < nvariables; ++i)
				{
					m_fields[i]->EvaluateBoundaryConditions(time);
				}
			}

			// Waters Inflow Boundary Condition
			if (m_fields[m_stress[0]]->GetBndConditions()[n]->GetUserDefined() == SpatialDomains::eWaters)
			{
				SetBoundaryWatersTAU(n,cnt,time);

				// forward transform to fill the modal coeffs
				m_fields[m_stress[0]]->GetBndCondExpansions()[n]->FwdTrans_BndConstrained(m_fields[m_stress[0]]->GetBndCondExpansions()[n]->GetPhys(),m_fields[m_stress[0]]->GetBndCondExpansions()[n]->UpdateCoeffs());
				m_fields[m_stress[1]]->GetBndCondExpansions()[n]->FwdTrans_BndConstrained(m_fields[m_stress[1]]->GetBndCondExpansions()[n]->GetPhys(),m_fields[m_stress[1]]->GetBndCondExpansions()[n]->UpdateCoeffs());
			}

			cnt +=m_fields[m_stress[0]]->GetBndCondExpansions()[n]->GetExpSize();
		}
	}

	void Viscoelastic::SetStressWallBoundaryConditions(const Array<OneD, const Array<OneD, NekDouble> > &inarray,const NekDouble aii_Dt)
	{

	}


	void Viscoelastic::Determinante(const NekDouble phystot, const Array<OneD, const Array<OneD, NekDouble> > &M, Array<OneD, NekDouble> &Det)
	{
		//int phystot = M[0].num_elements();
		//int phystot = m_fields[0]->GetTotPoints();
		Array<OneD,NekDouble> tmp(phystot);

		 Vmath::Vmul(phystot, M[6], 1, M[3],1, tmp, 1);
		 Vmath::Vmul(phystot, M[5], 1, M[4],1, Det, 1);
		 Vmath::Vsub(phystot, tmp, 1, Det,1, tmp, 1);
		 Vmath::Vmul(phystot, M[0], 1, tmp,1, Det, 1);

		 Vmath::Vmul(phystot, M[2], 1, M[6],1, tmp, 1);
		 Vmath::Vmul(phystot, M[1], 1, tmp,1, tmp, 1);
		 Vmath::Vsub(phystot, Det, 1, tmp,1, Det, 1);
	}

	void Viscoelastic::StressMatrixInverseOldB(const NekDouble phystot,const Array<OneD, const Array<OneD, NekDouble> > &M, Array<OneD, Array<OneD, NekDouble> > &Mout)
	{
		//int phystot = M[0].num_elements();
		//int phystot = m_fields[0]->GetTotPoints();
		Array<OneD,NekDouble> invDet(phystot);
		Array<OneD,NekDouble> tmp(phystot);

		Determinante(phystot, M, invDet);
		Vmath::Sdiv(phystot, 1.0, invDet, 1, invDet, 1);

		Vmath::Vmul(phystot, M[6], 1, M[3],1, tmp, 1);
		Vmath::Vmul(phystot, M[5], 1, M[4],1, Mout[0], 1);
		Vmath::Vsub(phystot, tmp, 1, Mout[0],1, Mout[0], 1);

		Vmath::Vmul(phystot, M[6], 1, M[1],1, Mout[1], 1);
		Vmath::Smul(phystot, -1.0, Mout[1], 1,Mout[1], 1);

		Vmath::Vmul(phystot, M[4], 1, M[1],1, Mout[2], 1);

		Vmath::Vmul(phystot, M[6], 1, M[2],1, Mout[3], 1);
		Vmath::Smul(phystot, -1.0, Mout[3], 1,Mout[3], 1);

		Vmath::Vmul(phystot, M[6], 1, M[0],1, Mout[4], 1);

		Vmath::Vmul(phystot, M[4], 1, M[0],1, Mout[5], 1);
		Vmath::Smul(phystot, -1.0, Mout[5], 1,Mout[5], 1);

		Vmath::Vmul(phystot, M[5], 1, M[2],1, Mout[6], 1);

		Vmath::Vmul(phystot, M[5], 1, M[0],1, Mout[7], 1);
		Vmath::Smul(phystot, -1.0, Mout[7], 1,Mout[7], 1);

		Vmath::Vmul(phystot, M[3], 1, M[0],1, tmp, 1);
		Vmath::Vmul(phystot, M[2], 1, M[1],1, Mout[8], 1);
		Vmath::Vsub(phystot, tmp, 1, Mout[8],1, Mout[8], 1);

		for(int i=0;i<9;i++)
		{
			Vmath::Vmul(phystot, invDet, 1,  Mout[i],1,  Mout[i], 1);
		}
	}

	void Viscoelastic::StressMatrixInversetimesRHS(const NekDouble phystot, Array<OneD, Array<OneD, NekDouble> > &M, Array<OneD, Array<OneD,NekDouble> > &RHS, Array<OneD, Array<OneD, NekDouble> > &tau)
	{
		//int phystot = M[0].num_elements();
		//int phystot = m_fields[0]->GetTotPoints();
		int matrixdim = 9;
		int i;
		Array<OneD,NekDouble> tmp(phystot);
		Array<OneD, Array<OneD, NekDouble> > Mout(matrixdim);

        Mout[0] = Array<OneD, NekDouble>(phystot*matrixdim);
        for(i = 1; i < matrixdim; ++i)
        {
        	Mout[i] = Mout[i-1] + phystot;
        }

        StressMatrixInverseOldB(phystot,M,Mout);

        Vmath::Vmul(phystot, Mout[0], 1, RHS[0],1, tmp, 1);
        Vmath::Vmul(phystot, Mout[1], 1, RHS[1],1, tau[0], 1);
        Vmath::Vadd(phystot, tmp, 1, tau[0],1, tau[0], 1);
        Vmath::Vmul(phystot, Mout[2], 1, RHS[2],1, tmp, 1);
        Vmath::Vadd(phystot, tmp, 1, tau[0],1, tau[0], 1);

        Vmath::Vmul(phystot, Mout[3], 1, RHS[0],1, tmp, 1);
        Vmath::Vmul(phystot, Mout[4], 1, RHS[1],1, tau[1], 1);
        Vmath::Vadd(phystot, tmp, 1, tau[1],1, tau[1], 1);
        Vmath::Vmul(phystot, Mout[5], 1, RHS[2],1, tmp, 1);
        Vmath::Vadd(phystot, tmp, 1, tau[1],1, tau[1], 1);

        Vmath::Vmul(phystot, Mout[6], 1, RHS[0],1, tmp, 1);
        Vmath::Vmul(phystot, Mout[7], 1, RHS[1],1, tau[2], 1);
        Vmath::Vadd(phystot, tmp, 1, tau[2],1, tau[2], 1);
        Vmath::Vmul(phystot, Mout[8], 1, RHS[2],1, tmp, 1);
        Vmath::Vadd(phystot, tmp, 1, tau[2],1, tau[2], 1);
	}

	void Viscoelastic::SetUpMatrixandRHSOldB(const NekDouble phystot, const Array<OneD, const Array<OneD, NekDouble> > &inarray, Array<OneD, Array<OneD, NekDouble> > &L,
		                Array<OneD, Array<OneD, NekDouble> > &M,Array<OneD, Array<OneD, NekDouble> > &RHS,const NekDouble aii_Dt)
	{
		NekDouble aiiDtinv = 1.0/aii_Dt;
		Array<OneD, NekDouble> tmp(phystot);
		// Compute 3x3 Matrix of LHS
		Vmath::Smul(phystot,-2.0,L[0],1,tmp,1);
		Vmath::Sadd(phystot,(aiiDtinv+1.0/m_weissenberg), tmp,1, M[0], 1);

		Vmath::Smul(phystot,-2.0,L[1],1,M[1],1);

		Vmath::Smul(phystot,-1.0,L[2],1,M[2],1);

		M[3] = Array<OneD, NekDouble>(phystot,(aiiDtinv+1.0/m_weissenberg));

		if(m_equationType == eViscoelasticDEVSSG || m_equationType == eViscoelasticDEVSSGALE)
		{
			Vmath::Vadd(phystot, L[0], 1, L[3],1, tmp, 1);
			Vmath::Vsub(phystot, M[3], 1, tmp,1, M[3], 1);
		}

		Vmath::Smul(phystot,-1.0,L[1],1,M[4],1);

		Vmath::Smul(phystot,-2.0,L[2],1,M[5],1);

		Vmath::Smul(phystot,-2.0,L[3],1,tmp,1);
		Vmath::Sadd(phystot,(aiiDtinv+1.0/m_weissenberg), tmp,1, M[6], 1);

		//Compute RHS

		// Calculate tyy
		Vmath::Smul(phystot,(2*(1.0-m_beta)/m_weissenberg),L[0],1,RHS[0],1);
		Vmath::Smul(phystot,aiiDtinv,inarray[0],1,tmp,1);

		Vmath::Vadd(phystot, tmp, 1, RHS[0],1, RHS[0], 1);

		Vmath::Smul(phystot,(2*(1.0-m_beta)/m_weissenberg),L[3],1,RHS[2],1);
		Vmath::Smul(phystot,aiiDtinv,inarray[2],1,tmp,1);

		Vmath::Vadd(phystot, tmp, 1, RHS[2],1, RHS[2], 1);

		Vmath::Vadd(phystot, L[1], 1, L[2],1, tmp, 1);
		Vmath::Smul(phystot,((1.0-m_beta)/m_weissenberg),tmp,1,tmp,1);
		Vmath::Smul(phystot,aiiDtinv,inarray[1],1,RHS[1],1);

		Vmath::Vadd(phystot, tmp, 1, RHS[1],1, RHS[1], 1);

	}

	void Viscoelastic::VelocityGradientTensor(Array<OneD, Array<OneD, NekDouble> > &outarray, Array<OneD, Array<OneD, NekDouble> > &L)
	{
		// Calculate velocity gradient tensor components of velocities of timestep n+1
		m_fields[m_velocity[0]]->PhysDeriv(0,outarray[m_velocity[0]], L[0]);
		m_fields[m_velocity[0]]->PhysDeriv(1,outarray[m_velocity[0]], L[1]);
		m_fields[m_velocity[1]]->PhysDeriv(0,outarray[m_velocity[1]], L[2]);
		m_fields[m_velocity[1]]->PhysDeriv(1,outarray[m_velocity[1]], L[3]);
	}

	void Viscoelastic::AddGiesekusTermtoRHS(const NekDouble phystot, const Array<OneD, const Array<OneD, NekDouble> > &stressin,
			Array<OneD, Array<OneD, NekDouble> > &RHS,Array<OneD, Array<OneD, NekDouble> > &RHSpN)
		{
			Array<OneD, NekDouble> tmp(phystot);

			//Add Giesekus-Terms to RHS

			// tmp = tauxx*tauxx+tauxy*tauxy
			// z = v*w + x*y
			Vmath::Vvtvvtp (phystot,stressin[0],1,stressin[0],1,stressin[1],1, stressin[1],1,tmp, 1);
			// txx += -m_alpha/(1-m_beta)*tmp
			Vmath::Svtvp(phystot,(-m_alpha/(1-m_beta)),tmp,1,RHS[0],1,RHSpN[0],1);

			// txy
			Vmath::Vvtvvtp (phystot,stressin[1],1,stressin[0],1,stressin[2],1, stressin[1],1,tmp, 1);
			// txx += -m_alpha/(1-m_beta)*tmp
			Vmath::Svtvp(phystot,(-m_alpha/(1-m_beta)),tmp,1,RHS[1],1,RHSpN[1],1);

			// tyy
			Vmath::Vvtvvtp (phystot,stressin[1],1,stressin[1],1,stressin[2],1, stressin[2],1,tmp, 1);
			// txx += -m_alpha/(1-m_beta)*tmp
			Vmath::Svtvp(phystot,(-m_alpha/(1-m_beta)),tmp,1,RHS[2],1,RHSpN[2],1);

		}

	// New routine with matrix inversion
	void Viscoelastic::CalcSourceTermImplicitGiesekus(const Array<OneD, const Array<OneD, NekDouble> > &inarray,
		                Array<OneD, Array<OneD, NekDouble> > &outarray,const NekDouble aii_Dt)
		{
			int phystot = m_fields[0]->GetTotPoints();
			int nvariables = m_fields.num_elements();
			int StressDim  = m_stress.num_elements();
			int VelDim  = m_velocity.num_elements();

			int nvar = VelDim*VelDim;
			int matrixdim = 7;
			int i;
			int it=1;
			static int outputcalls = 1;
			NekDouble maximum, max1;
			Array<OneD, NekDouble> Diff(phystot);
			maximum = 1.0;
			Array<OneD, NekDouble> max(StressDim);
			for(i = 0; i<StressDim; i++)
			{
				max[i]=1.0;
			}

			Array<OneD, Array<OneD, NekDouble> > M(matrixdim);

	        M[0] = Array<OneD, NekDouble>(phystot*matrixdim);
	        for(i = 1; i < matrixdim; ++i)
	        {
	        	M[i] = M[i-1] + phystot;
	        }

			Array<OneD, Array<OneD, NekDouble> > tau(StressDim);

	        tau[0] = Array<OneD, NekDouble>(phystot*StressDim);
	        for(i = 1; i < StressDim; ++i)
	        {
	        	tau[i] = tau[i-1] + phystot;
	        }

			Array<OneD, Array<OneD, NekDouble> > RHS(StressDim);

	        RHS[0] = Array<OneD, NekDouble>(phystot*StressDim);
	        for(i = 1; i < StressDim; ++i)
	        {
	        	RHS[i] = RHS[i-1] + phystot;
	        }

			// RHS plus nonlinear terms
	        Array<OneD, Array<OneD, NekDouble> > RHSpN(StressDim);

	        RHSpN[0] = Array<OneD, NekDouble>(phystot*StressDim);
	        for(i = 1; i < StressDim; ++i)
	        {
	        	RHSpN[i] = RHSpN[i-1] + phystot;
	        }

	        Array<OneD, Array<OneD, NekDouble> > L(VelDim*VelDim);

	        L[0] = Array<OneD, NekDouble>(phystot*VelDim*VelDim);
	        for(i = 1; i < (VelDim*VelDim); ++i)
	        {
	        	L[i] = L[i-1] + phystot;
	        }

			Array<OneD, Array<OneD, NekDouble> > stressin(StressDim);

	        stressin[0] = Array<OneD, NekDouble>(phystot*StressDim);
	        for(i = 1; i < StressDim; ++i)
	        {
	        	stressin[i] = stressin[i-1] + phystot;
	        }

	        for(i = 0; i < StressDim; ++i)
	       	{
	        	Vmath::Vcopy(phystot,inarray[m_stress[i]],1,stressin[i],1);
	       	}

	        //VelocityGradientTensor(outarray,L);

	        switch(m_equationType)
	        {
				 case eUnsteadyViscoelastic:
						 // Computes the Deformation Tensor and saves D[0]=dUdx, D[1]=dUdy, D[2]=dVdx, D[3]=dVdy
						 VelocityGradientTensor(outarray,L);
						 break;
				 case eViscoelasticDEVSSG:
				 case eViscoelasticDEVSSGALE:
						 // Computes the Deformation Tensor and saves D[0]=dUdx, D[1]=dUdy, D[2]=dVdx, D[3]=dVdy
							 for(i = 0; i < (VelDim*VelDim); ++i)
								{
									 m_G[i]->BwdTrans(m_G[i]->GetCoeffs(),L[i]);
								}
						 break;
	        }

			 // Loop over inverting matrix till converged.
	         SetUpMatrixandRHSOldB(phystot, stressin, L, M, RHS,aii_Dt);

			 // Start with tau(it=1) = tau^n
			 for(i = 0; i < StressDim; ++i)
			 {
				 Vmath::Vcopy(phystot,m_stressold[i],1,stressin[i],1);
			 }

		//while(maximum>1e-8 || it<=5)
			while(maximum>1e-9)
			{

				 /* calculates Giesekus Term and adds it to RHS */
				 AddGiesekusTermtoRHS(phystot, stressin,RHS,RHSpN);

				 /* Computes the inverse of M and multiplies it with RHS to obtain stress at new values and save them in tau */
				 StressMatrixInversetimesRHS(phystot, M,RHSpN,tau);

				for(i = 0; i < StressDim; ++i)
				{
					m_fields[m_stress[i]]->FwdTrans(tau[i],m_fields[m_stress[i]]->UpdateCoeffs());
					m_fields[m_stress[i]]->BwdTrans(m_fields[m_stress[i]]->GetCoeffs(),tau[i]);

					 /* Calculate maximum error */
					 Vmath::Vsub(phystot,tau[i],1,stressin[i],1,Diff,1);
					 Vmath::Vabs(phystot,Diff,1,Diff,1);
					 max[i] = Vmath::Vmax(phystot,Diff,1);
					 //max[i] = m_fields[m_stress[i]]->L2(stressinarray[m_stress[i]]);

					 /* form inarray to contain physical values from old time step */
					 /* use the updated stress for input at next iteration */
					 Vmath::Vcopy(phystot,tau[i],1,stressin[i],1);
					 //stressin[i]=tau[i];
				}

				max1 = (max[0]<max[1]) ? max[1] : max[0];
				maximum = (max1 <max[2]) ? max[2] : max1;

				if(!((outputcalls+1)%m_infosteps))
				{
					cout << "it:" << it << " max:" << maximum << endl;
				}


				it++;
			}

			for(i=0;i<StressDim;i++)
			{
				outarray[m_stress[i]]=tau[i];
			}

			outputcalls++;
		}

	void Viscoelastic::CalcSourceTermImplicitOldB(const Array<OneD, const Array<OneD, NekDouble> > &inarray,
		                Array<OneD, Array<OneD, NekDouble> > &outarray,const NekDouble aii_Dt)
		{
			int phystot = m_fields[0]->GetTotPoints();
			int nvariables = m_fields.num_elements();
			int StressDim  = m_stress.num_elements();
			int VelDim  = m_velocity.num_elements();

			int nvar = VelDim*VelDim;
			int matrixdim = 7;
			int i;

			Array<OneD, Array<OneD, NekDouble> > M(matrixdim);

	        M[0] = Array<OneD, NekDouble>(phystot*matrixdim);
	        for(i = 1; i < matrixdim; ++i)
	        {
	        	M[i] = M[i-1] + phystot;
	        }

			Array<OneD, Array<OneD, NekDouble> > tau(StressDim);

	        tau[0] = Array<OneD, NekDouble>(phystot*StressDim);
	        for(i = 1; i < StressDim; ++i)
	        {
	        	tau[i] = tau[i-1] + phystot;
	        }

			Array<OneD, Array<OneD, NekDouble> > RHS(StressDim);

	        RHS[0] = Array<OneD, NekDouble>(phystot*StressDim);
	        for(i = 1; i < StressDim; ++i)
	        {
	        	RHS[i] = RHS[i-1] + phystot;
	        }

	        Array<OneD, Array<OneD, NekDouble> > L(VelDim*VelDim);

	        L[0] = Array<OneD, NekDouble>(phystot*VelDim*VelDim);
	        for(i = 1; i < (VelDim*VelDim); ++i)
	        {
	        	L[i] = L[i-1] + phystot;
	        }

			Array<OneD, Array<OneD, NekDouble> > stressin(StressDim);

	        stressin[0] = Array<OneD, NekDouble>(phystot*StressDim);
	        for(i = 1; i < StressDim; ++i)
	        {
	        	stressin[i] = stressin[i-1] + phystot;
	        }

	        for(i = 0; i < StressDim; ++i)
	       	{
	       	    stressin[i] = inarray[m_stress[i]];
	       	}

	        //VelocityGradientTensor(outarray,L);

	        switch(m_equationType)
	        {
				 case eUnsteadyViscoelastic:
						 // Computes the Deformation Tensor and saves D[0]=dUdx, D[1]=dUdy, D[2]=dVdx, D[3]=dVdy
						 VelocityGradientTensor(outarray,L);
						 break;
				 case eViscoelasticDEVSSG:
				 case eViscoelasticDEVSSGALE:
						 // Computes the Deformation Tensor and saves D[0]=dUdx, D[1]=dUdy, D[2]=dVdx, D[3]=dVdy
							 for(i = 0; i < (VelDim*VelDim); ++i)
								{	// TODO: change this back !!
								 	 Vmath::Vcopy(m_G[0]->GetTotPoints(),m_G[i]->GetPhys(),1,L[i],1);
									// m_G[i]->BwdTrans(m_G[i]->GetCoeffs(),L[i]);
								}
						 break;
	        }

			 /* Fills 3x3 matrix M and RHS for the Oldroyd-B model*/
			 SetUpMatrixandRHSOldB(phystot, stressin, L, M, RHS,aii_Dt);

			 /* Computes the inverse of M and multiplies it with RHS to obtain stress at new values and save them in tau */
			 StressMatrixInversetimesRHS(phystot, M,RHS,tau);

			for(i=0;i<StressDim;i++)
			{
				outarray[m_stress[i]]=tau[i];
			}
		}

	void Viscoelastic::SetInflowBoundaryConditionsGiesekus(const Array<OneD, const Array<OneD, NekDouble> > &outarray,const NekDouble aii_Dt)
		{			// VelDim*VelDim*phystot = Number of components of rate of deformation tensor

/*			    int  i,el,n, nquad_e,id1;
			    int StressDim=m_stress.num_elements();
			    int VelDim=m_velocity.num_elements();
			    int    elmtid,nq,offset, edge;

			    StdRegions::StdExpansionSharedPtr elmt;
			   // NekDouble sumdrag=0,sumlift=0;
			    Array<OneD, const SpatialDomains::BoundaryConditionShPtr > txxBndConds;
			    Array<OneD, MultiRegions::ExpListSharedPtr> txxBndExp;
			    Array<OneD, MultiRegions::ExpListSharedPtr>  txyBndExp;
			    Array<OneD, MultiRegions::ExpListSharedPtr>  tyyBndExp;

			    txxBndConds = m_fields[m_stress[0]]->GetBndConditions();
			    txxBndExp   = m_fields[m_stress[0]]->GetBndCondExpansions();
			    txyBndExp   = m_fields[m_stress[1]]->GetBndCondExpansions();
			    tyyBndExp   = m_fields[m_stress[2]]->GetBndCondExpansions();
			    Array<OneD, int> ElmtID,EdgeID;

			    //Fills ElmtID and EdgeID with global ids of elements id number and edges id number of the boundary expansion
			    m_fields[m_stress[0]]->GetBoundaryToElmtMap(ElmtID,EdgeID);

			    StdRegions::StdExpansion1DSharedPtr txxEdgeExp, txyEdgeExp, tyyEdgeExp;
			    // Pointer to incoming information
			    Array<OneD, const NekDouble> U,V;

			    int maxpts = 0,cnt;

			    // find the maximum values of points in an element along boundary
			    for(cnt = n = 0; n < txxBndConds.num_elements(); ++n)
			    {
			//Loop over all elements in boundary region
			    	//cnt starts with element 0 then going through elements
			        for(i = 0; i < txxBndExp[n]->GetExpSize(); ++i)
			        {
			        	nq = m_fields[m_stress[0]]->GetExp(ElmtID[cnt++])->GetTotPoints();
			        	maxpts = max(maxpts, nq);
			        }
			    }

			    ASSERTL0(m_expdim == 2,"Not set up for 3D expansions");

			    int nvar = VelDim*VelDim;
			    int matrixdim = 7;

				Array<OneD, Array<OneD, NekDouble> > M(matrixdim);

				M[0] = Array<OneD, NekDouble>(maxpts*matrixdim);
				for(i = 1; i < matrixdim; ++i)
				{
					M[i] = M[i-1] + maxpts;
				}

				Array<OneD, Array<OneD, NekDouble> > tau(StressDim);

				tau[0] = Array<OneD, NekDouble>(maxpts*StressDim);
				for(i = 1; i < StressDim; ++i)
				{
					tau[i] = tau[i-1] + maxpts;
				}

				Array<OneD, Array<OneD, NekDouble> > RHS(StressDim);

				RHS[0] = Array<OneD, NekDouble>(maxpts*StressDim);
				for(i = 1; i < StressDim; ++i)
				{
					RHS[i] = RHS[i-1] + maxpts;
				}

				Array<OneD, Array<OneD, NekDouble> > D(VelDim*VelDim);

				D[0] = Array<OneD, NekDouble>(maxpts*VelDim*VelDim);
				for(i = 1; i < (VelDim*VelDim); ++i)
				{
					D[i] = D[i-1] + maxpts;
				}

			    Array<OneD, NekDouble> txx(maxpts),txy(maxpts),tyy(maxpts);

				Array<OneD, Array<OneD, NekDouble> > stressin(StressDim);*/

			/*	stressin[0] = Array<OneD, NekDouble>(maxpts*StressDim);
				for(i = 1; i < StressDim; ++i)
				{
					stressin[i] = stressin[i-1] + maxpts;
				} */

			// loop over all boundary regions
	/*		    for(cnt = n = 0; n < txxBndConds.num_elements(); ++n)
			    {
			        string type = txxBndConds[n]->GetUserDefined().GetEquation();
			//for the ones marked as Drag calculate drag along boundary edge
			        if(type == "Wall")
			        {
			            if(m_expdim == 2)
			            {
			                //loop over all elements along the boundary region
			            	for(el = 0; el < txxBndExp[n]->GetExpSize(); ++el,cnt++)
			                {
			                    // find element and edge of this expansion.
			                    // calculate drag;
			                    elmtid = ElmtID[cnt];
			                    elmt   = m_fields[m_stress[0]]->GetExp(elmtid);
			                    nq     = elmt->GetTotPoints();
			                    //offset to point to values of element
			                    offset = m_fields[m_stress[0]]->GetPhys_Offset(elmtid); */

			                    /*
			                     * Let U,V, stressin point to the values of the element
			                     */
	/*		                    U = outarray[m_velocity[0]] + offset;
			                    V = outarray[m_velocity[1]] + offset;
			                    if(m_equationType==eUnsteadyViscoelastic)
			                    {
			                    	for(i=0;i<StressDim;i++)
			                    	{
			                    		// should be m_stressold as a start
										stressin[i] = m_stressold[i] + offset;
			                    	}
			                    }

			                    // Calculating deformation tensor
			                    elmt->PhysDeriv(0,U,D[0]);
			                    elmt->PhysDeriv(1,U,D[1]);
			                    elmt->PhysDeriv(0,V,D[2]);
			                    elmt->PhysDeriv(1,V,D[3]); */

			        			/* Fills 3x3 matrix M and RHS for the Oldroyd-B model*/
			     //   			SetUpMatrixandRHSOldB(nq, stressin, D, M, RHS,aii_Dt);

			        			/* Computes the inverse of M and multiplies it with RHS to obtain stress at new values and save them in tau */
			    //    			StressMatrixInversetimesRHS(nq, M,RHS,tau);

			                    /*
			                     * From here we perform the operations along the edge
			                     */

			                    // Gets boundary edge expansion do i use i or edge ID???? as we want to copy values into
			                    // boundary expansion for the pressure
			                    // we want to access the information of the boundary expansion???
			                    // so we might need to exchange elmt->GetEdgeExp by the following line
			                    // should give actual edge in boundary expansion
			    /*                txxEdgeExp =  boost::dynamic_pointer_cast<StdRegions::StdExpansion1D> (txxBndExp[n]->GetExp(el));
			                    txyEdgeExp =  boost::dynamic_pointer_cast<StdRegions::StdExpansion1D> (txyBndExp[n]->GetExp(el));
			                    tyyEdgeExp =  boost::dynamic_pointer_cast<StdRegions::StdExpansion1D> (tyyBndExp[n]->GetExp(el));


			                    // grab edge on the boundary
			                    edge = EdgeID[cnt];

			                    // Get edge values and put into txx, txy, tyy (elements of cauchy stress)
			                    elmt->GetEdgePhysVals(edge,txxEdgeExp,tau[0],txx);
			                    elmt->GetEdgePhysVals(edge,txyEdgeExp,tau[1],txy);
			                    elmt->GetEdgePhysVals(edge,tyyEdgeExp,tau[2],tyy);

			                    nquad_e = txxEdgeExp->GetNumPoints(0);
			                    // id1: id of vertex in boundary condition expansion
			                    id1  = txxBndExp[n]->GetPhys_Offset(el);

			        		    Vmath::Vcopy(nquad_e,&txx[0], 1,&(txxBndExp[n]->UpdatePhys())[id1],1);
			        		    Vmath::Vcopy(nquad_e,&txy[0], 1,&(txyBndExp[n]->UpdatePhys())[id1],1);
			        		    Vmath::Vcopy(nquad_e,&tyy[0], 1,&(tyyBndExp[n]->UpdatePhys())[id1],1); */

			                 /*   cout << "txxonbc: " ;
			                    for(int j=0;j<nquad_e;j++)
			                    {
			                    	cout << txx[j] << ", ";
			                    }
			                    cout<<endl;

			                    cout << "txyonbc: " ;
			                    for(int j=0;j<nquad_e;j++)
			                    {
			                    	cout << txy[j] << ", ";
			                    }
			                    cout<<endl;

			                    cout << "tyyonbc: " ;
			                    for(int j=0;j<nquad_e;j++)
			                    {
			                    	cout << tyy[j] << ", ";
			                    }
			                    cout<<endl; */


			                    // calcuate (phi, dp/dn = [N-kinvis curl x curl v].n)
			                   // Pvals = PBndExp[n]->UpdateCoeffs()+PBndExp[n]->GetCoeff_Offset(i);
			                    // Decide if normals facing outwards
			                   // NegateNormals = (elmt->GetEorient(edge) == StdRegions::eForwards)? false:true;

			                   // Pbc->NormVectorIProductWRTBase(Uy,Vx,Pvals,NegateNormals);
			 /*               }
			            	txxBndExp[n]->FwdTrans_BndConstrained(txxBndExp[n]->GetPhys(),txxBndExp[n]->UpdateCoeffs());
			            	txyBndExp[n]->FwdTrans_BndConstrained(txyBndExp[n]->GetPhys(),txyBndExp[n]->UpdateCoeffs());
			            	tyyBndExp[n]->FwdTrans_BndConstrained(tyyBndExp[n]->GetPhys(),tyyBndExp[n]->UpdateCoeffs());
			            }
			        }
			        else if(type == "" || type == "TimeDependent")  // for all other boundary conditions

			        {
			            cnt += txxBndExp[n]->GetExpSize();
			        }
			            else
			            {
			                ASSERTL1(false,"Problems setting wall Stress BC");
			            }

			        } */
			}

	void Viscoelastic::SetStressWallBoundaryConditionsOldB(const Array<OneD, const Array<OneD, NekDouble> > &outarray,const NekDouble aii_Dt)
		{			// VelDim*VelDim*phystot = Number of components of rate of deformation tensor

			}

	void Viscoelastic::SetStressSymmetryBoundaryConditionsOldB(const Array<OneD, const Array<OneD, NekDouble> > &inarray,const Array<OneD, const Array<OneD, NekDouble> > &outarray,const NekDouble aii_Dt)
			{			// VelDim*VelDim*phystot = Number of components of rate of deformation tensor

				}

	// Waters analytical solution for velocity component u and tauxx, tauxy
	// physarray will be tau^n+1 (maybe positioned in SolveUnsteadyStokes after computation of tau^n+1 the except for Initial time)
	void Viscoelastic::SetBoundaryWatersTAU(int bcRegion, int cnt, NekDouble time)
	{
		int nvariables      = 2; //only stress components tauxx and tauxy have to be computed
	    int nTraceNumPoints = m_fields[m_stress[0]]->GetTrace()->GetNpoints();
	    int i;
	    Array<OneD, Array<OneD, NekDouble> > Fwd(nvariables);
        Fwd[0] = Array<OneD, NekDouble>(nTraceNumPoints*nvariables);
        Fwd[1] = Fwd[0] + nTraceNumPoints;

	    // get physical values of the forward trace (from exp to phys)
        // to still have other bcs
	//    for (int i = 0; i < nvariables; ++i)
	  //  {
	   /* 	m_fields[m_stress[0]]->ExtractTracePhys(physarray[0],Fwd[0]);
	    	m_fields[m_stress[0]]->ExtractTracePhys(physarray[m_stress[0]],Fwd[1]);
	    	m_fields[m_stress[0]]->ExtractTracePhys(physarray[m_stress[1]],Fwd[2]); */
	 //   }

	    for(int e = 0; e < m_fields[m_stress[0]]->GetBndCondExpansions()[bcRegion]->GetExpSize(); ++e)
	    {
	    	int npoints = m_fields[m_stress[0]]->GetBndCondExpansions()[bcRegion]->GetExp(e)->GetNumPoints(0);
	    	// id1: id of vertex in boundary condition expansion
	    	int id1  = m_fields[m_stress[0]]->GetBndCondExpansions()[bcRegion]->GetPhys_Offset(e) ;
	    	// id2: id of vertex in trace
	    	int id2  = m_fields[m_stress[0]]->GetTrace()->GetPhys_Offset(m_fields[m_stress[0]]->GetTraceMap()->GetBndCondTraceToGlobalTraceMap(cnt+e));

	    	Array<OneD,NekDouble> x0(npoints,0.0);
	    	Array<OneD,NekDouble> x1(npoints,0.0);
	    	Array<OneD,NekDouble> x2(npoints,0.0);

	    	m_fields[m_stress[0]]->GetBndCondExpansions()[bcRegion]->GetExp(e)->GetCoords(x0,x1,x2);

	    	NekDouble y;

			// Loop on all the points of that edge
			for(int j = 0; j < npoints; j++)
			{
// Fwd[1][kk] and Fwd[2][kk] and copy them separately into into boundary expansion
				y=x1[j];
				int kk = id2+j;
				CalcDirichletWatersTAU(y, Fwd[0][kk], Fwd[1][kk],time);
				// Test it first with steady state solution!! :)
				// SetDirichletWatersSteadyState(y, Fwd[0][kk], Fwd[1][kk]);

				//cout << "(x,y)=(" << x0[j] << "," << y << ")" << endl;
				//cout << "tauxx=" << Fwd[0][kk] << " ,tauxy=" << Fwd[1][kk] << endl;
		    }

		    Vmath::Vcopy(npoints,&Fwd[0][id2], 1,&(m_fields[m_stress[0]]->GetBndCondExpansions()[bcRegion]->UpdatePhys())[id1],1);
		    Vmath::Vcopy(npoints,&Fwd[1][id2], 1,&(m_fields[m_stress[1]]->GetBndCondExpansions()[bcRegion]->UpdatePhys())[id1],1);
/*			   m_fields[0]->GetBndCondExpansions()[bcRegion]->FwdTrans_BndConstrained(m_fields[0]->GetBndCondExpansions()[bcRegion]->GetPhys(),
			    		    		m_fields[0]->GetBndCondExpansions()[bcRegion]->UpdateCoeffs());
			    m_fields[m_stress[0]]->GetBndCondExpansions()[bcRegion]->FwdTrans_BndConstrained(m_fields[m_stress[0]]->GetBndCondExpansions()[bcRegion]->GetPhys(),
			    		    		m_fields[m_stress[0]]->GetBndCondExpansions()[bcRegion]->UpdateCoeffs());
			    m_fields[m_stress[1]]->GetBndCondExpansions()[bcRegion]->FwdTrans_BndConstrained(m_fields[m_stress[1]]->GetBndCondExpansions()[bcRegion]->GetPhys(),
			    		    		m_fields[m_stress[1]]->GetBndCondExpansions()[bcRegion]->UpdateCoeffs()); */

	    }

	  }

	void Viscoelastic::SetBoundaryWatersU2(int bcRegion, int cnt, NekDouble time)
	{
		NekDouble y;
		int npoints = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetNpoints();
		Array<OneD,NekDouble> x0(npoints,0.0);
		Array<OneD,NekDouble> x1(npoints,0.0);
		Array<OneD,NekDouble> x2(npoints,0.0);

		m_fields[0]->GetBndCondExpansions()[bcRegion]->GetCoords(x0,x1,x2);

		for(int j = 0; j < npoints; j++)
		{
			y=x1[j];
			(m_fields[0]->GetBndCondExpansions()[bcRegion]->UpdatePhys())[j]
			= CalcDirichletWatersU(y,time);
					//4.0*y*(1-y);
			//cout << "(x,y)=(" << x0[j] << "," << y << ")" << endl;
			//cout << "u=" << (m_fields[0]->GetBndCondExpansions()[bcRegion]->GetPhys())[j];
		}
	}


	/// Calculate time dependent Dirichlet Inflow/ Outflow Boundary Condition for the Velocity after Waters
	// method = 1 calculate velocity, method=2 TAUxy , method=3 TAUxx
	NekDouble Viscoelastic::CalcDirichletWatersU(NekDouble y, NekDouble time)
	{
		NekDouble t,Ay,DAy,El,sumu,u;
		NekDouble pi=3.1415926535;
		NekDouble N,alphaN,betaN2,betaN,alphaNstar,betaNstar,gammaN,aN,bN,pN,qN;
		NekDouble fac,a,b,GNt;
		NekDouble Umax = 1.0;

		// number of summands
		int ntimes = 20;

		// y = y / height;

		t = time*m_kinvis;	  // Non dimensionalized time in analytic solution

		Ay = 4*y*(1-y);  // parabolic profile
		DAy = -4*(2*y-1); // derivative of parabolic profile
		sumu = 0;

		if(m_equationType==eUnsteadyViscoelastic)
		{
		El = m_weissenberg*m_kinvis; // Elasticity number

		// if(dimension!=1.0)
		// {
		// 	El = El / ((dimension)*(dimension));
		// }

		for (int i = 1; i <= ntimes; ++i) {
			N = (2*i - 1)*pi;
			alphaN = 1 + m_beta * El * pow(N,2);
			betaN2 = pow(alphaN,2) - 4 * El * pow(N,2);
			betaN  = sqrt(fabs(betaN2));
			alphaNstar = alphaN / (2*El);
			betaNstar = betaN / (2*El);
			gammaN = 1 + pow(N,2) * El * (m_beta -2);
			aN = 1 + (gammaN/betaN);
			bN = 1 - (gammaN/betaN);
			pN = -alphaNstar + betaNstar;
			qN = -alphaNstar - betaNstar;

			if(betaN2>=0)
			{
				fac = 0.5;
				GNt =  GetGNexp(t,fac,aN,pN,bN,qN);
				// GetGNexp(t,fac,a,p,b,q)=fac*( (a*exp(p*t)) + (b*exp(q*t)))
			}
			else
			{
				fac = exp(-alphaNstar*t);
				a = 1.;
				b =  gammaN/betaN;
				GNt =  GetGNcossin(t,fac,a,betaNstar,b,betaNstar);
				//GetGNcossin(t,fac,a,p,b,q)=fac*(a*cos(p*t) + b*sin(q*t))
			}

			sumu += ((sin(N*y)/pow(N,3)) * GNt);
		}
		}
		// Newtonian solution for time-dependent Poiseuille Flow
		else{
		for (int i = 1; i <= ntimes; ++i)
		{
			N = (2*i - 1)*pi;
			sumu += ((sin(N*y)/pow(N,3)) * exp(-t*N*N));
		}
		}

		u = Ay - 32*sumu;

		u = Umax*u;

		//printf("(%f,%f)\n",u,y);fflush(stdout);

		return u;
	}

	// Calculate time dependent Dirichlet Inflow/ Outflow Boundary Condition for the Velocity after Waters
	// method = 1 calculate velocity, method=2 TAUxy , method=3 TAUxx
	void Viscoelastic::CalcDirichletWatersTAU(NekDouble y, NekDouble &tauxx, NekDouble &tauxy, NekDouble time)
		{
			NekDouble t,Ay,DAy,El,sumtauxy,Reeff;
			NekDouble pi=3.1415926535;
			NekDouble N,alphaN,betaN2,betaN,alphaNstar,betaNstar,gammaN,aN,bN,pN,qN;
			NekDouble cN2,cN4,hN;
			NekDouble fac,a,b,p,q,HNt,INt,JNt,KNMt1,KNMt2,KNMt;
			NekDouble sumxx1=0.0,sumxx2=0.0,sumxx3=0.0,sumxx4=0.0;
			NekDouble tauxx1=0.0,tauxx2=0.0,tauxx3=0.0,tauxx4=0.0;
			NekDouble Cxy=0.0, Cxx=0.0;

			NekDouble M,alphaM,betaM2,betaM,alphaMstar,betaMstar,gammaM,aM,bM,pM,qM;
			NekDouble cM2,hM,hNM,betaNMplus,betaNMminus,cNM2plus,cNM2minus,ANM,BNM,DNM,ENM;

			// number of summands
			int ntimes = 10;
			NekDouble Umax=1.0;

	//		y = y / height;
			t = time*m_kinvis;	  // Non dimensionalized time in analytic solution
			El = m_weissenberg*m_kinvis; // Elasticity number

			// Compute constants for y
			CalcDirichletWatersTAUCompConstants(y, Cxx, Cxy);

			//std::cout << Cxx << Cxy;

			// Reeff = umax*m_reynolds;

			// if(dimension!=1.0)
			// {
			// 	y = y / dimension;
			//	t = t / ((dimension)*(dimension));
			//	El = El / ((dimension)*(dimension));
			//	Reeff = umax*dimension*m_reynolds;
			// }

			Ay = 4*y*(1-y);  // parabolic profile
			DAy = -4*(2*y-1); // derivative of parabolic profile

		    sumtauxy=0;

			for (int i = 1; i <= ntimes; ++i) {
				N = (2*i - 1)*pi;
				alphaN = 1 + m_beta * El * pow(N,2);
				betaN2 = pow(alphaN,2) - 4 * El * pow(N,2);
				betaN  = sqrt(fabs(betaN2));
				alphaNstar = alphaN / (2*El);
				betaNstar = betaN / (2*El);
				gammaN = 1 + pow(N,2) * El * (m_beta -2);
				aN = 1 + (gammaN/betaN);
				bN = 1 - (gammaN/betaN);
				pN = -alphaNstar + betaNstar;
				qN = -alphaNstar - betaNstar;

				if(betaN2>=0)
				{
					fac = 0.5;

					//Values for Tauxy, tauxx2
					a = aN/(pN+(1/El));
					b = bN/(qN+(1/El));
					HNt =  GetGNexp(t,fac,a,pN,b,qN);
					// GetGNexp(t,fac,a,p,b,q)=fac*( (a*exp(p*t)) + (b*exp(q*t)))

					//tauxx1
					a = aN/pN;
					p = pN - 1/El;
					b = bN/qN;
					q = qN - 1/El;
					INt = GetGNexp(t,fac,a,p,b,q);

					//tauxx3
					a = aN/((pN+1/El)*(pN+1/El));
					b = bN/((qN+1/El)*(qN+1/El));
					JNt =  GetGNexp(t,fac,a,pN,b,qN);

				}
				else
				{
					hN  = -alphaNstar + (1/El);
					cN2 = pow(hN,2) + pow(betaNstar,2);

					//Values for Tauxy, tauxx2
					fac = (exp(-alphaNstar*t))/cN2;
					b = betaNstar + (hN*gammaN)/betaN;
					a =  hN - (betaNstar*gammaN)/betaN;
					HNt =  GetGNcossin(t,fac,a,betaNstar,b,betaNstar);
					//GetGNcossin(t,fac,a,p,b,q)=fac*(a*cos(p*t) + b*sin(q*t))

					//tauxx1
					fac = ((exp((-alphaNstar-1/El)*t))/(pow(alphaNstar,2)+pow(betaNstar,2)));
					a = -alphaNstar - betaNstar*(gammaN/betaN);
					b = (betaNstar-alphaNstar*(gammaN/betaN));
					INt =  GetGNcossin(t,fac,a,betaNstar,b,betaNstar);

					//tauxx3
					cN4 = cN2*cN2; //TODO: Was that meant from Carew?
					fac = exp(-alphaNstar*t)/cN4;
					a = hN*(hN - (betaNstar*gammaN)/betaN) - betaNstar*(betaNstar + (hN*gammaN)/betaN);
					b = hN*(betaNstar + (hN*gammaN)/betaN) + betaNstar*(hN - (betaNstar*gammaN)/betaN);
					JNt =  GetGNcossin(t,fac,a,betaNstar,b,betaNstar);
				}

				sumtauxy += ((cos(N*y)/pow(N,2)) * HNt); //and sumxx2
				sumxx1 += ((cos(N*y)/pow(N,2)) * INt);
				sumxx3 += ((cos(N*y)/pow(N,2)) * JNt);

				// mixed sum!!! KNM
				for (int j = 1; j <= ntimes; ++j) {
					M = (2*j - 1)*pi;
					alphaM = 1 + m_beta * El * pow(M,2);
					betaM2 = pow(alphaM,2) - 4 * El * pow(M,2);
					betaM  = sqrt(fabs(betaM2));
					alphaMstar = alphaM / (2*El);
					betaMstar = betaM / (2*El);
					gammaM = 1 + pow(M,2) * El * (m_beta -2);
					aM = 1 + (gammaM/betaM);
					bM = 1 - (gammaM/betaM);
					pM = -alphaMstar + betaMstar;
					qM = -alphaMstar - betaMstar;


					if(betaN2>=0 && betaM2>=0)
					{
						fac = 0.25;

						//tauxx4-1
						a = (aN*aM)/((pM+(1/El))*(pN+pM+(1/El)));
						p = pN + pM;
						b = (aN*bM)/((qM+(1/El))*(pN+qM+(1/El)));
						q = pN + qM;
						KNMt1 =  GetGNexp(t,fac,a,p,b,q);
						// GetGNexp(t,fac,a,p,b,q)=fac*( (a*exp(p*t)) + (b*exp(q*t)))

						//tauxx4-2
						a = (bN*aM)/((pM+(1/El))*(qN+pM+(1/El)));
						p = qN + pM;
						b = (bN*bM)/((qM+(1/El))*(qN+qM+(1/El)));
						q = qN + qM;
						KNMt2 =  GetGNexp(t,fac,a,p,b,q);

						KNMt = KNMt1 + KNMt2;
					}
					else if(betaN2<0 && betaM2<0)
					{
						hM  = -alphaMstar + (1/El);
						cM2 = pow(hM,2) + pow(betaMstar,2);

						fac = (exp(-(alphaNstar+alphaMstar)*t))/cM2;

						//mixed terms
						hNM = -alphaNstar -alphaMstar + 1/El;
						betaNMplus = betaNstar + betaMstar;
						betaNMminus = betaNstar - betaMstar;
						cNM2plus = pow(hNM,2) + pow(betaNMplus,2);
						cNM2minus = pow(hNM,2) + pow(betaNMminus,2);

						ANM = 0.5*(betaMstar + hM*(gammaM/betaM) + hM*(gammaN/betaN) - betaMstar*(gammaN/betaN)*(gammaM/betaM));
						BNM = 0.5*(-betaMstar - hM*(gammaM/betaM) + hM*(gammaN/betaN) - betaMstar*(gammaN/betaN)*(gammaM/betaM));
						DNM = 0.5*(hM - betaMstar*(gammaM/betaM) - betaMstar*(gammaN/betaN) - hM*(gammaN/betaN)*(gammaM/betaM));
						ENM = 0.5*(hM - betaMstar*(gammaM/betaM) + betaMstar*(gammaN/betaN) + hM*(gammaN/betaN)*(gammaM/betaM));


						//tauxx4-1
						b = (ANM*hNM + DNM*betaNMplus)/cNM2plus;
						a = (DNM*hNM - ANM*betaNMplus)/cNM2plus;
						KNMt1 =  GetGNcossin(t,fac,a,betaNMplus,b,betaNMplus);
						//GetGNcossin(t,fac,a,p,b,q)=fac*(a*cos(p*t) + b*sin(q*t))

						//tauxx4-2
						b = (BNM*hNM + ENM*betaNMminus)/cNM2minus;
						a = (ENM*hNM - BNM*betaNMminus)/cNM2minus;
						KNMt2 =  GetGNcossin(t,fac,a,betaNMminus,b,betaNMminus);

						KNMt = KNMt1 + KNMt2;
					}
					else
					{
						//printf("Warning: Mixed term for N=%f and M=%f is not considered\n",N,M);
						//fflush(stdout);
					}

					sumxx4 += ((cos(N*y)/pow(N,2))*(cos(M*y)/pow(M,2)) * KNMt);
				//end m	loop
				}
		    //end n loop
			}

			sumxx2 = sumtauxy;

			tauxy = ((1-m_beta)/El)*(El*DAy - 32*sumtauxy) + Cxy*exp(-t/El);

			tauxx2 = 2*m_reynolds*DAy*(1-m_beta)*(El*DAy - 32*sumxx2);
			tauxx1 = 2*m_reynolds*Cxy*(DAy*exp(-t/El)*t - 32*sumxx1);
			tauxx3 = - ((64*m_reynolds*DAy*(1-m_beta))/El)*sumxx3;
			tauxx4 = ((2048*m_reynolds*(1-m_beta))/El)*sumxx4;

			tauxx = tauxx1 + tauxx2 + tauxx3 + tauxx4 + Cxx*exp(-t/El);

			tauxx = pow(Umax,2)*tauxx;
			tauxy = Umax*tauxy;

			//cout << "tauxx in function" << tauxx << endl;
			//cout << "tauxy:" << tauxy << endl;

			// printf("Inflow BC: Cxy[%i]=%f, Cxx[%i]=%f for y=(%f)\n",index,Cxy[index],index,Cxx[index],y);fflush(stdout);
			// TODO: Is this necessary or not?
			// if(dimension!=1.0)
			// {
			// 	tauxx = tauxx / dimension;
			//	tauxy = tauxy / dimension;
			// }

			//printf("(%f,%f)\n",u,y);fflush(stdout);
		}

	// Calculate time dependent Dirichlet Inflow/ Outflow Boundary Condition for the Velocity after Waters
	// method = 1 calculate velocity, method=2 TAUxy , method=3 TAUxx
	void Viscoelastic::CalcDirichletWatersTAUCompConstants(NekDouble y, NekDouble &Cxx, NekDouble &Cxy)
	{
		NekDouble Ay,DAy,El,sumCxy;
		NekDouble pi=3.1415926535;
		NekDouble N,alphaN,betaN2,betaN,alphaNstar,betaNstar,gammaN,aN,bN,pN,qN;
		NekDouble cN2,cN4,hN;
		NekDouble fac,a,b,HN0,IN0,JN0,KNM01,KNM02,KNM0;
		NekDouble sumCxx1=0.0,sumCxx2=0.0,sumCxx3=0.0,sumCxx4=0.0;
		NekDouble Cxx1=0.0,Cxx2=0.0,Cxx3=0.0,Cxx4=0.0;

		NekDouble M,alphaM,betaM2,betaM,alphaMstar,betaMstar,gammaM,aM,bM,pM,qM;
		NekDouble cM2,hM,hNM,betaNMplus,betaNMminus,cNM2plus,cNM2minus,ANM,BNM,DNM,ENM;

		Ay = 4*y*(1-y);  // parabolic profile
		DAy = -4*(2*y-1); // derivative of parabolic profile
		El = m_weissenberg*m_kinvis; // Elasticity number
	    sumCxy=0;

	    int ntimes = 20;

		for (int i = 1; i <= ntimes; ++i) {
			N = (2*i - 1)*pi;
			alphaN = 1 + m_beta * El * pow(N,2);
			betaN2 = pow(alphaN,2) - 4 * El * pow(N,2);
			betaN  = sqrt(fabs(betaN2));
			alphaNstar = alphaN / (2*El);
			betaNstar = betaN / (2*El);
			gammaN = 1 + pow(N,2) * El * (m_beta -2);
			aN = 1 + (gammaN/betaN);
			bN = 1 - (gammaN/betaN);
			pN = -alphaNstar + betaNstar;
			qN = -alphaNstar - betaNstar;

			if(betaN2>=0)
			{
				fac = 0.5;

				//Values for Tauxy, tauxx2
				a = aN/(pN+(1/El));
				b = bN/(qN+(1/El));
				HN0 =  GetGN0exp(fac,a,b);
				// GetGN0exp(NekDouble fac, NekDouble a, NekDouble b)=fac*(a+b)

				//tauxx1
				a = aN/pN;
				b = bN/qN;
				IN0 = GetGN0exp(fac,a,b);

				//tauxx3
				a = aN/((pN+1/El)*(pN+1/El));
				b = bN/((qN+1/El)*(qN+1/El));
				JN0 =  GetGN0exp(fac,a,b);

			}
			else
			{
				hN  = -alphaNstar + (1/El);
				cN2 = pow(hN,2) + pow(betaNstar,2);

				//Values for Tauxy, tauxx2
				fac = 1./cN2;
				a =  hN - (betaNstar*gammaN)/betaN;
				HN0 =  GetGN0cossin(fac,a);
				//GetGN0cossin(NekDouble fac, NekDouble a)=fac*a

				//tauxx1
				fac = 1./(pow(alphaNstar,2)+pow(betaNstar,2));
				a = -alphaNstar - betaNstar*(gammaN/betaN);
				IN0 =  GetGN0cossin(fac,a);

				//tauxx3
				cN4 = cN2*cN2; //TODO: Was that meant from Carew?
				fac = 1/cN4;
				a = hN*(hN - (betaNstar*gammaN)/betaN) - betaNstar*(betaNstar + (hN*gammaN)/betaN);
				JN0 =  GetGN0cossin(fac,a);
			}

			sumCxy += ((cos(N*y)/pow(N,2)) * HN0);

			sumCxx1 += ((cos(N*y)/pow(N,2)) * IN0);
			sumCxx3 += ((cos(N*y)/pow(N,2)) * JN0);

			// mixed sum!!! KNM
			for (int j = 1; j <= ntimes; ++j) {
				M = (2*j - 1)*pi;
				alphaM = 1 + m_beta * El * pow(M,2);
				betaM2 = pow(alphaM,2) - 4 * El * pow(M,2);
				betaM  = sqrt(fabs(betaM2));
				alphaMstar = alphaM / (2*El);
				betaMstar = betaM / (2*El);
				gammaM = 1 + pow(M,2) * El * (m_beta -2);
				aM = 1 + (gammaM/betaM);
				bM = 1 - (gammaM/betaM);
				pM = -alphaMstar + betaMstar;
				qM = -alphaMstar - betaMstar;


				if(betaN2>=0 && betaM2>=0)
				{
					fac = 0.25;

					//tauxx4-1
					a = (aN*aM)/((pM+(1/El))*(pN+pM+(1/El)));
					b = (aN*bM)/((qM+(1/El))*(pN+qM+(1/El)));
					KNM01 =  GetGN0exp(fac,a,b);
					// GetGN0exp(NekDouble fac, NekDouble a, NekDouble b)=fac*(a+b)

					//tauxx4-2
					a = (bN*aM)/((pM+(1/El))*(qN+pM+(1/El)));
					b = (bN*bM)/((qM+(1/El))*(qN+qM+(1/El)));
					KNM02 =  GetGN0exp(fac,a,b);

					KNM0 = KNM01 + KNM02;
				}
				else if(betaN2<0 && betaM2<0)
				{
					hM  = -alphaMstar + (1/El);
					cM2 = pow(hM,2) + pow(betaMstar,2);

					fac = 1./cM2;

					//mixed terms
					hNM = -alphaNstar -alphaMstar + 1/El;
					betaNMplus = betaNstar + betaMstar;
					betaNMminus = betaNstar - betaMstar;
					cNM2plus = pow(hNM,2) + pow(betaNMplus,2);
					cNM2minus = pow(hNM,2) + pow(betaNMminus,2);

					ANM = 0.5*(betaMstar + hM*(gammaM/betaM) + hM*(gammaN/betaN) - betaMstar*(gammaN/betaN)*(gammaM/betaM));
					BNM = 0.5*(-betaMstar - hM*(gammaM/betaM) + hM*(gammaN/betaN) - betaMstar*(gammaN/betaN)*(gammaM/betaM));
					DNM = 0.5*(hM - betaMstar*(gammaM/betaM) - betaMstar*(gammaN/betaN) - hM*(gammaN/betaN)*(gammaM/betaM));
					ENM = 0.5*(hM - betaMstar*(gammaM/betaM) + betaMstar*(gammaN/betaN) + hM*(gammaN/betaN)*(gammaM/betaM));


					//tauxx4-1
					a = (DNM*hNM - ANM*betaNMplus)/cNM2plus;
					KNM01 =  GetGN0cossin(fac,a);
					//GetGN0cossin(NekDouble fac, NekDouble a)=fac*a

					//tauxx4-2
					a = (ENM*hNM - BNM*betaNMminus)/cNM2minus;
					KNM02 =  GetGN0cossin(fac,a);

					KNM0 = KNM01 + KNM02;
				}
				else
				{
					//printf("Warning: Mixed term for N=%f and M=%f is not considered\n",N,M);
					//fflush(stdout);
				}

				sumCxx4 += ((cos(N*y)/pow(N,2))*(cos(M*y)/pow(M,2)) * KNM0);
			//end m	loop
			}
	    //end n loop
		}

		sumCxx2 = sumCxy;

		Cxy = - (((1-m_beta)/El)*(El*DAy - 32*sumCxy));

		Cxx2 = 2*m_reynolds*DAy*(1-m_beta)*(El*DAy - 32*sumCxx2);
		Cxx1 = 2*m_reynolds*Cxy*(- 32*sumCxx1);
		Cxx3 = - ((64*m_reynolds*DAy*(1-m_beta))/El)*sumCxx3;
		Cxx4 = ((2048*m_reynolds*(1-m_beta))/El)*sumCxx4;

		Cxx = -(Cxx1 + Cxx2 + Cxx3 + Cxx4);
		// printf("Cxy[%i]=%f, Cxx[%i]=%f for y=(%f)\n",index,Cxy[index],index,Cxx[index],y);fflush(stdout);
	}

	void Viscoelastic::SetDirichletWatersSteadyState(NekDouble y, NekDouble &tauxx, NekDouble &tauxy)
	{
	   NekDouble Umax=1.0;
	   NekDouble u;
	   u     = 4.0*y*(1-y);
	   u = Umax*u;


	   NekDouble DUDY = 4*(1-2*y);
	 //  printf("y=%f,u=%f, DUDY=%f\n",y,u,DUDY); fflush(stdout);
	   tauxy = (1-m_beta)*DUDY;
	   tauxy = Umax*tauxy;
	   tauxx = 2*m_weissenberg*(1-m_beta)*DUDY*DUDY;
	   tauxx = Umax*Umax*tauxx;
	  // printf("y=%f,tauxx=%f, tauxy=%f\n",y,tauxx,tauxy); fflush(stdout);
	}

	void Viscoelastic::WriteWaterstoFile(int nchk,NekDouble x, NekDouble y, NekDouble xtol, NekDouble ytol)
	{
        int phystot = m_fields[0]->GetTotPoints();
        //static int count = 0;
        Array<OneD,NekDouble> x0(phystot,0.0);
        Array<OneD,NekDouble> x1(phystot,0.0);
        Array<OneD,NekDouble> x2(phystot,0.0);

        m_fields[0]->GetCoords(x0,x1,x2);

    	stringstream convert;

    	convert << "Waters" << m_sessionName << "x" << x << "y" << y << "Wi" << m_weissenberg << ".txt";
    	//fname = fname.substr(0,fname.find_last_of('.')) + ".txt";
    	//convert << "DragWi" << m_weissenberg << ".txt";
    	string fname = convert.str();



        if(nchk==0)
        {
        	ofstream outfile(fname.c_str());
            for(int i = 0; i < phystot; ++i)
            {
            	if(x0[i]<x+xtol && x0[i]>x-xtol)
            	{
            		if(x1[i]<y+ytol && x1[i]>=y-ytol)
            		{
            			outfile << "(x,y) = (" << x0[i] << "," << x1[i] << ")" << "\n";
            		}
            	}
            }
			// Writes header of Matlab file
			outfile.width(10);
			outfile << "t" << "\t";
			outfile.width(10);
			outfile << "u" << "\t";
			outfile.width(10);
			outfile << "tauxx" << "\t";
			outfile.width(10);
			outfile << "tauxy" << "\t"  << "\n";
			outfile.close();
        }

        //count++;

        ofstream outfile(fname.c_str(),ios::app);


        for(int i = 0; i < phystot; ++i)
        {
        	if(x0[i]<x+xtol && x0[i]>x-xtol)
        	{
        		if(x1[i]<y+ytol && x1[i]>=y-ytol)
        		{
					cout << "Writing outfile: " << fname <<"at " << "(x,y) = (" << x0[i] << "," << x1[i] << ")" << endl;
					outfile.width(10);
					outfile.precision(8);
					outfile << m_time << "\t";
					outfile.width(10);
					outfile.precision(8);
					outfile << m_fields[0]->GetPhys()[i] << "\t";
					outfile.width(10);
					outfile.precision(8);
					outfile << m_fields[m_stress[0]]->GetPhys()[i] << "\t";
					outfile.width(10);
					outfile.precision(8);
					outfile << m_fields[m_stress[1]]->GetPhys()[i] << "\t";
					outfile << std::endl;
        		}
        	}
        }

	}

void Viscoelastic::WriteMinMaxtoFileHeader(int vstart, int nvariables, string fname)
{
    ofstream outfile(fname.c_str());

    outfile.width(8);
    outfile << "t" << "\t";
    for(int i=0;i<nvariables;i++)
    {
    	outfile.width(10);
    	outfile << "Max" << GetVariable(vstart+i) << "\t";
    	outfile.width(10);
    	outfile << "Min" << GetVariable(vstart+i) << "\t";
    }

    outfile << "\n";
}

void Viscoelastic::WriteMinMaxtoFile(const Array<OneD, const Array<OneD, NekDouble> > &physfield, int vstart, int nvariables, string fname, bool IsFromMFields)
{
	int phystot = m_fields[0]->GetTotPoints();
	int i;
	Array<OneD,NekDouble> x0(phystot,0.0);
	Array<OneD,NekDouble> x1(phystot,0.0);
	Array<OneD,NekDouble> x2(phystot,0.0);

	m_fields[0]->GetCoords(x0,x1,x2);

	ofstream outfile(fname.c_str(),ios::app);

	cout << "Writing" << fname << endl;

	outfile.width(10);
	outfile.precision(8);
	outfile << m_time << "\t";

	if(IsFromMFields)
	{
    	for(i = 0; i < nvariables; ++i)
		{
			outfile.width(10);
			outfile.precision(8);
			outfile << Vmath::Vmax(phystot,m_fields[vstart+i]->GetPhys(),1) << "\t";
			outfile.width(10);
			outfile.precision(8);
			outfile << Vmath::Vmin(phystot,m_fields[vstart+i]->GetPhys(),1) << "\t";
		}
		outfile << endl;

	}
	else
	{
		for(i = 0; i < nvariables; ++i)
		{
			outfile.width(10);
			outfile.precision(8);
			outfile << Vmath::Vmax(phystot,physfield[vstart+i],1) << "\t";
			outfile.width(10);
			outfile.precision(8);
			outfile << Vmath::Vmin(phystot,physfield[vstart+i],1) << "\t";
		}
		outfile << endl;
	}

}

// inarray: values at time step t^n
void Viscoelastic::CalcDrag(const Array<OneD, const Array<OneD, NekDouble> > &inarray, NekDouble &sumdrag, NekDouble &sumlift)
{
    int  i,el,n;
   // NekDouble sumdrag=0,sumlift=0;
    Array<OneD, const SpatialDomains::BoundaryConditionShPtr > PBndConds;
    Array<OneD, MultiRegions::ExpListSharedPtr>  PBndExp;

    PBndConds = m_fields[0]->GetBndConditions();
    PBndExp   = m_fields[0]->GetBndCondExpansions();
    Array<OneD, int> ElmtID,EdgeID;

    //Fills ElmtID and EdgeID with global ids of elements id number and edges id number of the boundary expansion
    m_fields[0]->GetBoundaryToElmtMap(ElmtID,EdgeID);

    StdRegions::StdExpansion1DSharedPtr EdgeExp;
    Array<OneD, const NekDouble> U,V,txx,txy,tyy,P;
    Array<OneD, Array<OneD, NekDouble> > normals;
    	Array<OneD, const NekDouble> G11,G12,G21,G22;


    int maxpts = 0,cnt;

    int ncoeffs = m_fields[0]->GetNcoeffs();
    int phystot = m_fields[m_velocity[0]]->GetTotPoints();
    Array<OneD, NekDouble> pressurecoeff(ncoeffs,0.0);
    Array<OneD, NekDouble> pressurephys(phystot,0.0);

    // Project pressure onto velocity space
    m_pressure->BwdTrans(m_pressure->GetCoeffs(), m_pressure->UpdatePhys());

    // project onto velocity space without enforcing BCs
    m_fields[0]->FwdTrans_IterPerExp(m_pressure->GetPhys(),pressurecoeff);
    m_fields[0]->BwdTrans(pressurecoeff,pressurephys);


    // find the maximum values of points
    for(cnt = n = 0; n < PBndConds.num_elements(); ++n)
    {
//Loop over all elements in boundary region
    	//cnt starts with element 0 then going through elements
        for(i = 0; i < PBndExp[n]->GetExpSize(); ++i)
        {
            maxpts = max(maxpts, m_fields[0]->GetExp(ElmtID[cnt++])->GetTotPoints());
        }
    }



    ASSERTL0(m_expdim == 2,"Not set up for 3D expansions");

    Array<OneD, NekDouble> dUdx(10*maxpts);
    Array<OneD, NekDouble> dUdy = dUdx + maxpts;
    Array<OneD, NekDouble> dVdx  = dUdy + maxpts;
    Array<OneD, NekDouble> dVdy = dVdx  + maxpts;
    Array<OneD, NekDouble> temp = dVdy + maxpts;
    Array<OneD, NekDouble> cauchyxx = temp + maxpts;
    Array<OneD, NekDouble> cauchyxy = cauchyxx + maxpts;
    Array<OneD, NekDouble> cauchyyy = cauchyxy + maxpts;
    Array<OneD, NekDouble> drag = cauchyyy + maxpts;
    Array<OneD, NekDouble> lift = drag + maxpts;

    int    elmtid,nq,offset, edge;
    StdRegions::StdExpansionSharedPtr elmt;
    bool NegateNormals;

// loop over all boundary regions
    for(cnt = n = 0; n < m_fields[0]->GetBndConditions().num_elements(); ++n)
    {
//for the ones marked as Drag calculate drag along boundary edge
        if(m_fields[0]->GetBndConditions()[n]->GetUserDefined() == SpatialDomains::eDrag)
        {
            if(m_expdim == 2)
            {
                //loop over all elements along the boundary region
            	for(el = 0; el < PBndExp[n]->GetExpSize(); ++el,cnt++)
                {
                    // find element and edge of this expansion.
                    // calculate drag;
                    elmtid = ElmtID[cnt];
                    elmt   = m_fields[0]->GetExp(elmtid);
                    nq     = elmt->GetTotPoints();
                    //offset to point to values of element
                    offset = m_fields[0]->GetPhys_Offset(elmtid);

                    /*
                     * Calculate the Cauchy stress tensor in the elements along the boundary
                     */
                    U = inarray[m_velocity[0]] + offset;
                    V = inarray[m_velocity[1]] + offset;
                    P = pressurephys + offset;
                    if(m_equationType==eUnsteadyViscoelastic || m_equationType==eViscoelasticDEVSSG || m_equationType==eViscoelasticDEVSSGALE)
                    {
                        txx = inarray[m_stress[0]] + offset;
                        txy = inarray[m_stress[1]] + offset;
                        tyy = inarray[m_stress[2]] + offset;
                    }
                    if(m_equationType==eViscoelasticDEVSSG || m_equationType==eViscoelasticDEVSSGALE)
				   {
                        G11 = (m_G[0]->GetPhys())+ m_G[0]->GetPhys_Offset(elmtid);
                        G12 = (m_G[1]->GetPhys())+ m_G[1]->GetPhys_Offset(elmtid);
                        G21 = (m_G[2]->GetPhys())+ m_G[2]->GetPhys_Offset(elmtid);
                        G22 = (m_G[3]->GetPhys())+ m_G[3]->GetPhys_Offset(elmtid);
				   }
                    else if(m_equationType==eViscoelasticLogc)
                    {
                        txx = (m_fields[m_stress[0]]->GetPhys())+ m_fields[m_stress[0]]->GetPhys_Offset(elmtid);
                        txy = (m_fields[m_stress[1]]->GetPhys())+ m_fields[m_stress[1]]->GetPhys_Offset(elmtid);
                        tyy = (m_fields[m_stress[2]]->GetPhys())+ m_fields[m_stress[2]]->GetPhys_Offset(elmtid);
                    }

                    // Calculating deformation tensor
                    elmt->PhysDeriv(0,U,dUdx);
                    elmt->PhysDeriv(1,U,dUdy);
                    elmt->PhysDeriv(0,V,dVdx);
                    elmt->PhysDeriv(1,V,dVdy);

                    if(m_equationType==eUnsteadyViscoelastic)
                    {
                        // Calculate Cauchy stress tensor along boundary elements
                        Vmath::Vsub(nq,txx,1,P,1,cauchyxx,1);
                        Vmath::Svtvp(nq,2*m_beta,dUdx,1,cauchyxx,1,cauchyxx,1);

                        Vmath::Vadd(nq,dUdy,1,dVdx,1,temp,1);
                        Vmath::Svtvp(nq,m_beta,temp,1,txy,1,cauchyxy,1);

                        Vmath::Vsub(nq,tyy,1,P,1,cauchyyy,1);
                        Vmath::Svtvp(nq,2*m_beta,dVdy,1,cauchyyy,1,cauchyyy,1);
                    }
                    if(m_equationType==eViscoelasticDEVSSG || m_equationType==eViscoelasticDEVSSGALE)
					  {
						  // Calculate Cauchy stress tensor along boundary elements
						  Vmath::Vsub(nq,txx,1,P,1,cauchyxx,1);
						  Vmath::Svtvp(nq,2*(m_beta+m_theta),dUdx,1,cauchyxx,1,cauchyxx,1);
						  Vmath::Svtvp(nq,-2*(m_theta),G11,1,cauchyxx,1,cauchyxx,1);

						  Vmath::Vadd(nq,dUdy,1,dVdx,1,temp,1);
						  Vmath::Svtvp(nq,(m_beta+m_theta),temp,1,txy,1,cauchyxy,1);
						  Vmath::Vadd(nq,G12,1,G21,1,temp,1);
						  Vmath::Svtvp(nq,(-m_theta),temp,1,cauchyxy,1,cauchyxy,1);

						  Vmath::Vsub(nq,tyy,1,P,1,cauchyyy,1);
						  Vmath::Svtvp(nq,2*(m_beta+m_theta),dVdy,1,cauchyyy,1,cauchyyy,1);
						  Vmath::Svtvp(nq,-2*(m_theta),G22,1,cauchyyy,1,cauchyyy,1);
					  }



                    /*
                     * From here we perform the operations along the edge
                     */

                    // Gets boundary edge expansion do i use i or edge ID???? as we want to copy values into
                    // boundary expansion for the pressure
                    // we want to access the information of the boundary expansion???
                    // so we might need to exchange elmt->GetEdgeExp by the following line
                    // should give actual edge in boundary expansion
                    EdgeExp =  boost::dynamic_pointer_cast<StdRegions::StdExpansion1D> (PBndExp[n]->GetExp(el));

                    // grab edge on the boundary
                    edge = EdgeID[cnt];

                    //or ?????????
                    //EdgeExp = elmt->GetEdgeExp(edge)
                    //EdgeExp = m_fields[0]->GetExp(el)->GetEdgeExp(i);

                    //PBndExp[n]->SetUpPhysNormals(EdgeExp,edge);
                    //PBndExp[n]->GetNormals(normals);
                    //EdgeExp->SetUpPhysNormals(elmt,edge);

                    //normals	=  EdgeExp->GetMetricInfo()->GetNormal();
                    elmt->ComputeEdgeNormal(edge);
                    normals = elmt->GetEdgeNormal(edge);
                    int nquad_e = EdgeExp->GetNumPoints(0);


                    // Get edge values and put into txx, txy, tyy (elements of cauchy stress)
                    elmt->GetEdgePhysVals(edge,EdgeExp,cauchyxx,dUdx);
                    elmt->GetEdgePhysVals(edge,EdgeExp,cauchyxy,dUdy);
                    elmt->GetEdgePhysVals(edge,EdgeExp,cauchyyy,dVdy);

               /*     cout << "Cauchystressxx: " ;
                    for(int j=0;j<nquad_e;j++)
                    {
                    	cout << dUdx[j] << ", ";
                    }
                    cout<<endl;

                    cout << "Cauchystressxy: " ;
                    for(int j=0;j<nquad_e;j++)
                    {
                    	cout << dUdy[j] << ", ";
                    }
                    cout<<endl; */


                    // Does that give me the normals along the edge?
                    // NegateNormals = (elmt->GetEorient(edge) == StdRegions::eForwards)? false:true;
                    // Drag computed with inward normals????
                    NegateNormals = (elmt->GetEorient(edge) == StdRegions::eBackwards)? false:true;


                    Vmath::Vmul(nquad_e,normals[0],1,dUdx,1, drag,1);
                    Vmath::Vvtvp(nquad_e,normals[1],1,dUdy,1,drag,1,drag,1);

                   /* cout << "Normals: (n1,n2)= " ;
                    for(int j=0;j<nquad_e;j++)
                    {
                    	cout << "(" <<normals[0][j] << ", "<< normals[1][j] <<")";
                    }
                    cout<<endl; */

                    Vmath::Vmul(nquad_e,normals[0],1,dUdy,1, lift,1);
                    Vmath::Vvtvp(nquad_e,normals[1],1,dVdy,1,lift,1,lift,1);

                    if(NegateNormals == true)
                    {
                        Vmath::Neg(nquad_e,drag,1);
                        Vmath::Neg(nquad_e,lift,1);
                    }
                    // does this have to be the edge id or i? Integrate along edge
                    sumdrag += EdgeExp->Integral(drag);
                    sumlift += EdgeExp->Integral(lift);
                }
            }
        }
        else  // for all other boundary conditions
        {
            cnt += PBndExp[n]->GetExpSize();
        }

        }
/*        int n_element  = m_fields[0]->GetExpSize();                // number of element in the mesh
        for(int el = 0; el < n_element; ++el)
          {
    	int n_edges = m_fields[0]->GetExp(el)->GetNedges();     // Nedges for that element
    	int n_points = m_fields[0]->GetExp(el)->GetTotPoints(); // Nquadrature points for that element
    	Array<OneD, NekDouble> one2D(n_points, 1.0);
    	NekDouble Area = m_fields[0]->GetExp(el)->Integral(one2D);
    	Array<OneD, NekDouble> Lengths(n_edges);
    	for(int edge = 0; edge< n_edges; ++edge)
    	  {
    	    int n_pointsEdge = m_fields[0]->GetExp(el)->GetEdgeExp(edge,false)->GetTotPoints(); // Number of Quadrature Points on each edge
    	    Array<OneD, NekDouble> one1D(n_pointsEdge, 1.0);

    	    NekDouble L2 = 2.0*Area/L1;
    	    Lengths[edge] = (L1<L2) ? L1 : L2;
    	  }
    	MinLength[el] = Vmath::Vmin(n_edges,&Lengths[0],1);
          } */
   // sumdrag = 2*sumdrag;
    cout << "drag " << sumdrag << endl;
    cout << "lift " << sumlift << endl;
}
void Viscoelastic::WriteDragToFile(int count, const Array<OneD, const Array<OneD, NekDouble> > &inarray)
{
	NekDouble sumdrag=0,sumlift=0;
    // Maybe include session name!
	stringstream convert;

	//convert << "Drag" << m_sessionName << ".txt";
	convert << "Drag" << m_sessionName << ".txt";
	//fname = fname.substr(0,fname.find_last_of('.')) + ".txt";
	//convert << "DragWi" << m_weissenberg << ".txt";
	string fname = convert.str();
	fname = fname.substr(fname.find_last_of('/')+1);

	if(count==0)
	{
		WriteDragtoFileHeader(fname);
	}
    if(!((count+1)%m_infosteps))
    {
	ofstream outfile(fname.c_str(),ios::app);

	cout << "Writing " << fname << endl;

	CalcDrag(inarray,sumdrag,sumlift);

	outfile.width(10);
	outfile.precision(8);
	outfile << m_time << "\t";
	outfile.width(10);
	outfile.precision(8);
    outfile << sumdrag << "\t";
	outfile.width(10);
	outfile.precision(8);
    outfile << sumlift << "\t";

    outfile << "\n";
    }
}

void Viscoelastic::WriteDragtoFileHeader(string fname)
{
    ofstream outfile(fname.c_str());

    outfile.width(10);
    outfile << "t" << "\t";

	outfile.width(10);
	outfile << "Drag"  << "\t";
	outfile.width(10);
	outfile << "Lift" << "\t";

    outfile << "\n";
}
// inarray: values at time step t^n
void Viscoelastic::PrintStressAlongBoundary(const Array<OneD, const Array<OneD, NekDouble> > &inarray, string fname, int count)
{
    int  i,el,n;
   // NekDouble sumdrag=0,sumlift=0;
    Array<OneD, const SpatialDomains::BoundaryConditionShPtr > PBndConds;
    Array<OneD, MultiRegions::ExpListSharedPtr>  PBndExp;
int nvariables = inarray.num_elements()+1;
int phystot = inarray[0].num_elements();

    PBndConds = m_fields[m_stress[0]]->GetBndConditions();
    PBndExp   = m_fields[m_stress[0]]->GetBndCondExpansions();
    Array<OneD, int> ElmtID,EdgeID;

    //Fills ElmtID and EdgeID with global ids of elements id number and edges id number of the boundary expansion
    m_fields[m_stress[0]]->GetBoundaryToElmtMap(ElmtID,EdgeID);

    Array<OneD, StdRegions::StdExpansionSharedPtr > elmtVel(m_velocity.num_elements());
    StdRegions::StdExpansion1DSharedPtr EdgeExp;
    Array<OneD, const NekDouble> U,V,txx,txy,tyy,P;
    Array<OneD, Array<OneD, NekDouble> > normals;
    int maxpts = 0,cnt;

    Array<OneD, Array<OneD, NekDouble> > field(nvariables);
    field[0] = Array<OneD, NekDouble>(phystot*nvariables);
    for(i = 1; i < nvariables; ++i)
    {
    	field[i] = field[i-1] + phystot;
    }


    // find the maximum values of points
    for(cnt = n = 0; n < PBndConds.num_elements(); ++n)
    {
//Loop over all elements in boundary region
    	//cnt starts with element 0 then going through elements
        for(i = 0; i < PBndExp[n]->GetExpSize(); ++i)
        {
            maxpts = max(maxpts, m_fields[m_stress[0]]->GetExp(ElmtID[cnt++])->GetTotPoints());
        }
    }

    ASSERTL0(m_expdim == 2,"Not set up for 3D expansions");

    Array<OneD, NekDouble> dUdx(10*maxpts);
    Array<OneD, NekDouble> dUdy = dUdx + maxpts;
    Array<OneD, NekDouble> dVdx  = dUdy + maxpts;
    Array<OneD, NekDouble> dVdy = dVdx  + maxpts;
    Array<OneD, NekDouble> temp = dVdy + maxpts;
    Array<OneD, NekDouble> cauchyxx = temp + maxpts;
    Array<OneD, NekDouble> cauchyxy = cauchyxx + maxpts;
    Array<OneD, NekDouble> cauchyyy = cauchyxy + maxpts;
    Array<OneD, NekDouble> drag = cauchyyy + maxpts;
    Array<OneD, NekDouble> lift = drag + maxpts;

    int    elmtid,nq,offset, edge;
    StdRegions::StdExpansionSharedPtr elmt;
    bool NegateNormals;

// loop over all boundary regions
    for(cnt = n = 0; n < PBndConds.num_elements(); ++n)
    {
//for the ones marked as Drag calculate drag along boundary edge
        if(PBndConds[n]->GetUserDefined()== SpatialDomains::eDrag)
        {
            if(m_expdim == 2)
            {
            	int first=0;
                //loop over all elements along the boundary region
            	for(el = 0; el < PBndExp[n]->GetExpSize(); ++el,cnt++)
                {

                    // find element and edge of this expansion.
                    // calculate drag;
                    elmtid = ElmtID[cnt];
                    elmt   = m_fields[m_stress[0]]->GetExp(elmtid);
                    nq     = elmt->GetTotPoints();
                    //offset to point to values of element
                    offset = m_fields[m_stress[0]]->GetPhys_Offset(elmtid);

                    /*
                     * Calculate the Cauchy stress tensor in the elements along the boundary
                     */
                    U = inarray[m_velocity[0]] + offset;
                    V = inarray[m_velocity[1]] + offset;
                    P = (m_pressure->GetPhys())+ m_pressure->GetPhys_Offset(elmtid);
                    if(m_equationType==eUnsteadyViscoelastic)
                    {
                        txx = inarray[m_stress[0]] + offset;
                        txy = inarray[m_stress[1]] + offset;
                        tyy = inarray[m_stress[2]] + offset;
                    }
                    /*else if(m_equationType==eViscoelasticLogc)
                    {
                        txx = (m_fields[m_stress[0]]->GetPhys())+ m_fields[m_stress[0]]->GetPhys_Offset(elmtid);
                        txy = (m_fields[m_stress[1]]->GetPhys())+ m_fields[m_stress[1]]->GetPhys_Offset(elmtid);
                        tyy = (m_fields[m_stress[2]]->GetPhys())+ m_fields[m_stress[2]]->GetPhys_Offset(elmtid);
                    }*/

                    /*
                     * From here we perform the operations along the edge
                     */

                    // Gets boundary edge expansion do i use i or edge ID???? as we want to copy values into
                    // boundary expansion for the pressure
                    // we want to access the information of the boundary expansion???
                    // so we might need to exchange elmt->GetEdgeExp by the following line
                    // should give actual edge in boundary expansion
                    EdgeExp =  boost::dynamic_pointer_cast<StdRegions::StdExpansion1D> (PBndExp[n]->GetExp(el));

                    // grab edge on the boundary
                    edge = EdgeID[cnt];

                    //normals	=  EdgeExp->GetMetricInfo()->GetNormal();
                    normals	= elmt->GetEdgeNormal(edge);
                    int nquad_e = EdgeExp->GetNumPoints(0);

                    Array<OneD, NekDouble> x0(nquad_e);
                    Array<OneD, NekDouble> x1(nquad_e);
                    EdgeExp->GetCoords(x0,x1);

                    elmt->GetEdgePhysVals(edge,EdgeExp,U,field[0]);
                    elmt->GetEdgePhysVals(edge,EdgeExp,V,field[1]);
                    elmt->GetEdgePhysVals(edge,EdgeExp,txx,field[2]);
                    elmt->GetEdgePhysVals(edge,EdgeExp,txy,field[3]);
                    elmt->GetEdgePhysVals(edge,EdgeExp,tyy,field[4]);
                    elmt->GetEdgePhysVals(edge,EdgeExp,P,field[5]);

                    // Maybe include session name!
                    /*	stringstream convert;
                    	static int count=0;
                    	//convert << "Drag" << m_sessionName << ".txt";
                    	convert << "Values around cylinder" << m_sessionName << ".txt";
                    	//fname = fname.substr(0,fname.find_last_of('.')) + ".txt";
                    	//convert << "DragWi" << m_weissenberg << ".txt";
                    	string fname = convert.str();
                    	//fname = fname.substr(fname.find_last_of('/')+1); */

                   if(!((count+1)%m_infosteps))
                    {
                    	if(first==0)
                    	{
                    		//cout << "Writing " << fname << endl;
                    	    ofstream outfile(fname.c_str(),ios::app);

                    	    outfile.width(10);
                    	    outfile << "t" << "\t";

                    		outfile.width(10);
                    		outfile << "x"  << "\t";
                    		outfile.width(10);
                    		outfile << "y" << "\t";
                       		outfile.width(10);
							outfile << "u"  << "\t";
							outfile.width(10);
							outfile << "v" << "\t";
							outfile.width(10);
							outfile << "txx"  << "\t";
							outfile.width(10);
							outfile << "txy" << "\t";
							outfile.width(10);
							outfile << "tyy"  << "\t";
							outfile.width(10);
							outfile << "p" << "\t";

                    	    outfile << "\n";
                    	    first++;
                    	}
                       //if(!((count+1)%m_infosteps))
                       {
                    	ofstream outfile(fname.c_str(),ios::app);

                    	//cout << "Writing " << fname << endl;

                    	  for(int j=0;j<nquad_e;j++)
                    	                    {
                    	outfile.width(10);
                    	outfile.precision(8);
                    	outfile << m_time << "\t";
                    	outfile.width(10);
                    	outfile.precision(8);
                       outfile << x0[j] << "\t";
                    	outfile.width(10);
                    	outfile.precision(8);
                       outfile << x1[j] << "\t";
                       for(i = 0; i < nvariables; ++i)
                         {
                       	outfile.width(10);
                          	outfile.precision(8);
                          	outfile << field[i][j] << "\t";
                         }

                       outfile << "\n";
                    	                    }
                     //  count++;
                       }
                    }
                    // Does that give me the normals along the edge?
                    // NegateNormals = (elmt->GetEorient(edge) == StdRegions::eForwards)? false:true;
                    // Drag computed with inward normals????
                    NegateNormals = (elmt->GetEorient(edge) == StdRegions::eBackwards)? false:true;

                }
            }
        }
        else
        {
            cnt += PBndExp[n]->GetExpSize();
        }

        }
}

void Viscoelastic::AddGGTDivergenceWeakForm(Array<OneD, Array<OneD, NekDouble> > &Weakoutarray)
{
	int nCoeffs = m_fields[m_velocity[0]]->GetNcoeffs();
	int phystot = m_fields[m_velocity[0]]->GetTotPoints();
	Array<OneD, NekDouble> twoG11(phystot);
	Array<OneD, NekDouble> G12pG21(phystot);
	Array<OneD, NekDouble> twoG22(phystot);

	Array<OneD, NekDouble> tmp1(nCoeffs);
	Array<OneD, NekDouble> tmp2(nCoeffs);

	Vmath::Smul(phystot, 2.0, m_G[0]->GetPhys(), 1, twoG11, 1);
	Vmath::Vadd(phystot, m_G[1]->GetPhys(), 1, m_G[2]->GetPhys(), 1, G12pG21, 1);
	Vmath::Smul(phystot, 2.0, m_G[3]->GetPhys(), 1, twoG22, 1);

    m_fields[m_velocity[0]]->IProductWRTDerivBase(0, twoG11, tmp1);
    m_fields[m_velocity[0]]->IProductWRTDerivBase(1, G12pG21, tmp2);
    Vmath::Vadd(nCoeffs, tmp1, 1, tmp2, 1, tmp1, 1);
    Vmath::Smul(nCoeffs,(1.-m_beta), tmp1, 1, tmp1, 1);

    Vmath::Vadd(nCoeffs, Weakoutarray[0], 1, tmp1, 1, Weakoutarray[0], 1);

    m_fields[m_velocity[1]]->IProductWRTDerivBase(0, G12pG21, tmp1);
    m_fields[m_velocity[1]]->IProductWRTDerivBase(1, twoG22, tmp2);
    Vmath::Vadd(nCoeffs, tmp1, 1, tmp2, 1, tmp1, 1);
    Vmath::Smul(nCoeffs,(1.-m_beta), tmp1, 1, tmp1, 1);

    Vmath::Vadd(nCoeffs, Weakoutarray[1], 1, tmp1, 1, Weakoutarray[1], 1);
}

void Viscoelastic::AddGGTTimesNormalToVelocityNeumannBC(Array<OneD, Array<OneD, NekDouble> > &Weakoutarray)
{
	   int i,el,n,cnt, offset, phys_offset;
	    int edge, elmtid, nq;
	    Array<OneD, NekDouble> e_outarray;
	    Array<OneD, const NekDouble> G11,G12,G21,G22;
	    StdRegions::StdExpansionSharedPtr EdgeExp;
	    StdRegions::StdExpansionSharedPtr elmt;

	    Array<OneD, int> ElmtID,EdgeID;

	    //Fills ElmtID and EdgeID with global ids of elements id number and edges id number of the boundary expansion
	    m_fields[m_velocity[0]]->GetBoundaryToElmtMap(ElmtID,EdgeID);

	    int maxpts = 0;

	     // find the maximum values of points
	     for(cnt = n = 0; n < m_fields[m_velocity[0]]->GetBndConditions().num_elements(); ++n)
	     {
	 //Loop over all elements in boundary region
	     	//cnt starts with element 0 then going through elements
	         for(i = 0; i < m_fields[m_velocity[0]]->GetBndCondExpansions()[n]->GetExpSize(); ++i)
	         {
	             maxpts = max(maxpts, m_fields[m_velocity[0]]->GetExp(ElmtID[cnt++])->GetTotPoints());
	         }
	     }

	     Array<OneD, NekDouble> twoG11edge(5*maxpts);
	     Array<OneD, NekDouble> G12pG21edge = twoG11edge + maxpts;
	     Array<OneD, NekDouble> twoG22edge = G12pG21edge + maxpts;
	     Array<OneD, NekDouble> tmp1 = twoG22edge + maxpts;
	     Array<OneD, NekDouble> tmp2 = tmp1 + maxpts;

	    // Neumann BC for velocity component u
		for(cnt = n = 0; n < m_fields[m_velocity[0]]->GetBndConditions().num_elements(); ++n)
		{
			// Waters solution
			if((m_fields[m_velocity[0]]->GetBndConditions()[n])->GetBoundaryConditionType() == SpatialDomains::eNeumann)
			{
				// loop over elements along boundary
				for(el = 0; el < m_fields[m_velocity[0]]->GetBndCondExpansions()[n]->GetExpSize(); ++el,cnt++)
				{
					elmtid = ElmtID[cnt];
					elmt   = m_fields[m_velocity[0]]->GetExp(elmtid);
					nq     = elmt->GetTotPoints();
					offset = m_fields[m_velocity[0]]->GetPhys_Offset(elmtid);

					G11 = m_G[0]->GetPhys() + offset;
					G12 = m_G[1]->GetPhys() + offset;
					G21 = m_G[2]->GetPhys() + offset;

					Vmath::Smul(nq,2.0*(-m_theta),G11,1,tmp1,1);

					Vmath::Vadd(nq,G12,1,G21,1,tmp2,1);
					Vmath::Smul(nq,(-m_theta),tmp2,1,tmp2,1);

					EdgeExp =  boost::dynamic_pointer_cast<StdRegions::StdExpansion1D> (m_fields[m_velocity[0]]->GetBndCondExpansions()[n]->GetExp(el));

					// grab edge on the boundary
					edge = EdgeID[cnt];

					elmt->GetEdgePhysVals(edge,EdgeExp,tmp1,twoG11edge);
					elmt->GetEdgePhysVals(edge,EdgeExp,tmp2,G12pG21edge);

					/* int nquad_e = EdgeExp->GetNumPoints(0);

					for(int j=0;j<nquad_e;j++)
					{
						cout << txxedge[j] << ", ";
					}
					cout << endl; */

					e_outarray = Weakoutarray[0] + m_fields[m_velocity[0]]->GetCoeff_Offset(elmtid);

					// in: physical values of txx and txy, out:weak values of forcing function
					elmt->AddEdgeNormBoundaryInt(edge,EdgeExp,twoG11edge,
																G12pG21edge,
				                                                e_outarray);
				}
			}

			cnt +=m_fields[m_velocity[0]]->GetBndCondExpansions()[n]->GetExpSize();
		}

		for(cnt = n = 0; n < m_fields[m_velocity[1]]->GetBndConditions().num_elements(); ++n)
		{
			// Waters solution
			if((m_fields[m_velocity[1]]->GetBndConditions()[n])->GetBoundaryConditionType() == SpatialDomains::eNeumann)
			{
				// loop over elements along boundary
				for(el = 0; el < m_fields[m_velocity[1]]->GetBndCondExpansions()[n]->GetExpSize(); ++el,cnt++)
				{
					elmtid = ElmtID[cnt];
					elmt   = m_fields[m_velocity[1]]->GetExp(elmtid);
					nq     = elmt->GetTotPoints();
					offset = m_fields[m_velocity[1]]->GetPhys_Offset(elmtid);

					G12 = m_G[1]->GetPhys() + offset;
					G21 = m_G[2]->GetPhys() + offset;
					G22 = m_G[3]->GetPhys() + offset;

					Vmath::Vadd(nq,G12,1,G21,1,tmp1,1);
					Vmath::Smul(nq,(-m_theta),tmp1,1,tmp1,1);

					Vmath::Smul(nq,2.0*(-m_theta),G22,1,tmp2,1);

					EdgeExp =  boost::dynamic_pointer_cast<StdRegions::StdExpansion1D> (m_fields[m_velocity[1]]->GetBndCondExpansions()[n]->GetExp(el));

					// grab edge on the boundary
					edge = EdgeID[cnt];

					elmt->GetEdgePhysVals(edge,EdgeExp,tmp1,G12pG21edge);
					elmt->GetEdgePhysVals(edge,EdgeExp,tmp2,twoG22edge);

					/* int nquad_e = EdgeExp->GetNumPoints(0);

					for(int j=0;j<nquad_e;j++)
					{
						cout << txxedge[j] << ", ";
					}
					cout << endl; */

					e_outarray = Weakoutarray[1] + m_fields[m_velocity[1]]->GetCoeff_Offset(elmtid);

					// in: physical values of txx and txy, out:weak values of forcing function
					elmt->AddEdgeNormBoundaryInt(edge,EdgeExp,G12pG21edge,
																twoG22edge,
				                                                e_outarray);
				}
			}

			cnt +=m_fields[m_velocity[1]]->GetBndCondExpansions()[n]->GetExpSize();
		}
}

void Viscoelastic::AddGDivergenceWeakForm(Array<OneD, Array<OneD, NekDouble> > &Weakoutarray)
{
	int nCoeffs = m_fields[m_velocity[0]]->GetNcoeffs();
	int phystot = m_fields[m_velocity[0]]->GetTotPoints();

	Array<OneD, NekDouble> tmp1(nCoeffs);
	Array<OneD, NekDouble> tmp2(nCoeffs);

    m_fields[m_velocity[0]]->IProductWRTDerivBase(0, m_G[0]->GetPhys(), tmp1);
    m_fields[m_velocity[0]]->IProductWRTDerivBase(1, m_G[1]->GetPhys(), tmp2);
    Vmath::Vadd(nCoeffs, tmp1, 1, tmp2, 1, tmp1, 1);
    Vmath::Smul(nCoeffs,(1.-m_beta), tmp1, 1, tmp1, 1);

    Vmath::Vadd(nCoeffs, Weakoutarray[0], 1, tmp1, 1, Weakoutarray[0], 1);

    m_fields[m_velocity[1]]->IProductWRTDerivBase(0, m_G[2]->GetPhys(), tmp1);
    m_fields[m_velocity[1]]->IProductWRTDerivBase(1, m_G[3]->GetPhys(), tmp2);
    Vmath::Vadd(nCoeffs, tmp1, 1, tmp2, 1, tmp1, 1);
    Vmath::Smul(nCoeffs,(1.-m_beta), tmp1, 1, tmp1, 1);

    Vmath::Vadd(nCoeffs, Weakoutarray[1], 1, tmp1, 1, Weakoutarray[1], 1);
}

void Viscoelastic::AddGTimesNormalToVelocityNeumannBC(Array<OneD, Array<OneD, NekDouble> > &Weakoutarray)
{
	   int i,el,n,cnt, offset, phys_offset;
	    int edge, elmtid, nq;
	    Array<OneD, NekDouble> e_outarray;
	    Array<OneD, const NekDouble> G11,G12,G21,G22;
	    StdRegions::StdExpansionSharedPtr EdgeExp;
	    StdRegions::StdExpansionSharedPtr elmt;

	    Array<OneD, int> ElmtID,EdgeID;

	    //Fills ElmtID and EdgeID with global ids of elements id number and edges id number of the boundary expansion
	    m_fields[m_velocity[0]]->GetBoundaryToElmtMap(ElmtID,EdgeID);

	    int maxpts = 0;

	     // find the maximum values of points
	     for(cnt = n = 0; n < m_fields[m_velocity[0]]->GetBndConditions().num_elements(); ++n)
	     {
	 //Loop over all elements in boundary region
	     	//cnt starts with element 0 then going through elements
	         for(i = 0; i < m_fields[m_velocity[0]]->GetBndCondExpansions()[n]->GetExpSize(); ++i)
	         {
	             maxpts = max(maxpts, m_fields[m_velocity[0]]->GetExp(ElmtID[cnt++])->GetTotPoints());
	         }
	     }

	     Array<OneD, NekDouble> Giiedge(5*maxpts);
	     Array<OneD, NekDouble> Gijedge = Giiedge + maxpts;
	     Array<OneD, NekDouble> tmp1 = Gijedge + maxpts;
	     Array<OneD, NekDouble> tmp2 = tmp1 + maxpts;

	    // Neumann BC for velocity component u
		for(cnt = n = 0; n < m_fields[m_velocity[0]]->GetBndConditions().num_elements(); ++n)
		{
			// Waters solution
			if((m_fields[m_velocity[0]]->GetBndConditions()[n])->GetBoundaryConditionType() == SpatialDomains::eNeumann)
			{
				// loop over elements along boundary
				for(el = 0; el < m_fields[m_velocity[0]]->GetBndCondExpansions()[n]->GetExpSize(); ++el,cnt++)
				{
					elmtid = ElmtID[cnt];
					elmt   = m_fields[m_velocity[0]]->GetExp(elmtid);
					nq     = elmt->GetTotPoints();
					offset = m_fields[m_velocity[0]]->GetPhys_Offset(elmtid);

					G11 = m_G[0]->GetPhys() + offset;
					G12 = m_G[1]->GetPhys() + offset;

					Vmath::Smul(nq,(-m_theta),G11,1,tmp1,1);
					Vmath::Smul(nq,(-m_theta),G12,1,tmp2,1);

					EdgeExp =  boost::dynamic_pointer_cast<StdRegions::StdExpansion1D> (m_fields[m_velocity[0]]->GetBndCondExpansions()[n]->GetExp(el));

					// grab edge on the boundary
					edge = EdgeID[cnt];

					elmt->GetEdgePhysVals(edge,EdgeExp,tmp1,Giiedge);
					elmt->GetEdgePhysVals(edge,EdgeExp,tmp2,Gijedge);

					/* int nquad_e = EdgeExp->GetNumPoints(0);

					for(int j=0;j<nquad_e;j++)
					{
						cout << txxedge[j] << ", ";
					}
					cout << endl; */

					e_outarray = Weakoutarray[0] + m_fields[m_velocity[0]]->GetCoeff_Offset(elmtid);

					// in: physical values of txx and txy, out:weak values of forcing function
					elmt->AddEdgeNormBoundaryInt(edge,EdgeExp,Giiedge,
																Gijedge,
				                                                e_outarray);
				}
			}

			cnt +=m_fields[m_velocity[0]]->GetBndCondExpansions()[n]->GetExpSize();
		}

		for(cnt = n = 0; n < m_fields[m_velocity[1]]->GetBndConditions().num_elements(); ++n)
		{
			// Waters solution
			if((m_fields[m_velocity[1]]->GetBndConditions()[n])->GetBoundaryConditionType() == SpatialDomains::eNeumann)
			{
				// loop over elements along boundary
				for(el = 0; el < m_fields[m_velocity[1]]->GetBndCondExpansions()[n]->GetExpSize(); ++el,cnt++)
				{
					elmtid = ElmtID[cnt];
					elmt   = m_fields[m_velocity[1]]->GetExp(elmtid);
					nq     = elmt->GetTotPoints();
					offset = m_fields[m_velocity[1]]->GetPhys_Offset(elmtid);

					G21 = m_G[2]->GetPhys() + offset;
					G22 = m_G[3]->GetPhys() + offset;

					Vmath::Smul(nq,(-m_theta),G21,1,tmp1,1);

					Vmath::Smul(nq,(-m_theta),G22,1,tmp2,1);

					EdgeExp =  boost::dynamic_pointer_cast<StdRegions::StdExpansion1D> (m_fields[m_velocity[1]]->GetBndCondExpansions()[n]->GetExp(el));

					// grab edge on the boundary
					edge = EdgeID[cnt];

					elmt->GetEdgePhysVals(edge,EdgeExp,tmp1,Gijedge);
					elmt->GetEdgePhysVals(edge,EdgeExp,tmp2,Giiedge);

					/* int nquad_e = EdgeExp->GetNumPoints(0);

					for(int j=0;j<nquad_e;j++)
					{
						cout << txxedge[j] << ", ";
					}
					cout << endl; */

					e_outarray = Weakoutarray[1] + m_fields[m_velocity[1]]->GetCoeff_Offset(elmtid);

					// in: physical values of txx and txy, out:weak values of forcing function
					elmt->AddEdgeNormBoundaryInt(edge,EdgeExp,Gijedge,
																Giiedge,
				                                                e_outarray);
				}
			}

			cnt +=m_fields[m_velocity[1]]->GetBndCondExpansions()[n]->GetExpSize();
		}
}

void Viscoelastic::CalcGL2Projection(Array<OneD, Array<OneD, NekDouble> > &outarray)
{
	int VelDim = m_velocity.num_elements();
	int LDim = VelDim*VelDim;
	int phystot = m_G[0]->GetTotPoints();
	int i;

	Array<OneD, Array<OneD, NekDouble> > L(LDim);

	L[0] = Array<OneD, NekDouble>(phystot*LDim);
	for(i = 1; i < (LDim); ++i)
	{
		L[i] = L[i-1] + phystot;
	}

	VelocityGradientTensor(outarray,L);

	for(i = 0; i < (LDim); ++i)
	{
		m_G[i]->FwdTrans(L[i],m_G[i]->UpdateCoeffs());
		m_G[i]->BwdTrans(m_G[i]->GetCoeffs(),m_G[i]->UpdatePhys());
	}

}

void Viscoelastic::UpdateG()
{
	int VelDim = m_velocity.num_elements();
	int LDim = VelDim*VelDim;
	int phystot = m_G[0]->GetTotPoints();
	int i;
	MultiRegions::ExpList2DSharedPtr ExpField;

	for(i = 0; i < (LDim); ++i)
	{
		if(ExpField = boost::dynamic_pointer_cast<
	        			MultiRegions::ExpList2D>(m_G[i]))
	  	{
		//TODO: UNDO THIS!!!	ExpField->UpdateExpList2D(m_session,m_graph,"txx");
	  	}
		else
		{
			cout << "Warning G is not updated" << endl;
		}
	}
}

// Result is saved in m_stressold as this is what is used for rhs of momentum equation
void Viscoelastic::SolveOlroydBPseudoSpectral(Array<OneD, Array<OneD, NekDouble> > &fieldsit,
        Array<OneD, Array<OneD, NekDouble> > &fieldsn,
        const NekDouble aii_Dt)
{
	int vstart = m_stress[0];
	int ncoeffs = m_fields[m_stress[0]]->GetNcoeffs();
	int phystot = m_fields[m_stress[0]]->GetTotPoints();
	int StressDim  = m_stress.num_elements();
	int VelDim = m_velocity.num_elements();
	int nvariables = m_fields.num_elements();
	int i;
	Array <OneD, Array<OneD, NekDouble> > rhs(m_nConvectiveFields);
	rhs[0] = Array<OneD, NekDouble>(phystot*m_nConvectiveFields,0.0);
	for(i = 1; i < m_nConvectiveFields; ++i)
	{
		rhs[i] = rhs[i-1] + phystot;
	}

	//Compute Advection with values from tau_it
	Array<OneD, Array<OneD, NekDouble> > Advfield(VelDim);
    for(i = 0; i < VelDim; ++i)
    {
	   Advfield[i] = Array<OneD, NekDouble>(phystot);
	   Vmath::Vcopy(phystot,fieldsit[m_velocity[i]],1,Advfield[i],1);
    }

    if(m_ALE)
    {
	//  ComputeVelocityMinusMeshVelocity(Advfield);
    }


	if(m_projectionType == MultiRegions::eGalerkin)
	{
		Array<OneD, NekDouble>  Deriv = Array<OneD, NekDouble> (phystot*VelDim);
        for(i = 0; i < StressDim; ++i)
        {
        	//AdvectionNonConservativeForm(Advfield,m_stressold,rhs,Deriv);
        	//Vmath::Neg(nqtot,rhs[m_stress[i]],1);
        }
	}
	else
	{
		Array<OneD, Array<OneD, NekDouble> > Weakoutarray(StressDim);
		Weakoutarray[0] = Array<OneD, NekDouble>(ncoeffs*StressDim);
		for(i = 1; i < StressDim; ++i)
		{
		   Weakoutarray[i] = Weakoutarray[i-1] + ncoeffs;
		}

		WeakDGAdvection(Advfield,m_stressold, Weakoutarray,true,true,StressDim,vstart);

		// Compute -N(tau)
		for(i = 0; i < StressDim; ++i)
		{
			Vmath::Neg(ncoeffs,Weakoutarray[i],1);
		}


		for(i = 0; i < StressDim; ++i)
		{
			 m_fields[m_stress[i]]->MultiplyByElmtInvMass(Weakoutarray[i],
												Weakoutarray[i]);
			 m_fields[m_stress[i]]->BwdTrans(Weakoutarray[i],rhs[m_stress[i]]);
		}

	}

	for(i = 0; i < m_stress.num_elements(); ++i)
	{
		Vmath::Smul(phystot,aii_Dt,rhs[m_stress[i]],1,rhs[m_stress[i]],1);
		Vmath::Vadd(phystot,rhs[m_stress[i]],1,fieldsn[m_stress[i]],1,rhs[m_stress[i]],1);
		//cout << "Rhs[" << i << "]=" << Vmath::Vmax(phystot,rhs[i],1) << endl;
	}

	CalcSourceTermImplicitOldB(rhs,fieldsit,aii_Dt);
	for(int i = 0; i < m_stress.num_elements(); ++i)
	{
	   if(m_projectionType == MultiRegions::eGalerkin)
	   {
		   m_fields[m_stress[i]]->FwdTrans(fieldsit[m_stress[i]],m_fields[m_stress[i]]->UpdateCoeffs());
	   }
	   else
	   {
		   m_fields[m_stress[i]]->FwdTrans(fieldsit[m_stress[i]],m_fields[m_stress[i]]->UpdateCoeffs());
	   }

		 m_fields[m_stress[i]]->BwdTrans(m_fields[m_stress[i]]->GetCoeffs(),m_stressold[i]);
		 m_fields[m_stress[i]]->BwdTrans(m_fields[m_stress[i]]->GetCoeffs(),m_fields[m_stress[i]]->UpdatePhys());
		 //m_fields[m_stress[i]]->SetPhys(outarray[m_stress[i]]);
		 //m_fields[m_stress[i]]->SetPhysState(true);
	}
}



} //end of namespace
