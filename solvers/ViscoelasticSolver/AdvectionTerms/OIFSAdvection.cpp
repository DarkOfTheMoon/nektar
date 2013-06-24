///////////////////////////////////////////////////////////////////////////////
//
// File OIFSAdvection.cpp
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
// Description: Evaluation of the Navier Stokes advective term
//
///////////////////////////////////////////////////////////////////////////////

#include <ViscoelasticSolver/AdvectionTerms/OIFSAdvection.h>
#include <cstdio>
#include <cstdlib>

namespace Nektar
{
	string OIFSAdvection::className  = GetAdvectionTermFactory().RegisterCreatorFunction("OIFS", OIFSAdvection::create);

	/**
	 * Constructor. Creates ...
	 *
	 * \param
	 * \param
	 */

	OIFSAdvection::OIFSAdvection(
			const LibUtilities::SessionReaderSharedPtr&        pSession,
			const SpatialDomains::MeshGraphSharedPtr&          pGraph):
		AdvectionTerm(pSession, pGraph)

	{

	   	int i,j;
		std::string velids[] = {"u","v","w"};

		m_velocityid = Array<OneD,int>(m_spacedim);
		int numfields = m_spacedim;

		for(i = 0; i <m_nConvectiveFields; ++i)
		{
			for(j = 0; j < numfields; ++j)
			{
				std::string var = m_session->GetVariable(j);
				if(m_session->NoCaseStringCompare(velids[i],var) == 0)
				{
					m_velocityid[i] = j;
					break;
				}

				if(j == numfields)
				{
					std::string error = "Failed to find field: " + var;
					ASSERTL0(false,error.c_str());
				}
			}
		}

		int VelDim     = m_velocityid.num_elements();
		// Assume all fields but last to be convected by velocity.
		m_nConvectiveFields=VelDim;
		m_fields   = Array<OneD, MultiRegions::ExpListSharedPtr>(m_nConvectiveFields);
		m_RK4itmax = 1;

	}

	OIFSAdvection::~OIFSAdvection()
	{
	}

    // Initialise OIFS time integration with starting field and according to Timeintegrationmethod for Navier-Stokes equation
    void OIFSAdvection::DoInitialise(Array<OneD, MultiRegions::ExpListSharedPtr > &pFields,
    								LibUtilities::TimeIntegrationMethod IntMethod,
    								 const Array<OneD, const Array<OneD, NekDouble> > &y_0,
    								 const NekDouble &timestep,
    								 const NekDouble &time)
    {
		m_ode.DefineOdeRhs     (&OIFSAdvection::DoOdeRhs,        this);
		m_ode.DefineProjection (&OIFSAdvection::DoOdeProjection, this);

		m_timeIntMethod = LibUtilities::eClassicalRungeKutta4;
		int i,j;
	    for(i=0;i<m_nConvectiveFields;i++)
	    {
	    	m_fields[i] = pFields[i];
	    }

		int phystot = m_fields[0]->GetTotPoints();
		int VelDim = m_velocityid.num_elements();

	    switch(IntMethod)
		{
			case LibUtilities::eIMEXOrder1:
			{
			    m_numMultiSteps = 1;

			}
			break;
			case LibUtilities::eIMEXOrder2:
			{
                m_numMultiSteps = 2;
			}
			break;
			default:
			ASSERTL0(0,"Integration method not setup: Options include ImexOrder1, ImexOrder2");
			break;
		}

	    m_velocity = Array<OneD, Array<OneD, Array<OneD, NekDouble> > > (m_numMultiSteps);
	    for(i=0;i<m_numMultiSteps;i++)
		{
			m_velocity[i] = Array<OneD, Array<OneD, NekDouble> >(VelDim);
			for(j=0;j<VelDim;j++)
			{
				m_velocity[i][j] = Array<OneD, NekDouble>(phystot);
			}
		}


	    for(i=0;i<VelDim;i++)
		{
			m_velocity[0][i] = y_0[i];
		}
        m_intScheme = Array<OneD, LibUtilities::TimeIntegrationSchemeSharedPtr>(m_numMultiSteps);
		m_utilde = Array<OneD, LibUtilities::TimeIntegrationSolutionSharedPtr>(m_numMultiSteps);


		for(i=0;i<m_numMultiSteps;i++)
		{
			// Used in the first time step to initalize the scheme
			LibUtilities::TimeIntegrationSchemeKey IntKey(m_timeIntMethod);
			m_intScheme[i] = LibUtilities::TimeIntegrationSchemeManager()[IntKey];
			// Initialise the scheme for the actual time integration scheme
			m_utilde[i] = m_intScheme[i]->InitializeScheme(timestep,y_0,time,m_ode);
		}
	    m_time =  Array<OneD, NekDouble>(m_numMultiSteps,0.0);
	    m_time[0] = time;
	    m_timestep[0] = timestep;
	    m_ncalls = 1;
    }


	//Advection function
	// Inarray should contain solution y_n and y_n-1??
	// Outarray contains y~_n, y~_n-1 to be set into time integration for navier stokes solver
	void OIFSAdvection:: v_DoAdvection(	Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
											    const Array<OneD, const Array<OneD, NekDouble> > &pInarray,
												Array<OneD, Array<OneD, NekDouble> > &pOutarray,
											    Array<OneD, NekDouble> &pWk)
	{
		ASSERTL0(false,"This advection form is not defined in this class");
	}
	//Advection function
	// Takes the solutionvector of the time integration of the navier stokes equation , and uses m_solvector[0] .. m_solvector[steps]
	// to compute RK4 solution of advection equation
	// the result of the RK4 computation is then saved in solvector[0] ... solvection[steps]
	void OIFSAdvection::DoAdvection(	Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
										LibUtilities::TimeIntegrationSolutionSharedPtr &solvector,
										const NekDouble timestep, const NekDouble time)
	{
		int i,j,k;
		int VelDim = m_velocityid.num_elements();
		int numfields = m_nConvectiveFields;
		int phystot = m_fields[0]->GetTotPoints();
		m_nsteps = 1;
		if(m_ncalls>3)
		{
			m_nsteps = m_numMultiSteps;
		}
		NekDouble dt;
		if(m_numMultiSteps>1)
		{
			m_timestep[1] = m_timestep[0];
			m_time[1] = m_time[0];
		}
		m_timestep[0] = timestep/m_RK4itmax;
		m_time[0] = time;

        /// Array holding all dependent variables.
		// TODO: Maybe this is only necessary at inititialisation as its pointers
        for(i=0;i<numfields;i++)
        {
        	m_fields[i] = pFields[i];
        }

        // Set up wrapper to fields data storage.
		Array<OneD, Array<OneD, NekDouble> >   fields(numfields);
		Array<OneD, Array<OneD, NekDouble> >  y_n;
		NekDouble    t_n;

	   /* for(i=0;i<m_nsteps;i++)
		{
	    	y_n = solvector->GetValue(i);
			for(j=0;j<VelDim;j++)
			{
				m_velocity[i][j] = y_n[j];
			}
		} */

		if(m_numMultiSteps>1)
		{
			for(j=0;j<VelDim;j++)
			{
				Vmath::Vcopy(phystot,m_velocity[0][j],1,m_velocity[1][j],1);
			}
		}

		for(j=0;j<VelDim;j++)
		{
			m_velocity[0][j] = pFields[j]->GetPhys();
		}

		for(i=0;i<m_nsteps;i++)
		{
			// Copy solution y_n/y_n-1 into m_utilde solution
			//	m_utilde[i]->UpdateSolution() = solvector->GetValue(i);
			y_n = solvector->GetValue(i);
			//m_time[i] = solvector->GetValueTime(i);

			// Set the required value in the input solution
			// vector of the current scheme
			m_utilde[i]->SetValue(0,y_n,m_time[i]);
		}


		for(k=0;k<m_RK4itmax;k++)
		{
			for(i=0;i<m_nsteps;i++)
			{
				//cout << "tn= " << m_time[0] << ", tn-1= " << m_time[1] << ", dt= " << dt << endl;
				for(j=0;j<=i;j++)
				{
					m_dt = m_timestep[j];
					fields = m_intScheme[i]->TimeIntegrate(m_timestep[j],m_utilde[i],m_ode);
					cout << "Max=(" << Vmath::Vmax(phystot,fields[0],1) << "," << Vmath::Vmax(phystot,fields[1],1) << endl;
				}
			}
		}
		// Copy results into solvector for Navier Stokes equation
		for(i=0;i<m_nsteps;i++)
		{
			solvector->SetValue(i,m_utilde[i]->GetSolution(),solvector->GetValueTime(i));
		}

		m_ncalls++;
	 }


    //Evaluation of the advective terms
    void OIFSAdvection::ComputeAdvectionTerm(
            Array<OneD, MultiRegions::ExpListSharedPtr > &pFields,
            const Array<OneD, Array<OneD, NekDouble> > &pV,
            const Array<OneD, const NekDouble> &pU,
            Array<OneD, NekDouble> &pOutarray,
            int pVelocityComponent,
            Array<OneD, NekDouble> &pWk)
    {
        // use dimension of Velocity vector to dictate dimension of operation
        int ndim       = pV.num_elements();

        // ToDo: here we should add a check that V has right dimension

        int nPointsTot = pFields[0]->GetNpoints();
        Array<OneD, NekDouble> grad0,grad1,grad2;

        grad0 = Array<OneD, NekDouble> (nPointsTot);

        // Evaluate V\cdot Grad(u)
        switch(ndim)
        {
        case 1:
            pFields[0]->PhysDeriv(pU,grad0);
            Vmath::Vmul(nPointsTot,grad0,1,pV[0],1,pOutarray,1);
            break;
        case 2:
            grad1 = Array<OneD, NekDouble> (nPointsTot);
            pFields[0]->PhysDeriv(pU,grad0,grad1);
            Vmath::Vmul (nPointsTot,grad0,1,pV[0],1,pOutarray,1);
            Vmath::Vvtvp(nPointsTot,grad1,1,pV[1],1,pOutarray,1,pOutarray,1);
            break;
        case 3:
            grad1 = Array<OneD, NekDouble> (nPointsTot);
            grad2 = Array<OneD, NekDouble> (nPointsTot);
            pFields[0]->PhysDeriv(pU,grad0,grad1,grad2);
            Vmath::Vmul( nPointsTot,grad0,1,pV[0],1,pOutarray,1);
            Vmath::Vvtvp(nPointsTot,grad1,1,pV[1],1,pOutarray,1,pOutarray,1);
            Vmath::Vvtvp(nPointsTot,grad2,1,pV[2],1,pOutarray,1,pOutarray,1);
            break;
        default:
            ASSERTL0(false,"dimension unknown");
        }
    }

void OIFSAdvection::WeakDGAdvection(
                  const Array<OneD, Array<OneD, NekDouble> >& InField,
                  Array<OneD, Array<OneD, NekDouble> >& OutField,
                  bool NumericalFluxIncludesNormal,
                  bool InFieldIsInPhysSpace)
      {
          int i;
          int nVelDim         = m_spacedim;
          int nPointsTot      = m_fields[0]->GetNpoints();
          int ncoeffs         = m_fields[0]->GetNcoeffs();
          int nTracePointsTot = m_fields[0]->GetTrace()->GetNpoints();
          int nvariables      = m_fields.num_elements();


          Array<OneD, Array<OneD, NekDouble> > fluxvector(nVelDim);
          Array<OneD, Array<OneD, NekDouble> > physfield (nvariables);

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
                  m_fields[i]->BwdTrans(InField[i],physfield[i]);
              }
          }

          // Get the advection part (without numerical flux)
          for(i = 0; i < nvariables; ++i)
          {
              // Get the ith component of the  flux vector in (physical space)
              v_GetFluxVector(i, physfield, fluxvector);

              // Calculate the i^th value of (\grad_i \phi, F)
              WeakAdvectionGreensDivergenceForm(fluxvector,OutField[i]);
          }

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
              v_NumericalFlux(physfield, numflux);

              // Evaulate  <\phi, \hat{F}\cdot n> - OutField[i]
              for(i = 0; i < nvariables; ++i)
              {
                  Vmath::Neg(ncoeffs,OutField[i],1);
                  m_fields[i]->AddTraceIntegral(numflux[i],OutField[i]);
                  m_fields[i]->SetPhysState(false);
              }
          }
          // if the NumericalFlux function does not include the
          // normal in the output
          else
          {
        	  cout << "Warning: Type of flux should still be implemented" << endl;
            /*  Array<OneD, Array<OneD, NekDouble> > numfluxX   (nvariables);
              Array<OneD, Array<OneD, NekDouble> > numfluxY   (nvariables);

              for(i = 0; i < nvariables; ++i)
              {
                  numfluxX[i]   = Array<OneD, NekDouble>(nTracePointsTot);
                  numfluxY[i]   = Array<OneD, NekDouble>(nTracePointsTot);
              }

              // Evaluate numerical flux in physical space which may in
              // general couple all component of vectors
              v_NumericalFlux(physfield, numfluxX, numfluxY);

              // Evaulate  <\phi, \hat{F}\cdot n> - OutField[i]
              for(i = 0; i < nvariables; ++i)
              {
                  Vmath::Neg(ncoeffs,OutField[i],1);
                  m_fields[i]->AddTraceIntegral(numfluxX[i], numfluxY[i],
                                                OutField[i]);
                  m_fields[i]->SetPhysState(false);
              } */
          }
      }

void OIFSAdvection::WeakAdvectionGreensDivergenceForm(
                  const Array<OneD, Array<OneD, NekDouble> > &F,
                  Array<OneD, NekDouble> &outarray)
      {
          // use dimension of Velocity vector to dictate dimension of operation
          int ndim    = F.num_elements();
          int nCoeffs = m_fields[0]->GetNcoeffs();

          Array<OneD, NekDouble> iprod(nCoeffs);
          Vmath::Zero(nCoeffs, outarray, 1);

          for (int i = 0; i < ndim; ++i)
          {
              m_fields[0]->IProductWRTDerivBase(i, F[i], iprod);
              Vmath::Vadd(nCoeffs, iprod, 1, outarray, 1, outarray, 1);
          }
      }





void OIFSAdvection::v_GetFluxVector(const int i, Array<OneD, Array<OneD, NekDouble> > &physfield,
                           Array<OneD, Array<OneD, NekDouble> > &flux)
    {
        ASSERTL1(flux.num_elements() == m_velocity.num_elements(),"Dimension of flux array and velocity array do not match");

        for(int j = 0; j < flux.num_elements(); ++j)
        {
            Vmath::Vmul(m_fields[0]->GetNpoints(),physfield[i],1,
                m_velocity[0][j],1,flux[j],1);
        }
    }

void OIFSAdvection::v_NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &physfield, Array<OneD, Array<OneD, NekDouble> > &numflux)
    {
        int i;

        int nTraceNumPoints = m_fields[0]->GetTrace()->GetNpoints();
        int nvel = m_spacedim; //m_velocity.num_elements();
        Array<OneD, Array<OneD, NekDouble> > traceNormals = Array<OneD, Array<OneD, NekDouble> >(m_spacedim);
        for(i = 0; i < m_spacedim; ++i)
		{
        	traceNormals[i] = Array<OneD, NekDouble> (m_fields[0]->GetTrace()->GetNpoints());
		}
        m_fields[0]->GetTrace()->GetNormals(traceNormals);

        Array<OneD, NekDouble > Fwd(nTraceNumPoints);
        Array<OneD, NekDouble > Bwd(nTraceNumPoints);
        Array<OneD, NekDouble > Vn (nTraceNumPoints,0.0);

        // Get Edge Velocity - Could be stored if time independent
        for(i = 0; i < nvel; ++i)
        {
            m_fields[0]->ExtractTracePhys(m_velocity[0][i], Fwd);
            Vmath::Vvtvp(nTraceNumPoints,traceNormals[i],1,Fwd,1,Vn,1,Vn,1);
        }

        for(i = 0; i < numflux.num_elements(); ++i)
        {
            m_fields[i]->GetFwdBwdTracePhys(physfield[i],Fwd,Bwd);
            //evaulate upwinded m_fields[i]
            m_fields[i]->GetTrace()->Upwind(Vn,Fwd,Bwd,numflux[i]);
            // calculate m_fields[i]*Vn
            Vmath::Vmul(nTraceNumPoints,numflux[i],1,Vn,1,numflux[i],1);
        }
    }

    void OIFSAdvection::DoOdeRhs(	const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
                                     Array<OneD,        Array<OneD, NekDouble> >&outarray,
                                     const NekDouble time)
    {
        int i,j;
        int nvariables = inarray.num_elements();
        int npoints =  m_fields[0]->GetNpoints();

        switch (m_projectionType)
        {
        case MultiRegions::eDiscontinuousGalerkin:
            {
                int ncoeffs    = inarray[0].num_elements();
                Array<OneD, Array<OneD, NekDouble> > WeakAdv(nvariables);

                WeakAdv[0] = Array<OneD, NekDouble>(ncoeffs*nvariables);
                for(i = 1; i < nvariables; ++i)
                {
                    WeakAdv[i] = WeakAdv[i-1] + ncoeffs;
                }

                WeakDGAdvection(inarray, WeakAdv,true,true);

                for(i = 0; i < nvariables; ++i)
                {
                    m_fields[i]->MultiplyByElmtInvMass(WeakAdv[i],
                                                       WeakAdv[i]);
                    m_fields[i]->BwdTrans(WeakAdv[i],outarray[i]);
                    Vmath::Neg(npoints,outarray[i],1);
                }

                break;
            }
            case MultiRegions::eGalerkin:
            {
            	int nqtot      = m_fields[0]->GetTotPoints();
                int VelDim     = m_velocityid.num_elements();
        		Array<OneD, Array<OneD, NekDouble> > velocity(VelDim);
        		for(i = 0; i < VelDim; ++i)
				{
        			velocity[i] = Array<OneD, NekDouble>(nqtot,0.0);
				}

        		Array<OneD, NekDouble > Deriv;
        		int k = m_nsteps-1;
                for(i = 0; i < VelDim; ++i)
                {
                	if(k==0)
                	{
                		Vmath::Vcopy(nqtot,m_velocity[0][i],1,velocity[i],1);
                	}
                	else
                	{
                		for(j=0;j<nqtot;j++)
                		{
                			velocity[i][j] = ((time - m_time[1])/m_dt)*m_velocity[0][i][j] + (1-((time - m_time[1])/m_dt))*m_velocity[1][i][j];
                		}
                	}
                }

              //  cout << "t= " << time << ", (u,v)=( " << Vmath::Vmax(nqtot,velocity[0],1) << ", " << Vmath::Vmax(nqtot,velocity[1],1) << ")" << endl;

                Deriv = Array<OneD, NekDouble> (nqtot*VelDim);

                // Calculate -V\cdot Grad(u);
                for(i = 0; i < nvariables; ++i)
                {
                	ComputeAdvectionTerm(m_fields,velocity,inarray[i],outarray[i],i,Deriv);
                    Vmath::Neg(npoints,outarray[i],1);
                }
                break;
            }
        }
    }



    /**
     *
     */
    void OIFSAdvection::DoOdeProjection(const Array<OneD,
                                            const Array<OneD, NekDouble> >&inarray,
                                            Array<OneD,       Array<OneD, NekDouble> >&outarray,
                                            const NekDouble time)
    {
        int i;
        int nvariables = inarray.num_elements();
       // SetBoundaryConditions(time);

        switch(m_projectionType)
        {
        case MultiRegions::eDiscontinuousGalerkin:
            {
                // Just copy over array
                int npoints =  m_fields[0]->GetNpoints();

                for(i = 0; i < nvariables; ++i)
                {
                    Vmath::Vcopy(npoints,inarray[i],1,outarray[i],1);
                }
            }
            break;
        case MultiRegions::eGalerkin:
            {
                Array<OneD, NekDouble> coeffs(m_fields[0]->GetNcoeffs());
                int npoints =  m_fields[0]->GetNpoints();

                for(i = 0; i < nvariables; ++i)
                {
                    //m_fields[i]->FwdTrans(inarray[i],coeffs,false);
                    //m_fields[i]->BwdTrans_IterPerExp(coeffs,outarray[i]);
                	Vmath::Vcopy(npoints,inarray[i],1,outarray[i],1);
                }
                break;
            }
        default:
            ASSERTL0(false,"Unknown projection scheme");
            break;
        }
    }

    /*
       * Calculate weak DG advection in the form
       * \f$ \langle\phi, \hat{F}\cdot n\rangle - (\nabla \phi \cdot F) \f$
       * @param   InField         Fields.
       * @param   OutField        Storage for result.
       * @param   NumericalFluxIncludesNormal     Default: true.
       * @param   InFieldIsPhysSpace              Default: false.
       * @param   nvariables      Number of fields.
       */
            /**
       * Computes the weak Green form of advection terms (without boundary
       * integral), i.e. \f$ (\nabla \phi \cdot F) \f$ where for example
       * \f$ F=uV \f$.
       * @param   F           Fields.
       * @param   outarray    Storage for result.
       *
       * \note Assuming all fields are of the same expansion and order so that we
       * can use the parameters of m_fields[0].
       */
              } //end of namespace




