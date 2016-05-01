///////////////////////////////////////////////////////////////////////////////
//
// File: ForcingMovingBody.cpp
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
// Description: Moving Body m_forcing (movement of a body in a domain is achieved
// via a m_forcing term which is the results of a coordinate system motion)
//
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/Forcing/ForcingMovingBody.h>
#include <MultiRegions/ExpList.h>
#include <IncNavierStokesSolver/Forcing/spline.h>

namespace Nektar
{
std::string ForcingMovingBody::className = SolverUtils::GetForcingFactory().
            RegisterCreatorFunction("MovingBody",
                                    ForcingMovingBody::create,
                                    "Moving Body Forcing");

ForcingMovingBody::ForcingMovingBody(
        const LibUtilities::SessionReaderSharedPtr& pSession)
    : Forcing(pSession)
{
}

void ForcingMovingBody::v_InitObject(
        const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
        const unsigned int& pNumForcingFields,
        const TiXmlElement* pForce)
{
    // Just 3D homogenous 1D problems can use this techinque
    ASSERTL0(pFields[0]->GetExpType()==MultiRegions::e3DH1D,
             "Moving body implemented just for 3D Homogenous 1D expansions.");

    // At this point we know in the xml file where those quantities
    // are declared (equation or file) - via a function name which is now
    // stored in funcNameD etc. We need now to fill in with this info the
    // m_zta and m_eta vectors (actuallythey are matrices) Array to control if
    // the motion is determined by an equation or is from a file.(not Nektar++)
    // check if we need to load a file or we have an equation
    CheckIsFromFile(pForce);

    // Initialise movingbody filter
    InitialiseFilter(pFields, pForce);

    // Initialise the cable model
    InitialiseCableModel(pFields);

    // Load mapping
    m_mapping = GlobalMapping::Mapping::Load(m_session, pFields);
    m_mapping->SetTimeDependent( true );

    if(m_vdim > 0)
    {
        m_mapping->SetFromFunction( false );
    }

    m_zta = Array<OneD, Array< OneD, NekDouble> > (3);
    m_eta = Array<OneD, Array< OneD, NekDouble> > (3);
    // What are this bi-dimensional vectors ------------------------------------------
    // m_zta[0] = m_zta                     |  m_eta[0] = m_eta                      |
    // m_zta[1] = d(m_zta)/dt               |  m_eta[1] = d(m_eta)/dt                |
    // m_zta[2] = dd(m_zta)/ddtt            |  m_eta[2] = dd(m_eta)/ddtt             |
    //--------------------------------------------------------------------------------
    int phystot = pFields[0]->GetTotPoints();
    for(int i = 0; i < m_zta.num_elements(); i++)
    {
        m_zta[i] = Array<OneD, NekDouble>(phystot,0.0);
        m_eta[i] = Array<OneD, NekDouble>(phystot,0.0);
    }
}

void ForcingMovingBody::v_Apply(
        const Array<OneD, MultiRegions::ExpListSharedPtr>&  pFields,
        const Array<OneD, Array<OneD, NekDouble> >&         inarray,
              Array<OneD, Array<OneD, NekDouble> >&         outarray,
        const NekDouble&                                    time)
{
    m_movingBodyCalls++;

    Array<OneD, Array<OneD, NekDouble> > HydroForces (2);
	HydroForces[0] = Array<OneD, NekDouble>(m_np,0.0);
	HydroForces[1] = Array<OneD, NekDouble>(m_np,0.0); 

	// Update teh fluid forces
    m_MovBodyfilter->UpdateForce(pFields, HydroForces, time);
	
    // for "free" type (m_vdim = 2), the cable vibrates both in streamwise and crossflow
    // directions, for "constrained" type (m_vdim = 1), the cable only vibrates in crossflow
    // direction, and for "forced" type (m_vdim = 0), the calbe vibrates specifically along
    // a given function or file

    if(m_vdim == 1 || m_vdim == 2)
    {
        // For free vibration case, displacements, velocities and acceleartions
        // are obtained through solving structure dynamic model
        StructDynModel(pFields, HydroForces, time);
        // Convert result to format required by mapping
        int physTot = pFields[0]->GetTotPoints();
        Array< OneD, Array< OneD, NekDouble> >  coords(3);
        Array< OneD, Array< OneD, NekDouble> >  coordsVel(3);

        for(int i = 0; i < 3; i++)
        {
            coords[i] = Array< OneD, NekDouble> (physTot, 0.0);
            coordsVel[i] = Array< OneD, NekDouble> (physTot, 0.0);
        }

        // Get original coordinates
        pFields[0]->GetCoords(coords[0], coords[1], coords[2]);

        // Add displacement to coordinates
        Vmath::Vadd(physTot, coords[0], 1, m_zta[0], 1, coords[0], 1);
        Vmath::Vadd(physTot, coords[1], 1, m_eta[0], 1, coords[1], 1);
        // Copy velocities
        Vmath::Vcopy(physTot, m_zta[1], 1, coordsVel[0], 1);
        Vmath::Vcopy(physTot, m_eta[1], 1, coordsVel[1], 1);
        
        // Update mapping
        m_mapping->UpdateMapping(time, coords, coordsVel);
    }
    else if(m_vdim == 0)
    {
        // For forced vibration case, load from specified file or function
        int cnt = 0;
        for(int j = 0; j < m_funcName.num_elements(); j++)
        {
            if(m_IsFromFile[cnt] && m_IsFromFile[cnt+1])
            {
                ASSERTL0(false, "Motion loading from file needs specific "
                                "implementation: Work in Progress!");
            }
            else
            {
                EvaluateFunction(pFields, m_session, m_motion[0], m_zta[j],
                                 m_funcName[j], time);
                EvaluateFunction(pFields, m_session, m_motion[1], m_eta[j],
                                 m_funcName[j], time);
                cnt = cnt + 2;
            }
        }
        
        // Update mapping
        m_mapping->UpdateMapping(time);

        // Convert result from mapping       
        int physTot = pFields[0]->GetTotPoints();
        Array< OneD, Array< OneD, NekDouble> >  coords(3);
        Array< OneD, Array< OneD, NekDouble> >  coordsVel(3);
        for(int i =0; i<3; i++)
        {
            coords[i] = Array< OneD, NekDouble> (physTot, 0.0);
            coordsVel[i] = Array< OneD, NekDouble> (physTot, 0.0);
        }
        // Get original coordinates
        pFields[0]->GetCoords(coords[0], coords[1], coords[2]);

        // Get Coordinates and coord velocity from mapping
        m_mapping->GetCartesianCoordinates(m_zta[0], m_eta[0], coords[2]);
        m_mapping->GetCoordVelocity(coordsVel);

        // Calculate displacement
        Vmath::Vsub(physTot, m_zta[0], 1, coords[0], 1, m_zta[0], 1);
        Vmath::Vsub(physTot, m_eta[0], 1, coords[1], 1, m_eta[0], 1);

        Vmath::Vcopy(physTot, coordsVel[0], 1, m_zta[1], 1);
        Vmath::Vcopy(physTot, coordsVel[1], 1, m_eta[1], 1);
    }
    else
    {
        ASSERTL0(false, 
                 "Unrecogized vibration type for cables");
    }

	// write the loading and responses into the files
    m_MovBodyfilter->WriteForces(pFields, m_q  , time);
    m_MovBodyfilter->WriteMotion(pFields, m_vib, time);
}

/**
 *
 */
void ForcingMovingBody::StructDynModel(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        	  Array<OneD, Array<OneD, NekDouble> > &HydroForces,
              NekDouble time)
{
    LibUtilities::CommSharedPtr vcomm = pFields[0]->GetComm();
    int colrank = vcomm->GetColumnComm()->GetRank();
    int nproc   = vcomm->GetColumnComm()->GetSize();

    int npts, nstrips, ntotpoints;
    m_session->LoadParameter("HomModesZ", npts);
    m_session->LoadParameter("Strip_Z", nstrips, 1);

	ntotpoints = nstrips*npts;

	Array<OneD, Array<OneD, NekDouble> > fces(2);
	fces[0] =  Array<OneD, NekDouble>(ntotpoints,0.0);
	fces[1] =  Array<OneD, NekDouble>(ntotpoints,0.0);

    Array<OneD, Array<OneD, NekDouble> > bnddis(2);
    Array<OneD, Array<OneD, NekDouble> > bndvel(2);
	for (int i = 0; i < 2; i++)
	{
		bnddis[i] = Array<OneD, NekDouble>(m_np);
		bndvel[i] = Array<OneD, NekDouble>(m_np);
	}
	
    //fill the force vectors
	if(colrank != 0)
	{
		vcomm->GetColumnComm()->Send(0, HydroForces[0]);
		vcomm->GetColumnComm()->Send(0, HydroForces[1]);
	}
	else
	{
		for(int n = 0; n < nstrips; n++)
		{
			for(int i = 0; i < nproc/nstrips; i++)
			{
                int id = n+i*nstrips;
                if(id != 0 && colrank == id)
                {
                    int offset = (i+n*nproc/nstrips)*m_np;
                    for(int j = 0; j < 2; j++)
                    {
                        Array<OneD, NekDouble> tmp (m_np);
                        Array<OneD, NekDouble> tmp0(m_np);
                        Array<OneD, NekDouble> tmp1(m_np);
                        vcomm->GetColumnComm()->Recv(id, tmp0);
                        vcomm->GetColumnComm()->Recv(id, tmp1);
                		Vmath::Vcopy(m_np, tmp0, 1, tmp = fces[0]+offset, 1);
                		Vmath::Vcopy(m_np, tmp1, 1, tmp = fces[1]+offset, 1);
                    }
                }
			}	
		}	
		for(int j = 0; j < 2; j++)
		{
			Vmath::Vcopy(m_np, HydroForces[j], 1, fces[j], 1);	
		}
	}	
	
    if(colrank == 0)
    {
		int ns;
		m_session->LoadParameter("HomStructModesZ", ns);


    	// stores the vibration varialbes
    	Array<OneD, Array<OneD, Array<OneD, NekDouble> > > vibvars (2);
    	for(int i = 0; i < 2; i++)
    	{
        	vibvars[i] = Array<OneD, Array<OneD, NekDouble> > (3);
        	for(int j = 0; j < 3; j++)
        	{
            	vibvars[i][j] = Array<OneD, NekDouble>(ntotpoints,0.0);
        	}	
    	}

        for(int n = 0, cn = 1; n < m_vdim; n++, cn--)
        {  
			if(nstrips > 1)
			{
				int cnt = m_np*nproc/nstrips;
				for(int i = 0; i < ns-1; i++)
				{
					for(int j = 0; j < cnt; j++)
					{
						m_q[cn][i] += fces[cn][i*cnt+j]/cnt;
					}
				} 
				m_q[cn][ns-1] = m_q[cn][0];
			}
			else
			{
				int cnt = npts/(ns-1); 
				if(cnt > 1)
				{
					for(int i = 1; i < ns-1; i++)
					{
						for(int j = 0; j < cnt+1; j++)
						{
							m_q[cn][i] += fces[cn][i*cnt+j-cnt/2]/(cnt+1);
						}
					}

					for(int j = 0; j < cnt/2+1; j++)
					{
						m_q[cn][0] += fces[cn][j]/(cnt/2+1);
					}
				}
				else
				{
					Vmath::Vcopy(ns-1,fces[cn],1,m_q[cn],1);
					m_q[cn][ns-1] = m_q[cn][0];
				}
			}
				
        	// Fictitious mass method used to 
			// stablize the explicit coupling at
        	// relatively lower mass ratio
        	bool fictmass;
        	m_session->MatchSolverInfo(
				"FictitiousMassMethod", "True", fictmass, false);
        	if(fictmass)
        	{
            	NekDouble fictrho, fictdamp;
            	m_session->LoadParameter("FictMass", fictrho);
            	m_session->LoadParameter("FictDamp", fictdamp);

            	static NekDouble Betaq_Coeffs[2][2] = 
                                {{1.0,  0.0},{2.0, -1.0}};

            	// only consider second order approximation 
				// for fictitious variables
            	int  intOrder= 2;
            	int  nint    = min(m_movingBodyCalls,intOrder);
            	int  nlevels = m_fV[0].num_elements();

                RollOver(m_fV[cn]);
                RollOver(m_fA[cn]);

                Vmath::Vcopy(ns-1, m_vib[cn][1], 1, m_fV[cn][0], 1);
                Vmath::Vcopy(ns-1, m_vib[cn][2], 1, m_fA[cn][0], 1);

                // Extrapolate to n+1
                Vmath::Smul(ns-1, 
                            Betaq_Coeffs[nint-1][nint-1],
                            m_fV[cn][nint-1],    1,
                            m_fV[cn][nlevels-1], 1);
                Vmath::Smul(ns-1, 
                            Betaq_Coeffs[nint-1][nint-1],
                            m_fA[cn][nint-1],    1,
                            m_fA[cn][nlevels-1], 1);

                for(int n = 0; n < nint-1; ++n)
                {
                    Vmath::Svtvp(ns-1, 
                                 Betaq_Coeffs[nint-1][n],
                                 m_fV[cn][n],1,m_fV[cn][nlevels-1],1,
                                 m_fV[cn][nlevels-1],1);
                    Vmath::Svtvp(ns-1, 
                                 Betaq_Coeffs[nint-1][n],
                                 m_fA[cn][n],1,m_fA[cn][nlevels-1],1,
                                 m_fA[cn][nlevels-1],1);
                }

                // Add the fictitious forces on the RHS of the equation
                Vmath::Svtvp(ns-1, fictdamp,m_fV[cn][nlevels-1],1,
                             m_q[cn],1,m_q[cn],1);
                Vmath::Svtvp(ns-1, fictrho, m_vib[cn][2],1,
                             m_q[cn],1,m_q[cn],1);
            }
        
			std::string solver = m_session->GetSolverInfo("StructDynSolver");
			if(boost::iequals(solver, "ModalDecompostion"))
			{
				ModalDecompositionMethod(m_q[cn], m_vib[cn]);
            }
			else if(boost::iequals(solver, "FiniteElementMethod"))
			{
				FiniteElementMethod(m_q[cn], m_vib[cn], m_rot[cn]);
			}

			//cubic spline interpolation
			std::vector<double> Xp(ns);
			std::vector<double> tmp(ns);
			for(int i = 0; i < ns; i++)
    		{
        		Xp[i] = m_length*i/(ns-1);
    		}
	
            tk::spline s;
    		for(int n = 0; n < 3; n++)
			{
				for(int i = 0; i < ns; i++)
    			{
					tmp[i] = m_vib[cn][n][i];
    			}
	
            	s.set_points(Xp,tmp);
	
				for(int i = 0; i < npts; i++)
				{
					if(nstrips > 1)
					{
						for(int j = 0; j < nstrips; j++)
						{
							NekDouble x = m_length*j/nstrips;
							vibvars[cn][n][i+j*npts] = s(x);
						}
					}
					else
					{
						NekDouble x = m_length*i/npts;
						vibvars[cn][n][i] = s(x);
					}
				}
            }
        }

        for(int n = 0; n < nstrips; n++)
        {
            for(int i = 0; i < nproc/nstrips; i++)
            {
                int id = n+i*nstrips;
                if(id != 0)
                {
                    int offset = (i+n*nproc/nstrips)*m_np;
                    for(int j = 0; j < 2; j++)
                    {
                        Array<OneD, NekDouble> tmp0(m_np);
                        Array<OneD, NekDouble> tmp1(m_np);
                        Vmath::Vcopy(m_np, vibvars[j][0]+offset, 1, tmp0, 1);
                        Vmath::Vcopy(m_np, vibvars[j][1]+offset, 1, tmp1, 1);
                        vcomm->GetColumnComm()->Send(id, tmp0);
                        vcomm->GetColumnComm()->Send(id, tmp1);
                    }
                }
            }
        }
		
        for(int j = 0; j < 2; j++)
        {
            Vmath::Vcopy(m_np, vibvars[j][0], 1, bnddis[j], 1);
            Vmath::Vcopy(m_np, vibvars[j][1], 1, bndvel[j], 1);
        }
    }


    // get the boundary displacements and velocities
	// for fluid solver
	if(colrank != 0)
	{
		for(int j = 0; j < 2; j++)
		{
			vcomm->GetColumnComm()->Recv(0, bnddis[j]);
			vcomm->GetColumnComm()->Recv(0, bndvel[j]);
		}
	}

    // Set the m_forcing term
    for(int p = 0; p < m_np; p++)
    {
        Array<OneD, NekDouble> tmp;
    	int n = pFields[0]->GetPlane(p)->GetTotPoints();
        int offset  = p * n;
        Vmath::Fill(n, bnddis[0][p], tmp = m_zta[0] + offset, 1);
        Vmath::Fill(n, bnddis[1][p], tmp = m_eta[0] + offset, 1);
        Vmath::Fill(n, bndvel[0][p], tmp = m_zta[1] + offset, 1);
        Vmath::Fill(n, bndvel[1][p], tmp = m_eta[1] + offset, 1);
    }
}
 
/**
 *
 */
void ForcingMovingBody::ModalDecompositionMethod(
        const Array<OneD, NekDouble> &HydroForces,
              Array<OneD, Array<OneD, NekDouble> > &vibrations)
{  
    std::string supptype = m_session->GetSolverInfo("SupportType");

    int npts = HydroForces.num_elements()-1;
        
    Array<OneD, Array<OneD, NekDouble> > fft_i(4);
    Array<OneD, Array<OneD, NekDouble> > fft_o(4);

    for(int i = 0; i < 4; i++)
    {
        fft_i[i] = Array<OneD, NekDouble>(npts, 0.0);
        fft_o[i] = Array<OneD, NekDouble>(npts, 0.0);
    }

    Vmath::Vcopy(npts, HydroForces, 1, fft_i[0], 1);
	for(int i = 0; i < 3; i++)
	{
    	Vmath::Vcopy(npts, vibrations[i], 1, fft_i[i+1], 1);
    }

    // Implement Fourier transformation of the motion variables
    if(boost::iequals(supptype, "Free-Free"))
    {
        for(int j = 0 ; j < 4; ++j)
        {
            m_FFT->FFTFwdTrans(fft_i[j], fft_o[j]);
        }
    }
    else if(boost::iequals(supptype, "Pinned-Pinned"))
    {
        //TODO:
        int N = fft_i[0].num_elements();

        for(int var = 0; var < 4; var++)
        {
            for(int i = 0; i < N; i++)
            {
                fft_o[var][i] = 0;

                for(int k = 0; k < N; k++)
                {
                    fft_o[var][i] +=
                        fft_i[var][k]*
						sin(M_PI/N*k*(i+1));
                }
            }
        }
    }
    else
    {
        ASSERTL0(false,
                    "Unrecognized support type for cable's motion");
    }

    // solve the ODE in the wave space
    for(int i = 0; i < npts; ++i)
    {
        int nrows = 3;

        Array<OneD, NekDouble> tmp0,tmp1,tmp2;
        tmp0 = Array<OneD, NekDouble> (3,0.0);
        tmp1 = Array<OneD, NekDouble> (3,0.0);
        tmp2 = Array<OneD, NekDouble> (3,0.0);

        for(int var = 0; var < 3; var++)
        {
            tmp0[var] = fft_o[var+1][i];
        }

        tmp2[0] = fft_o[0][i];

        Blas::Dgemv('N', nrows, nrows, 1.0,
                    &(m_CoeffMat_B[i]->GetPtr())[0],
                    nrows, &tmp0[0], 1,
                    0.0,   &tmp1[0], 1);
        Blas::Dgemv('N', nrows, nrows, 1.0/m_structrho,
                    &(m_CoeffMat_A[i]->GetPtr())[0],
                    nrows, &tmp2[0], 1,
                    1.0,   &tmp1[0], 1);

        for(int var = 0; var < 3; var++)
        {
            fft_i[var][i] = tmp1[var];
        }
    }

    // get physical coeffients via Backward fourier transformation of wave
    // coefficients
    if(boost::iequals(supptype, "Free-Free"))
    {
        for(int var = 0; var < 3; var++)
        {
            m_FFT->FFTBwdTrans(fft_i[var], fft_o[var]);
        }
    }
    else if(boost::iequals(supptype, "Pinned-Pinned"))
    {
        //TODO:
        int N = fft_i[0].num_elements();

        for(int var = 0; var < 3; var++)
        {
            for(int i = 0; i < N; i++)
            {
                fft_o[var][i] = 0;

                for(int k = 0; k < N; k++)
                {
                    fft_o[var][i] +=
                        fft_i[var][k]*
						sin(M_PI/N*i*(k+1))*2/N;
                }
            }
        }
    }
    else
    {
        ASSERTL0(false,
        	"Unrecognized support type for cable's motion");
    }
    
	for(int i = 0; i < 3; i++)
	{
    	Vmath::Vcopy(npts, fft_o[i], 1, vibrations[i], 1);
	}
}

/**
 *
 */
void ForcingMovingBody::FiniteElementMethod(
			  Array<OneD, NekDouble> &HydroForces,
              Array<OneD, Array<OneD, NekDouble> > &vibs,
              Array<OneD, Array<OneD, NekDouble> > &rots)
{
	NekDouble delta, alfa;
    NekDouble c0,c1,c2,c3,c4,c5,c6,c7;
    delta = 0.5;
    alfa  = 0.25;
    c0 = 1.0/(alfa*m_timestep*m_timestep);
    c1 = delta/(alfa*m_timestep);
    c2 = 1.0/(alfa*m_timestep);
    c3 = 1.0/(2.0*alfa)-1.0;
    c4 = delta/alfa-1.0;
    c5 = m_timestep/2.0*(delta/alfa-2.0);
    c6 = m_timestep*(1.0-delta);
    c7 = delta*m_timestep;

	int npts;
	m_session->LoadParameter("HomStructModesZ", npts);

	NekDouble gauss[2],eg[2];
	NekDouble N0[2],N1[2],N2[2],N3[2];	
	NekDouble Hyf[2];
	NekDouble eq[4];
	NekDouble elength;

    Array<OneD, NekDouble>    Q(2*npts,0.0);
    Array<OneD, NekDouble> tmp (2*npts,0.0);
    Array<OneD, NekDouble> tmp0(2*npts,0.0);
    Array<OneD, NekDouble> tmp1(2*npts,0.0);
    Array<OneD, NekDouble> tmp2(2*npts,0.0);

	std::vector<double> Xp(npts), Yp(npts);

	elength = m_length/(npts-1);

	for(int i = 0; i < npts; i++)
	{
		Xp[i] = i*elength;
		Yp[i] = HydroForces[i];
	}

    tk::spline s;
    s.set_points(Xp,Yp); 
   
	gauss[0] =-0.57735026919;
	gauss[1] = 0.57735026919;

    eg[0] = gauss[0]*0.5+0.5;
    eg[1] = gauss[1]*0.5+0.5;

	for(int i = 0; i < 2; i++)
	{
		N0[i] = 1.0-3.0*eg[i]*eg[i]+2.0*eg[i]*eg[i]*eg[i];
		N1[i] = (eg[i]-2.0*eg[i]*eg[i]+eg[i]*eg[i]*eg[i])*elength;
		N2[i] = 3.0*eg[i]*eg[i]-2.0*eg[i]*eg[i]*eg[i];
		N3[i] = (-eg[i]*eg[i]+eg[i]*eg[i]*eg[i])*elength;
	}

	//assemble forces in right hand side
	for(int i = 0; i < npts-1; i++)
	{
		NekDouble x0,x1;
		x0 = (eg[0]+i)*elength;
		x1 = (eg[1]+i)*elength;

		Hyf[0] = s(x0);
		Hyf[1] = s(x1);

		eq[0] = 0.5*elength*(Hyf[0]*N0[0]+Hyf[1]*N0[1]);
    	eq[1] = 0.5*elength*(Hyf[0]*N1[0]+Hyf[1]*N1[1]);
    	eq[2] = 0.5*elength*(Hyf[0]*N2[0]+Hyf[1]*N2[1]);
    	eq[3] = 0.5*elength*(Hyf[0]*N3[0]+Hyf[1]*N3[1]);

		for(int j = 0; j < 4; j++)
		{
			Q[2*i+j] = Q[2*i+j] + eq[j];
		}
	}
	
	Array<OneD, NekDouble>    motions (2*npts);
	Array<OneD, NekDouble>  d_motions (2*npts);
	Array<OneD, NekDouble> dd_motions (2*npts);
	for(int i = 0; i < npts; i++)
	{
		   motions[2*i]   = vibs[0][i];
		   motions[2*i+1] = rots[0][i];
		 d_motions[2*i]   = vibs[1][i];
		 d_motions[2*i+1] = rots[1][i];
		dd_motions[2*i]   = vibs[2][i];
		dd_motions[2*i+1] = rots[2][i];
	}	

    Vmath::Smul (2*npts, c0, motions, 1, tmp, 1);
    Vmath::Svtvp(2*npts, c2, d_motions, 1, tmp, 1, tmp, 1);
    Vmath::Svtvp(2*npts, c3, dd_motions, 1, tmp, 1, tmp, 1);

    Blas::Dgemv('N', 2*npts, 2*npts,
                1.0, &(m_CoeffMat_M->GetPtr())[0],
                2*npts, &tmp [0],1,
                0.0,    &tmp0[0],1);

    Vmath::Vadd(2*npts, Q, 1, tmp0, 1, Q, 1);

    Vmath::Smul (2*npts, c1, motions, 1, tmp, 1);
    Vmath::Svtvp(2*npts, c4, d_motions, 1, tmp, 1, tmp, 1);
    Vmath::Svtvp(2*npts, c5, dd_motions, 1, tmp, 1, tmp, 1);

    Blas::Dgemv('N', 2*npts, 2*npts,
                1.0, &(m_CoeffMat_D->GetPtr())[0],
                2*npts, &tmp [0],1,
                0.0,    &tmp0[0],1);

    Vmath::Vadd(2*npts, Q, 1, tmp0, 1, Q, 1);

	//boundary condition
	Q[0] = 0.0;
	Q[2*npts-2] = 0.0;

	//tmp0 = a_(t+dt);
    Blas::Dgemv('N', 2*npts, 2*npts,
                1.0, &(m_CoeffMat_K->GetPtr())[0],
                2*npts, &Q[0], 1,
                0.0,  &tmp0[0],1);

	//tmp1 = dda_(t+dt);
	Vmath::Vsub(2*npts, tmp0, 1, motions, 1, tmp1, 1);
	Vmath::Smul(2*npts, c0, tmp1, 1, tmp1, 1);
	Vmath::Svtvp(2*npts, -c2, d_motions, 1, tmp1, 1, tmp1, 1);
	Vmath::Svtvp(2*npts, -c3, dd_motions, 1, tmp1, 1, tmp1, 1);

	//tmp2 = da_(t+dt); 
	Vmath::Smul(2*npts, c6, dd_motions, 1, tmp2, 1);
	Vmath::Svtvp(2*npts, c7, tmp1, 1, tmp2, 1, tmp2, 1);
    Vmath::Vadd(2*npts, d_motions, 1, tmp2, 1, tmp2, 1);

	Vmath::Vcopy(2*npts, tmp0, 1, motions, 1);
	Vmath::Vcopy(2*npts, tmp2, 1, d_motions, 1);
	Vmath::Vcopy(2*npts, tmp1, 1, dd_motions, 1);

	for(int i = 0; i < npts; i++)
	{
		vibs[0][i] = tmp0[2*i];
		rots[0][i] = tmp0[2*i+1];
		vibs[1][i] = tmp2[2*i];
		rots[1][i] = tmp2[2*i+1];
		vibs[2][i] = tmp1[2*i];
		rots[2][i] = tmp1[2*i+1];
	}	
}
  
/**
 *
 */
void ForcingMovingBody::InitialiseCableModel(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields)
{
    m_movingBodyCalls = 0;

    LibUtilities::CommSharedPtr vcomm = pFields[0]->GetComm();
    int colrank = vcomm->GetColumnComm()->GetRank();

    // number of vibration modes
    int ns;
	// number of fourier modes in fluid solver
	int npts;
    m_session->LoadParameter("HomStructModesZ", ns);
	m_session->LoadParameter("HomModesZ", npts);
    m_session->LoadParameter("LC", m_length);
    m_session->LoadParameter("TimeStep", m_timestep, 0.01);
    m_FFT = LibUtilities::GetNektarFFTFactory().CreateInstance(
                                            "NekFFTW", npts);
  
	// distributed loading on the cable
	m_q = Array<OneD, Array<OneD, NekDouble> >(2);
	m_q[0] = Array<OneD, NekDouble>(ns,0.0);
	m_q[1] = Array<OneD, NekDouble>(ns,0.0);
  
	// responeses of the cable
	m_vib = Array<OneD, Array<OneD, Array<OneD, NekDouble> > >(2);
	for (int i = 0; i < m_motion.num_elements(); i++)
	{
		m_vib[i] = Array<OneD, Array<OneD, NekDouble> > (3);
		for (int j = 0; j < 3; j++)
		{
			m_vib[i][j] = Array<OneD, NekDouble>(ns, 0.0);	
		}
	}

    Array<OneD, unsigned int> ZIDs;
    ZIDs = pFields[0]->GetZIDs();
    m_np = ZIDs.num_elements();

    std::string vibtype = m_session->GetSolverInfo("VibrationType");

    if(boost::iequals(vibtype, "Constrained"))
    {
        m_vdim = 1;
    }
    else if (boost::iequals(vibtype, "Free"))
    {
        m_vdim = 2;
    }
    else if (boost::iequals(vibtype, "Forced"))
    {
        m_vdim = 0;
    }

    // load the structural dynamic parameters from xml file
    m_session->LoadParameter("StructRho",  m_structrho);
    m_session->LoadParameter("StructDamp", m_structdamp, 0.0);

    // Identify whether the fictitious mass method is active for explicit
    // coupling of fluid solver and structural dynamics solver
    bool fictmass;
    m_session->MatchSolverInfo("FictitiousMassMethod", "True",
                                fictmass, false);
    if(fictmass)
    {
        NekDouble fictrho, fictdamp;
        m_session->LoadParameter("FictMass", fictrho);
        m_session->LoadParameter("FictDamp", fictdamp);
        m_structrho  += fictrho;
        m_structdamp += fictdamp;

        // Storage array of Struct Velocity and Acceleration used for
        // extrapolation of fictitious force
        m_fV = Array<OneD, Array<OneD, Array<OneD, NekDouble> > > (2);
        m_fA = Array<OneD, Array<OneD, Array<OneD, NekDouble> > > (2);
        for (int i = 0; i < m_motion.num_elements(); ++i)
        {
            m_fV[i] = Array<OneD, Array<OneD, NekDouble> > (2);
            m_fA[i] = Array<OneD, Array<OneD, NekDouble> > (2);

            for(int n = 0; n < 2; ++n)
            {
                m_fV[i][n] = Array<OneD, NekDouble>(ns, 0.0);
                m_fA[i][n] = Array<OneD, NekDouble>(ns, 0.0);
            }
        }
    }

    // Setting the coefficient matrices for solving structural dynamic ODEs
	std::string solver = m_session->GetSolverInfo("StructDynSolver");
	if(boost::iequals(solver, "ModalDecompostion"))
	{
    	SetDynEqCoeffMatrix();
    }
	else if(boost::iequals(solver, "FiniteElementMethod"))
	{
    	SetStiffnessMatrix();
	}

    // Set initial condition for cable's motion
    int cnt = 0;
    
    for(int j = 0; j < m_funcName.num_elements(); j++)
    {
        // loading from the specified files through inputstream
        if(m_IsFromFile[cnt] && m_IsFromFile[cnt+1])
        {
            std::ifstream inputStream;

            int nzpoints;

            m_session->LoadParameter("HomStructModesZ", nzpoints);

            if (vcomm->GetRank() == 0)
            {
                std::string filename = m_session->GetFunctionFilename(
                    m_funcName[j], m_motion[0]);
                
                // Open intputstream for cable motions
                inputStream.open(filename.c_str());

                // Import the head string from the file
                Array<OneD, std::string> tmp(9);
                for(int n = 0; n < tmp.num_elements(); n++)
                {
                    inputStream >> tmp[n];
                }

                NekDouble time, z_cds;
                // Import the motion variables from the file
                for (int n = 0; n < nzpoints; n++)
                {
                    inputStream >> setprecision(6) >> time;
                    inputStream >> setprecision(6) >> z_cds;
					for(int i = 0; i < 2; i++)
					{
						for(int j = 0; j < 3; j++)
						{
                    		inputStream >> setprecision(8) >> m_vib[i][j][n];
						}
					}
                }
                // Close inputstream for cable motions
                inputStream.close();
            }
            cnt = cnt + 2;
        }
        else //Evaluate from the functions specified in xml file
        {
            if(colrank == 0)
            {	
                Array<OneD, NekDouble> x0(ns, 0.0);
                Array<OneD, NekDouble> x1(ns, 0.0);
                Array<OneD, NekDouble> x2(ns, 0.0);

                for (int p = 0; p < ns; p++)
                {
                    x2[p] = m_length/(ns-1)*p;
                }

                Array<OneD, LibUtilities::EquationSharedPtr> ffunc(2);
				for(int i = 0; i < 2; i++)
				{
                	ffunc[i] = m_session->GetFunction(m_funcName[j], m_motion[i]);
                	ffunc[i]->Evaluate(x0, x1, x2, 0.0, m_vib[i][j]);
            	}
			}
            cnt = cnt + 2;
        }
    }
}

/**
 *
 */
void ForcingMovingBody::SetDynEqCoeffMatrix()
{
    int nplanes;

    nplanes = m_session->GetParameter("HomStructModesZ")-1; 

    m_CoeffMat_A = Array<OneD, DNekMatSharedPtr>(nplanes);
    m_CoeffMat_B = Array<OneD, DNekMatSharedPtr>(nplanes);

    NekDouble tmp1, tmp2, tmp3;
    NekDouble tmp4, tmp5, tmp6, tmp7;

    // load the structural dynamic parameters from xml file
    NekDouble cabletension;
    NekDouble bendingstiff;
    NekDouble structstiff;
    m_session->LoadParameter("StructStiff",  structstiff,  0.0);
    m_session->LoadParameter("CableTension", cabletension, 0.0);
    m_session->LoadParameter("BendingStiff", bendingstiff, 0.0);

    tmp1 =   m_timestep * m_timestep;
    tmp2 =  structstiff / m_structrho;
    tmp3 = m_structdamp / m_structrho;
    tmp4 = cabletension / m_structrho;
    tmp5 = bendingstiff / m_structrho;

    // solve the ODE in the wave space for cable motion to obtain disp, vel and
    // accel

    std::string supptype = m_session->GetSolverInfo("SupportType");

    for(int plane = 0; plane < nplanes; plane++)
    {
        int nel = 3;
        m_CoeffMat_A[plane]
                = MemoryManager<DNekMat>::AllocateSharedPtr(nel,nel);
        m_CoeffMat_B[plane]
                = MemoryManager<DNekMat>::AllocateSharedPtr(nel,nel);

        // Initialised to avoid compiler warnings.
        unsigned int K = 0;
        NekDouble beta = 0.0;

        if (boost::iequals(supptype, "Free-Free"))
        {
            K = plane/2;
            beta = 2.0 * M_PI/m_length;
        }
        else if(boost::iequals(supptype, "Pinned-Pinned"))
        {
            K = plane+1;
            beta = M_PI/m_length;
        }
        else
        {
            ASSERTL0(false,
                        "Unrecognized support type for cable's motion");
        }

        tmp6 = beta * K;
        tmp6 = tmp6 * tmp6;
        tmp7 = tmp6 * tmp6;
        tmp7 = tmp2 + tmp4 * tmp6 + tmp5 * tmp7;

        (*m_CoeffMat_A[plane])(0,0) = tmp7;
        (*m_CoeffMat_A[plane])(0,1) = tmp3;
        (*m_CoeffMat_A[plane])(0,2) = 1.0;
        (*m_CoeffMat_A[plane])(1,0) = 0.0;
        (*m_CoeffMat_A[plane])(1,1) = 1.0;
        (*m_CoeffMat_A[plane])(1,2) =-m_timestep/2.0;
        (*m_CoeffMat_A[plane])(2,0) = 1.0;
        (*m_CoeffMat_A[plane])(2,1) = 0.0;
        (*m_CoeffMat_A[plane])(2,2) =-tmp1/4.0;
        (*m_CoeffMat_B[plane])(0,0) = 0.0;
        (*m_CoeffMat_B[plane])(0,1) = 0.0;
        (*m_CoeffMat_B[plane])(0,2) = 0.0;
        (*m_CoeffMat_B[plane])(1,0) = 0.0;
        (*m_CoeffMat_B[plane])(1,1) = 1.0;
        (*m_CoeffMat_B[plane])(1,2) = m_timestep/2.0;
        (*m_CoeffMat_B[plane])(2,0) = 1.0;
        (*m_CoeffMat_B[plane])(2,1) = m_timestep;
        (*m_CoeffMat_B[plane])(2,2) = tmp1/4.0;

        m_CoeffMat_A[plane]->Invert();
        (*m_CoeffMat_B[plane]) =
            (*m_CoeffMat_A[plane]) * (*m_CoeffMat_B[plane]);
    }
}

/**
 *
 */
void ForcingMovingBody::SetStiffnessMatrix()
{
    // load the structural dynamic parameters
    NekDouble cabletension;
	NekDouble bendingstiff;
		
    m_session->LoadParameter("CableTension", cabletension, 0.0);
    m_session->LoadParameter("BendingStiff", bendingstiff, 0.0);
    int npts;
    m_session->LoadParameter("HomStructModesZ", npts);

    m_CoeffMat_K = MemoryManager<DNekMat>::AllocateSharedPtr(
                            2*npts,2*npts);
    m_CoeffMat_M = MemoryManager<DNekMat>::AllocateSharedPtr(
                            2*npts,2*npts);
    m_CoeffMat_D = MemoryManager<DNekMat>::AllocateSharedPtr(
                            2*npts,2*npts);


	NekDouble elength = m_length/(npts-1);
	NekDouble ei = bendingstiff/(elength*elength*elength);
	NekDouble ma = elength/420.0;
	NekDouble k_b[4][4],k_t[4][4],m[4][4],d[4][4];

    k_b[0][0] = 12.0;
    k_b[0][1] = 6.0*elength;
    k_b[0][2] = -12.0;
    k_b[0][3] = 6.0*elength;

    k_b[1][0] = k_b[0][1];
    k_b[1][1] = 4.0*elength*elength;
    k_b[1][2] = -6.0*elength;
    k_b[1][3] = 2.0*elength*elength;

    k_b[2][0] = k_b[0][2];
    k_b[2][1] = k_b[1][2];
    k_b[2][2] = 12.0;
    k_b[2][3] = -6.0*elength;

    k_b[3][0] = k_b[0][3];
    k_b[3][1] = k_b[1][3];
    k_b[3][2] = k_b[2][3];
    k_b[3][3] = 4.0*elength*elength;

    k_t[0][0] = 6.0/5.0;
    k_t[0][1] = elength/10.0;
    k_t[0][2] = -6.0/5.0;
    k_t[0][3] = elength/10.0;

    k_t[1][0] = k_t[0][1];
    k_t[1][1] = 2.0/15.0*elength*elength;
    k_t[1][2] = -elength/10.0;
    k_t[1][3] = -elength*elength/30.0;

    k_t[2][0] = k_t[0][2];
    k_t[2][1] = k_t[1][2];
    k_t[2][2] = 6.0/5.0;
    k_t[2][3] = -elength/10.0;

    k_t[3][0] = k_t[0][3];
    k_t[3][1] = k_t[1][3];
    k_t[3][2] = k_t[2][3];
    k_t[3][3] = 2.0/15.0*elength*elength;

    m[0][0] = 156.0;
    m[0][1] = 22.0*elength;
    m[0][2] = 54.0;
    m[0][3] =-13.0*elength;

    m[1][0] = m[0][1];
    m[1][1] = 4.0*elength*elength;
    m[1][2] = 13.0*elength;
    m[1][3] =-3.0*elength*elength;

    m[2][0] = m[0][2];
    m[2][1] = m[1][2];
    m[2][2] = 156.0;
    m[2][3] =-22.0*elength;

    m[3][0] = m[0][3];
    m[3][1] = m[1][3];
    m[3][2] = m[2][3];
    m[3][3] = 4.0*elength*elength;

    for(int i = 0; i < 4; i++)
    {
    	for(int j = 0; j < 4; j++ )
        {
        	k_b[i][j] = k_b[i][j] * ei;
            k_t[i][j] = k_t[i][j] * cabletension/elength;
			d[i][j]   = m[i][j] * ma * m_structdamp;
            m[i][j]   = m[i][j] * ma * m_structrho;
        }
    }

    for(int i = 0; i < 2*npts; i++)
    {
        for(int j = 0; j < 2*npts; j++)
        {
        	(*m_CoeffMat_K)(i,j) = 0.0;
            (*m_CoeffMat_M)(i,j) = 0.0;
            (*m_CoeffMat_D)(i,j) = 0.0;

        }
    }

    for(int n = 0; n < npts-1; n++)
    {
    	for(int i = 0; i < 4; i++)
        {
        	for(int j = 0; j < 4; j++)
            {
            	(*m_CoeffMat_K)(2*n+i,2*n+j) += k_b[i][j]+ k_t[i][j];
                (*m_CoeffMat_M)(2*n+i,2*n+j) += m[i][j];
                (*m_CoeffMat_D)(2*n+i,2*n+j) += d[i][j];

            }
        }
		
	}

	NekDouble alfa,delta;
	NekDouble c0,c1;

	alfa  = 0.25;
	delta = 0.5;
	c0 = 1.0/(alfa*m_timestep*m_timestep);
	c1 = delta/(alfa*m_timestep);

    for(int i = 0; i < 2*npts; i++)
    {
        for(int j = 0; j < 2*npts; j++)
        {
            (*m_CoeffMat_K)(i,j) += 
            	c0*(*m_CoeffMat_M)(i,j)+c1*(*m_CoeffMat_D)(i,j);
        }
    }
	
	//boundary condition is imposed here;
	(*m_CoeffMat_K)(0,0) = 1.0e20;
	(*m_CoeffMat_K)(2*npts-2,2*npts-2) = 1.0e20;

	//calculate the inverse of the stiffness matrix
	m_CoeffMat_K->Invert();
}

/**
 * Function to roll time-level storages to the next step layout.
 * The stored data associated with the oldest time-level
 * (not required anymore) are moved to the top, where they will
 * be overwritten as the solution process progresses.
 */
void ForcingMovingBody::RollOver(Array<OneD, Array<OneD, NekDouble> > &input)
{
    int nlevels = input.num_elements();

    Array<OneD, NekDouble> tmp;
    tmp = input[nlevels-1];

    for(int n = nlevels-1; n > 0; --n)
    {
        input[n] = input[n-1];
    }

    input[0] = tmp;
}

/**
 *
 */
void ForcingMovingBody::CheckIsFromFile(const TiXmlElement* pForce)
{

    m_funcName = Array<OneD, std::string> (3);
    m_motion = Array<OneD, std::string> (2);
    m_motion[0] = "x";
    m_motion[1] = "y";

    m_IsFromFile = Array<OneD, bool> (6);
    // Loading the x-dispalcement (m_zta) and the y-displacement (m_eta)
    // Those two variables are both functions of z and t and the may come
    // from an equation (forced vibration) or from another solver which, given
    // the aerodynamic forces at the previous step, calculates the 
    // displacements.

    //Get the body displacement: m_zta and m_eta
    const TiXmlElement* funcNameElmt_D
                    = pForce->FirstChildElement("DISPLACEMENTS");
    ASSERTL0(funcNameElmt_D,
             "MOVINGBODYFORCE tag has to specify a function name which "
             "prescribes the body displacement as d(z,t).");

    m_funcName[0] = funcNameElmt_D->GetText();
    ASSERTL0(m_session->DefinesFunction(m_funcName[0]),
             "Function '" + m_funcName[0] + "' not defined.");

    //Get the body velocity of movement: d(m_zta)/dt and d(m_eta)/dt
    const TiXmlElement* funcNameElmt_V
                    = pForce->FirstChildElement("VELOCITIES");
    ASSERTL0(funcNameElmt_D,
             "MOVINGBODYFORCE tag has to specify a function name which "
             "prescribes the body velocity of movement as v(z,t).");

    m_funcName[1] = funcNameElmt_V->GetText();
    ASSERTL0(m_session->DefinesFunction(m_funcName[1]),
             "Function '" + m_funcName[1] + "' not defined.");


    //Get the body acceleration: dd(m_zta)/ddt and dd(m_eta)/ddt
    const TiXmlElement* funcNameElmt_A
                    = pForce->FirstChildElement("ACCELERATIONS");
    ASSERTL0(funcNameElmt_A,
             "MOVINGBODYFORCE tag has to specify a function name which "
             "prescribes the body acceleration as a(z,t).");

    m_funcName[2] = funcNameElmt_A->GetText();
    ASSERTL0(m_session->DefinesFunction(m_funcName[2]),
             "Function '" + m_funcName[2] + "' not defined.");

    LibUtilities::FunctionType vType;

    // Check Displacement x
    vType = m_session->GetFunctionType(m_funcName[0],m_motion[0]);
    if(vType == LibUtilities::eFunctionTypeExpression)
    {
        m_IsFromFile[0] = false;
    }
    else if (vType == LibUtilities::eFunctionTypeFile)
    {
        m_IsFromFile[0] = true;
    }
    else
    {
        ASSERTL0(false, "The displacements in x must be specified via an "
                        "equation <E> or a file <F>");
    }

    // Check Displacement y
    vType = m_session->GetFunctionType(m_funcName[0],m_motion[1]);
    if(vType == LibUtilities::eFunctionTypeExpression)
    {
        m_IsFromFile[1] = false;
    }
    else if (vType == LibUtilities::eFunctionTypeFile)
    {
        m_IsFromFile[1] = true;
    }
    else
    {
        ASSERTL0(false, "The displacements in y must be specified via an "
                        "equation <E> or a file <F>");
    }

    // Check Velocity x
    vType = m_session->GetFunctionType(m_funcName[1],m_motion[0]);
    if (vType == LibUtilities::eFunctionTypeExpression)
    {
        m_IsFromFile[2] = false;
    }
    else if (vType == LibUtilities::eFunctionTypeFile)
    {
        m_IsFromFile[2] = true;
    }
    else
    {
        ASSERTL0(false, "The velocities in x must be specified via an equation "
                        "<E> or a file <F>");
    }

    // Check Velocity y
    vType = m_session->GetFunctionType(m_funcName[1],m_motion[1]);
    if (vType == LibUtilities::eFunctionTypeExpression)
    {
        m_IsFromFile[3] = false;
    }
    else if (vType == LibUtilities::eFunctionTypeFile)
    {
        m_IsFromFile[3] = true;
    }
    else
    {
        ASSERTL0(false, "The velocities in y must be specified via an equation "
                        "<E> or a file <F>");
    }

    // Check Acceleration x
    vType = m_session->GetFunctionType(m_funcName[2],m_motion[0]);
    if (vType == LibUtilities::eFunctionTypeExpression)
    {
        m_IsFromFile[4] = false;
    }
    else if (vType == LibUtilities::eFunctionTypeFile)
    {
        m_IsFromFile[4] = true;
    }
    else
    {
        ASSERTL0(false, "The accelerations in x must be specified via an "
                        "equation <E> or a file <F>");
    }

    // Check Acceleration y
    vType = m_session->GetFunctionType(m_funcName[2],m_motion[1]);
    if (vType == LibUtilities::eFunctionTypeExpression)
    {
        m_IsFromFile[5] = false;
    }
    else if (vType == LibUtilities::eFunctionTypeFile)
    {
        m_IsFromFile[5] = true;
    }
    else
    {
        ASSERTL0(false, "The accelerations in y must be specified via an "
                        "equation <E> or a file <F>");
    }
}

/**
 *
 */
void ForcingMovingBody::InitialiseFilter(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        const TiXmlElement* pForce)
{
    // Get the outputfile name, output frequency and 
    // the boundary's ID for the cable's wall
    std::string typeStr = pForce->Attribute("TYPE");
    std::map<std::string, std::string> vParams;

    const TiXmlElement *param = pForce->FirstChildElement("PARAM");
    while (param)
    {
        ASSERTL0(param->Attribute("NAME"),
                 "Missing attribute 'NAME' for parameter in filter "
                 + typeStr + "'.");
        std::string nameStr = param->Attribute("NAME");

        ASSERTL0(param->GetText(), "Empty value string for param.");
        std::string valueStr = param->GetText();

        vParams[nameStr] = valueStr;

        param = param->NextSiblingElement("PARAM");
    }

    // Creat the filter for MovingBody, where we performed the calculation of
    // fluid forces and write both the aerodynamic forces and motion variables
    // into the output files
    m_MovBodyfilter = MemoryManager<FilterMovingBody>::
                                    AllocateSharedPtr(m_session, vParams);

    // Initialise the object of MovingBody filter
    m_MovBodyfilter->Initialise(pFields, 0.0);

}

}
