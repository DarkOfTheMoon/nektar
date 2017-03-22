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
#include <LibUtilities/LinearAlgebra/Lapack.hpp>
using namespace std;

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
    InitialiseVibrationModel(pFields);
	
	LibUtilities::CommSharedPtr vcomm = pFields[0]->GetComm();
    int colrank = vcomm->GetColumnComm()->GetRank();
	if(colrank == 0)
	{
		 m_q = Array<OneD, Array <OneD, NekDouble> > (2);
         m_q[0] = Array<OneD, NekDouble>(m_nv+1, 0.0);
         m_q[1] = Array<OneD, NekDouble>(m_nv+1, 0.0);
	}
    // Load mapping
    m_mapping = GlobalMapping::Mapping::Load(m_session, pFields);
    m_mapping->SetTimeDependent( true );

    m_mapping->SetFromFunction( false );

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

    // for "free" type (m_vdim = 2), the cable vibrates both in streamwise and crossflow
    // directions, for "constrained" type (m_vdim = 1), the cable only vibrates in crossflow
    // direction, and for "forced" type (m_vdim = 0), the calbe vibrates specifically along
    // a given function or file
    if(m_vdim == 1 || m_vdim == 2)
    {
        // For free vibration case, displacements, velocities and acceleartions
        // are obtained through solving structure dynamic model
        EvaluateVibrationModel(pFields, time);
        
        // Convert result to format required by mapping
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
	
		//Get initial mapping coordinates
        if(m_IsMapFromFile[0] && m_IsMapFromFile[1])
        {
            ASSERTL0(false, "Motion loading from file needs specific "
                                "implementation: Work in Progress!");
        }
        else
        {
            EvaluateFunction(pFields, m_session, m_motion[0], coords[0],
                                 m_mapfuncName[0], time);
            EvaluateFunction(pFields, m_session, m_motion[1], coords[1],
                                 m_mapfuncName[0], time);
        }
	
        // Add displacement to coordinates
        Vmath::Vadd(physTot, coords[0], 1, m_zta[0], 1, coords[0], 1);
        Vmath::Vadd(physTot, coords[1], 1, m_eta[0], 1, coords[1], 1);
        // Copy velocities
        Vmath::Vcopy(physTot, m_zta[1], 1, coordsVel[0], 1);
        Vmath::Vcopy(physTot, m_eta[1], 1, coordsVel[1], 1);
        

        // Update mapping
        m_mapping->UpdateMapping(time, coords, coordsVel);

        LibUtilities::CommSharedPtr vcomm = pFields[0]->GetComm();

        // Write the motions into the file
        m_MovBodyfilter->WriteMotion(m_vib, time);
    }
    else if(m_vdim == 0)
    {
        ASSERTL0(!m_IsHomostrip,
             "Force vibration implemented just for full resolution.");
    	Array<OneD, Array<OneD, NekDouble> > Hydroforces (2);
    	Hydroforces[0] = Array<OneD, NekDouble> (m_np, 0.0);
    	Hydroforces[1] = Array<OneD, NekDouble> (m_np, 0.0);
		m_MovBodyfilter->UpdateForce(pFields, Hydroforces, time);
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

        // Add displacement to coordinates
        Vmath::Vadd(physTot, coords[0], 1, m_zta[0], 1, coords[0], 1);
        Vmath::Vadd(physTot, coords[1], 1, m_eta[0], 1, coords[1], 1);
        // Copy velocities
        Vmath::Vcopy(physTot, m_zta[1], 1, coordsVel[0], 1);
        Vmath::Vcopy(physTot, m_eta[1], 1, coordsVel[1], 1);
        // Update mapping
        m_mapping->UpdateMapping(time, coords, coordsVel);
    }
    else
    {
        ASSERTL0(false, 
                 "Unrecogized vibration type for cable's dynamic model");
    }
}

/**
 *
 */
void ForcingMovingBody::EvaluateVibrationModel(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
		const NekDouble                                   &time)
{
    Array<OneD, Array<OneD, NekDouble> > Hydroforces (2);
    Hydroforces[0] = Array<OneD, NekDouble> (m_np, 0.0);
    Hydroforces[1] = Array<OneD, NekDouble> (m_np, 0.0);

	// Update hydroforces	
    m_MovBodyfilter->UpdateForce(pFields, Hydroforces, time);
	
    LibUtilities::CommSharedPtr vcomm = pFields[0]->GetComm();
    int colrank = vcomm->GetColumnComm()->GetRank();
    int nproc   = vcomm->GetColumnComm()->GetSize();

    int ntot;
	ntot =m_nstrips*m_nz;

    //the hydrodynamic forces
    Array<OneD, Array <OneD, NekDouble> > fces(2);
    //forces in x-direction
    fces[0] = Array <OneD, NekDouble> (ntot,0.0);
    //forces in y-direction
    fces[1] = Array <OneD, NekDouble> (ntot,0.0);

    //reshuffle the forces
    if(colrank == 0)
    {
        Vmath::Vcopy(m_np, Hydroforces[0], 1, fces[0], 1);
        Vmath::Vcopy(m_np, Hydroforces[1], 1, fces[1], 1);
        Array<OneD, NekDouble> tmp (m_np);
        Array<OneD, NekDouble> tmp0(m_np);
        Array<OneD, NekDouble> tmp1(m_np);
		for(int n = 0; n < m_nstrips; n++)
		{
        	for (int i = 0; i < nproc/m_nstrips; i++)
        	{
				int idrank = n+i*m_nstrips;
				if(idrank != 0)
				{
					int offset = (i+n*nproc/m_nstrips)*m_np;
                	vcomm->GetColumnComm()->Recv(idrank, tmp0);
                	vcomm->GetColumnComm()->Recv(idrank, tmp1);
					Vmath::Vcopy(m_np, tmp0, 1, tmp = fces[0]+offset, 1);		
					Vmath::Vcopy(m_np, tmp1, 1, tmp = fces[1]+offset, 1);		
				}
            }
        }
    }
    else
    {
    	vcomm->GetColumnComm()->Send(0, Hydroforces[0]);
        vcomm->GetColumnComm()->Send(0, Hydroforces[1]);
    }
	
    if(colrank == 0)
    {
        // Fictitious mass method used to stablize the explicit coupling at
        // relatively lower mass ratio
        if(m_IsFictmass)
        {
            NekDouble fictrho, fictdamp;
            m_session->LoadParameter("FictMass", fictrho);
            m_session->LoadParameter("FictDamp", fictdamp);

            static NekDouble Betaq_Coeffs[2][2] = 
                                {{1.0,  0.0},{2.0, -1.0}};

            // only consider second order approximation for fictitious variables
            int  intOrder= 2;
            int  nint    = min(m_movingBodyCalls,intOrder);
            int  nlevels = m_fV[0].num_elements();

            for(int i = 0; i < m_motion.num_elements(); ++i)
            {
                RollOver(m_fV[i]);
                RollOver(m_fA[i]);

                Vmath::Vcopy(ntot, m_spv[i][1], 1, m_fV[i][0], 1);
                Vmath::Vcopy(ntot, m_spv[i][2], 1, m_fA[i][0], 1);

                // Extrapolate to n+1
                Vmath::Smul(ntot, 
                            Betaq_Coeffs[nint-1][nint-1],
                            m_fV[i][nint-1],    1,
                            m_fV[i][nlevels-1], 1);
                Vmath::Smul(ntot, 
                            Betaq_Coeffs[nint-1][nint-1],
                            m_fA[i][nint-1],    1,
                            m_fA[i][nlevels-1], 1);

                for(int n = 0; n < nint-1; ++n)
                {
                    Vmath::Svtvp(ntot, 
                                 Betaq_Coeffs[nint-1][n],
                                 m_fV[i][n],1,m_fV[i][nlevels-1],1,
                                 m_fV[i][nlevels-1],1);
                    Vmath::Svtvp(ntot, 
                                 Betaq_Coeffs[nint-1][n],
                                 m_fA[i][n],1,m_fA[i][nlevels-1],1,
                                 m_fA[i][nlevels-1],1);
                }

                // Add the fictitious forces on the RHS of the equation
                Vmath::Svtvp(ntot, fictdamp,m_fV[i][nlevels-1],1,
                             fces[i],1,fces[i],1);
                Vmath::Svtvp(ntot, fictrho, m_fA[i][nlevels-1],1,
                             fces[i],1,fces[i],1);
            }
        }
    }

    //structural solver runs on the root process
    if(colrank == 0)
    {
        // Tensioned cable model is evaluated in wave space
        for(int n = 0, cn = 1; n < m_vdim; n++, cn--)
        {  
			Vmath::Zero(m_nv+1,m_q[0],1);
			Vmath::Zero(m_nv+1,m_q[1],1);

			int n_ave = m_nz/m_nv; 
			if(m_nstrips == 1)
			{
				if (n_ave > 1)
				{
					for(int i = 1; i < m_nv; i++)
					{
						for(int j = 0; j < n_ave+1; j++)
						{
							m_q[cn][i] += 
								fces[cn][i*n_ave+j-n_ave/2]/(n_ave+1);
						}
					}
					for(int j = 0; j < n_ave/2+1; j++)
					{
						m_q[cn][0] += fces[cn][j]/(n_ave/2+1);
					}
				}
				else
				{
					Vmath::Vcopy(m_nv,fces[cn],1,m_q[cn],1);
				}
			}
			else
			{
				for(int i = 0; i < m_nstrips; i++)
				{
					for(int j = 0; j < m_nz; j++)
					{
						m_q[cn][i] += fces[cn][i*m_nz+j]/m_nz;
					}
				}
			}
			m_q[cn][m_nv]=m_q[cn][0];
		}

    	std::string StructDynSolver = 
		m_session->GetSolverInfo("StructDynSolver");
   		if(boost::iequals(StructDynSolver, "ModalDecomposition"))
   		{
        	for(int n = 0, cn = 1; n < m_vdim; n++, cn--)
        	{
				ModalDecompositionMethod(m_q[cn], m_vib[cn]);
			}
		}
        else if (boost::iequals(StructDynSolver, "SHARPy"))
        {
        	EvaluateSHARPy(m_q,m_vib);
        }

        // Tensioned cable model is evaluated in wave space
        for(int n = 0, cn = 1; n < m_vdim; n++, cn--)
        {
			std::vector<double> Xp(m_nv+1);
			std::vector<double> temp(m_nv+1);
	
			for(int j = 0; j < 3; j++)
			{
				//cubic spline interpolation
    			for(int i = 0; i < m_nv+1; i++)
    			{
        			Xp[i] = m_length*i/m_nv;
					temp[i] = m_vib[cn][j][i];
    			}

            	tk::spline s;
            	s.set_points(Xp,temp);
				NekDouble x;
				for(int n = 0; n < m_nstrips; n++)
				{
					for(int i = 0; i < m_nz; i++)
					{
						if(m_nstrips == 1)
						{
							x = m_length*i/m_nz;
						}
						else
						{
							x = m_length*n/m_nv;
						}
						m_spv[cn][j][n*m_nz+i] = s(x);
					}
				}
			}
        }
    }

    Array<OneD, Array<OneD, NekDouble> > dis(2);
    Array<OneD, Array<OneD, NekDouble> > vel(2);
	for (int i = 0; i < 2; i++)
	{
		dis[i] = Array<OneD, NekDouble>(m_np);
		vel[i] = Array<OneD, NekDouble>(m_np);
	}

    // send physical coeffients to all planes of each processor
	if(colrank != 0)
	{
		for(int j = 0; j < 2; j++)
		{
			vcomm->GetColumnComm()->Recv(0, dis[j]);
			vcomm->GetColumnComm()->Recv(0, vel[j]);
		}
	}
	else
	{
		for(int n = 0; n < m_nstrips; n++)
		{
            for (int i = 0; i < nproc/m_nstrips; ++i)
            {
				int idrank = n+i*m_nstrips;
				if(idrank != 0)
				{
					int offset = (i+n*nproc/m_nstrips)*m_np;	
					for(int j = 0; j < 2; j++)
					{
						Array<OneD, NekDouble> tmp0(m_np);
						Array<OneD, NekDouble> tmp1(m_np);
						Vmath::Vcopy(m_np, m_spv[j][0]+offset, 1, tmp0, 1);	
						Vmath::Vcopy(m_np, m_spv[j][1]+offset, 1, tmp1, 1);	
						vcomm->GetColumnComm()->Send(idrank, tmp0);
						vcomm->GetColumnComm()->Send(idrank, tmp1);
					}
				}
			}
			
			for(int j = 0; j < 2; j++)
			{
				Vmath::Vcopy(m_np, m_spv[j][0], 1, dis[j], 1);	
				Vmath::Vcopy(m_np, m_spv[j][1], 1, vel[j], 1);	
			}
		}	
	}

    // Set the m_forcing term based on the motion of the cable
    for(int plane = 0; plane < m_np; plane++)
    {
    	int n = pFields[0]->GetPlane(plane)->GetTotPoints();

        Array<OneD, NekDouble> tmp(n,0.0);
		
        int offset  = plane * n;
        Vmath::Fill(n, dis[0][plane], tmp = m_zta[0] + offset, 1);
        Vmath::Fill(n, dis[1][plane], tmp = m_eta[0] + offset, 1);
        Vmath::Fill(n, vel[0][plane], tmp = m_zta[1] + offset, 1);
        Vmath::Fill(n, vel[1][plane], tmp = m_eta[1] + offset, 1);
    }
}
 
/**
 *
 */
void ForcingMovingBody::ModalDecompositionMethod(
     	const Array<OneD, NekDouble> &HydroForces,
              Array<OneD, Array<OneD, NekDouble> > &motions)
{  

    Array<OneD, Array<OneD, NekDouble> > fft_i(4);
    Array<OneD, Array<OneD, NekDouble> > fft_o(4);

    for(int i = 0; i < 4; i++)
    {
        fft_i[i] = Array<OneD, NekDouble>(m_nv, 0.0);
        fft_o[i] = Array<OneD, NekDouble>(m_nv, 0.0);
    }

    Vmath::Vcopy(m_nv, HydroForces, 1, fft_i[0], 1);
	for(int i = 0; i < 3; i++)
	{
		Vmath::Vcopy(m_nv, motions[i],1,fft_i[i+1],1);
	}
    
    std::string supptype = m_session->GetSolverInfo("SupportType");
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
						sin(M_PI/(N)*(k+0.5)*(i+1));
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
    for(int i = 0; i < m_nv; ++i)
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
						sin(M_PI/(N)*(k+0.5)*(i+1))*2/N;
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
		Vmath::Vcopy(m_nv,fft_o[i],1,motions[i],1);
		motions[i][m_nv] = fft_o[i][0];
	}
}

/**
 *
 */
void ForcingMovingBody::InitialiseSHARPy()
{
    //----------------------------Initial Displacement: ImpStart vs Static Solution
    // Initialise static beam data.
    // 1. Read input data.
	m_session->LoadParameter("BeamLength_SHARPy", m_BeamLength);
	m_session->LoadParameter("NumNodesElem_SHARPy", m_NumNodesElem,2);
	m_OutFile 				= Array<OneD, char>(25,"sol.txt");			 //
	m_session->LoadParameter("MaxElNod_SHARPy", m_MaxElNod,3);  		 // Max number of nodes per element.
	m_session->LoadParameter("NumElems_SHARPy", m_NumElems,20); 		 // Elements discreising the beam
	m_session->LoadParameter("ElemProj_SHARPy", m_ElemProj,0);  		 // Element info computed in (1) global frame (2) fixed element frame (3) moving element frame.
	m_session->LoadParameter("MaxIterations_SHARPy", m_MaxIterations,20);// Newton-Raphson iterations for nonlinear solution
	m_session->LoadParameter("NumLoadSteps_SHARPy", m_NumLoadSteps,10);  // Number of load increments
	m_session->LoadParameter("NumGauss_SHARPy", m_NumGauss,3);  		 // Gauss points in element integration
	m_session->LoadParameter("Solution_SHARPy", m_Solution,912);		 // 902/912: cbeam3 linear/nonlinear flexible-body dynamic 
	m_session->LoadParameter("DimMat_SHARPy", m_DimMat,24);     
	m_session->LoadParameter("DeltaCurved_SHARPy", m_DeltaCurved,1.e-4); // Min. angle for two vectors to be parallel
	m_session->LoadParameter("MinDelta_SHARPy", m_MinDelta,1.e-6);       // Convergence param for Newton-Rhaphson iterations
	m_session->LoadParameter("NewmarkDamp_SHARPy", m_NewmarkDamp,0.01);  // Numerical damping in Newmark integration scheme
	m_FollowerForce 		= true;								         // Forces follow local deflections
	m_FollowerForceRig 		= true;										 // Forces follow the body-fixed frame
	m_PrintInfo				= true;									 // 
	m_OutInBframe 			= true;									 // Print velocities in local-frame (if not, use body-attached frame)
	m_OutInaframe 			= true;									 // Print velocities in body-fixed-frame (if not, use inertial frame)

	m_iStep 				= 0;
	m_gamma					= 1.0/2.0+m_NewmarkDamp;
	m_beta					= 1.0/4.0*(m_gamma+0.5)*(m_gamma+0.5);

	m_NumNodes 				= Array<OneD, int> (m_NumElems);                             
	m_MemNo 				= Array<OneD, int> (m_NumElems);								
	m_Conn 					= Array<OneD, int> (m_MaxElNod*m_NumElems);						
	m_Master_Array 			= Array<OneD, int> (m_MaxElNod*m_NumElems*2);				
	m_Length 				= Array<OneD, NekDouble> (m_NumElems,0.0);                           
	m_PreCurv 				= Array<OneD, NekDouble> (3*m_NumElems,0.0);						
	m_Psi 					= Array<OneD, NekDouble> (3*m_NumElems,0.0);							
	m_Vector				= Array<OneD, NekDouble> (3*m_NumElems,0.0);                          
	m_Mass_Array 			= Array<OneD, NekDouble> (6*m_NumElems*6,0.0);                   
	m_Stiff_Array 			= Array<OneD, NekDouble> (6*m_NumElems*6,0.0);                  
	m_InvStiff_Array 		= Array<OneD, NekDouble> (6*m_NumElems*6,0.0);               
	m_RBMass_Array 			= Array<OneD, NekDouble> (m_MaxElNod*m_NumElems*6*6,0.0);		

	Input_elem_SHARPy();

	m_BoundConds 			= Array<OneD, int> (m_NumNodes_tot);
	m_Nod_Sflag 			= Array<OneD, int> (m_NumNodes_tot);		
	m_Master 				= Array<OneD, int> (2*m_NumNodes_tot);   
	m_Vdof 					= Array<OneD, int> (m_NumNodes_tot);       						
	m_Fdof 					= Array<OneD, int> (m_NumNodes_tot);       								
	m_ListIN 				= Array<OneD, int> (m_NumNodes_tot);

	m_PosIni_Array 			= Array<OneD, NekDouble> (m_NumNodes_tot*3,0.0);
	m_ForceStatic_Array 	= Array<OneD, NekDouble> (m_NumNodes_tot*6,0.0);
	m_PhiNodes 				= Array<OneD, NekDouble> (m_NumNodes_tot,0.0);
	m_PsiIni_Array 			= Array<OneD, NekDouble> (m_NumElems*m_MaxElNod*3,0.0);
	m_PosDefor_Array 		= Array<OneD, NekDouble> (m_NumNodes_tot*3,0.0);											
	m_PsiDefor_Array		= Array<OneD, NekDouble> (m_NumElems*m_MaxElNod*3,0.0);										
	m_PosDeforDot_Array		= Array<OneD, NekDouble> (m_NumNodes_tot*3,0.0);			
	m_PsiDeforDot_Array 	= Array<OneD, NekDouble> (m_NumElems*m_MaxElNod*3,0.0);		
	m_PosDeforDDot_Array 	= Array<OneD, NekDouble> (m_NumNodes_tot*3,0.0);			
	m_PsiDeforDDot_Array 	= Array<OneD, NekDouble> (m_NumElems*m_MaxElNod*3,0.0);	
	m_Quat_Array 			= Array<OneD, NekDouble> (4,0.0);
	m_Cao_Array 			= Array<OneD, NekDouble> (3*3,0.0);
	m_Quat_Array[0] 		= 1.0;
	m_Cao_Array[0]  		= 1.0;
	m_Cao_Array[4]  		= 1.0;
	m_Cao_Array[8]  		= 1.0;

	Input_node_SHARPy();

	// 2. Compute initial (undeformed) geometry.
	SHARPy::Wrap_xbeam_undef_geom(
		m_NumElems,m_NumNodes_tot,
		&m_NumNodes[0],&m_MemNo[0],
		&m_Conn[0],&m_Master_Array[0],
		&m_Length[0],&m_PreCurv[0],
		&m_Psi[0],&m_Vector[0],&m_Mass_Array[0],
		&m_Stiff_Array[0],&m_InvStiff_Array[0],
		&m_RBMass_Array[0],&m_PosIni_Array[0],
		&m_PhiNodes[0],&m_PsiIni_Array[0],
		m_FollowerForce, m_FollowerForceRig,          
		m_PrintInfo, m_OutInBframe, m_OutInaframe,    
		m_ElemProj, m_MaxIterations, m_NumLoadSteps,  
		m_NumGauss, m_Solution, m_DeltaCurved,       
		m_MinDelta, m_NewmarkDamp);         

	// 3. Identify nodal degrees of freedom.
	SHARPy::Wrap_xbeam_undef_dofs(
		m_NumElems,m_NumNodes_tot,
		&m_NumNodes[0],&m_MemNo[0],
		&m_Conn[0],&m_Master_Array[0],
		&m_Length[0],&m_PreCurv[0],
		&m_Psi[0],&m_Vector[0],&m_Mass_Array[0],
		&m_Stiff_Array[0],&m_InvStiff_Array[0],
		&m_RBMass_Array[0],&m_BoundConds[0],
		&m_Master[0],&m_Vdof[0],&m_Fdof[0],m_NumDof,
		&m_Nod_Sflag[0]);

	// Add RBMasses (Rigid-Body masses)
	////

	//Calculate initial displacements.
	Vmath::Vcopy(m_NumNodes_tot*3, m_PosIni_Array, 1, m_PosDefor_Array, 1);
	Vmath::Vcopy(m_NumElems*m_MaxElNod*3, m_PsiIni_Array, 1, m_PsiDefor_Array, 1);

	SHARPy::Wrap_cbeam3_solv_nlnstatic (
		m_NumDof,m_NumElems, 
		&m_NumNodes[0], &m_MemNo[0], &m_Conn[0],       
		&m_Master_Array[0],                           
		&m_Length[0], &m_PreCurv[0],               
		&m_Psi[0], &m_Vector[0], &m_Mass_Array[0],   
		&m_Stiff_Array[0],                        
		&m_InvStiff_Array[0], &m_RBMass_Array[0],  
		m_NumNodes_tot, 
		&m_Master[0], &m_Vdof[0], &m_Fdof[0],       
		&m_ForceStatic_Array[0],                    
		&m_PosIni_Array[0],&m_PsiIni_Array[0],        
		&m_PosDefor_Array[0],&m_PsiDefor_Array[0],  
		m_FollowerForce, m_FollowerForceRig,       
		m_PrintInfo, m_OutInBframe, m_OutInaframe,    
		m_ElemProj, m_MaxIterations, m_NumLoadSteps,  
		m_NumGauss, m_Solution, m_DeltaCurved,        
		m_MinDelta, m_NewmarkDamp);               

	m_NumDof2 		= m_NumDof*m_NumDof;
	m_X 			= Array<OneD, NekDouble> (m_NumDof,0.0);
	m_DX 			= Array<OneD, NekDouble> (m_NumDof,0.0);
	m_DXdt 			= Array<OneD, NekDouble> (m_NumDof,0.0);
	m_DXddt 		= Array<OneD, NekDouble> (m_NumDof,0.0);
	m_Force_Array 	= Array<OneD, NekDouble> (m_NumNodes_tot*6,0.0);      
	m_Vrel 			= Array<OneD, NekDouble> (6,0.0);				
	m_VrelDot 		= Array<OneD, NekDouble> (6,0.0);						
	m_Mglobal_Array = Array<OneD, NekDouble> (m_NumDof2,0.0);			
	m_Mvel_Array 	= Array<OneD, NekDouble> (m_NumDof*6,0.0);			
	m_Cglobal_Array = Array<OneD, NekDouble> (m_NumDof2,0.0);			
	m_Cvel_Array 	= Array<OneD, NekDouble> (m_NumDof*6,0.0);			
	m_Kglobal_Array = Array<OneD, NekDouble> (m_NumDof2,0.0);			
	m_Fglobal_Array = Array<OneD, NekDouble> (m_NumDof2,0.0);			
	m_Qglobal 		= Array<OneD, NekDouble> (m_NumDof,0.0);			
	m_Asys_Array 	= Array<OneD, NekDouble>(m_NumDof2,0.0);

	//---------------------------- Initialise Structural Variables
	// 1. Input data for transient dynamic solution.
	m_session->LoadParameter("t0_SHARPy", m_t0,0.0);    // 
	m_session->LoadParameter("tf_SHARPy", m_tf,100.0);  //
	m_session->LoadParameter("dt_SHARPy", m_dt,0.005);  //
	m_NumSteps = (m_tf - m_t0)/m_dt;

	m_Time = Array<OneD, NekDouble> (m_NumSteps);
	for(int i = 0; i < m_NumSteps; i++)
	{
		m_Time[i] = m_t0 + m_dt*NekDouble(i);
	} 

	//Initialize Local variables.
	for(int k = 0; k < m_NumNodes_tot; k++)
	{
		m_ListIN[k]=m_Vdof[k];
	}

	Vmath::Zero(m_NumDof2,m_Asys_Array,1);
	Vmath::Zero(m_NumDof2,m_Mglobal_Array,1);
	Vmath::Zero(m_NumDof2,m_Cglobal_Array,1);
	Vmath::Zero(m_NumDof2,m_Kglobal_Array,1);
	Vmath::Zero(m_NumDof2,m_Fglobal_Array,1);
	Vmath::Zero(m_NumDof,m_Qglobal,1);

	//Extract initial displacements and velocities.
	SHARPy::Wrap_cbeam3_solv_disp2state(
		m_NumNodes_tot,m_NumDof,m_NumElems,
		&m_Master[0],&m_Vdof[0],&m_Fdof[0],
		&m_PosDefor_Array[0],&m_PsiDefor_Array[0],
		&m_PosDeforDot_Array[0],&m_PsiDeforDot_Array[0],
		&m_X[0],&m_DXdt[0]);

	//Compute initial acceleration (we are neglecting qdotdot in Kmass).
	Vmath::Smul(m_NumNodes_tot*3,0.0,m_PosDeforDot_Array,1,m_PosDeforDDot_Array,1);
	Vmath::Smul(m_NumElems*m_MaxElNod*3,0.0,m_PsiDeforDot_Array,1,m_PsiDeforDDot_Array,1);
	Vmath::Vadd(m_NumNodes_tot*6,m_Force_Array,1,m_ForceStatic_Array,1,m_Force_Array,1);

	SHARPy::Wrap_cbeam3_asbly_dynamic(
		m_NumElems,m_NumNodes_tot,&m_NumNodes[0],
		&m_MemNo[0],&m_Conn[0],&m_Master_Array[0],&m_Length[0],
		&m_PreCurv[0],&m_Psi[0],&m_Vector[0],&m_Mass_Array[0],
		&m_Stiff_Array[0],&m_InvStiff_Array[0],&m_RBMass_Array[0],
		&m_Master[0],&m_Vdof[0],&m_Fdof[0],
		&m_PosIni_Array[0],&m_PsiIni_Array[0],
		&m_PosDefor_Array[0],&m_PsiDefor_Array[0],
		&m_PosDeforDot_Array[0],&m_PsiDeforDot_Array[0],
		&m_PosDeforDDot_Array[0],&m_PsiDeforDDot_Array[0],
		&m_Force_Array[0],&m_Vrel[0],&m_VrelDot[0],
		m_NumDof,m_DimMat,
		m_ms,&m_Mglobal_Array[0],&m_Mvel_Array[0],
		m_cs,&m_Cglobal_Array[0],&m_Cvel_Array[0],
		m_ks,&m_Kglobal_Array[0],
		m_fs,&m_Fglobal_Array[0],&m_Qglobal[0],
		m_FollowerForce, m_FollowerForceRig,        
		m_PrintInfo, m_OutInBframe, m_OutInaframe,   
		m_ElemProj, m_MaxIterations, m_NumLoadSteps,  
		m_NumGauss, m_Solution, m_DeltaCurved,      
		m_MinDelta, m_NewmarkDamp,
		&m_Cao_Array[0]);

	Array<OneD, NekDouble> Temp0_Array(m_NumDof,0.0);
	Array<OneD, NekDouble> Temp1_Array(m_NumDof,0.0);
	
	SHARPy::Wrap_fem_m2v(
		m_NumNodes_tot,6,
		&m_Force_Array[0],m_NumDof,
		&Temp0_Array[0],&m_ListIN[0]);

	DNekMatSharedPtr fglobal 
		= MemoryManager<DNekMat>::AllocateSharedPtr(
			m_NumDof,m_NumDof);

	for(int i = 0; i < m_NumDof; i++) 
	{    
		for(int j = 0; j < m_NumDof; j++) 
		{
			(*fglobal)(i,j)= m_Fglobal_Array[i*m_NumDof+j];
		}
	}

	Blas::Dgemv('N', m_NumDof, m_NumDof,
				1.0, &(fglobal->GetPtr())[0],
				m_NumDof, &Temp0_Array[0], 1,
				0.0,  &Temp1_Array[0],1);

	Vmath::Vsub(m_NumDof,m_Qglobal,1,Temp1_Array,1,m_Qglobal,1);
	Vmath::Smul(m_NumDof,-1.0,m_Qglobal,1,m_Qglobal,1);
	Vmath::Vadd(m_NumDof2,m_Mglobal_Array,1,m_Asys_Array,1,m_Asys_Array,1);

	int N = m_NumDof;
	Array<OneD, int> ipivot (N); 
	int info =0;
	Lapack::Dgetrf(N, N, m_Asys_Array.get(), N, ipivot.get(), info);
	if( info < 0 )
	{    
		std::string message = "ERROR: The " + boost::lexical_cast<std::string>(-info) +
								"th parameter had an illegal parameter for dgetrf";
	}
	else if( info > 0 )
	{    
		std::string message = "ERROR: Element u_" + boost::lexical_cast<std::string>(info) +
		boost::lexical_cast<std::string>(info) + " is 0 from dgetrf";
		ASSERTL0(false, message.c_str());
	}    

	// N means no transponse (direct matrix)
	int ncolumns_b =1;
	Vmath::Vcopy(m_NumDof, m_Qglobal, 1, m_DXddt, 1);
	Lapack::Dgetrs('N', N, ncolumns_b, m_Asys_Array.get(), N, ipivot.get(), m_DXddt.get(), N, info);
	if( info < 0 )
	{    
		std::string message = "ERROR: The " + boost::lexical_cast<std::string>(-info) +
		"th parameter had an illegal parameter for dgetrf";
		ASSERTL0(false, message.c_str());
	}
	else if( info > 0 )
	{    
		std::string message = "ERROR: Element u_" + boost::lexical_cast<std::string>(info) +
		boost::lexical_cast<std::string>(info) + " is 0 from dgetrf";
		ASSERTL0(false, message.c_str());
	}
}

/**
 *
 **/
void ForcingMovingBody::EvaluateSHARPy(const Array<OneD, Array<OneD, NekDouble> > &HydroForces,
				Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &motions)
{
	Array<OneD, NekDouble> temp0(m_NumNodes_tot,0.0);
	//Vmath::Vcopy(m_NumNodes_tot,HydroForces[2],1,m_Force_Array,1);
	Vmath::Vcopy(m_NumNodes_tot,HydroForces[0],1,temp0 = m_Force_Array +   m_NumNodes_tot,1);
	Vmath::Vcopy(m_NumNodes_tot,HydroForces[1],1,temp0 = m_Force_Array + 2*m_NumNodes_tot,1);

	// dimensionalize the hydroforces applying on the cylinder.
	NekDouble D,U_inf,rho_f;
    m_session->LoadParameter("Diameter_CrossSection", D,1.0);
    m_session->LoadParameter("Velocity_Infinity", U_inf,1.0);
    m_session->LoadParameter("Density_Fluid", rho_f,1.0);

	NekDouble coeff = rho_f*U_inf*U_inf*D;
	coeff = m_BeamLength/(m_NumNodes_tot-1)*coeff;
	Vmath::Smul(6*m_NumNodes_tot,coeff,m_Force_Array,1,m_Force_Array,1);
	
	m_dt = m_Time[m_iStep+1]-m_Time[m_iStep];

	int size = 4;
	Array<OneD, NekDouble> Temp0_Array(size*size,0.0);
	Array<OneD, NekDouble> Temp1_Array(size*size,0.0);
	Array<OneD, NekDouble> Temp2_Array(size*size,0.0);
	Array<OneD, NekDouble> Temp3_Array(size*size,0.0);
	Array<OneD, NekDouble> Unit4_Array(size*size,0.0);

	for(int i = 0; i < 4; i++)
	{
		Unit4_Array[i*4+i] = 1.0;
	}

	NekDouble tmp = 0.25*(m_Time[m_iStep+1]-m_Time[m_iStep]);
	
	// Update transformation matrix for given angular velocity
	SHARPy::Wrap_xbeam_QuadSkew(4, &m_Vrel[3], &Temp0_Array[0]);
	Vmath::Svtvp(size*size, tmp, Temp0_Array, 1, Unit4_Array, 1, Temp1_Array, 1);
	SHARPy::Wrap_lu_invers(size, &Temp1_Array[0], &Temp2_Array[0]);
	Vmath::Svtvp(size, -1.0*tmp, Temp0_Array, 1, Unit4_Array, 1,  Temp1_Array, 1);
	SHARPy::Wrap_matmul(size, size, size, &Temp1_Array[0], &m_Quat_Array[0], &Temp3_Array[0]);
	SHARPy::Wrap_matmul(size, size, size, &Temp2_Array[0], &Temp3_Array[0], &m_Quat_Array[0]);
	SHARPy::Wrap_xbeam_Rot(&m_Quat_Array[0], &m_Cao_Array[0]);
	
	// Predictor step.
	Array<OneD,NekDouble> Temp_Array(m_NumDof, 0.0);
	Vmath::Smul(m_NumDof, (0.5-m_beta)*m_dt*m_dt, m_DXddt, 1, Temp_Array, 1);
	Vmath::Svtvp(m_NumDof, m_dt, m_DXdt, 1, Temp_Array, 1, Temp_Array, 1);
	Vmath::Vadd(m_NumDof, m_X, 1, Temp_Array, 1, m_X, 1);
	Vmath::Svtvp(m_NumDof, (1.0-m_gamma)*m_dt, m_DXddt, 1, m_DXdt, 1, m_DXdt, 1);
	Vmath::Zero(m_NumDof, m_DXddt, 1);

	int Iter;
	// Iteration until convergence.
	for(Iter = 0; Iter < m_MaxIterations; Iter++)
	{
		if (Iter == m_MaxIterations)
		{
			ASSERTL0(false,"Solution did not converge in SHARPy");
		}

		// Update nodal positions and velocities .
		SHARPy::Wrap_cbeam3_solv_state2disp(
			m_NumElems,m_NumNodes_tot,&m_NumNodes[0],
			&m_MemNo[0],&m_Conn[0],&m_Master_Array[0],
			&m_Length[0],&m_PreCurv[0],&m_Psi[0],
			&m_Vector[0],&m_Mass_Array[0],&m_Stiff_Array[0],
			&m_InvStiff_Array[0],&m_RBMass_Array[0],
			&m_Master[0],&m_Vdof[0],&m_Fdof[0],
			&m_PosIni_Array[0],&m_PsiIni_Array[0],
			m_NumDof,&m_X[0],&m_DXdt[0],&m_PosDefor_Array[0],
			&m_PsiDefor_Array[0],&m_PosDeforDot_Array[0],
			&m_PsiDeforDot_Array[0]);

		// Compute system functionals and matrices. (Use initial accelerations for Kgyr).
		Vmath::Zero(m_NumDof, m_Qglobal, 1);
		Vmath::Zero(m_NumDof*6, m_Mvel_Array, 1);
		Vmath::Zero(m_NumDof*6, m_Cvel_Array, 1);
		Vmath::Zero(m_NumDof2,m_Mglobal_Array,1);
		Vmath::Zero(m_NumDof2,m_Cglobal_Array,1);
		Vmath::Zero(m_NumDof2,m_Kglobal_Array,1);
		Vmath::Zero(m_NumDof2,m_Fglobal_Array,1);
		Vmath::Smul(m_NumNodes_tot*3,0.0,m_PosDefor_Array,1,m_PosDeforDDot_Array,1);
		Vmath::Smul(m_NumElems*m_MaxElNod*3,0.0,m_PsiDefor_Array,1,m_PsiDeforDDot_Array,1);
		Vmath::Vadd(m_NumNodes_tot*6,m_Force_Array,1,m_ForceStatic_Array,1,m_Force_Array,1);

		SHARPy::Wrap_cbeam3_asbly_dynamic(
			m_NumElems,m_NumNodes_tot,&m_NumNodes[0],
			&m_MemNo[0],&m_Conn[0],&m_Master_Array[0],&m_Length[0],
			&m_PreCurv[0],&m_Psi[0],&m_Vector[0],&m_Mass_Array[0],
			&m_Stiff_Array[0],&m_InvStiff_Array[0],&m_RBMass_Array[0],
			&m_Master[0],&m_Vdof[0],&m_Fdof[0],
			&m_PosIni_Array[0],&m_PsiIni_Array[0],
			&m_PosDefor_Array[0],&m_PsiDefor_Array[0],
			&m_PosDeforDot_Array[0],&m_PsiDeforDot_Array[0],
			&m_PosDeforDDot_Array[0],&m_PsiDeforDDot_Array[0],
			&m_Force_Array[0],&m_Vrel[0],&m_VrelDot[0],m_NumDof,m_DimMat,
			m_ms,&m_Mglobal_Array[0],&m_Mvel_Array[0],
			m_cs,&m_Cglobal_Array[0],&m_Cvel_Array[0],
			m_ks,&m_Kglobal_Array[0],
			m_fs,&m_Fglobal_Array[0],&m_Qglobal[0],
			m_FollowerForce, m_FollowerForceRig,        
			m_PrintInfo, m_OutInBframe, m_OutInaframe,    
			m_ElemProj, m_MaxIterations, m_NumLoadSteps,  
			m_NumGauss, m_Solution, m_DeltaCurved,        
			m_MinDelta, m_NewmarkDamp,        
			&m_Cao_Array[0]);

		// Compute admissible error.
		Array<OneD, NekDouble>tmp0(m_NumDof,0.0);
		Array<OneD, NekDouble>tmp1(m_NumDof,0.0);
		Vmath::Vabs(m_NumDof, m_Qglobal, 1, tmp0, 1);
		NekDouble MinDelta;
		MinDelta = m_MinDelta*std::max(1.0, Vmath::Vmax(m_NumDof,tmp0,1));

		// Compute the residual.
		Array<OneD,NekDouble> Temp4_Array(m_NumDof,0.0);
		Array<OneD,NekDouble> Temp5_Array(m_NumDof,0.0);
		Array<OneD,NekDouble> Temp6_Array(m_NumDof,0.0);
		Array<OneD,NekDouble> Temp7_Array(m_NumDof,0.0);

		DNekMatSharedPtr mglobal = 
			MemoryManager<DNekMat>::AllocateSharedPtr(
				m_NumDof,m_NumDof);

		for(int i = 0; i < m_NumDof; i++)
		{
			for(int j = 0; j < m_NumDof; j++)
			{
				(*mglobal)(i,j) = m_Mglobal_Array[i*m_NumDof+j];
			}
		}

		Blas::Dgemv('N', m_NumDof, m_NumDof,
					1.0, &(mglobal->GetPtr())[0],
					m_NumDof, &m_DXddt[0], 1,
					0.0, &Temp4_Array[0], 1);
		
		SHARPy::Wrap_matmul(
			m_NumDof,6,1,
			&m_Mvel_Array[0], &m_VrelDot[0], 
			&Temp5_Array[0]);
		
		Vmath::Vadd(m_NumDof, Temp4_Array, 1, Temp5_Array, 1, Temp5_Array, 1);
		Vmath::Vadd(m_NumDof, Temp5_Array, 1, m_Qglobal,   1, m_Qglobal,   1);

		SHARPy::Wrap_fem_m2v(
			m_NumNodes_tot,6,
			&m_Force_Array[0], 
			m_NumDof, &Temp6_Array[0], 
			&m_ListIN[0]);

		DNekMatSharedPtr fglobal = 
			MemoryManager<DNekMat>::AllocateSharedPtr(
				m_NumDof,m_NumDof);

		for(int i = 0; i < m_NumDof; i++)
		{
			for(int j = 0; j < m_NumDof; j++)
			{
				(*fglobal)(i,j)= m_Fglobal_Array[i*m_NumDof+j];
			}
		}

		Blas::Dgemv('N', m_NumDof, m_NumDof,
					1.0, &(fglobal->GetPtr())[0],
					m_NumDof, &Temp6_Array[0], 1,
					0.0, &Temp7_Array[0], 1);

		Vmath::Vsub(m_NumDof, m_Qglobal, 1, Temp7_Array, 1, m_Qglobal, 1);

		// Check convergence.
		Vmath::Vabs(m_NumDof, m_Qglobal, 1, tmp0, 1);
		Vmath::Vabs(m_NumDof, m_DX,      1, tmp1, 1);
		Vmath::Vadd(m_NumDof, tmp0,      1, tmp1, 1, tmp1, 1);

		if (Vmath::Vmax(m_NumDof,tmp1,1) < MinDelta)
		{
			if (m_PrintInfo)
			{
				cout<<"Subiteration = "<<Iter<<"; Delta = "
					<<Vmath::Vmax(m_NumDof,tmp0,1)<<endl;
			}
			break;
		}

		// Calculate Jacobian
		NekDouble factor1, factor2;
		factor1 = m_gamma/(m_beta*m_dt);
		factor2 = 1.0/(m_beta*m_dt*m_dt);
		
		Vmath::Zero(m_NumDof2,m_Asys_Array,1);
		Vmath::Vadd(m_NumDof2,m_Kglobal_Array,1,m_Asys_Array,1,m_Asys_Array,1);
		Vmath::Svtvp(m_NumDof2,factor1,m_Cglobal_Array,1,m_Asys_Array,1,m_Asys_Array,1);
		Vmath::Svtvp(m_NumDof2,factor2,m_Mglobal_Array,1,m_Asys_Array,1,m_Asys_Array,1);

		// Calculation of the correction.
		Vmath::Smul(m_NumDof, -1.0, m_Qglobal, 1, m_Qglobal, 1);

		Array<OneD, int> ipivot (m_NumDof);
		int info = 0;
		Lapack::Dgetrf(m_NumDof, m_NumDof, m_Asys_Array.get(), m_NumDof, ipivot.get(), info);
		if( info < 0 )
		{
			std::string message = "ERROR: The " + boost::lexical_cast<std::string>(-info) +
						"th parameter had an illegal parameter for dgetrf";
		}
		else if( info > 0 )
		{
			std::string message = "ERROR: Element u_" + boost::lexical_cast<std::string>(info) +
				boost::lexical_cast<std::string>(info) + " is 0 from dgetrf";
				ASSERTL0(false, message.c_str());
		}

		// N means no transponse (direct matrix)
		int ncolumns_b = 1;
		Vmath::Vcopy(m_NumDof, m_Qglobal,1,m_DX,1);
		Lapack::Dgetrs( 'N', m_NumDof, ncolumns_b, m_Asys_Array.get(), 
				m_NumDof, ipivot.get(), m_DX.get(), m_NumDof, info);
		if( info < 0 )
		{
			std::string message = "ERROR: The " + boost::lexical_cast<std::string>(-info) +
				"th parameter had an illegal parameter for dgetrf";
				ASSERTL0(false, message.c_str());
		}
		else if( info > 0 )
		{
			std::string message = "ERROR: Element u_" + boost::lexical_cast<std::string>(info) +
				boost::lexical_cast<std::string>(info) + " is 0 from dgetrf";
				ASSERTL0(false, message.c_str());
		}

		//Corrector step
		Vmath::Vadd(m_NumDof,  m_X, 1,  m_DX, 1, m_X, 1);
		Vmath::Svtvp(m_NumDof, factor1, m_DX, 1, m_DXdt,  1, m_DXdt,  1);
		Vmath::Svtvp(m_NumDof, factor2, m_DX, 1, m_DXddt, 1, m_DXddt, 1);
	}

	// Update nodal positions and velocities on the current converged time step.
	SHARPy::Wrap_cbeam3_solv_state2disp(
		m_NumElems,m_NumNodes_tot,&m_NumNodes[0],
		&m_MemNo[0],&m_Conn[0],&m_Master_Array[0],
		&m_Length[0],&m_PreCurv[0],&m_Psi[0],
		&m_Vector[0],&m_Mass_Array[0],
		&m_Stiff_Array[0],&m_InvStiff_Array[0],
		&m_RBMass_Array[0],&m_Master[0],&m_Vdof[0],&m_Fdof[0],
		&m_PosIni_Array[0],&m_PsiIni_Array[0],
		m_NumDof,&m_X[0],&m_DXdt[0],
		&m_PosDefor_Array[0],&m_PsiDefor_Array[0],
		&m_PosDeforDot_Array[0],&m_PsiDeforDot_Array[0]);

	//displacement and velocity for current step
	//Write output to export to main program (obsolete!).
	Array<OneD, NekDouble> temp(3*m_NumNodes_tot,0.0);
	Vmath::Vsub(3*m_NumNodes_tot,  m_PosDefor_Array, 1,m_PosIni_Array,1,temp,1);
	Vmath::Vcopy(m_NumNodes_tot, temp +  m_NumNodes_tot,1,motions[0][0],1);
	//Vmath::Vcopy(m_NumNodes_tot, temp +1*m_NumNodes_tot,1,motions[0][0],1);
	Vmath::Vcopy(m_NumNodes_tot, temp +2*m_NumNodes_tot,1,motions[1][0],1);
	Vmath::Vcopy(m_NumNodes_tot, m_PosDeforDot_Array +  m_NumNodes_tot,1,motions[0][1],1);
	//Vmath::Vcopy(m_NumNodes_tot, m_PosDeforDot_Array +1*m_NumNodes_tot,1,motions[1][1],1);
	Vmath::Vcopy(m_NumNodes_tot, m_PosDeforDot_Array +2*m_NumNodes_tot,1,motions[1][1],1);
	
	//nondimensionalize the displacement and velocity
	Vmath::Smul(m_NumNodes_tot,1.0/D,motions[0][0],1,motions[0][0],1);
	Vmath::Smul(m_NumNodes_tot,1.0/D,motions[1][0],1,motions[1][0],1);
	Vmath::Smul(m_NumNodes_tot,1.0/U_inf,motions[0][1],1,motions[0][1],1);
	Vmath::Smul(m_NumNodes_tot,1.0/U_inf,motions[1][1],1,motions[1][1],1);

	m_iStep++;
}
 
/**
 *
 */
void ForcingMovingBody::InitialiseVibrationModel(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields)
{
    m_session->LoadParameter("Strip_Z", m_nstrips,1);
    m_session->LoadParameter("VibModesZ", m_nv);
	m_session->LoadParameter("HomModesZ", m_nz);
    m_session->LoadParameter("LC", m_length);
    m_session->LoadParameter("Omega", m_omega);				//
    m_session->LoadParameter("Amp", m_A);      				//
    m_session->LoadParameter("WaveNum", m_wave_number); 	//
    m_session->LoadParameter("StructRho",  m_structrho);
    m_session->LoadParameter("StructDamp", m_structdamp, 0.0);
    m_session->LoadParameter("TimeStep", m_timestep, 0.01);
  
	m_spv = Array<OneD, Array<OneD, Array<OneD, NekDouble> > > (2);
	m_vib = Array<OneD, Array<OneD, Array<OneD, NekDouble> > > (2);
	for (int i = 0; i < m_motion.num_elements(); i++)
	{
		m_spv[i] = Array<OneD, Array<OneD, NekDouble> > (3);
		m_vib[i] = Array<OneD, Array<OneD, NekDouble> > (3);
		for (int j = 0; j < 3; j++)
		{
			m_spv[i][j] = Array<OneD, NekDouble>(m_nstrips*m_nz,0.0);
			m_vib[i][j] = Array<OneD, NekDouble> (m_nv+1, 0.0);
		}
	}

    Array<OneD, unsigned int> ZIDs;
    ZIDs = pFields[0]->GetZIDs();
    m_np = ZIDs.num_elements();

	// Set the DOF for vibration
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

    m_session->MatchSolverInfo("HomoStrip","True",m_IsHomostrip,false);
    m_session->MatchSolverInfo("FictitiousMassMethod", "True",
                                m_IsFictmass, false);
    if(m_IsFictmass)
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
                m_fV[i][n] = Array<OneD, NekDouble>(m_nstrips*m_nz, 0.0);
                m_fA[i][n] = Array<OneD, NekDouble>(m_nstrips*m_nz, 0.0);
            }
        }
    }

	// Get the matrices for struct solvers
    std::string StructDynSolver = 
            m_session->GetSolverInfo("StructDynSolver");
    if(boost::iequals(StructDynSolver, "ModalDecomposition"))
    {
    	SetModeMatrix();
		//
    	m_FFT = LibUtilities::GetNektarFFTFactory().CreateInstance(
                                               "NekFFTW", m_nv);
	}
	else if(boost::iequals(StructDynSolver, "SHARPy"))
	{
		InitialiseSHARPy();
	}

    m_movingBodyCalls = 0;

    // Set initial condition for cable's motion
    int cnt = 0;
    LibUtilities::CommSharedPtr vcomm = pFields[0]->GetComm();
    int colrank = vcomm->GetColumnComm()->GetRank();
    for(int j = 0; j < m_funcName.num_elements(); j++)
    {
        // loading from the specified files through inputstream
        if(m_IsFromFile[cnt] && m_IsFromFile[cnt+1])
        {
            std::ifstream inputStream;

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
                for (int n = 0; n < m_nv+1; n++)
                {
                    inputStream >> setprecision(6) >> time;
                    inputStream >> setprecision(6) >> z_cds;
					for(int i = 0; i < 2; i++)
					{
						for (int j = 0; j < 3; j++)
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
                Array<OneD, NekDouble> x0(m_nv+1, 0.0);
                Array<OneD, NekDouble> x1(m_nv+1, 0.0);
                Array<OneD, NekDouble> x2(m_nv+1, 0.0);
                Array<OneD, NekDouble> tmp(m_nv+1,0.0);

                for (int p = 0; p < m_nv+1; p++)
                {
                    x2[p] = m_length/m_nv*p;
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
void ForcingMovingBody::SetModeMatrix()
{
    int nplanes;

    nplanes = m_nv; 

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
	m_mapfuncName= Array<OneD, std::string> (2);
    m_motion = Array<OneD, std::string> (2);
    m_motion[0] = "x";
    m_motion[1] = "y";

    m_IsFromFile = Array<OneD, bool> (6);
    m_IsMapFromFile = Array<OneD, bool> (2);
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

    //Get the Mapping:
    const TiXmlElement* funcNameElmt_M
                    = pForce->FirstChildElement("MAPPING");
    ASSERTL0(funcNameElmt_M,
             "MOVINGBODYFORCE tag has to specify a function name which "
             "prescribes the body displacement as m(z,t).");

    m_mapfuncName[0] = funcNameElmt_M->GetText();
    ASSERTL0(m_session->DefinesFunction(m_mapfuncName[0]),
             "Function '" + m_mapfuncName[0] + "' not defined.");
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


    // Check Mappling x
    vType = m_session->GetFunctionType(m_mapfuncName[0],m_motion[0]);
    if (vType == LibUtilities::eFunctionTypeExpression)
    {
        m_IsMapFromFile[0] = false;
    }
    else if (vType == LibUtilities::eFunctionTypeFile)
    {
        m_IsMapFromFile[0] = true;
    }
    else
    {
        ASSERTL0(false, "The mapping in x must be specified via an "
                        "equation <E> or a file <F>");
    }

    // Check Mappling y
    vType = m_session->GetFunctionType(m_mapfuncName[0],m_motion[1]);
    if (vType == LibUtilities::eFunctionTypeExpression)
    {
        m_IsMapFromFile[1] = false;
    }
    else if (vType == LibUtilities::eFunctionTypeFile)
    {
        m_IsMapFromFile[1] = true;
    }
    else
    {
        ASSERTL0(false, "The mapping in y must be specified via an "
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

/**
 *
 */
void ForcingMovingBody::Input_node_SHARPy()
{	
	//Initial position vector of grid points.
	for(int i = 0; i < m_NumNodes_tot; i++)
	{
		m_PosIni_Array[i] = m_BeamLength*i/(m_NumNodes_tot-1);
	}
	//Initial pretwist angle.
	Vmath::Zero(m_NumNodes_tot,m_PhiNodes,1);

	//Static point forces.
	Vmath::Zero(6*m_NumNodes_tot,m_ForceStatic_Array,1);

	//Boundary conditions.
	Vmath::Zero(m_NumNodes_tot,m_BoundConds,1);
	m_BoundConds[0] = 1; //Clamp
	m_BoundConds[m_NumNodes_tot-1] = 1; //Clamp
}

/**
 *
 **/
void ForcingMovingBody::Input_elem_SHARPy()
{
	// Connectivies.
	Vmath::Zero(m_MaxElNod*m_NumElems,m_Conn,1);
	switch(m_NumNodesElem)  
    {   
    	case 2:
			{  
    			for(int i = 0; i < m_NumElems; i++)
				{
            		m_Conn[m_MaxElNod*i]   = i+1;
            		m_Conn[m_MaxElNod*i+1] = i+2;
            		m_NumNodes[i] = 2;
				}
				m_NumNodes_tot    = m_NumElems+1;
			}	
    		break;  
    	case 3:  
			{  
    			for(int i = 0; i < m_NumElems; i++)
				{
        		    m_Conn[m_MaxElNod*i]   = 2*i+1;
            		m_Conn[m_MaxElNod*i+1] = 2*i+3;
            		m_Conn[m_MaxElNod*i+2] = 2*i+2;
            		m_NumNodes[i] = 3;
				}
				m_NumNodes_tot    = 2*m_NumElems+1;
			}
    		break;  
        default:
            ASSERTL0(false,"Number of Nodes in each element is not correct in SHARPy");
    } 
	
	// Store element stiffness/mass (constant)

	NekDouble I,A,mu,E,G;

    m_session->LoadParameter("MomentofInertia_SHARPy",    I);
    m_session->LoadParameter("AreaofCrossSection_SHARPy", A);
    m_session->LoadParameter("YoungModulus_SHARPy",       E);
    m_session->LoadParameter("PoissonRatio_SHARPy",       mu);

    G = E/(2.0*(1.0+mu));

	NekDouble BeamMass[6][6];
	NekDouble BeamMassMatrix[6*m_NumElems][6];
	NekDouble BeamStiffness[6][6];
	NekDouble BeamStiffnessMatrix[6*m_NumElems][6];
	NekDouble BeamInvStiffnessMatrix[6*m_NumElems][6];

	for(int i = 0; i < 6; i++)
	{
		for(int j = 0; j < 6; j++)
		{
			BeamStiffness[i][j] = 0.0;
			BeamMass[i][j]      = 0.0;
		}
	}	
	
	for(int i = 0; i < 6; i++)
	{
		for(int j = 0; j < 6; j++)
		{
			BeamMass[i][j] 		= 0.0;
			BeamStiffness[i][j] = 0.0;
		}
	}

    BeamMass[0][0]  	= m_structrho*A;
    BeamMass[1][1]  	= BeamMass[0][0];
    BeamMass[2][2]  	= BeamMass[0][0];
    BeamMass[3][3]  	= G*I*2.0;
    BeamMass[4][4]  	= G*I;
    BeamMass[5][5]  	= G*I;
    BeamStiffness[0][0] = E*A;
    BeamStiffness[1][1] = G*A;
    BeamStiffness[2][2] = G*A;
    BeamStiffness[3][3] = 2.0*G*I;
    BeamStiffness[4][4] = E*I;
    BeamStiffness[5][5] = E*I;	
    
	int size = 6;
	Array<OneD, NekDouble> temp1(size*size,0.0);
	Array<OneD, NekDouble> temp2(size*size,0.0);
  
	for(int i = 0; i < size; i++)
	{
		for(int j = 0; j < size; j++)
		{
			temp1[size*i+j] = BeamStiffness[i][j];
		}
	}	
		
	SHARPy::Wrap_lu_invers(size, &temp1[0], &temp2[0]);

  	for (int i = 0; i < m_NumElems; i++)
	{
		for(int j = 0; j < size; j++)
		{
			for(int k = 0; k < size; k++)
			{
    			BeamStiffnessMatrix[i*size+j][k]     = BeamStiffness[j][k];
    			BeamMassMatrix[i*size+j][k]          = BeamMass[j][k];
    			BeamInvStiffnessMatrix[i*size+j][k]  = temp2[size*j+k];
			}
		}
  	}
	for(int i = 0; i < m_NumElems*size; i++)
	{
		for(int j = 0; j < size; j++)
		{
			int cnt = i+m_NumElems*size*j;
			m_Stiff_Array[cnt]    = BeamStiffnessMatrix[i][j];
			m_Mass_Array[cnt]     = BeamMassMatrix[i][j];
			m_InvStiff_Array[cnt] = BeamInvStiffnessMatrix[i][j];
		}
	}
	// Define lumped masses at element nodes.
	for (int i = 0; i < m_NumElems; i++)
	{
		 m_RBMass_Array[i] = 0.0;		
	}
	// Element orientation.
	for (int i = 0; i < m_NumElems; i++)
	{
		m_Vector[3*i+1] = 1.0;
	}
	// Define element types.
	//for (int i = 0; i < m_NumElems; i++)
	//{
		//m_MemNo[i] = 0;
	//}
}
}
