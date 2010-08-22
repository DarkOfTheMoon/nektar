///////////////////////////////////////////////////////////////////////////////
//
// File UnsteadySystem.cpp
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
// Description: Generic timestepping for Unsteady solvers
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>

#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>
#include <ADRSolver/EquationSystems/UnsteadySystem.h>

namespace Nektar
{
    /**
     * @class UnsteadySystem
     *
     * Provides the underlying timestepping framework for unsteady solvers
     * including the general timestepping routines. This class is not intended
     * to be directly instantiated, but rather is a base class on which to
     * define unsteady solvers.
     *
     * For details on implementing unsteady solvers see
     * \ref sectionADRSolverModuleImplementation here
     */

    /**
     * Processes SolverInfo parameters from the session file and sets up
     * timestepping-specific code.
     * @param   pSession        Session object to read parameters from.
     */
    UnsteadySystem::UnsteadySystem(SessionReaderSharedPtr& pSession)
        : EquationSystem(pSession)
    {
        // Load SolverInfo parameters
        pSession->MatchSolverInfo("DIFFUSIONADVANCEMENT","Explicit",
                                    m_explicitDiffusion,true);
        pSession->MatchSolverInfo("ADVECTIONADVANCEMENT","Explicit",
                                    m_explicitAdvection,true);
        pSession->MatchSolverInfo("REACTIONADVANCEMENT", "Explicit",
                                    m_explicitReaction, true);

        // Determine TimeIntegrationMethod to use.
        ASSERTL0(pSession->DefinesSolverInfo("TIMEINTEGRATIONMETHOD"),
                "No TIMEINTEGRATIONMETHOD defined in session.");
        int i;
        for (i = 0; i < (int)LibUtilities::SIZE_TimeIntegrationMethod; ++i)
        {
            bool match;
            pSession->MatchSolverInfo("TIMEINTEGRATIONMETHOD",
                    LibUtilities::TimeIntegrationMethodMap[i], match, false);
            if (match)
            {
                m_timeIntMethod = (LibUtilities::TimeIntegrationMethod) i;
                break;
            }
        }
        ASSERTL0(i != (int) LibUtilities::SIZE_TimeIntegrationMethod,
                                        "Invalid time integration type.");

        // Load generic input parameters
        pSession->LoadParameter("IO_InfoSteps", m_infosteps, 0);
    }


    /**
     *
     */
    UnsteadySystem::~UnsteadySystem()
    {

    }


    /**
     * Initialises the time integration scheme (as specified in the session
     * file), and perform the time integration.
     */
    void UnsteadySystem::v_DoSolve()
    {
        int i,n,nchk = 0;
        int ncoeffs = m_fields[0]->GetNcoeffs();
        int nvariables = m_fields.num_elements();

        // Set up wrapper to fields data storage.
        Array<OneD, Array<OneD, NekDouble> >   fields(nvariables);
        Array<OneD, Array<OneD, NekDouble> >   tmp(nvariables);

        for(i = 0; i < nvariables; ++i)
        {
            m_fields[i]->SetPhysState(false);
            fields[i]  = m_fields[i]->UpdatePhys();
        }

        // Declare an array of TimeIntegrationSchemes For multi-stage
        // methods, this array will have just one entry containing the
        // actual multi-stage method...
        // For multi-steps method, this can have multiple entries
        //  - the first scheme will used for the first timestep (this
        //    is an initialization scheme)
        //  - the second scheme will used for the second timestep
        //    (this is an initialization scheme)
        //  - ...
        //  - the last scheme will be used for all other time-steps
        //    (this will be the actual scheme)

        Array<OneD, LibUtilities::TimeIntegrationSchemeSharedPtr> IntScheme;
        LibUtilities::TimeIntegrationSolutionSharedPtr u;
        int numMultiSteps;

        switch(m_timeIntMethod)
        {
        case LibUtilities::eIMEXdirk_2_3_2:
        case LibUtilities::eIMEXdirk_3_4_3:
        case LibUtilities::eDIRKOrder2:
        case LibUtilities::eDIRKOrder3:
        case LibUtilities::eBackwardEuler:
        case LibUtilities::eForwardEuler:
        case LibUtilities::eClassicalRungeKutta4:
        case LibUtilities::eIMEXOrder1:
            {
                numMultiSteps = 1;

                IntScheme = Array<OneD, LibUtilities::TimeIntegrationSchemeSharedPtr>(numMultiSteps);

                LibUtilities::TimeIntegrationSchemeKey IntKey(m_timeIntMethod);
                IntScheme[0] = LibUtilities::TimeIntegrationSchemeManager()[IntKey];

                u = IntScheme[0]->InitializeScheme(m_timestep,fields,m_time,m_ode);
                break;
            }
        case LibUtilities::eAdamsBashforthOrder2:
            {
                numMultiSteps = 2;

                IntScheme = Array<OneD, LibUtilities::TimeIntegrationSchemeSharedPtr>(numMultiSteps);

                // Used in the first time step to initalize the scheme
                LibUtilities::TimeIntegrationSchemeKey IntKey0(LibUtilities::eForwardEuler);

                // Used for all other time steps
                LibUtilities::TimeIntegrationSchemeKey IntKey1(m_timeIntMethod);
                IntScheme[0] = LibUtilities::TimeIntegrationSchemeManager()[IntKey0];
                IntScheme[1] = LibUtilities::TimeIntegrationSchemeManager()[IntKey1];

                // Initialise the scheme for the actual time integration scheme
                u = IntScheme[1]->InitializeScheme(m_timestep,fields,m_time,m_ode);
                break;
            }
        case LibUtilities::eIMEXOrder2:
            {
                numMultiSteps = 2;

                IntScheme = Array<OneD, LibUtilities::TimeIntegrationSchemeSharedPtr>(numMultiSteps);

                // Used in the first time step to initalize the scheme
                LibUtilities::TimeIntegrationSchemeKey IntKey0(LibUtilities::eIMEXOrder1);

                // Used for all other time steps
                LibUtilities::TimeIntegrationSchemeKey IntKey1(m_timeIntMethod);
                IntScheme[0] = LibUtilities::TimeIntegrationSchemeManager()[IntKey0];
                IntScheme[1] = LibUtilities::TimeIntegrationSchemeManager()[IntKey1];

                // Initialise the scheme for the actual time integration scheme
                u = IntScheme[1]->InitializeScheme(m_timestep,fields,m_time,m_ode);
                break;
            }

        case LibUtilities::eIMEXOrder3:
            {
                numMultiSteps = 3;

                IntScheme = Array<OneD, LibUtilities::TimeIntegrationSchemeSharedPtr>(numMultiSteps);

                // Used in the first time step to initalize the scheme
                LibUtilities::TimeIntegrationSchemeKey IntKey0(LibUtilities::eIMEXOrder1);
                LibUtilities::TimeIntegrationSchemeKey IntKey1(LibUtilities::eIMEXOrder2);

                // Used for all other time steps
                LibUtilities::TimeIntegrationSchemeKey IntKey2(m_timeIntMethod);
                IntScheme[0] = LibUtilities::TimeIntegrationSchemeManager()[IntKey0];
                IntScheme[1] = LibUtilities::TimeIntegrationSchemeManager()[IntKey1];
                IntScheme[2] = LibUtilities::TimeIntegrationSchemeManager()[IntKey2];

                // Initialise the scheme for the actual time integration scheme
                u = IntScheme[2]->InitializeScheme(m_timestep,fields,m_time,m_ode);
                break;
            }
            default:
            {
                ASSERTL0(false,"populate switch statement for integration scheme");
            }
        }

        std::string outname = m_session->GetFilename() + ".his";
        std::ofstream hisFile (outname.c_str());

        // Perform integration in time.
        for(n = 0; n < m_steps; ++n)
        {
            // Integrate over timestep.
            if( n < numMultiSteps-1)
            {
                // Use initialisation schemes if time step is less than the
                // number of steps in the scheme.
                fields = IntScheme[n]->TimeIntegrate(m_timestep,u,m_ode);
            }
            else
            {
                fields = IntScheme[numMultiSteps-1]->TimeIntegrate(m_timestep,u,m_ode);
            }

            // Increment time.
            m_time += m_timestep;

            // Write out status information.
            if(!((n+1)%m_infosteps))
            {
                cout << "Steps: " << n+1 << "\t Time: " << m_time << "\t " << endl;
//                cout << "\r" << setw(3) << int((NekDouble)n/m_steps*100) << "%:\t"
//                     << setw(7) << "Steps: " << setw(6) << n+1
//                     << setw(7) << "Time: "
//                     << setw(8) << fixed << setprecision(2) << m_time
//                     << "\t " << flush;
            }

            // Write out checkpoint files.
            if(n&&(!((n+1)%m_checksteps)))
            {
                for(i = 0; i < nvariables; ++i)
                {
                    m_fields[i]->FwdTrans(fields[i],m_fields[i]->UpdateCoeffs());
                }
                Checkpoint_Output(nchk++);
                WriteHistoryData(hisFile);
            }
        }

        // At the end of the time integration, store final solution.
        for(i = 0; i < nvariables; ++i)
        {
            m_fields[i]->UpdatePhys() = fields[i];
        }
    }


    /**
     *
     */
    void UnsteadySystem::v_DoInitialise()
    {
        SetInitialConditions();
    }


    /**
     *
     */
    void UnsteadySystem::v_PrintSummary(std::ostream &out)
    {
        EquationSystem::v_PrintSummary(out);
        out << "\tAdvection       : " << (m_explicitAdvection ? "explicit" : "implicit") << endl;
        out << "\tDiffusion       : " << (m_explicitDiffusion ? "explicit" : "implicit") << endl;
        out << "\tReaction        : " << (m_explicitReaction  ? "explicit" : "implicit") << endl;
        out << "\tIntegration Type: " << LibUtilities::TimeIntegrationMethodMap[m_timeIntMethod] << endl;
        out << "\tTime Step       : " << m_timestep << endl;
        out << "\tNo. of Steps    : " << m_steps << endl;
        out << "\tCheckpoints     : " << m_checksteps << " steps" << endl;
    }


    /**
     *
     */
    void UnsteadySystem::v_NumericalFlux(
                Array<OneD, Array<OneD, NekDouble> > &physfield,
                Array<OneD, Array<OneD, NekDouble> > &numflux)
    {
        ASSERTL0(false, "This function is not implemented for this equation.");
    }


    /**
     *
     */
    void UnsteadySystem::v_NumericalFlux(
                Array<OneD, Array<OneD, NekDouble> > &physfield,
                Array<OneD, Array<OneD, NekDouble> > &numfluxX,
                Array<OneD, Array<OneD, NekDouble> > &numfluxY )
    {
        ASSERTL0(false, "This function is not implemented for this equation.");
    }


    /**
     *
     */
    void UnsteadySystem::v_NumFluxforScalar(
                Array<OneD, Array<OneD, NekDouble>  > &ufield,
                Array<OneD, Array<OneD, Array<OneD,NekDouble> > > &uflux)
    {
        int i,j;
        int nTraceNumPoints = GetTraceNpoints();
        int nvariables = m_fields.num_elements();
        int nqvar = uflux.num_elements();

        Array<OneD, NekDouble > Fwd(nTraceNumPoints);
        Array<OneD, NekDouble > Bwd(nTraceNumPoints);
        Array<OneD, NekDouble > Vn (nTraceNumPoints,0.0);
        Array<OneD, NekDouble > fluxtemp (nTraceNumPoints,0.0);

        // Get the sign of (v \cdot n), v = an arbitrary vector

        //  Evaulate upwind flux of uflux = \hat{u} \phi \cdot u = u^{(+,-)} n
        for (j = 0; j < nqvar; ++j)
        {
            for(i = 0; i < nvariables ; ++i)
            {
                //  Compute Forward and Backward value of ufield of i direction

                m_fields[i]->GetFwdBwdTracePhys(ufield[i],Fwd,Bwd);

                // if Vn >= 0, flux = uFwd, i.e.,
                //  edge::eForward, if V*n>=0 <=> V*n_F>=0, pick uflux = uFwd
                //  edge::eBackward, if V*n>=0 <=> V*n_B<0, pick uflux = uFwd

                // else if Vn < 0, flux = uBwd, i.e.,
                //  edge::eForward, if V*n<0 <=> V*n_F<0, pick uflux = uBwd
                //  edge::eBackward, if V*n<0 <=> V*n_B>=0, pick uflux = uBwd

                m_fields[i]->GetTrace()->Upwind(m_traceNormals[j],Fwd,Bwd,fluxtemp);

                // Imposing weak boundary condition with flux
                // if Vn >= 0, uflux = uBwd at Neumann, i.e.,
                //  edge::eForward, if V*n>=0 <=> V*n_F>=0, pick uflux = uBwd
                //  edge::eBackward, if V*n>=0 <=> V*n_B<0, pick uflux = uBwd

                // if Vn >= 0, uflux = uFwd at Neumann, i.e.,
                //  edge::eForward, if V*n<0 <=> V*n_F<0, pick uflux = uFwd
                //  edge::eBackward, if V*n<0 <=> V*n_B>=0, pick uflux = uFwd

                if(m_fields[0]->GetBndCondExpansions().num_elements())
                {
                    WeakPenaltyforScalar(i,ufield[i],fluxtemp);
                }

                // if Vn >= 0, flux = uFwd*(tan_{\xi}^- \cdot \vec{n} ), i.e,
                // edge::eForward, uFwd \(\tan_{\xi}^Fwd \cdot \vec{n} )
                // edge::eBackward, uFwd \(\tan_{\xi}^Bwd \cdot \vec{n} )

                // else if Vn < 0, flux = uBwd*(tan_{\xi}^- \cdot \vec{n} ), i. e,
                // edge::eForward, uBwd \(\tan_{\xi}^Fwd \cdot \vec{n} )
                // edge::eBackward, uBwd \(\tan_{\xi}^Bwd \cdot \vec{n} )

                Vmath::Vmul(nTraceNumPoints,m_traceNormals[j],1,fluxtemp,1,uflux[j][i],1);

            }
        }
    }


    /**
     *
     */
    void UnsteadySystem::v_NumFluxforVector(
                Array<OneD, Array<OneD, NekDouble> > &ufield,
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &qfield,
                Array<OneD, Array<OneD, NekDouble> > &qflux)
    {
        int nTraceNumPoints = GetTraceNpoints();
        int nvariables = m_fields.num_elements();
        int nqvar = qfield.num_elements();

        NekDouble C11 = 1.0;
        Array<OneD, NekDouble > Fwd(nTraceNumPoints);
        Array<OneD, NekDouble > Bwd(nTraceNumPoints);
        Array<OneD, NekDouble > Vn (nTraceNumPoints,0.0);

        Array<OneD, NekDouble > qFwd(nTraceNumPoints);
        Array<OneD, NekDouble > qBwd(nTraceNumPoints);
        Array<OneD, NekDouble > qfluxtemp(nTraceNumPoints,0.0);

        Array<OneD, NekDouble > uterm(nTraceNumPoints);

        // Evaulate upwind flux of qflux = \hat{q} \cdot u = q \cdot n - C_(11)*(u^+ - u^-)
        for(int i = 0; i < nvariables; ++i)
        {
            qflux[i] = Array<OneD, NekDouble> (nTraceNumPoints,0.0);
            for(int j = 0; j < nqvar; ++j)
            {
                //  Compute Forward and Backward value of ufield of jth         direction
                m_fields[i]->GetFwdBwdTracePhys(qfield[j][i],qFwd,qBwd);

                // if Vn >= 0, flux = uFwd, i.e.,
                //  edge::eForward, if V*n>=0 <=> V*n_F>=0, pick qflux = qBwd = q+
                //  edge::eBackward, if V*n>=0 <=> V*n_B<0, pick qflux = qBwd = q-

                // else if Vn < 0, flux = uBwd, i.e.,
                //  edge::eForward, if V*n<0 <=> V*n_F<0, pick qflux = qFwd = q-
                //  edge::eBackward, if V*n<0 <=> V*n_B>=0, pick qflux = qFwd   =q+

                m_fields[i]->GetTrace()->Upwind(m_traceNormals[j],qBwd,qFwd,qfluxtemp);
                Vmath::Vmul(nTraceNumPoints,m_traceNormals[j],1,qfluxtemp,1,qfluxtemp,1);

                // Generate Stability term = - C11 ( u- - u+ )
                m_fields[i]->GetFwdBwdTracePhys(ufield[i],Fwd,Bwd);
                Vmath::Vsub(nTraceNumPoints,Fwd,1,Bwd,1,uterm,1);
                Vmath::Smul(nTraceNumPoints,-1.0*C11,uterm,1,uterm,1);

                //  Flux = {Fwd,Bwd}*(nx,ny,nz) + uterm*(nx,ny)
                Vmath::Vadd(nTraceNumPoints,uterm,1,qfluxtemp,1,qfluxtemp,1);

                // Imposing weak boundary condition with flux
                if(m_fields[0]->GetBndCondExpansions().num_elements())
                {
                    WeakPenaltyforVector(i,j,qfield[j][i],qfluxtemp,C11);
                }

                // q_hat \cdot n = (q_xi \cdot n_xi) or (q_eta \cdot n_eta)
                // n_xi = n_x*tan_xi_x + n_y*tan_xi_y + n_z*tan_xi_z
                // n_xi = n_x*tan_eta_x + n_y*tan_eta_y + n_z*tan_eta_z
                Vmath::Vadd(nTraceNumPoints,qfluxtemp,1,qflux[i],1,qflux[i],1);
            }
        }
    }


    /**
     *
     */
    void UnsteadySystem::v_GetFluxVector(const int i, const int j,
                Array<OneD, Array<OneD, NekDouble> > &physfield,
                Array<OneD, Array<OneD, NekDouble> > &flux)
    {
        for(int k = 0; k < flux.num_elements(); ++k)
        {
            Vmath::Zero(GetNpoints(),flux[k],1);
        }
        Vmath::Vcopy(GetNpoints(),physfield[i],1,flux[j],1);
    }


    /**
     *
     */
    void UnsteadySystem::WeakPenaltyforScalar(const int var,
                const Array<OneD, const NekDouble> &physfield,
                      Array<OneD, NekDouble> &penaltyflux,
                      NekDouble time)
    {
        int i, j, e, npoints, id1, id2;
        // Number of boundary regions
        int nbnd = m_fields[var]->GetBndCondExpansions().num_elements();
        int Nfps, numBDEdge;
        int nTraceNumPoints = GetTraceNpoints();
        int cnt = 0;

        Array<OneD, NekDouble > uplus(nTraceNumPoints);

        m_fields[var]->ExtractTracePhys(physfield,uplus);
        for(i = 0; i < nbnd; ++i)
        {
            // Number of boundary expansion related to that region
            numBDEdge = m_fields[var]->GetBndCondExpansions()[i]->GetExpSize();
            // Evaluate boundary values g_D or g_N from input files
            SpatialDomains::ConstInitialConditionShPtr ifunc =  m_boundaryConditions->GetInitialCondition(0);
            npoints = m_fields[var]->GetBndCondExpansions()[i]->GetNpoints();
            Array<OneD,NekDouble> BDphysics(npoints);
            Array<OneD,NekDouble> x0(npoints,0.0);
            Array<OneD,NekDouble> x1(npoints,0.0);
            Array<OneD,NekDouble> x2(npoints,0.0);

            m_fields[var]->GetBndCondExpansions()[i]->GetCoords(x0,x1,x2);
            for(j = 0; j < npoints; j++)
            {
                BDphysics[j] = ifunc->Evaluate(x0[j],x1[j],x2[j],time);
            }

            // Weakly impose boundary conditions by modifying flux values
            for (e = 0; e < numBDEdge ; ++e)
            {
                // Number of points on the expansion
                Nfps = m_fields[var]->GetBndCondExpansions()[i]->GetExp(e)->GetNumPoints(0) ;
                id1 = m_fields[var]->GetBndCondExpansions()[i]->GetPhys_Offset(e);
                id2 = m_fields[0]->GetTrace()->GetPhys_Offset(m_fields[0]->GetTraceMap()->GetBndCondTraceToGlobalTraceMap(cnt++));

                // For Dirichlet boundary condition: uflux = g_D
                if(m_fields[var]->GetBndConditions()[i]->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                {
                    Vmath::Vcopy(Nfps,&BDphysics[id1],1,&penaltyflux[id2],1);
                }
                // For Neumann boundary condition: uflux = u+
                else if((m_fields[var]->GetBndConditions()[i])->GetBoundaryConditionType() == SpatialDomains::eNeumann)
                {
                    Vmath::Vcopy(Nfps,&uplus[id2],1,&penaltyflux[id2],1);
                }
            }
        }
    }


    /**
     * Diffusion: Imposing weak boundary condition for q with flux
     *  uflux = g_D  on Dirichlet boundary condition
     *  uflux = u_Fwd  on Neumann boundary condition
     */
    void UnsteadySystem::WeakPenaltyforVector(
                const int var,
                const int dir,
                const Array<OneD, const NekDouble> &physfield,
                      Array<OneD, NekDouble> &penaltyflux,
                      NekDouble C11,
                      NekDouble time)
    {
        int i, j, e, npoints, id1, id2;
        int nbnd = m_fields[var]->GetBndCondExpansions().num_elements();
        int numBDEdge, Nfps;
        int nTraceNumPoints = GetTraceNpoints();
        Array<OneD, NekDouble > uterm(nTraceNumPoints);
        Array<OneD, NekDouble > qtemp(nTraceNumPoints);
        int cnt = 0;

        m_fields[var]->ExtractTracePhys(physfield,qtemp);

        for(i = 0; i < nbnd; ++i)
        {
            numBDEdge = m_fields[var]->GetBndCondExpansions()[i]->GetExpSize();
            // Evaluate boundary values g_D or g_N from input files
            SpatialDomains::ConstInitialConditionShPtr ifunc = m_boundaryConditions->GetInitialCondition(0);
            npoints = m_fields[var]->GetBndCondExpansions()[i]->GetNpoints();

            Array<OneD,NekDouble> BDphysics(npoints);
            Array<OneD,NekDouble> x0(npoints,0.0);
            Array<OneD,NekDouble> x1(npoints,0.0);
            Array<OneD,NekDouble> x2(npoints,0.0);

            m_fields[var]->GetBndCondExpansions()[i]->GetCoords(x0,x1,x2);
            for(j = 0; j < npoints; j++)
            {
                BDphysics[j] = ifunc->Evaluate(x0[j],x1[j],x2[j],time);
            }

            // Weakly impose boundary conditions by modifying flux values
            for (e = 0; e < numBDEdge ; ++e)
            {
                Nfps = m_fields[var]->GetBndCondExpansions()[i]->GetExp(e)->GetNumPoints(0);

                id1 = m_fields[var]->GetBndCondExpansions()[i]->GetPhys_Offset(e);
                id2 = m_fields[0]->GetTrace()->GetPhys_Offset(m_fields[0]->GetTraceMap()->GetBndCondTraceToGlobalTraceMap(cnt++));

                // For Dirichlet boundary condition: qflux = q+ - C_11 (u+ -    g_D) (nx, ny)
                if(m_fields[var]->GetBndConditions()[i]->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                {
                    Vmath::Vmul(Nfps,&m_traceNormals[dir][id2],1,&qtemp[id2],1, &penaltyflux[id2],1);
                }
                // For Neumann boundary condition: qflux = g_N
                else if((m_fields[var]->GetBndConditions()[i])->GetBoundaryConditionType() == SpatialDomains::eNeumann)
                {
                    Vmath::Vmul(Nfps,&m_traceNormals[dir][id2],1,&BDphysics[id1],1,&penaltyflux[id2],1);
                }
            }
        }
    }
}