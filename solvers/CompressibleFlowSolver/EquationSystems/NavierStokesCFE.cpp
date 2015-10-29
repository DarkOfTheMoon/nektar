///////////////////////////////////////////////////////////////////////////////
//
// File NavierStokesCFE.cpp
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
// Description: Navier Stokes equations in conservative variables
//
///////////////////////////////////////////////////////////////////////////////

#include <CompressibleFlowSolver/EquationSystems/NavierStokesCFE.h>

namespace Nektar
{
    string NavierStokesCFE::className =
        SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
            "NavierStokesCFE", NavierStokesCFE::create,
            "NavierStokes equations in conservative variables.");

    NavierStokesCFE::NavierStokesCFE(
            const LibUtilities::SessionReaderSharedPtr& pSession)
        : CompressibleFlowSystem(pSession)
    {
    }

    void NavierStokesCFE::v_InitObject()
    {
        CompressibleFlowSystem::v_InitObject();

        if(m_session->DefinesSolverInfo("PROBLEMTYPE"))
        {

            std::string ProblemTypeStr = m_session->
                                         GetSolverInfo("PROBLEMTYPE");
            for(int i = 0; i < (int) SIZE_ProblemType; ++i)
            {
                if (NoCaseStringCompare(ProblemTypeMap[i], ProblemTypeStr) == 0)
                {
                    m_problemType = (ProblemType)i;
                    break;
                }
            }
        }
        else
        {
            m_problemType = (ProblemType)0;
        }

        if ((m_explicitAdvection == true) && (m_explicitDiffusion == true))
        {
            m_ode.DefineOdeRhs     (&NavierStokesCFE::DoOdeRhs,        this);
            m_ode.DefineProjection (&NavierStokesCFE::DoOdeProjection, this);
        }
        else if ((m_explicitAdvection == true) &&
                 (m_explicitDiffusion == false))
        {
            m_ode.DefineImplicitSolve(&NavierStokesCFE::DoImplicitSolve, this);
            m_ode.DefineOdeRhs       (&NavierStokesCFE::DoOdeRhs,        this);
        }
        else
        {
            ASSERTL0(false, "Fully implicit compressible Navier-Stokes "
                     "equations not implemented.");
        }
    }

    NavierStokesCFE::~NavierStokesCFE()
    {

    }

    void NavierStokesCFE::v_GenerateSummary(SolverUtils::SummaryList& s)
    {
        CompressibleFlowSystem::v_GenerateSummary(s);
        SolverUtils::AddSummaryItem(s, "Problem Type",
                                    ProblemTypeMap[m_problemType]);
    }

    void NavierStokesCFE::v_SetInitialConditions(
        NekDouble initialtime,
        bool dumpInitialConditions,
        const int domain)
    {
        EquationSystem::v_SetInitialConditions(initialtime, false);

        // insert white noise in initial condition
        NekDouble Noise;
        int phystot = m_fields[0]->GetTotPoints();
        Array<OneD, NekDouble> noise(phystot);

        m_session->LoadParameter("Noise", Noise,0.0);
        int m_nConvectiveFields =  m_fields.num_elements();

        if (Noise > 0.0)
        {
            for (int i = 0; i < m_nConvectiveFields; i++)
            {
                Vmath::FillWhiteNoise(phystot, Noise, noise, 1,
                                      m_comm->GetColumnComm()->GetRank()+1);
                Vmath::Vadd(phystot, m_fields[i]->GetPhys(), 1,
                            noise, 1, m_fields[i]->UpdatePhys(), 1);
                m_fields[i]->FwdTrans_IterPerExp(m_fields[i]->GetPhys(),
                                                 m_fields[i]->UpdateCoeffs());
            }
        }

        CompressibleFlowSystem::v_SetInitialConditions();

        if (dumpInitialConditions)
        {
            // Dump initial conditions to file
            Checkpoint_Output(0);
        }
    }

    void NavierStokesCFE::DoOdeRhs(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
              Array<OneD,       Array<OneD, NekDouble> > &outarray,
        const NekDouble                                   time)
    {
        int i;
        int nvariables = inarray.num_elements();
        int npoints    = GetNpoints();


        Array<OneD, Array<OneD, NekDouble> > advVel(m_spacedim);
        
        Array<OneD, Array<OneD, NekDouble> > outarrayAdv(nvariables);
        Array<OneD, Array<OneD, NekDouble> > outarrayDiff(nvariables);
        for (i = 0; i < nvariables; ++i)
        {
            outarrayAdv[i] = Array<OneD, NekDouble>(npoints, 0.0);
            outarrayDiff[i] = Array<OneD, NekDouble>(npoints, 0.0);
        }

        // Advection term in physical rhs form
        m_advection->Advect(nvariables, m_fields, advVel, inarray,
                            outarrayAdv, time);
        
        // Extract pressure and temperature
        Array<OneD, NekDouble > pressure   (npoints, 0.0);
        Array<OneD, NekDouble > temperature(npoints, 0.0);
        GetPressure(inarray, pressure);
        GetTemperature(inarray, pressure, temperature);

        if (m_explicitDiffusion == true)
        {
            Array<OneD, Array<OneD, NekDouble> > inarrayDiff(nvariables-1);
            for (i = 0; i < nvariables-1; ++i)
            {
                inarrayDiff[i] = Array<OneD, NekDouble>(npoints, 0.0);
            }
            
            // Extract and copy velocities
            for (i = 1; i < nvariables-1; ++i)
            {
                Vmath::Vdiv(npoints, inarray[i], 1, inarray[0], 1,
                            inarrayDiff[i-1], 1);
            }

            // Copy temperature
            Vmath::Vcopy(npoints, temperature, 1, inarrayDiff[nvariables-2], 1);
            
            // Diffusion term in physical RHS form
            m_diffusion->Diffuse(nvariables, m_fields, inarrayDiff,
                                 outarrayDiff);
        }
        else if (m_explicitDiffusion == false)
        {
            Array<OneD, Array<OneD, NekDouble> > inarrayDiff(nvariables);
            for (i = 0; i < nvariables; ++i)
            {
                inarrayDiff[i] = Array<OneD, NekDouble>(npoints, 0.0);
            }
            
            // Extract and copy velocities
            for (i = 1; i < nvariables-1; ++i)
            {
                Vmath::Vdiv(npoints, inarray[i], 1, inarray[0], 1,
                            inarrayDiff[i-1], 1);
            }
            
            // Copy temperature
            Vmath::Vcopy(npoints, temperature, 1, inarrayDiff[nvariables-2], 1);
            
            // Copy density
            Vmath::Vcopy(npoints, inarray[0], 1, inarrayDiff[nvariables-1], 1);

            // Diffusion term in physical RHS form
            m_diffusion->Diffuse(nvariables, m_fields, inarrayDiff,
                                 outarrayDiff);
        }
        
        // Calculate the overall RHS
        for (i = 0; i < nvariables; ++i)
        {
            Vmath::Vsub(npoints, outarrayDiff[i], 1, outarrayAdv[i], 1,
                        outarray[i], 1);
        }

        // Add sponge layer if defined in the session file
        std::vector<SolverUtils::ForcingSharedPtr>::const_iterator x;
        for (x = m_forcing.begin(); x != m_forcing.end(); ++x)
        {
            (*x)->Apply(m_fields, inarray, outarray, time);
        }
    }

    void NavierStokesCFE::DoOdeProjection(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
              Array<OneD,       Array<OneD, NekDouble> > &outarray,
        const NekDouble                                   time)
    {
        int i;
        int nvariables = inarray.num_elements();

        switch(m_projectionType)
        {
            case MultiRegions::eDiscontinuous:
            {
                // Just copy over array
                int npoints = GetNpoints();

                for(i = 0; i < nvariables; ++i)
                {
                    Vmath::Vcopy(npoints, inarray[i], 1, outarray[i], 1);
                }
                
                SetBoundaryConditions(outarray, time);
                break;
            }
            case MultiRegions::eGalerkin:
            case MultiRegions::eMixed_CG_Discontinuous:
            {
                ASSERTL0(false, "No Continuous Galerkin for full compressible "
                                "Navier-Stokes equations");
                break;
            }
            default:
                ASSERTL0(false, "Unknown projection scheme");
                break;
        }
    }
    
    
    /* @brief Compute the diffusion term implicitly.
     *
     * @param inarray    Given fields.
     * @param outarray   Calculated solution.
     * @param time       Time.
     * @param lambda     Diffusion coefficient.
     */
    void NavierStokesCFE::DoImplicitSolve(
        const Array<OneD, const Array<OneD, NekDouble> >&inarray,
              Array<OneD,       Array<OneD, NekDouble> >&outarray,
        const NekDouble time,
        const NekDouble dt)
    {
        int nq = m_fields[0]->GetNpoints();
        int nvariables = inarray.num_elements();
        StdRegions::ConstFactorMap factors;
        
        // Forcing term for the Helmholtz problem (previous time-level)
        Array<OneD, Array< OneD, NekDouble> > F(nvariables);
        F[0] = Array<OneD, NekDouble> (nq*nvariables);
        for (int n = 1; n < nvariables; ++n)
        {
            F[n] = F[n-1] + nq;
        }
        
        // ---------------------------------------------------------------------
        // Correct spatially varying coefficients - need to comment them out
        // because HelmSolve call accepts only Stdregions::ConstFactorMap input
        // ---------------------------------------------------------------------
        /*
        Array<OneD, NekDouble> factorLambdaMomentum(nq, 0.0);
        Array<OneD, NekDouble> factorLambdaEnergy(nq, 0.0);
        Vmath::Smul(nq, -1.0 / (dt * m_mu),
                    inarray[0], 1, factorLambdaMomentum, 1);
        Vmath::Smul(nq, (m_gamma / m_Prandtl),
                    factorLambdaMomentum, 1,
                    factorLambdaEnergy, 1);
        
        // Defining forcing term as the previous time-level solution scaled
        // by nabla operator coefficient and timestep (see chapter 5 of
        // Karnidakis and Sherwin book)
        for (int i = 1; i < nvariables-1; ++i)
        {
            Vmath::Vmul(nq, factorLambdaMomentum, 1, inarray[i], 1, F[i], 1);
        }
        Vmath::Vmul(nq, factorLambdaEnergy, 1, inarray[nvariables-1], 1,
                    F[nvariables-1], 1);
        */
        // ---------------------------------------------------------------------
        
        // Setting boundary conditions
        SetBoundaryConditions(outarray, time);
        
        // Defining the scalar factors needed by the HelmSolve call
        factors[StdRegions::eFactorLambda] = m_rhoInf / (dt * m_mu);
        
        // Solve the momentum equations with Helmholtz solver
        for (int i = 1; i < nvariables-2; ++i)
        {
            m_fields[i]->HelmSolve(F[i], m_fields[i]->UpdateCoeffs(),
                                   NullFlagList, factors);
            
            Vmath::Vcopy(outarray[i].num_elements(),
                         m_fields[i]->GetCoeffs(), 1,
                         outarray[i], 1);
            
            m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(), outarray[i]);
        }
        
        // Defining the scalar factors needed by the HelmSolve call
        factors[StdRegions::eFactorLambda] = (m_Prandtl * m_rhoInf) /
                                                (m_gamma * dt * m_mu);
        
        // Solve the energy equation with Helmholtsz solver
        m_fields[nvariables-1]->HelmSolve(
            F[nvariables-1],
            m_fields[nvariables-1]->UpdateCoeffs(),
            NullFlagList, factors);
        
        Vmath::Vcopy(outarray[nvariables-1].num_elements(),
                     m_fields[nvariables-1]->GetCoeffs(), 1,
                     outarray[nvariables-1], 1);
        
        m_fields[nvariables-1]->BwdTrans(m_fields[nvariables-1]->GetCoeffs(),
                                         outarray[nvariables-1]);
    }

    
    
    /* @brief Set the boundary conditions.
     *
     * @param inarray    Given fields.
     * @param time       Time.
     */
    void NavierStokesCFE::SetBoundaryConditions(
        Array<OneD, Array<OneD, NekDouble> > &inarray,
        NekDouble                             time)
    {
        int cnt = 0;
        std::string varName;

        // Loop over Boundary Regions
        for (int n = 0; n < m_fields[0]->GetBndConditions().num_elements(); ++n)
        {
            std::string type = m_fields[0]->
            GetBndConditions()[n]->GetUserDefined();
            SetCommonBC(type, n, time, cnt, inarray);
        }
    }
}
