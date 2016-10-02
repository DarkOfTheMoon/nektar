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
#include <MultiRegions/GlobalLinSys.h>

using namespace std;

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

    NavierStokesCFE::~NavierStokesCFE()
    {

    }

    /**
     * @brief Initialization object for CompressibleFlowSystem class.
     */
    void NavierStokesCFE::v_InitObject()
    {
        CompressibleFlowSystem::v_InitObject();

        string diffName;
        m_session->LoadSolverInfo("DiffusionType", diffName, "LDGNS");

        m_diffusion = SolverUtils::GetDiffusionFactory()
                                    .CreateInstance(diffName, diffName);

        if (m_specHP_dealiasing)
        {
            if (m_explicitDiffusion)
            {
                m_diffusion->SetFluxVectorNS(
                    &NavierStokesCFE::GetViscousFluxVectorDeAlias,
                    this);
            }
            else
            {
                ASSERTL0(false,
                    "Imex viscous flux not implemented with dealiasing");
            }
        }
        else
        {
            if (m_explicitDiffusion)
            {
                m_diffusion->SetFluxVectorNS(&NavierStokesCFE::
                        GetViscousFluxVector, this);
            }
            else
            {
                m_diffusion->SetFluxVectorNS(&NavierStokesCFE::
                        GetViscousFluxVectorSemiImplicit, this);
            }
        }

        // Concluding initialisation of diffusion operator
        m_diffusion->InitObject         (m_session, m_fields);
    }

    void NavierStokesCFE::v_DoDiffusion(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
              Array<OneD,       Array<OneD, NekDouble> > &outarray,
            const Array<OneD, Array<OneD, NekDouble> >   &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >   &pBwd)
    {
        int i;
        int nvariables = inarray.num_elements();
        int nInputs    = m_explicitDiffusion ? (nvariables-1) : nvariables;
        int npoints    = GetNpoints();
        int nTracePts  = GetTraceTotPoints();

        Array<OneD, Array<OneD, NekDouble> > outarrayDiff(nvariables);

        Array<OneD, Array<OneD, NekDouble> > inarrayDiff(nInputs);
        Array<OneD, Array<OneD, NekDouble> > inFwd(nInputs);
        Array<OneD, Array<OneD, NekDouble> > inBwd(nInputs);

        for (i = 0; i < nvariables; ++i)
        {
            outarrayDiff[i] = Array<OneD, NekDouble>(npoints);
        }

        for (i = 0; i < nInputs; ++i)
        {
            inarrayDiff[i] = Array<OneD, NekDouble>(npoints);
            inFwd[i]       = Array<OneD, NekDouble>(nTracePts);
            inBwd[i]       = Array<OneD, NekDouble>(nTracePts);
        }

        // Extract pressure
        //    (use inarrayDiff[0] as a temporary storage for the pressure)
        m_varConv->GetPressure(inarray, inarrayDiff[0]);

        // Extract temperature
        m_varConv->GetTemperature(inarray, inarrayDiff[0],
                inarrayDiff[nvariables-2]);

        // Extract velocities
        m_varConv->GetVelocityVector(inarray, inarrayDiff);

        // Density for semi-implicit
        if( !m_explicitDiffusion )
        {
            Vmath::Vcopy(npoints, inarray[0], 1, inarrayDiff[nvariables-1], 1);
        }

        // Repeat calculation for trace space
        if (pFwd == NullNekDoubleArrayofArray || 
            pBwd == NullNekDoubleArrayofArray)
        {
            inFwd = NullNekDoubleArrayofArray;
            inBwd = NullNekDoubleArrayofArray;
        }
        else
        {
            m_varConv->GetPressure(pFwd,    inFwd[0]);
            m_varConv->GetPressure(pBwd,    inBwd[0]);

            m_varConv->GetTemperature(pFwd,    inFwd[0],
                inFwd[nvariables-2]);
            m_varConv->GetTemperature(pBwd,    inBwd[0],
                inBwd[nvariables-2]);

            m_varConv->GetVelocityVector(pFwd, inFwd);
            m_varConv->GetVelocityVector(pBwd, inBwd);

            if( !m_explicitDiffusion )
            {
                Vmath::Vcopy(nTracePts, pFwd[0], 1, inFwd[nvariables-1], 1);
                Vmath::Vcopy(nTracePts, pBwd[0], 1, inBwd[nvariables-1], 1);
            }
        }

        // Diffusion term in physical rhs form
        m_diffusion->Diffuse(nvariables, m_fields, inarrayDiff, outarrayDiff,
                             inFwd, inBwd);

        for (i = 0; i < nvariables; ++i)
        {
            Vmath::Vadd(npoints,
                        outarrayDiff[i], 1,
                        outarray[i], 1,
                        outarray[i], 1);
        }
    }

    /**
     * @brief Perform implicit solution when using IMEX
     */
    void NavierStokesCFE::v_DoImplicitSolve(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
              Array<OneD,       Array<OneD, NekDouble> > &outarray,
        const NekDouble                                   time,
        const NekDouble                                   aii_Dt)
    {
        int nq = m_fields[0]->GetNpoints();
        StdRegions::ConstFactorMap factors;
        // For variable factors case
        StdRegions::VarCoeffMap    varcoeff;
        Array< OneD, NekDouble>    kinvis;

        // Forcing term for the Helmholtz problem (previous time-level)
        Array< OneD, NekDouble> F(nq);

        // Setting boundary conditions
        SetBoundaryConditions(outarray, time);

        // Density has no implicit part -> outarray = inarray
        Vmath::Vcopy(nq, inarray[0], 1, outarray[0], 1);

        // Variable factors type
        StdRegions::VarCoeffType varcoefftypes[3]
            = {StdRegions::eVarCoeffD00,
               StdRegions::eVarCoeffD11,
               StdRegions::eVarCoeffD22};

        if (m_variableCoeffs)
        {
            // Initialise variable factors
            for (int i = 0; i < m_spacedim; ++i)
            {
                varcoeff[varcoefftypes[i]] =
                        Array<OneD, NekDouble>(nq, 0.0);
            }
            // Calculate kinematic viscosity
            kinvis = Array<OneD, NekDouble>     (nq);
            Array<OneD, NekDouble>          tmp1(nq);
            Array<OneD, NekDouble>          tmp2(nq);
            if (m_ViscosityType == "Variable")
            {
                // Calculate pressure
                m_varConv->GetPressure(inarray, tmp1);
                // Extract temperature
                m_varConv->GetTemperature(inarray, tmp1, tmp2);
                // Get viscosity
                m_varConv->GetDynamicViscosity(tmp2, tmp1);
            }
            else
            {
                Vmath::Fill(nq, m_mu, tmp1, 1);
            }
            Vmath::Vdiv( nq, tmp1, 1, m_fields[0]->GetPhys(), 1, kinvis, 1);
        }
        else
        {
            varcoeff = StdRegions::NullVarCoeffMap;
        }

        // Defining the scalar factors for momentum Helmholtz equation
        factors[StdRegions::eFactorTau]    = 1.0;
        factors[StdRegions::eFactorLambda] = m_rhoInf / (aii_Dt * m_mu);

        // Solve the momentum equations with Helmholtz solver
        for (int i = 1; i < m_spacedim+1; ++i)
        {
            if (m_variableCoeffs)
            {
                NekDouble fac;
                // Set diagonal terms
                for (int j = 0; j < m_spacedim; ++j)
                {
                    if (j == (i-1))
                    {
                        fac = 4.0/3.0;
                    }
                    else
                    {
                        fac = 1.0;
                    }
                    Vmath::Smul( nq, fac*m_rhoInf/m_mu, kinvis, 1,
                                    varcoeff[varcoefftypes[j]], 1);
                }
            }

            Vmath::Smul(nq, -factors[StdRegions::eFactorLambda],
                        inarray[i], 1, F, 1);

            m_fields[i]->HelmSolve(F, m_fields[i]->UpdateCoeffs(),
                                   NullFlagList, factors, varcoeff);
            m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(), outarray[i]);

            if (m_variableCoeffs)
            {
                if (LibUtilities::NekManager<MultiRegions::GlobalLinSysKey,
                                    MultiRegions::GlobalLinSys>::
                                    PoolCreated(std::string("GlobalLinSys")))
                {
                    LibUtilities::NekManager<MultiRegions::GlobalLinSysKey,
                                        MultiRegions::GlobalLinSys>::
                                        ClearManager(std::string("GlobalLinSys"));
                }
            }
        }

        // Defining the scalar factors for energy Helmholtz equation
        factors[StdRegions::eFactorLambda] = (m_Prandtl * m_rhoInf) /
                                                (m_gamma * aii_Dt * m_mu);

        // Solve the energy equation with Helmholtz solver
        if (m_variableCoeffs)
        {
            // Set diagonal terms
            for (int j = 0; j < m_spacedim; ++j)
            {
                Vmath::Smul( nq, 1.0*m_rhoInf/m_mu    , kinvis, 1,
                            varcoeff[varcoefftypes[j]], 1);
            }
        }

        Vmath::Smul(nq, -factors[StdRegions::eFactorLambda],
                        inarray[m_spacedim+1], 1, F, 1);
        m_fields[m_spacedim+1]->HelmSolve(
            F,
            m_fields[m_spacedim+1]->UpdateCoeffs(),
            NullFlagList, factors, varcoeff);

        m_fields[m_spacedim+1]->BwdTrans(m_fields[m_spacedim+1]->GetCoeffs(),
                                         outarray[m_spacedim+1]);

        // Clear managers to avoid using too much memory
        if (m_variableCoeffs)
        {
            if (LibUtilities::NekManager<MultiRegions::GlobalLinSysKey,
                                    MultiRegions::GlobalLinSys>::
                                    PoolCreated(std::string("GlobalLinSys")))
            {
                LibUtilities::NekManager<MultiRegions::GlobalLinSysKey,
                                    MultiRegions::GlobalLinSys>::
                                    ClearManager(std::string("GlobalLinSys"));
            }
            LibUtilities::NekManager<LocalRegions::MatrixKey,
                DNekScalMat, LocalRegions::MatrixKey::opLess>::ClearManager();
            LibUtilities::NekManager<LocalRegions::MatrixKey,
                DNekScalBlkMat, LocalRegions::MatrixKey::opLess>::ClearManager();
        }
    }

    /* @brief Return the flux vector for the LDG diffusion problem.
     * \todo Complete the viscous flux vector
     */
    void NavierStokesCFE::GetViscousFluxVector(
        const Array<OneD, Array<OneD, NekDouble> >               &physfield,
              Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &derivativesO1,
              Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &viscousTensor)
    {
        int i, j;
        int nVariables = m_fields.num_elements();
        int nPts       = physfield[0].num_elements();

        // Stokes hypothesis
        const NekDouble lambda = -2.0/3.0;

        // Auxiliary variables
        Array<OneD, NekDouble > mu                 (nPts, 0.0);
        Array<OneD, NekDouble > thermalConductivity(nPts, 0.0);
        Array<OneD, NekDouble > divVel             (nPts, 0.0);

        // Variable viscosity through the Sutherland's law
        if (m_ViscosityType == "Variable")
        {
            m_varConv->GetDynamicViscosity(physfield[nVariables-2], mu);
            NekDouble tRa = m_Cp / m_Prandtl;
            Vmath::Smul(nPts, tRa, mu, 1, thermalConductivity, 1);
        }
        else
        {
            Vmath::Fill(nPts, m_mu, mu, 1);
            Vmath::Fill(nPts, m_thermalConductivity,
                        thermalConductivity, 1);
        }

        // Velocity divergence
        for (j = 0; j < m_spacedim; ++j)
        {
            Vmath::Vadd(nPts, divVel, 1, derivativesO1[j][j], 1,
                        divVel, 1);
        }

        // Velocity divergence scaled by lambda * mu
        Vmath::Smul(nPts, lambda, divVel, 1, divVel, 1);
        Vmath::Vmul(nPts, mu,  1, divVel, 1, divVel, 1);

        // Viscous flux vector for the rho equation = 0
        for (i = 0; i < m_spacedim; ++i)
        {
            Vmath::Zero(nPts, viscousTensor[i][0], 1);
        }

        // Viscous stress tensor (for the momentum equations)
        for (i = 0; i < m_spacedim; ++i)
        {
            for (j = i; j < m_spacedim; ++j)
            {
                Vmath::Vadd(nPts, derivativesO1[i][j], 1,
                                  derivativesO1[j][i], 1,
                                  viscousTensor[i][j+1], 1);

                Vmath::Vmul(nPts, mu, 1,
                                  viscousTensor[i][j+1], 1,
                                  viscousTensor[i][j+1], 1);

                if (i == j)
                {
                    // Add divergence term to diagonal
                    Vmath::Vadd(nPts, viscousTensor[i][j+1], 1,
                                  divVel, 1,
                                  viscousTensor[i][j+1], 1);
                }
                else
                {
                    // Copy to make symmetric
                    Vmath::Vcopy(nPts, viscousTensor[i][j+1], 1,
                                       viscousTensor[j][i+1], 1);
                }
            }
        }

        // Terms for the energy equation
        for (i = 0; i < m_spacedim; ++i)
        {
            Vmath::Zero(nPts, viscousTensor[i][m_spacedim+1], 1);
            // u_j * tau_ij
            for (j = 0; j < m_spacedim; ++j)
            {
                Vmath::Vvtvp(nPts, physfield[j], 1,
                               viscousTensor[i][j+1], 1,
                               viscousTensor[i][m_spacedim+1], 1,
                               viscousTensor[i][m_spacedim+1], 1);
            }
            // Add k*T_i
            Vmath::Vvtvp(nPts, thermalConductivity, 1,
                               derivativesO1[i][m_spacedim], 1,
                               viscousTensor[i][m_spacedim+1], 1,
                               viscousTensor[i][m_spacedim+1], 1);
        }
    }

    /**
     * @brief Return the flux vector for the LDG diffusion problem.
     * \todo Complete the viscous flux vector
     */
    void NavierStokesCFE::GetViscousFluxVectorDeAlias(
        const Array<OneD, Array<OneD, NekDouble> >               &physfield,
              Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &derivativesO1,
              Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &viscousTensor)
    {
        int i, j;
        int nVariables = m_fields.num_elements();
        // Factor to rescale 1d points in dealiasing.
        NekDouble OneDptscale = 2;
        // Get number of points to dealias a cubic non-linearity
        int nPts      = m_fields[0]->Get1DScaledTotPoints(OneDptscale);
        int nPts_orig = physfield[0].num_elements();

        // Stokes hypothesis
        const NekDouble lambda = -2.0/3.0;

        // Auxiliary variables
        Array<OneD, NekDouble > mu                 (nPts, 0.0);
        Array<OneD, NekDouble > thermalConductivity(nPts, 0.0);
        Array<OneD, NekDouble > divVel             (nPts, 0.0);

        // Variable viscosity through the Sutherland's law
        if (m_ViscosityType == "Variable")
        {
            m_varConv->GetDynamicViscosity(physfield[nVariables-2], mu);
            NekDouble tRa = m_Cp / m_Prandtl;
            Vmath::Smul(nPts, tRa, mu, 1, thermalConductivity, 1);
        }
        else
        {
            Vmath::Fill(nPts, m_mu, mu, 1);
            Vmath::Fill(nPts, m_thermalConductivity,
                        thermalConductivity, 1);
        }

        // Interpolate inputs and initialise interpolated output
        Array<OneD, Array<OneD, NekDouble> > vel_interp(m_spacedim);
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > >
                                             deriv_interp(m_spacedim);
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > >
                                             out_interp(m_spacedim);
        for (i = 0; i < m_spacedim; ++i)
        {
            // Interpolate velocity
            vel_interp[i]   = Array<OneD, NekDouble> (nPts);
            m_fields[0]->PhysInterp1DScaled(
                OneDptscale, physfield[i], vel_interp[i]);

            // Interpolate derivatives
            deriv_interp[i] = Array<OneD,Array<OneD,NekDouble> > (m_spacedim+1);
            for (j = 0; j < m_spacedim+1; ++j)
            {
                deriv_interp[i][j] = Array<OneD, NekDouble> (nPts);
                m_fields[0]->PhysInterp1DScaled(
                    OneDptscale, derivativesO1[i][j], deriv_interp[i][j]);
            }

            // Output (start from j=1 since flux is zero for rho)
            out_interp[i] = Array<OneD,Array<OneD,NekDouble> > (m_spacedim+2);
            for (j = 1; j < m_spacedim+2; ++j)
            {
                out_interp[i][j] = Array<OneD, NekDouble> (nPts);
            }
        }

        // Velocity divergence
        for (j = 0; j < m_spacedim; ++j)
        {
            Vmath::Vadd(nPts, divVel, 1, deriv_interp[j][j], 1,
                        divVel, 1);
        }

        // Velocity divergence scaled by lambda * mu
        Vmath::Smul(nPts, lambda, divVel, 1, divVel, 1);
        Vmath::Vmul(nPts, mu,  1, divVel, 1, divVel, 1);

        // Viscous flux vector for the rho equation = 0 (no need to dealias)
        for (i = 0; i < m_spacedim; ++i)
        {
            Vmath::Zero(nPts_orig, viscousTensor[i][0], 1);
        }

        // Viscous stress tensor (for the momentum equations)
        for (i = 0; i < m_spacedim; ++i)
        {
            for (j = i; j < m_spacedim; ++j)
            {
                Vmath::Vadd(nPts, deriv_interp[i][j], 1,
                                  deriv_interp[j][i], 1,
                                  out_interp[i][j+1], 1);

                Vmath::Vmul(nPts, mu, 1,
                                  out_interp[i][j+1], 1,
                                  out_interp[i][j+1], 1);

                if (i == j)
                {
                    // Add divergence term to diagonal
                    Vmath::Vadd(nPts, out_interp[i][j+1], 1,
                                  divVel, 1,
                                  out_interp[i][j+1], 1);
                }
                else
                {
                    // Make symmetric
                    out_interp[j][i+1] = out_interp[i][j+1];
                }
            }
        }

        // Terms for the energy equation
        for (i = 0; i < m_spacedim; ++i)
        {
            Vmath::Zero(nPts, out_interp[i][m_spacedim+1], 1);
            // u_j * tau_ij
            for (j = 0; j < m_spacedim; ++j)
            {
                Vmath::Vvtvp(nPts, vel_interp[j], 1,
                               out_interp[i][j+1], 1,
                               out_interp[i][m_spacedim+1], 1,
                               out_interp[i][m_spacedim+1], 1);
            }
            // Add k*T_i
            Vmath::Vvtvp(nPts, thermalConductivity, 1,
                               deriv_interp[i][m_spacedim], 1,
                               out_interp[i][m_spacedim+1], 1,
                               out_interp[i][m_spacedim+1], 1);
        }

        // Project to original space
        for (i = 0; i < m_spacedim; ++i)
        {
            for (j = 1; j < m_spacedim+2; ++j)
            {
                m_fields[0]->PhysGalerkinProjection1DScaled(
                    OneDptscale,
                    out_interp[i][j],
                    viscousTensor[i][j]);
            }
        }
    }

    /**
     * @brief Return the flux vector for the LDG diffusion problem using
     *        semi-implicit time integration.
     */
    void NavierStokesCFE::GetViscousFluxVectorSemiImplicit(
        const Array<OneD, Array<OneD, NekDouble> >               &physfield,
              Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &derivativesO1,
              Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &viscousTensor)
    {
        int i, j;
        int nVariables = m_fields.num_elements();
        int nPts       = physfield[0].num_elements();

        // Stokes hypothesis
        const NekDouble lambda = -2.0/3.0;

        // Auxiliary variables
        Array<OneD, NekDouble > mu                 (nPts, 0.0);
        Array<OneD, NekDouble > rhoRef             (nPts, 0.0);
        Array<OneD, NekDouble > thermalConductivity(nPts, 0.0);
        Array<OneD, NekDouble > divVel             (nPts, 0.0);
        Array<OneD, NekDouble > tmp                (nPts, 0.0);
        Array<OneD, NekDouble > tmp2               (nPts, 0.0);

        // Variable viscosity through the Sutherland's law
        if (m_ViscosityType == "Variable")
        {
            m_varConv->GetDynamicViscosity(physfield[nVariables-2], mu);
            NekDouble tRa = m_Cp / m_Prandtl;
            Vmath::Smul(nPts, tRa, mu, 1, thermalConductivity, 1);
        }
        else
        {
            Vmath::Fill(nPts, m_mu, mu, 1);
            Vmath::Fill(nPts, m_thermalConductivity,
                        thermalConductivity, 1);
        }

        // Density used in implicit part:
        //    - rhoInf when using constant coefficients
        //    - rho    when using variable coefficients
        if ( m_variableCoeffs)
        {
            Vmath::Vcopy(nPts, physfield[nVariables-1], 1, rhoRef, 1);
        }
        else
        {
            Vmath::Fill(nPts, m_rhoInf, rhoRef, 1);
        }

        // Velocity divergence
        for (j = 0; j < m_spacedim; ++j)
        {
            Vmath::Vadd(nPts, divVel, 1, derivativesO1[j][j], 1,
                        divVel, 1);
        }

        // Velocity divergence scaled by lambda * mu
        Vmath::Smul(nPts, lambda, divVel, 1, divVel, 1);
        Vmath::Vmul(nPts, mu,  1, divVel, 1, divVel, 1);

        // Viscous flux vector for the rho equation = 0
        for (i = 0; i < m_spacedim; ++i)
        {
            Vmath::Zero(nPts, viscousTensor[i][0], 1);
        }

        // Viscous stress tensor (for the momentum equations)
        for (i = 0; i < m_spacedim; ++i)
        {
            for (j = i; j < m_spacedim; ++j)
            {
                Vmath::Vadd(nPts, derivativesO1[i][j], 1,
                                  derivativesO1[j][i], 1,
                                  viscousTensor[i][j+1], 1);

                Vmath::Vmul(nPts, mu, 1,
                                  viscousTensor[i][j+1], 1,
                                  viscousTensor[i][j+1], 1);

                if (i == j)
                {
                    // Add divergence term to diagonal
                    Vmath::Vadd(nPts, viscousTensor[i][j+1], 1,
                                  divVel, 1,
                                  viscousTensor[i][j+1], 1);
                }
                else
                {
                    // Copy to make symmetric
                    Vmath::Vcopy(nPts, viscousTensor[i][j+1], 1,
                                       viscousTensor[j][i+1], 1);
                }
            }
        }

        // Terms for the energy equation

        // Calculate -k*T/rhoRef
        Vmath::Vdiv(nPts, thermalConductivity, 1,
                            rhoRef, 1,
                            tmp, 1);
        Vmath::Vmul(nPts, tmp, 1,
                            physfield[nVariables-2], 1,
                            tmp, 1);
        Vmath::Neg(nPts, tmp, 1);

        // Calculate flux for energy equation
        for (i = 0; i < m_spacedim; ++i)
        {
            // -(k*T/rhoRef)*rho_i
            Vmath::Vmul(nPts, tmp, 1,
                            derivativesO1[i][m_spacedim+1], 1,
                            viscousTensor[i][m_spacedim+1], 1);

            if ( !m_variableCoeffs)
            {
                // + k * T_i
                Vmath::Vvtvp(nPts, thermalConductivity, 1,
                                derivativesO1[i][m_spacedim], 1,
                                viscousTensor[i][m_spacedim+1], 1,
                                viscousTensor[i][m_spacedim+1], 1);

                // - k * rho/rhoRef * T_i
                Vmath::Vmul(nPts, thermalConductivity, 1,
                                physfield[nVariables-1], 1,
                                tmp2, 1);
                Vmath::Vdiv(nPts, tmp2, 1,
                                rhoRef, 1,
                                tmp2, 1);
                Vmath::Vmul(nPts, derivativesO1[i][m_spacedim], 1,
                                tmp2, 1,
                                tmp2, 1);
                Vmath::Vsub(nPts, viscousTensor[i][m_spacedim+1], 1,
                                tmp2, 1,
                                viscousTensor[i][m_spacedim+1], 1);
            }

            for (j = 0; j < m_spacedim; ++j)
            {
                // - mu*gamma/Pr * rho/rhoRef * u_j * u_j,i
                Vmath::Smul(nPts, -m_gamma/m_Prandtl,
                            mu, 1,
                            tmp2, 1);
                if ( !m_variableCoeffs)
                {
                    Vmath::Vmul(nPts, physfield[nVariables-1], 1,
                                tmp2, 1,
                                tmp2, 1);
                    Vmath::Vdiv(nPts, tmp2, 1,
                                rhoRef, 1,
                                tmp2, 1);
                }
                Vmath::Vmul(nPts, physfield[j], 1,
                            tmp2, 1,
                            tmp2, 1);
                Vmath::Vvtvp(nPts, tmp2, 1,
                               derivativesO1[i][j], 1,
                               viscousTensor[i][m_spacedim+1], 1,
                               viscousTensor[i][m_spacedim+1], 1);

                // - mu*gamma/(2*rhoRef*Pr) * u_j*u_j * rho_i
                Vmath::Smul(nPts, -m_gamma/(2*m_Prandtl),
                            mu, 1,
                            tmp2, 1);
                Vmath::Vdiv(nPts, tmp2, 1,
                            rhoRef, 1,
                            tmp2, 1);
                Vmath::Vmul(nPts, physfield[j], 1,
                            tmp2, 1,
                            tmp2, 1);
                Vmath::Vmul(nPts, physfield[j], 1,
                            tmp2, 1,
                            tmp2, 1);
                Vmath::Vvtvp(nPts, tmp2, 1,
                               derivativesO1[i][m_spacedim+1], 1,
                               viscousTensor[i][m_spacedim+1], 1,
                               viscousTensor[i][m_spacedim+1], 1);

                // + u_j * tau_ij
                Vmath::Vvtvp(nPts, physfield[j], 1,
                               viscousTensor[i][j+1], 1,
                               viscousTensor[i][m_spacedim+1], 1,
                               viscousTensor[i][m_spacedim+1], 1);
            }
        }
        // Correct fluxes for momentum equation using Imex
        NekDouble fac;
        for (i = 0; i < m_spacedim; ++i)
        {
            for (j = 0; j < m_spacedim; ++j)
            {
                if ( (i == j) && m_variableCoeffs)
                {
                    fac = 4.0/3.0;
                }
                else
                {
                    fac = 1.0;
                }

                // - fac*mu* rho/rhoRef * u_j,i
                Vmath::Vmul(nPts, mu, 1,
                                  derivativesO1[i][j], 1,
                                  tmp2, 1);
                Vmath::Smul(nPts, fac,
                                  tmp2, 1,
                                  tmp2, 1);
                if ( !m_variableCoeffs)
                {
                    Vmath::Vmul(nPts, physfield[nVariables-1], 1,
                                tmp2, 1,
                                tmp2, 1);
                    Vmath::Vdiv(nPts, tmp2, 1,
                                rhoRef, 1,
                                tmp2, 1);
                }
                Vmath::Vsub(nPts, viscousTensor[i][j+1], 1,
                                  tmp2, 1,
                                  viscousTensor[i][j+1], 1);

                // - fac*mu*u_j/rhoRef * rho_i
                Vmath::Vmul(nPts, mu, 1,
                                  physfield[j], 1,
                                  tmp2, 1);
                Vmath::Vdiv(nPts, tmp2, 1,
                                  rhoRef, 1,
                                  tmp2, 1);
                Vmath::Vmul(nPts, tmp2, 1,
                                  derivativesO1[i][m_spacedim+1], 1,
                                  tmp2, 1);
                Vmath::Smul(nPts, fac,
                                  tmp2, 1,
                                  tmp2, 1);
                Vmath::Vsub(nPts, viscousTensor[i][j+1], 1,
                                  tmp2, 1,
                                  viscousTensor[i][j+1], 1);
            }
        }
    }

}
