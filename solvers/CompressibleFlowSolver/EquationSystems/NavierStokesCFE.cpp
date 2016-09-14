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
                        fac = 1.0; //4.0/3.0;
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
}
