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
        int npoints    = GetNpoints();
        int nTracePts  = GetTraceTotPoints();

        Array<OneD, Array<OneD, NekDouble> > outarrayDiff(nvariables);

        Array<OneD, Array<OneD, NekDouble> > inarrayDiff(nvariables-1);
        Array<OneD, Array<OneD, NekDouble> > inFwd(nvariables-1);
        Array<OneD, Array<OneD, NekDouble> > inBwd(nvariables-1);

        for (i = 0; i < nvariables; ++i)
        {
            outarrayDiff[i] = Array<OneD, NekDouble>(npoints);
        }

        for (i = 0; i < nvariables-1; ++i)
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
        int nvariables = inarray.num_elements();
        StdRegions::ConstFactorMap factors;

        // Forcing term for the Helmholtz problem (previous time-level)
        Array< OneD, NekDouble> F(nq);

        // TODO: Variable factors

        // Setting boundary conditions
        SetBoundaryConditions(outarray, time);

        // Density has no implicit part -> outarray = inarray / aii_Dt
        Vmath::Smul(nq, 1.0/aii_Dt,
                        inarray[0], 1, outarray[0], 1);

        // Defining the scalar factors for momentum Helmholtz equation
        factors[StdRegions::eFactorTau]    = 1.0;
        factors[StdRegions::eFactorLambda] = m_rhoInf / (aii_Dt * m_mu);

        // Solve the momentum equations with Helmholtz solver
        for (int i = 1; i < nvariables-1; ++i)
        {
            Vmath::Smul(nq, -factors[StdRegions::eFactorLambda],
                        inarray[i], 1, F, 1);
            m_fields[i]->HelmSolve(F, m_fields[i]->UpdateCoeffs(),
                                   NullFlagList, factors);
            m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(), outarray[i]);
        }

        // Defining the scalar factors for energy Helmholtz equation
        factors[StdRegions::eFactorLambda] = (m_Prandtl * m_rhoInf) /
                                                (m_gamma * aii_Dt * m_mu);

        // Solve the energy equation with Helmholtz solver
        Vmath::Smul(nq, -factors[StdRegions::eFactorLambda],
                        inarray[nvariables-1], 1, F, 1);
        m_fields[nvariables-1]->HelmSolve(
            F[nvariables-1],
            m_fields[nvariables-1]->UpdateCoeffs(),
            NullFlagList, factors);

        m_fields[nvariables-1]->BwdTrans(m_fields[nvariables-1]->GetCoeffs(),
                                         outarray[nvariables-1]);
    }
}
