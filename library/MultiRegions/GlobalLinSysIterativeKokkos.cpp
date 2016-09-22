///////////////////////////////////////////////////////////////////////////////
//
// File: GlobalLinSysIterativeKokkos.cpp
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
// Description: GlobalLinSysIterativeKokkos definition
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/GlobalLinSysIterativeKokkos.h>

using namespace std;

namespace Nektar
{
namespace MultiRegions
{
/**
 * @class GlobalLinSysIterativeKokkos
 *
 * Solves a linear system using iterative methods.
 */

/// Constructor for full direct matrix solve.
GlobalLinSysIterativeKokkos::GlobalLinSysIterativeKokkos(
        const GlobalLinSysKey &pKey,
        const boost::weak_ptr<ExpList> &pExpList,
        const boost::shared_ptr<AssemblyMap>
        &pLocToGloMap)
        : GlobalLinSys(pKey, pExpList, pLocToGloMap),
          m_rhs_magnitude(NekConstants::kNekUnsetDouble),
          m_rhs_mag_sm(0.9),
          m_precon(NullPreconditionerSharedPtr),
          m_totalIterations(0),
          m_numPrevSols(0)
{
    LibUtilities::SessionReaderSharedPtr vSession
                                    = pExpList.lock()->GetSession();

    m_tolerance = pLocToGloMap->GetIterativeTolerance();
    m_maxiter   = pLocToGloMap->GetMaxIterations();

    LibUtilities::CommSharedPtr vComm = m_expList.lock()->GetComm()->GetRowComm();
    m_root    = (vComm->GetRank())? false : true;
    m_verbose = (vSession->DefinesCmdLineArgument("verbose"))? true :false;
}

GlobalLinSysIterativeKokkos::~GlobalLinSysIterativeKokkos()
{
}


/**
 *
 */
void GlobalLinSysIterativeKokkos::v_SolveLinearSystem(
            const int nGlobal,
            const Array<OneD,const NekDouble> &pInput,
                  Array<OneD,      NekDouble> &pOutput,
            const AssemblyMapSharedPtr &plocToGloMap,
            const int nDir)
{
    // applying plain Conjugate Gradient
    DoConjugateGradient(nGlobal, pInput, pOutput, plocToGloMap, nDir);
}


/**  
 * Solve a global linear system using the conjugate gradient method.  
 * We solve only for the non-Dirichlet modes. The operator is evaluated  
 * using an auxiliary function v_DoMatrixMultiply defined by the  
 * specific solver. Distributed math routines are used to support  
 * parallel execution of the solver.  
 *  
 * The implemented algorithm uses a reduced-communication reordering of  
 * the standard PCG method (Demmel, Heath and Vorst, 1993)  
 *  
 * @param       pInput      Input residual  of all DOFs.  
 * @param       pOutput     Solution vector of all DOFs.  
 */
void GlobalLinSysIterativeKokkos::DoConjugateGradient(
    const int                          nGlobal,
    const Array<OneD,const NekDouble> &pInput,
          Array<OneD,      NekDouble> &pOutput,
    const AssemblyMapSharedPtr        &plocToGloMap,
    const int                          nDir)
{
    if (!m_precon)
    {
        v_UniqueMap();
        m_precon = CreatePrecon(plocToGloMap);
        m_precon->BuildPreconditioner();
    }

    // Get the communicator for performing data exchanges
    LibUtilities::CommSharedPtr vComm
        = m_expList.lock()->GetComm()->GetRowComm();

    // Get vector sizes
    int nNonDir = nGlobal - nDir;

    // Allocate array storage
    Array<OneD, NekDouble> w_A    (nGlobal, 0.0);
    Array<OneD, NekDouble> s_A    (nGlobal, 0.0);
    Array<OneD, NekDouble> p_A    (nNonDir, 0.0);
    Array<OneD, NekDouble> r_A    (nNonDir, 0.0);
    Array<OneD, NekDouble> q_A    (nNonDir, 0.0);
    Array<OneD, NekDouble> tmp;

    // Create NekVector wrappers for linear algebra operations
    NekVector<NekDouble> in (nNonDir,pInput  + nDir,      eWrapper);
    NekVector<NekDouble> out(nNonDir,tmp = pOutput + nDir,eWrapper);
    NekVector<NekDouble> w  (nNonDir,tmp = w_A + nDir,    eWrapper);
    NekVector<NekDouble> s  (nNonDir,tmp = s_A + nDir,    eWrapper);
    NekVector<NekDouble> p  (nNonDir,p_A,                 eWrapper);
    NekVector<NekDouble> r  (nNonDir,r_A,                 eWrapper);
    NekVector<NekDouble> q  (nNonDir,q_A,                 eWrapper);

    int k;
    NekDouble alpha, beta, rho, rho_new, mu, eps,  min_resid;
    Array<OneD, NekDouble> vExchange(3,0.0);

    // Copy initial residual from input
    r = in;
    // zero homogeneous out array ready for solution updates
    // Should not be earlier in case input vector is same as
    // output and above copy has been peformed
    Vmath::Zero(nNonDir,tmp = pOutput + nDir,1);


    // evaluate initial residual error for exit check
    vExchange[2] = Vmath::Dot2(nNonDir,
                               r_A,
                               r_A,
                               m_map + nDir);

    vComm->AllReduce(vExchange, Nektar::LibUtilities::ReduceSum);

    eps       = vExchange[2];

    if(m_rhs_magnitude == NekConstants::kNekUnsetDouble)
    {
        NekVector<NekDouble> inGlob (nGlobal, pInput, eWrapper);
        Set_Rhs_Magnitude(inGlob);
    }

    // If input residual is less than tolerance skip solve.
    if (eps < m_tolerance * m_tolerance * m_rhs_magnitude)
    {
        if (m_verbose && m_root)
        {
            cout << "CG iterations made = " << m_totalIterations
                 << " using tolerance of "  << m_tolerance
                 << " (error = " << sqrt(eps/m_rhs_magnitude)
                 << ", rhs_mag = " << sqrt(m_rhs_magnitude) <<  ")"
                 << endl;
        }
        return;
    }

    m_totalIterations = 1;
    m_precon->DoPreconditioner(r_A, tmp = w_A + nDir);

    v_DoMatrixMultiply(w_A, s_A);

    k = 0;

    vExchange[0] = Vmath::Dot2(nNonDir,
                               r_A,
                               w_A + nDir,
                               m_map + nDir);

    vExchange[1] = Vmath::Dot2(nNonDir,
                               s_A + nDir,
                               w_A + nDir,
                               m_map + nDir);

    vComm->AllReduce(vExchange, Nektar::LibUtilities::ReduceSum);

    rho       = vExchange[0];
    mu        = vExchange[1];
    min_resid = m_rhs_magnitude;
    beta      = 0.0;
    alpha     = rho/mu;

    // Continue until convergence
    while (true)
    {
        if(k >= m_maxiter)
        {
            if (m_root)
            {
                cout << "CG iterations made = " << m_totalIterations
                     << " using tolerance of "  << m_tolerance
                     << " (error = " << sqrt(eps/m_rhs_magnitude)
                     << ", rhs_mag = " << sqrt(m_rhs_magnitude) <<  ")"
                     << endl;
            }
            ASSERTL0(false,
                     "Exceeded maximum number of iterations");
        }

//        Kokkos::parallel_for(nNonDir, [&] (const size_t i) {
//            // Compute new search direction p_k, q_k
//            p_A[i]            = beta * p_A[i] + w_A[nDir + i];
//            // Update solution x_{k+1}
//            q_A[i]            = beta * q_A[i] + s_A[nDir + i];
//            // Update solution x_{k+1}
//            pOutput[nDir + i] = alpha * p_A[i] + pOutput[nDir + i];
//            // Update residual vector r_{k+1}
//            r_A[i]            = -alpha * q_A[i] + r[i];
//        });

        // Compute new search direction p_k, q_k
        Vmath::Svtvp(nNonDir, beta, &p_A[0], 1, &w_A[nDir], 1, &p_A[0], 1);
        Vmath::Svtvp(nNonDir, beta, &q_A[0], 1, &s_A[nDir], 1, &q_A[0], 1);

        // Update solution x_{k+1}
        Vmath::Svtvp(nNonDir, alpha, &p_A[0], 1, &pOutput[nDir], 1, &pOutput[nDir], 1);

        // Update residual vector r_{k+1}
        Vmath::Svtvp(nNonDir, -alpha, &q_A[0], 1, &r_A[0], 1, &r_A[0], 1);

        // Apply preconditioner
        m_precon->DoPreconditioner(r_A, tmp = w_A + nDir);

        // Perform the method-specific matrix-vector multiply operation.
        v_DoMatrixMultiply(w_A, s_A);

//        // <r_{k+1}, w_{k+1}>
//        Kokkos::parallel_reduce(nNonDir, [&] (const size_t i, NekDouble & ret)
//        {
//            ret += m_map[nDir + i] == 1 ? r_A[i] * w_A[nDir + i] : 0;
//        }, vExchange[0]);
//        Kokkos::parallel_reduce(nNonDir, [&] (const size_t i, NekDouble & ret)
//        {
//            ret += m_map[nDir + i] == 1 ? s_A[nDir + i] * w_A[nDir + i] : 0;
//        }, vExchange[1]);
//        Kokkos::parallel_reduce(nNonDir, [&] (const size_t i, NekDouble & ret)
//        {
//            ret += m_map[nDir + i] == 1 ? r_A[i] * r_A[i] : 0;
//        }, vExchange[2]);

        // <r_{k+1}, w_{k+1}>
        vExchange[0] = Vmath::Dot2(nNonDir,
                                   r_A,
                                   w_A + nDir,
                                   m_map + nDir);
        // <s_{k+1}, w_{k+1}>
        vExchange[1] = Vmath::Dot2(nNonDir,
                                   s_A + nDir,
                                   w_A + nDir,
                                   m_map + nDir);

        // <r_{k+1}, r_{k+1}>
        vExchange[2] = Vmath::Dot2(nNonDir,
                                   r_A,
                                   r_A,
                                   m_map + nDir);

        // Perform inner-product exchanges
        vComm->AllReduce(vExchange, Nektar::LibUtilities::ReduceSum);

        rho_new = vExchange[0];
        mu      = vExchange[1];
        eps     = vExchange[2];

        m_totalIterations++;
        // test if norm is within tolerance
        if (eps < m_tolerance * m_tolerance * m_rhs_magnitude)
        {
            if (m_verbose && m_root)
            {
                cout << "CG iterations made = " << m_totalIterations
                     << " using tolerance of "  << m_tolerance
                     << " (error = " << sqrt(eps/m_rhs_magnitude)
                     << ", rhs_mag = " << sqrt(m_rhs_magnitude) <<  ")"
                     << endl;
            }
            break;
        }
        min_resid = min(min_resid, eps);

        // Compute search direction and solution coefficients
        beta  = rho_new/rho;
        alpha = rho_new/(mu - rho_new*beta/alpha);
        rho   = rho_new;
        k++;
    }
}

void GlobalLinSysIterativeKokkos::Set_Rhs_Magnitude(
    const NekVector<NekDouble> &pIn)
{
    Array<OneD, NekDouble> vExchange(1);
    vExchange[0] = Vmath::Dot(pIn.GetDimension(),&pIn[0],&pIn[0]);

    m_expList.lock()->GetComm()->GetRowComm()->AllReduce(
        vExchange, Nektar::LibUtilities::ReduceSum);

    // To ensure that very different rhs values are not being
    // used in subsequent solvers such as the velocit solve in
    // INC NS. If this works we then need to work out a better
    // way to control this.
    NekDouble new_rhs_mag = (vExchange[0] > 1e-6)? vExchange[0] : 1.0;

    if(m_rhs_magnitude == NekConstants::kNekUnsetDouble)
    {
        m_rhs_magnitude = new_rhs_mag;
    }
    else
    {
        m_rhs_magnitude = (m_rhs_mag_sm*(m_rhs_magnitude) +
                           (1.0-m_rhs_mag_sm)*new_rhs_mag);
    }
}

}
}
