///////////////////////////////////////////////////////////////////////////////
//
// File GlobalLinSys.cpp
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
// Description: GlobalLinSys definition
//
///////////////////////////////////////////////////////////////////////////////

#include <map>
#include <MultiRegions/GlobalLinSysIterativeKokkosFull.h>
#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>

using namespace std;

namespace Nektar
{
    namespace MultiRegions
    {
        /**
         * @class GlobalLinSysIterativeCG
         *
         * This class implements a conjugate gradient matrix solver.
         * Preconditioning is implemented using a Jacobi (diagonal)
         * preconditioner.
         */

        /**
         * Registers the class with the Factory.
         */
        string GlobalLinSysIterativeKokkosFull::className
                = GetGlobalLinSysFactory().RegisterCreatorFunction(
                    "IterativeKokkosFull",
                    GlobalLinSysIterativeKokkosFull::create,
                    "Iterative solver for full matrix system using Kokkos.");


        /**
         * Constructor for full direct matrix solve.
         * @param   pKey        Key specifying matrix to solve.
         * @param   pExp        Shared pointer to expansion list for applying
         *                      matrix evaluations.
         * @param   pLocToGloMap Local to global mapping.
         */
        GlobalLinSysIterativeKokkosFull::GlobalLinSysIterativeKokkosFull(
                    const GlobalLinSysKey &pKey,
                    const boost::weak_ptr<ExpList> &pExp,
                    const boost::shared_ptr<AssemblyMap> &pLocToGloMap)
            : GlobalLinSys         (pKey, pExp, pLocToGloMap),
              GlobalLinSysIterativeKokkos(pKey, pExp, pLocToGloMap)
        {
            ASSERTL1(m_linSysKey.GetGlobalSysSolnType()==eIterativeKokkosFull,
                     "This routine should only be used when using an Iterative "
                     "conjugate gradient matrix solve.");

            // Copy the matrix stuff over into Kokkos Views
            int nblks = GetNumBlocks();
            int data_size = 0;
            m_vecOffset = Kokkos::View<int*>("vecoffset", nblks);
            m_matOffset = Kokkos::View<int*>("matoffset", nblks);
            m_matsize   = Kokkos::View<int*>("size", nblks);
            m_rows      = Kokkos::View<int*>("rows", nblks);
            DNekScalMatSharedPtr loc_mat;

            for (int b = 0; b < GetNumBlocks(); ++b) {
                loc_mat = GetBlock(b);
                int r = loc_mat->GetRows();
                int c = loc_mat->GetColumns();
                m_matsize(b) = r*c;
                m_rows(b) = r;
                m_vecOffset(b) = (b == 0 ? 0 : m_vecOffset(b-1) + m_rows(b-1));
                m_matOffset(b) = (b == 0 ? 0 : m_matOffset(b-1) + m_matsize(b-1));
                data_size += m_matsize(b);
            }

            m_data = Kokkos::View<NekDouble*>("data", data_size);

            for (int b = 0; b < GetNumBlocks(); ++b) {
                loc_mat = GetBlock(b);
                int r = loc_mat->GetRows();
                int c = loc_mat->GetColumns();
                for (int i = 0; i < r; ++i) {
                    for (int j = 0; j < c; ++j) {
                        m_data(m_matOffset(b) + i*c + j) = (*loc_mat)(i,j);
                    }
                }
            }
        }


        /**
         *
         */
        GlobalLinSysIterativeKokkosFull::~GlobalLinSysIterativeKokkosFull()
        {

        }


        /**
         * Solve a global linear system with Dirichlet forcing using a
         * conjugate gradient method. This routine performs handling of the
         * Dirichlet forcing terms and wraps the underlying iterative solver
         * used for the remaining degrees of freedom.
         *
         * Consider solving for \f$x\f$, the matrix system \f$Ax=b\f$, where
         * \f$b\f$ is known. To enforce the Dirichlet terms we instead solve
         * \f[A(x-x_0) = b - Ax_0 \f]
         * where \f$x_0\f$ is the Dirichlet forcing.
         *
         * @param           pInput      RHS of linear system, \f$b\f$.
         * @param           pOutput     On input, values of dirichlet degrees
         *                              of freedom with initial guess on other values.
         *                              On output, the solution \f$x\f$.
         * @param           pLocToGloMap    Local to global mapping.
         * @param           pDirForcing Precalculated Dirichlet forcing.
         */
        void GlobalLinSysIterativeKokkosFull::v_Solve(
                    const Array<OneD, const NekDouble>  &pInput,
                          Array<OneD,       NekDouble>  &pOutput,
                    const AssemblyMapSharedPtr &pLocToGloMap,
                    const Array<OneD, const NekDouble>  &pDirForcing)
        {
            boost::shared_ptr<MultiRegions::ExpList> expList = m_expList.lock();
            bool vCG;
            if ((m_locToGloMap = boost::dynamic_pointer_cast<AssemblyMapCG>(
                     pLocToGloMap)))
            {
                vCG = true;
            }
            else if ((m_locToGloMap = boost::dynamic_pointer_cast<
                          AssemblyMapDG>(pLocToGloMap)))
            {
                vCG = false;
            }
            else
            {
                ASSERTL0(false, "Unknown map type");
            }

            bool dirForcCalculated = (bool) pDirForcing.num_elements();
            int nDirDofs  = pLocToGloMap->GetNumGlobalDirBndCoeffs();
            int nGlobDofs = pLocToGloMap->GetNumGlobalCoeffs();
            int nDirTotal = nDirDofs;
            
            expList->GetComm()->AllReduce(nDirTotal, LibUtilities::ReduceSum);
            
            Array<OneD, NekDouble> tmp(nGlobDofs), tmp2;

            if(nDirTotal)
            {
                // calculate the Dirichlet forcing
                if(dirForcCalculated)
                {
                    Vmath::Vsub(nGlobDofs, pInput.get(), 1,
                                pDirForcing.get(), 1,
                                tmp.get(), 1);
                }
                else
                {
                    // Calculate the dirichlet forcing B_b (== X_b) and
                    // substract it from the rhs
                    expList->GeneralMatrixOp(
                        m_linSysKey, pOutput, tmp, eGlobal);

                    Vmath::Vsub(nGlobDofs, pInput.get(), 1,
                                           tmp.get(),    1,
                                           tmp.get(),    1);
                }
                if (vCG)
                {
                    Array<OneD, NekDouble> out(nGlobDofs,0.0);

                    // solve for perturbation from intiial guess in pOutput
                    SolveLinearSystem(
                        nGlobDofs, tmp, out, pLocToGloMap, nDirDofs);
                    Vmath::Vadd(nGlobDofs-nDirDofs,    &out    [nDirDofs], 1,
                                &pOutput[nDirDofs], 1, &pOutput[nDirDofs], 1);
                }
                else
                {
                    ASSERTL0(false, "Need DG solve if using Dir BCs");
                }
            }
            else
            {
                Vmath::Vcopy(nGlobDofs, pInput, 1, tmp, 1);
                SolveLinearSystem(nGlobDofs, tmp, pOutput, pLocToGloMap);
            }
        }


        /**
         *
         */
        void GlobalLinSysIterativeKokkosFull::v_DoMatrixMultiply(
                const Array<OneD, NekDouble>& pInput,
                      Array<OneD, NekDouble>& pOutput)
        {
            boost::shared_ptr<MultiRegions::ExpList> expList = m_expList.lock();
            // Perform matrix-vector operation A*d_i
//            expList->GeneralMatrixOp(m_linSysKey,
//                                     pInput, pOutput, eGlobal);

            // Copy the matrix stuff over into Kokkos Views
            int nblks = GetNumBlocks();
//            int data_size = 0;
//            Kokkos::View<int*> vecOffset("vecoffset", nblks);
//            Kokkos::View<int*> matOffset("matoffset", nblks);
//            Kokkos::View<int*> matsize("size", nblks);
//            Kokkos::View<int*> rows("rows", nblks);
//            DNekScalMatSharedPtr loc_mat;
//
//            for (int b = 0; b < GetNumBlocks(); ++b) {
//                loc_mat = GetBlock(b);
//                int r = loc_mat->GetRows();
//                int c = loc_mat->GetColumns();
//                matsize(b) = r*c;
//                rows(b) = r;
//                vecOffset(b) = (b == 0 ? 0 : vecOffset(b-1) + rows(b-1));
//                matOffset(b) = (b == 0 ? 0 : matOffset(b-1) + matsize(b-1));
//                data_size += matsize(b);
//            }
//
//            Kokkos::View<NekDouble*> data("data", data_size);
//
//            for (int b = 0; b < GetNumBlocks(); ++b) {
//                loc_mat = GetBlock(b);
//                int r = loc_mat->GetRows();
//                int c = loc_mat->GetColumns();
//                for (int i = 0; i < r; ++i) {
//                    for (int j = 0; j < c; ++j) {
//                        data(matOffset(b) + i*c + j) = (*loc_mat)(i,j);
//                    }
//                }
//            }

            // Copy input
            int nlocdof = expList->GetNcoeffs();
            Array<OneD, NekDouble> tmp1(nlocdof);
            m_locToGloMap->GlobalToLocal(pInput, tmp1);

            Kokkos::View<NekDouble*> input("input", nlocdof);
            Kokkos::View<NekDouble*> output("output", nlocdof);
            for (int i = 0; i < nlocdof; ++i) {
                input(i) = tmp1[i];
            }

            // Do matrix-vector multiply
            Kokkos::parallel_for(nblks, [&] (const int& b) {
                const int elOffset = m_vecOffset(b), nRows = m_rows(b);
                const int matOff = m_matOffset(b);
                for (int i = 0; i < nRows; ++i)
                {
                    output(i + elOffset) = 0.0;
                    for (int j = 0; j < nRows; ++j)
                    {
                        output(i + elOffset) += m_data(matOff + i*nRows + j) * input(j + elOffset);
                    }
                }
            });

            for (int i = 0; i < nlocdof; ++i)
            {
                tmp1[i] = output(i);
            }

            m_locToGloMap->Assemble(tmp1, pOutput);

            // Apply robin boundary conditions to the solution.
            if(m_robinBCInfo.size() > 0)
            {
                ASSERTL0(false,
                        "Robin boundaries not set up in IterativeFull solver.");
                int nGlobal = m_locToGloMap->GetNumGlobalCoeffs();
                int nLocal  = m_locToGloMap->GetNumLocalCoeffs();
                int nDir    = m_locToGloMap->GetNumGlobalDirBndCoeffs();
                int nNonDir = nGlobal - nDir;
                Array<OneD, NekDouble> robin_A(nGlobal, 0.0);
                Array<OneD, NekDouble> robin_l(nLocal,  0.0);
                Array<OneD, NekDouble> tmp;
                NekVector<NekDouble> robin(nNonDir,
                                           tmp = robin_A + nDir, eWrapper);

                // Operation: p_A = A * d_A
                // First map d_A to local solution
                m_locToGloMap->GlobalToLocal(pInput, robin_l);

                // Iterate over all the elements computing Robin BCs where
                // necessary
                for (int n = 0; n < expList->GetNumElmts(); ++n)
                {
                    int nel = expList->GetOffset_Elmt_Id(n);
                    int offset = expList->GetCoeff_Offset(n);
                    int ncoeffs = expList->GetExp(nel)->GetNcoeffs();

                    if(m_robinBCInfo.count(nel) != 0) // add robin mass matrix
                    {
                        RobinBCInfoSharedPtr rBC;
                        Array<OneD, NekDouble> tmp;
                        StdRegions::StdExpansionSharedPtr vExp = expList->GetExp(nel);

                        // add local matrix contribution
                        for(rBC = m_robinBCInfo.find(nel)->second;rBC; rBC = rBC->next)
                        {
                            vExp->AddRobinEdgeContribution(rBC->m_robinID,rBC->m_robinPrimitiveCoeffs, tmp = robin_l + offset);
                        }
                    }
                    else
                    {
                        Vmath::Zero(ncoeffs, &robin_l[offset], 1);
                    }
                }

                // Map local Robin contribution back to global coefficients
                m_locToGloMap->LocalToGlobal(robin_l, robin_A);
                // Add them to the output of the GeneralMatrixOp
                Vmath::Vadd(nGlobal, pOutput, 1, robin_A, 1, pOutput, 1);
            }
        }

        /**
         *
         */
        void GlobalLinSysIterativeKokkosFull::v_UniqueMap()
        {
            m_map = m_locToGloMap->GetGlobalToUniversalMapUnique();
        }

    }
}
