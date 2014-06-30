///////////////////////////////////////////////////////////////////////////////
//
// File: AdvectionWeakDG.cpp
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
// Description: Weak DG advection class.
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Advection/AdvectionWeakDG.h>
#include <iostream>
#include <iomanip>

namespace Nektar
{
    namespace SolverUtils
    {
        std::string AdvectionWeakDG::type[] = {
            GetAdvectionFactory().RegisterCreatorFunction("WeakDG",
                                                    AdvectionWeakDG::create),
            GetAdvectionFactory().RegisterCreatorFunction("AdjointWeakDG",
                                                          AdvectionWeakDG::create),
            GetAdvectionFactory().RegisterCreatorFunction("AdjointWeakDGNS",
                                                          AdvectionWeakDG::create)};

        AdvectionWeakDG::AdvectionWeakDG(std::string advType):m_advType(advType)
        {
        }

        /**
         * @brief Initialise AdvectionWeakDG objects and store them before
         * starting the time-stepping.
         *
         * @param pSession  Pointer to session reader.
         * @param pFields   Pointer to fields.
         */
        void AdvectionWeakDG::v_InitObject(
            LibUtilities::SessionReaderSharedPtr        pSession,
            Array<OneD, MultiRegions::ExpListSharedPtr> pFields)
        {
            Advection::v_InitObject(pSession, pFields);
        }

        /**
         * @brief Compute the advection term at each time-step using the
         * Discontinuous Glaerkin approach (DG).
         *
         * @param nConvectiveFields   Number of fields.
         * @param fields              Pointer to fields.
         * @param advVel              Advection velocities.
         * @param inarray             Solution at the previous time-step.
         * @param outarray            Advection term to be passed at the
         *                            time integration class.
         */
        void AdvectionWeakDG::v_Advect(
            const int                                         nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> >        &advVel,
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                  Array<OneD, Array<OneD, NekDouble> >        &outarray)
        {
            int nDim            = fields[0]->GetCoordim(0);
            int nPointsTot      = fields[0]->GetTotPoints();
            int nCoeffs         = fields[0]->GetNcoeffs();
            int nTracePointsTot = fields[0]->GetTrace()->GetTotPoints();
            int i, j;

            Array<OneD, Array<OneD, NekDouble> >        outarray_tmp(nConvectiveFields);
            Array<OneD, Array<OneD, NekDouble> >        outarray_tmp2(nConvectiveFields);
            Array<OneD, Array<OneD, NekDouble> >        outarray_tmp3(nConvectiveFields);
            Array<OneD, Array<OneD, NekDouble> >        outarray_tmp4(nConvectiveFields);
            Array<OneD, Array<OneD, NekDouble> >        outarray_tmp5(nConvectiveFields);
            
            Array<OneD, Array<OneD, NekDouble> > tmp(nConvectiveFields);
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > fluxvector(
                                                        nConvectiveFields);
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > addfluxvector(
                                                        nConvectiveFields);
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > JacTransVec(
                                                        nConvectiveFields);
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > addJacTransVec(
                                                        nConvectiveFields);
            // Allocate storage for flux vector F(u).
            for (i = 0; i < nConvectiveFields; ++i)
            {
                fluxvector[i] =
                    Array<OneD, Array<OneD, NekDouble> >(m_spaceDim);
                
                for (j = 0; j < m_spaceDim; ++j)
                {
                    fluxvector[i][j] = Array<OneD, NekDouble>(nPointsTot);
                }
            }

            ASSERTL1(m_riemann,
                     "Riemann solver must be provided for AdvectionWeakDG.");
            
            // Retreive the direct solution for the forward problem
            
            if (m_advType == "WeakDG")
            {
                for (i = 0; i < nConvectiveFields; ++i)
                {
                    fluxvector[i] =
                    Array<OneD, Array<OneD, NekDouble> >(m_spaceDim);
                    
                    for (j = 0; j < m_spaceDim; ++j)
                    {
                        fluxvector[i][j] = Array<OneD, NekDouble>(nPointsTot);
                    }
                }
                
                m_fluxVector(inarray, fluxvector);
                // Get the advection part (without numerical flux)
                for(i = 0; i < nConvectiveFields; ++i)
                {
                    tmp[i] = Array<OneD, NekDouble>(nCoeffs, 0.0);
                    
                    outarray_tmp[i] = Array<OneD, NekDouble>(nCoeffs, 0.0);
                    
                    for (j = 0; j < nDim; ++j)
                    {
                        fields[i]->IProductWRTDerivBase(j,
                                                        fluxvector[i][j],
                                                        outarray[i]);
                        
                        
                        Vmath::Vadd(nCoeffs,
                                    outarray[i], 1,
                                    tmp[i], 1,
                                    tmp[i], 1);
                    }
                }
            }
            
            if (m_advType == "AdjointWeakDG")
            {
                for (i = 0; i < nConvectiveFields; ++i)
                {
                    fluxvector[i] =
                    Array<OneD, Array<OneD, NekDouble> >(m_spaceDim);
                    
                    JacTransVec[i] =
                    Array<OneD, Array<OneD, NekDouble> >(m_spaceDim);
                    
                    for (j = 0; j < m_spaceDim; ++j)
                    {
                        fluxvector[i][j] = Array<OneD, NekDouble>(nPointsTot);
                        JacTransVec[i][j] = Array<OneD, NekDouble>(nPointsTot);
                    }
                }
                
                m_fluxVector(inarray, fluxvector);
                
                m_JacTransposeDivVector(inarray, JacTransVec);
                
                // Get the advection part (without numerical flux)
                for(i = 0; i < nConvectiveFields; ++i)
                {
                    tmp[i] = Array<OneD, NekDouble>(nCoeffs, 0.0);
                    
                    outarray_tmp[i] = Array<OneD, NekDouble>(nCoeffs, 0.0);
                    outarray_tmp2[i] = Array<OneD, NekDouble>(nCoeffs, 0.0);
                    
                    
                    
                    for (j = 0; j < nDim; ++j)
                    {
                        fields[i]->IProductWRTDerivBase(j,
                                            fluxvector[i][j],
                                            outarray[i]);
                        
                        // since the adjoint for the convective part is NOT given in conservative form, an additional term where the derivatives of the transposed jacobians are multiplied with the basis hence,
                        
                        fields[i]->IProductWRTBase(JacTransVec[i][j],
                                                   outarray_tmp[i]);
                        
                        Vmath::Vadd(nCoeffs,
                                    outarray[i], 1,
                                    outarray_tmp[i], 1,
                                    outarray_tmp2[i], 1);
                        
                        Vmath::Vadd(nCoeffs,
                                    outarray_tmp2[i], 1,
                                    tmp[i], 1,
                                    tmp[i], 1);
                    }
                }
            }
            
            if (m_advType == "AdjointWeakDGNS")
            {
                for (i = 0; i < nConvectiveFields; ++i)
                {
                    fluxvector[i] =
                    Array<OneD, Array<OneD, NekDouble> >(m_spaceDim);
                    
                    addfluxvector[i] =
                    Array<OneD, Array<OneD, NekDouble> >(m_spaceDim);
                    
                    JacTransVec[i] =
                    Array<OneD, Array<OneD, NekDouble> >(m_spaceDim);
                    
                    addJacTransVec[i] =
                    Array<OneD, Array<OneD, NekDouble> >(m_spaceDim);
                    
                    for (j = 0; j < m_spaceDim; ++j)
                    {
                        fluxvector[i][j] =
                                        Array<OneD, NekDouble>(nPointsTot, 0.0);
                        addfluxvector[i][j] =
                                        Array<OneD, NekDouble>(nPointsTot, 0.0);
                        JacTransVec[i][j] =
                                        Array<OneD, NekDouble>(nPointsTot, 0.0);
                        addJacTransVec[i][j] =
                                        Array<OneD, NekDouble>(nPointsTot,0.0);
                    }
                }
                
                m_fluxVector(inarray, fluxvector);
                m_AddfluxVector(inarray, addfluxvector);
                
                m_JacTransposeDivVector(inarray, JacTransVec);
                m_AddJacTransposeDivVector(inarray, addJacTransVec);
                
                // Get the advection part (without numerical flux)
                for(i = 0; i < nConvectiveFields; ++i)
                {
                    tmp[i] = Array<OneD, NekDouble>(nCoeffs, 0.0);
                    
                    outarray_tmp[i] = Array<OneD, NekDouble>(nCoeffs, 0.0);
                    outarray_tmp2[i] = Array<OneD, NekDouble>(nCoeffs, 0.0);
                    outarray_tmp3[i] = Array<OneD, NekDouble>(nCoeffs, 0.0);
                    outarray_tmp4[i] = Array<OneD, NekDouble>(nCoeffs, 0.0);
                    outarray_tmp5[i] = Array<OneD, NekDouble>(nCoeffs, 0.0);
                    
                    for (j = 0; j < nDim; ++j)
                    {
                        fields[i]->IProductWRTDerivBase(j,
                                                        fluxvector[i][j],
                                                        outarray[i]);
                        
                        fields[i]->IProductWRTDerivBase(j,
                                                        addfluxvector[i][j],
                                                        outarray_tmp3[i]);
                        
                        // since the adjoint for the convective part is NOT given in conservative form, an additional term where the derivatives of the transposed jacobians are multiplied with the basis hence,
                        
                        fields[i]->IProductWRTBase(JacTransVec[i][j],
                                                   outarray_tmp[i]);
                        
                        fields[i]->IProductWRTBase(addJacTransVec[i][j],
                                                   outarray_tmp5[i]);
                        
                        Vmath::Vadd(nCoeffs,
                                    outarray[i], 1,
                                    outarray_tmp[i], 1,
                                    outarray_tmp2[i], 1);
                        
                        Vmath::Vadd(nCoeffs,
                                    outarray_tmp2[i], 1,
                                    outarray_tmp3[i], 1,
                                    outarray_tmp4[i], 1);
                        
                        Vmath::Vadd(nCoeffs,
                                    outarray_tmp4[i], 1,
                                    outarray_tmp5[i], 1,
                                    outarray[i], 1);
                        
                        Vmath::Vadd(nCoeffs,
                                    outarray[i], 1,
                                    tmp[i], 1,
                                    tmp[i], 1);
                    }
                }
            }
            
            // Store forwards/backwards space along trace space
            Array<OneD, Array<OneD, NekDouble> > Fwd    (nConvectiveFields);
            Array<OneD, Array<OneD, NekDouble> > Bwd    (nConvectiveFields);
            Array<OneD, Array<OneD, NekDouble> > numflux(nConvectiveFields);
            Array<OneD, Array<OneD, NekDouble> > numflux2(nConvectiveFields);

            for(i = 0; i < nConvectiveFields; ++i)
            {
                Fwd[i]      = Array<OneD, NekDouble>(nTracePointsTot, 0.0);
                Bwd[i]      = Array<OneD, NekDouble>(nTracePointsTot, 0.0);
                numflux[i]  = Array<OneD, NekDouble>(nTracePointsTot, 0.0);
                numflux2[i] = Array<OneD, NekDouble>(nTracePointsTot, 0.0);
                fields[i]->GetFwdBwdTracePhys(inarray[i], Fwd[i], Bwd[i]);
            }

            m_riemann->Solve(Fwd, Bwd, numflux);
            
            // Evaulate <\phi, \hat{F}\cdot n> - OutField[i]g
            for(i = 0; i < nConvectiveFields; ++i)
            {
                Vmath::Neg                      (nCoeffs, tmp[i], 1);
                fields[i]->AddTraceIntegral     (numflux[i], tmp[i]);
                fields[i]->MultiplyByElmtInvMass(tmp[i], tmp[i]);
                fields[i]->BwdTrans             (tmp[i], outarray[i]);
            }
        }
    }//end of namespace SolverUtils
}//end of namespace Nektar
