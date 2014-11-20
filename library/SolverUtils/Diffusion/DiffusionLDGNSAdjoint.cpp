///////////////////////////////////////////////////////////////////////////////
//
// File: DiffusionLDGNS.cpp
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
// Description: LDGNS diffusion class.
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Diffusion/DiffusionLDGNSAdjoint.h>
#include <iostream>
#include <iomanip>

namespace Nektar
{
    namespace SolverUtils
    {
        std::string DiffusionLDGNSAdjoint::type = GetDiffusionFactory().
        RegisterCreatorFunction("LDGNSAdjoint", DiffusionLDGNSAdjoint::create);
        
        DiffusionLDGNSAdjoint::DiffusionLDGNSAdjoint()
        {
        }
        
        void DiffusionLDGNSAdjoint::v_InitObject(
                                                 LibUtilities::SessionReaderSharedPtr        pSession,
                                                 Array<OneD, MultiRegions::ExpListSharedPtr> pFields)
        {
            m_session = pSession;
            m_session->LoadParameter ("Gamma",         m_gamma, 1.4);
            m_session->LoadParameter ("GasConstant",   m_gasConstant, 287.058);
            m_session->LoadParameter ("Twall",         m_Twall, 300.15);
            m_session->LoadSolverInfo("ViscosityType", m_ViscosityType,
                                      "Constant");
            m_session->LoadParameter ("mu",            m_mu, 1.78e-05);
            m_session->LoadParameter ("thermalConductivity",
                                      m_thermalConductivity, 0.0257);
            m_session->LoadParameter ("rhoInf",        m_rhoInf, 1.225);
            m_session->LoadParameter ("pInf",          m_pInf, 101325);
            m_session->LoadParameter ("rhoInfPrimal",   m_rhoInfPrimal, 0.0);
            m_session->LoadParameter ("uInfPrimal",     m_uInfPrimal, 1.0);
            m_session->LoadParameter ("vInfPrimal",     m_vInfPrimal, 1.0);
            m_session->LoadParameter ("pInfPrimal",     m_pInfPrimal, 1.0);
            m_session->LoadParameter ("alphaInfPrimal", m_alphaInfDir, 0.0);
            m_session->LoadParameter ("Lref",           m_Lref, 1.0);
            m_session->LoadSolverInfo("Target",         m_Target,  "NoN");
            
            // Setting up the normals
            int i;
            int nDim = pFields[0]->GetCoordim(0);
            int nTracePts = pFields[0]->GetTrace()->GetTotPoints();
            
            m_spaceDim = nDim;
            if (pSession->DefinesSolverInfo("HOMOGENEOUS"))
            {
                m_spaceDim = 3;
            }
            
            m_diffDim = m_spaceDim - nDim;
            
            m_traceVel = Array<OneD, Array<OneD, NekDouble> >(m_spaceDim);
            m_traceNormals = Array<OneD, Array<OneD, NekDouble> >(m_spaceDim);
            for(i = 0; i < m_spaceDim; ++i)
            {
                m_traceVel[i] = Array<OneD, NekDouble> (nTracePts, 0.0);
                m_traceNormals[i] = Array<OneD, NekDouble> (nTracePts);
            }
            pFields[0]->GetTrace()->GetNormals(m_traceNormals);
        }
        
        /**
         * @brief Calculate weak DG Diffusion in the LDG form for the
         * Navier-Stokes (NS) equations:
         *
         * \f$ \langle\psi, \hat{u}\cdot n\rangle
         *   - \langle\nabla\psi \cdot u\rangle
         *     \langle\phi, \hat{q}\cdot n\rangle -
         *     (\nabla \phi \cdot q) \rangle \f$
         *
         * The equations that need a diffusion operator are those related
         * with the velocities and with the energy.
         *
         */
        void DiffusionLDGNSAdjoint::v_Diffuse(
                    const int                                         nConvectiveFields,
                    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                    const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                          Array<OneD, Array<OneD, NekDouble> >        &outarray)
        {
            int i, j;
            int nDim      = fields[0]->GetCoordim(0);
            int nScalars  = inarray.num_elements();
            int nPts      = fields[0]->GetTotPoints();
            int nCoeffs   = fields[0]->GetNcoeffs();
            int nTracePts = fields[0]->GetTrace()->GetTotPoints();
            
            Array<OneD, NekDouble>               tmp1(nCoeffs);
            Array<OneD, Array<OneD, NekDouble> > tmp2(nConvectiveFields);
            
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > >
            numericalFluxO1(m_spaceDim);
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > >
            derivativesO1(m_spaceDim);
            
            Array<OneD, Array<OneD, NekDouble> > fluxvector(m_spaceDim);
            
            for (j = 0; j < m_spaceDim; ++j)
            {
                numericalFluxO1[j] = Array<OneD, Array<OneD, NekDouble> >(
                                                                nScalars);
                derivativesO1[j]   = Array<OneD, Array<OneD, NekDouble> >(
                                                                nScalars);
                
                for (i = 0; i < nScalars; ++i)
                {
                    numericalFluxO1[j][i] = Array<OneD, NekDouble>(
                                                                   nTracePts, 0.0);
                    derivativesO1[j][i]   = Array<OneD, NekDouble>(nPts, 0.0);
                }
            }
            
            // Compute the numerical fluxes for the first order derivatives
            v_NumericalFluxO1(fields, inarray, numericalFluxO1);
            
            for (j = 0; j < nDim; ++j)
            {
                for (i = 0; i < nScalars; ++i)
                {
                    fields[i]->IProductWRTDerivBase (j, inarray[i], tmp1);
                    Vmath::Neg                      (nCoeffs, tmp1, 1);
                    fields[i]->AddTraceIntegral     (numericalFluxO1[j][i],
                                                     tmp1);
                    fields[i]->SetPhysState         (false);
                    fields[i]->MultiplyByElmtInvMass(tmp1, tmp1);
                    fields[i]->BwdTrans             (tmp1, derivativesO1[j][i]);
                }
            }
            
            // For 3D Homogeneous 1D only take derivatives in 3rd direction
            if (m_diffDim == 1)
            {
                for (i = 0; i < nScalars; ++i)
                {
                    derivativesO1[2][i] = m_homoDerivs[i];
                }
            }
            
            // Initialisation viscous tensor
            m_viscTensor = Array<OneD, Array<OneD, Array<OneD, NekDouble> > >
            (m_spaceDim);
            Array<OneD, Array<OneD, NekDouble> > viscousFlux(nConvectiveFields);
            
            for (j = 0; j < m_spaceDim; ++j)
            {
                m_viscTensor[j] = Array<OneD, Array<OneD, NekDouble> >(
                                                                       nScalars);
                for (i = 0; i < nScalars; ++i)
                {
                    m_viscTensor[j][i] = Array<OneD, NekDouble>(nPts, 0.0);
                }
            }
            
            for (i = 0; i < nConvectiveFields; ++i)
            {
                viscousFlux[i] = Array<OneD, NekDouble>(nPts, 0.0);
            }
            
            m_fluxVectorNS(inarray, derivativesO1, m_viscTensor);
            
            // Compute u from q_{\eta} and q_{\xi}
            // Obtain numerical fluxes
            v_NumericalFluxO2(fields, inarray, m_viscTensor, viscousFlux);
            
            for (i = 0; i < nConvectiveFields; ++i)
            {
                tmp2[i] = Array<OneD, NekDouble>(nCoeffs, 0.0);
                
                for (j = 0; j < nDim; ++j)
                {
                    fields[i]->IProductWRTDerivBase(j, m_viscTensor[j][i], tmp1);
                    Vmath::Vadd(nCoeffs, tmp1, 1, tmp2[i], 1, tmp2[i], 1);
                }
                
                // Evaulate  <\phi, \hat{F}\cdot n> - outarray[i]
                Vmath::Neg                      (nCoeffs, tmp2[i], 1);
                fields[i]->AddTraceIntegral     (viscousFlux[i], tmp2[i]);
                fields[i]->SetPhysState         (false);
                fields[i]->MultiplyByElmtInvMass(tmp2[i], tmp2[i]);
                fields[i]->BwdTrans             (tmp2[i], outarray[i]);
            }
        }
        
        /**
         * @brief Builds the numerical flux for the 1st order derivatives
         *
         */
        void DiffusionLDGNSAdjoint::v_NumericalFluxO1(
                const Array<OneD, MultiRegions::ExpListSharedPtr>        &fields,
                const Array<OneD, Array<OneD, NekDouble> >               &inarray,
                      Array<OneD, Array<OneD, Array<OneD, NekDouble> > >
                                                &numericalFluxO1)
        {
            int i, j;
            int nTracePts  = fields[0]->GetTrace()->GetTotPoints();
            int nScalars   = inarray.num_elements();
            int nDim       = fields[0]->GetCoordim(0);
            
            Array<OneD, NekDouble > Vn      (nTracePts, 0.0);
            Array<OneD, NekDouble > fluxtemp(nTracePts, 0.0);
            
            // Get the normal velocity Vn
            for(i = 0; i < nDim; ++i)
            {
                
                fields[0]->AverageTracePhys(inarray[i], m_traceVel[i]);
                Vmath::Vvtvp(nTracePts, m_traceNormals[i], 1,
                             m_traceVel[i], 1, Vn, 1, Vn, 1);
                
                //fields[0]->ExtractTracePhys(inarray[i], m_traceVel[i]);
                //Vmath::Vvtvp(nTracePts, m_traceNormals[i], 1,
                //             m_traceVel[i], 1, Vn, 1, Vn, 1);
            }
            
            // Store forwards/backwards space along trace space
            Array<OneD, Array<OneD, NekDouble> > Fwd    (nScalars);
            Array<OneD, Array<OneD, NekDouble> > Bwd    (nScalars);
            Array<OneD, Array<OneD, NekDouble> > numflux(nScalars);
            
            for (i = 0; i < nScalars; ++i)
            {
                Fwd[i]     = Array<OneD, NekDouble>(nTracePts);
                Bwd[i]     = Array<OneD, NekDouble>(nTracePts);
                numflux[i] = Array<OneD, NekDouble>(nTracePts);
                fields[i]->GetFwdBwdTracePhys(inarray[i], Fwd[i], Bwd[i]);
                fields[0]->GetTrace()->Upwind(Vn, Fwd[i], Bwd[i], numflux[i]);
            }
            
            // Modify the values in case of boundary interfaces
            if (fields[0]->GetBndCondExpansions().num_elements())
            {
                v_WeakPenaltyO1(fields, inarray, numflux);
            }
            
            // Splitting the numerical flux into the dimensions
            for (j = 0; j < m_spaceDim; ++j)
            {
                for (i = 0; i < nScalars; ++i)
                {
                    Vmath::Vmul(nTracePts, m_traceNormals[j], 1,
                                numflux[i], 1, numericalFluxO1[j][i], 1);
                }
            }
        }
        
        /**
         * @brief Imposes appropriate bcs for the 1st order derivatives
         *
         */
        void DiffusionLDGNSAdjoint::v_WeakPenaltyO1(
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                      Array<OneD, Array<OneD, NekDouble> >        &penaltyfluxO1)
        {
            int cnt = 0;
            int i, j, e;
            int id1, id2;
            int eMax;
            
            int nBCEdgePts, nBndEdges, nBndRegions;
            
            int nTracePts = fields[0]->GetTrace()->GetTotPoints();
            int nScalars  = inarray.num_elements();
            const Array<OneD, const int> &traceBndMap
            = fields[0]->GetTraceBndMap();
            
            Array<OneD, NekDouble> tmp1(nTracePts, 0.0);
            Array<OneD, NekDouble> tmp2(nTracePts, 0.0);
            
            Array< OneD, Array<OneD, NekDouble > > scalarVariables(nScalars);
            Array< OneD, Array<OneD, NekDouble > > Fwd(nScalars);
            Array< OneD, Array<OneD, NekDouble > > Fwdnew(nScalars);
            // Extract internal values of the scalar variables for Neumann bcs
            for (i = 0; i < nScalars; ++i)
            {
                scalarVariables[i] = Array<OneD, NekDouble>(nTracePts, 0.0);
                
                Fwd[i] = Array<OneD, NekDouble>(nTracePts, 0.0);
                Fwdnew[i] = Array<OneD, NekDouble>(nTracePts);
                fields[i]->ExtractTracePhys(inarray[i], Fwd[i]);
                fields[i]->ExtractTracePhys(inarray[i], Fwdnew[i]);
            }
            
            nBndRegions = fields[0]->
            GetBndCondExpansions().num_elements();
            
            for (int n = 0; n < nBndRegions; ++n)
            {
                eMax = fields[0]->GetBndCondExpansions()[n]->GetExpSize();
                
                if (fields[0]->GetBndConditions()[n]->
                    GetBoundaryConditionType() ==
                    SpatialDomains::eDirichlet &&
                    fields[0]->GetBndConditions()[n]->
                    GetUserDefined() ==
                    SpatialDomains::eAdjointWall)
                {
                    for (e = 0; e < eMax; ++e)
                    {
                        nBCEdgePts = fields[0]->GetBndCondExpansions()[n]->
                        GetExp(e)->GetTotPoints();
                        id1 = fields[0]->GetBndCondExpansions()[n]->
                        GetPhys_Offset(e);
                        id2 = fields[0]->GetTrace()->GetPhys_Offset(traceBndMap[cnt+e]);
                    
                        Array<OneD, NekDouble> tmp(nBCEdgePts, 0.0);
                        Array<OneD, NekDouble> tmp1(nBCEdgePts, 0.0);
                        Array<OneD, NekDouble> tmp2(nBCEdgePts, 0.0);
                        Array<OneD, NekDouble> tmp3(nBCEdgePts, 0.0);
                        
                        /*
                         NekDouble Cinf = 0.5 * m_rhoInfPrimal
                         * (m_uInfPrimal*m_uInfPrimal+m_vInfPrimal*m_vInfPrimal) * m_Lref;
                         */
                        NekDouble norm_fac = 1.0;
                        
                        Array<OneD, Array<OneD, NekDouble> > DragDir(m_spaceDim);
                        Array<OneD, Array<OneD, NekDouble> > LiftDir(m_spaceDim);
                        
                        NekDouble Dx = norm_fac * cos (m_alphaInfDir);
                        NekDouble Dy = norm_fac * sin (m_alphaInfDir);
                            
                        NekDouble Lx = -norm_fac * sin (m_alphaInfDir);
                        NekDouble Ly =  norm_fac * cos (m_alphaInfDir);
                            
                        DragDir[0] = Array<OneD, NekDouble> (nBCEdgePts, Dx);
                        DragDir[1] = Array<OneD, NekDouble> (nBCEdgePts, Dy);
                            
                        LiftDir[0] = Array<OneD, NekDouble> (nBCEdgePts, Lx);
                        LiftDir[1] = Array<OneD, NekDouble> (nBCEdgePts, Ly);
                            
                            
                        Array<OneD, NekDouble> zeros(nBCEdgePts, 0.0);
                            
                        Vmath::Vcopy(nBCEdgePts,
                                        &zeros[0], 1,
                                        &penaltyfluxO1[0][id2], 1);
                            
                        Vmath::Vcopy(nBCEdgePts,
                                        &zeros[0], 1,
                                        &penaltyfluxO1[1][id2], 1);
                        
                        Vmath::Vcopy(nBCEdgePts,
                                        &zeros[0], 1,
                                        &penaltyfluxO1[2][id2], 1);
                        
                        Vmath::Vcopy(nBCEdgePts,
                                        &zeros[0], 1,
                                        &penaltyfluxO1[m_spaceDim+1][id2], 1);
                            
                        if (m_Target == "Lift")
                        {
                                
                            Vmath::Vcopy(nBCEdgePts,
                                            &LiftDir[0][0], 1,
                                            &penaltyfluxO1[1][id2], 1);
                            
                            Vmath::Vcopy(nBCEdgePts,
                                            &LiftDir[1][0], 1,
                                            &penaltyfluxO1[2][id2], 1);
                            
                        }
                        if (m_Target == "Drag")
                        {
                            Vmath::Vcopy(nBCEdgePts,
                                            &DragDir[0][0], 1,
                                            &penaltyfluxO1[1][id2], 1);
                            
                            Vmath::Vcopy(nBCEdgePts,
                                            &DragDir[1][0], 1,
                                            &penaltyfluxO1[2][id2], 1);
                        }
                    }
                    
                    cnt += fields[0]->GetBndCondExpansions()[n]->GetExpSize();
                }
        
                if (fields[0]->GetBndConditions()[n]->
                    GetBoundaryConditionType() ==
                    SpatialDomains::eDirichlet &&
                    fields[0]->GetBndConditions()[n]->
                    GetUserDefined() !=
                    SpatialDomains::eAdjointWall)
                {
                    //cout << n << "  " << eMax << endl;
                    
                    for (e = 0; e < eMax; ++e)
                    {

                        nBCEdgePts = fields[0]->GetBndCondExpansions()[n]->
                        GetExp(e)->GetTotPoints();
                        id1 = fields[0]->GetBndCondExpansions()[n]->
                        GetPhys_Offset(e);
                        id2 = fields[0]->GetTrace()->GetPhys_Offset(traceBndMap[cnt+e]);
                        
                        // Standard value 0 for farfield BCs
                        Array<OneD, NekDouble> zeros(nBCEdgePts, 0.0);
                        
                        for (i = 0;i < nScalars; ++i)
                        {
                            Vmath::Vcopy(nBCEdgePts,
                                         &(fields[i]->GetBndCondExpansions()[n]->
                                           UpdatePhys())[id1], 1,
                                         &penaltyfluxO1[i][id2], 1);
                        }
                        
                    }
                    
                    cnt += fields[0]->GetBndCondExpansions()[n]->GetExpSize();
                }
            }
        }
        /**
         * @brief Build the numerical flux for the 2nd order derivatives
         *
         */
        void DiffusionLDGNSAdjoint::v_NumericalFluxO2(
                const Array<OneD, MultiRegions::ExpListSharedPtr>        &fields,
                const Array<OneD, Array<OneD, NekDouble> >               &ufield,
                      Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &qfield,
                      Array<OneD, Array<OneD, NekDouble> >               &qflux)
        {
            int i, j;
            int nTracePts = fields[0]->GetTrace()->GetTotPoints();
            int nVariables   = fields.num_elements();
            int nDim         = fields[0]->GetCoordim(0);
            
            Array<OneD, NekDouble > Fwd(nTracePts);
            Array<OneD, NekDouble > Bwd(nTracePts);
            Array<OneD, NekDouble > Vn (nTracePts, 0.0);
            
            Array<OneD, NekDouble > qFwd     (nTracePts);
            Array<OneD, NekDouble > qBwd     (nTracePts);
            Array<OneD, NekDouble > qfluxtemp(nTracePts, 0.0);
            
            // Get the normal velocity Vn
            for(i = 0; i < nDim; ++i)
            {
                fields[0]->AverageTracePhys(ufield[i], m_traceVel[i]);
                Vmath::Vvtvp(nTracePts, m_traceNormals[i], 1,
                             m_traceVel[i], 1, Vn, 1, Vn, 1);
                
                //fields[0]->ExtractTracePhys(ufield[i], m_traceVel[i]);
                //Vmath::Vvtvp(nTracePts, m_traceNormals[i], 1,
                //             m_traceVel[i], 1, Vn, 1, Vn, 1);
            }
            
            // Evaulate Riemann flux
            // qflux = \hat{q} \cdot u = q \cdot n
            for (i = 0; i < nVariables; ++i)
            {
                qflux[i] = Array<OneD, NekDouble> (nTracePts, 0.0);
                for (j = 0; j < nDim; ++j)
                {
                    // Compute qFwd and qBwd value of qfield in position 'ji'
                    fields[i]->GetFwdBwdTracePhys(qfield[j][i], qFwd, qBwd);
                    
                    // Get Riemann flux of qflux --> LDG implies upwind
                    fields[i]->GetTrace()->Upwind(Vn, qBwd, qFwd, qfluxtemp);
                    
                    // Multiply the Riemann flux by the trace normals
                    Vmath::Vmul(nTracePts, m_traceNormals[j], 1, qfluxtemp, 1,
                                qfluxtemp, 1);
                    
                    // Impose weak boundary condition with flux
                    if (fields[0]->GetBndCondExpansions().num_elements())
                    {
                        v_WeakPenaltyO2(fields, i, j, qfield[j][i], qfluxtemp);
                    }
                    
                    // Store the final flux into qflux
                    Vmath::Vadd(nTracePts, qfluxtemp, 1, qflux[i], 1,
                                qflux[i], 1);
                }
            }
        }
        
        
        /**
         * @brief Imposes appropriate bcs for the 2nd order derivatives
         *
         */
        void DiffusionLDGNSAdjoint::v_WeakPenaltyO2(
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const int                                          var,
                const int                                          dir,
                const Array<OneD, const NekDouble>                &qfield,
                      Array<OneD,       NekDouble>                &penaltyflux)
        {
            int cnt = 0;
            int nBndEdges, nBndEdgePts;
            int i, e;
            int id2;
            
            int nTracePts   = fields[0]->GetTrace()->GetTotPoints();
            int nBndRegions = fields[var]->GetBndCondExpansions().num_elements();
            
            Array<OneD, NekDouble > uterm(nTracePts);
            Array<OneD, NekDouble > qtemp(nTracePts);
            
            // Extract the physical values of the solution at the boundaries
            fields[var]->ExtractTracePhys(qfield, qtemp);
            
            // Loop on the boundary regions to apply appropriate bcs
            for (i = 0; i < nBndRegions; ++i)
            {
                // Number of boundary regions related to region 'i'
                nBndEdges = fields[var]->
                GetBndCondExpansions()[i]->GetExpSize();
                
                // Weakly impose bcs by modifying flux values
                for (e = 0; e < nBndEdges; ++e)
                {
                    nBndEdgePts = fields[var]->
                    GetBndCondExpansions()[i]->GetExp(e)->GetNumPoints(0);
                    
                    id2 = fields[0]->GetTrace()->
                    GetPhys_Offset(fields[0]->GetTraceMap()->
                                   GetBndCondTraceToGlobalTraceMap(cnt++));
                    
                    // In case of Dirichlet bcs:
                    // uflux = gD
                    if(fields[var]->GetBndConditions()[i]->
                       GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                    {
                        Vmath::Vmul(nBndEdgePts, 
                                    &m_traceNormals[dir][id2], 1, 
                                    &qtemp[id2], 1, 
                                    &penaltyflux[id2], 1);
                    }
                    // 3.4) In case of Neumann bcs: 
                    // uflux = u+
                    else if((fields[var]->GetBndConditions()[i])->
                            GetBoundaryConditionType() == SpatialDomains::eNeumann)
                    {
                        ASSERTL0(false, 
                                 "Neumann bcs not implemented for LDGNS");
                        
                        /*
                         Vmath::Vmul(nBndEdgePts, 
                         &m_traceNormals[dir][id2], 1, 
                         &(fields[var]->
                         GetBndCondExpansions()[i]->
                         UpdatePhys())[id1], 1, 
                         &penaltyflux[id2], 1);
                         */
                    }
                }
            }
        }
    }//end of namespace SolverUtils
}//end of namespace Nektar
