///////////////////////////////////////////////////////////////////////////////
//
// File CompressibleFlowSystem.cpp
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
// Description: Auxiliary functions for the compressible flow system
//
///////////////////////////////////////////////////////////////////////////////
#include <iomanip>

#include <CompressibleFlowSolver/EquationSystems/CompressibleFlowSystem.h>
#include <LocalRegions/TriExp.h>
#include <LocalRegions/QuadExp.h>
#include <LocalRegions/HexExp.h>
#include <MultiRegions/ExpList.h>
#include <LibUtilities/Foundations/InterpCoeff.h>
#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>
#include <math.h>
#define PI 3.14159265
namespace Nektar
{
    string CompressibleFlowSystem::className =
    SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
                    "CompressibleFlowSystem",
                    CompressibleFlowSystem::create,
                    "Auxiliary functions for the compressible flow system.");
    
    CompressibleFlowSystem::CompressibleFlowSystem(
                    const LibUtilities::SessionReaderSharedPtr& pSession)
    : UnsteadySystem(pSession)
    {
    }
    
    /**
     * @brief Initialization object for CompressibleFlowSystem class.
     */
    void CompressibleFlowSystem::v_InitObject()
    {
        UnsteadySystem::v_InitObject();
        
        ASSERTL0(m_session->DefinesSolverInfo("UPWINDTYPE"),
                 "No UPWINDTYPE defined in session.");
        
        // Do not forwards transform initial condition
        m_homoInitialFwd = false;
        
        // Set up locations of velocity vector.
        m_vecLocs = Array<OneD, Array<OneD, NekDouble> >(1);
        m_vecLocs[0] = Array<OneD, NekDouble>(m_spacedim);
        for (int i = 0; i < m_spacedim; ++i)
        {
            m_vecLocs[0][i] = 1 + i;
        }
        
        // Get gamma parameter from session file.
        ASSERTL0(m_session->DefinesParameter("Gamma"),
                 "Compressible flow sessions must define a Gamma parameter.");
        m_session->LoadParameter("Gamma", m_gamma, 1.4);
        
        // Get E0 parameter from session file.
        ASSERTL0(m_session->DefinesParameter("pInf"),
                 "Compressible flow sessions must define a pInf parameter.");
        m_session->LoadParameter("pInf", m_pInf, 101325);
        
        // Get rhoInf parameter from session file.
        ASSERTL0(m_session->DefinesParameter("rhoInf"),
                 "Compressible flow sessions must define a rhoInf parameter.");
        m_session->LoadParameter("rhoInf", m_rhoInf, 1.225);
        
        // Get uInf parameter from session file.
        ASSERTL0(m_session->DefinesParameter("uInf"),
                 "Compressible flow sessions must define a uInf parameter.");
        m_session->LoadParameter("uInf", m_uInf, 0.1);
        
        // Get vInf parameter from session file.
        if (m_spacedim == 2 || m_spacedim == 3)
        {
            ASSERTL0(m_session->DefinesParameter("vInf"),
                     "Compressible flow sessions must define a vInf parameter"
                     "for 2D/3D problems.");
            m_session->LoadParameter("vInf", m_vInf, 0.0);
        }
        
        // Get wInf parameter from session file.
        if (m_spacedim == 3)
        {
            ASSERTL0(m_session->DefinesParameter("wInf"),
                     "Compressible flow sessions must define a wInf parameter"
                     "for 3D problems.");
            m_session->LoadParameter("wInf", m_wInf, 0.0);
        }
        m_session->LoadParameter ("GasConstant",    m_gasConstant,   287.058);
        m_session->LoadParameter ("Twall",          m_Twall,         300.15);
        m_session->LoadSolverInfo("ViscosityType",  m_ViscosityType, "Constant");
        m_session->LoadSolverInfo("Target",         m_Target,        "NoN");
        m_session->LoadParameter ("mu",             m_mu,            1.78e-05);
        m_session->LoadParameter ("thermalConductivity",
                                  m_thermalConductivity, 0.0257);
        m_session->LoadParameter ("adjointSwitch",  m_adjointSwitch, 0.0);
        m_session->LoadParameter ("rhoInfPrimal",   m_rhoInfPrimal,  0.0);
        m_session->LoadParameter ("uInfPrimal",     m_uInfPrimal,    1.0);
        m_session->LoadParameter ("vInfPrimal",     m_vInfPrimal,    1.0);
        m_session->LoadParameter ("pInfPrimal",     m_pInfPrimal,    1.0);
        m_session->LoadParameter ("alphaInfPrimal", m_alphaInfDir,   0.0);
        m_session->LoadParameter ("Lref",           m_Lref,          1.0);
        m_session->LoadParameter ("alpha",          m_alpha,         0.0);
        m_session->LoadParameter ("Skappa",         m_Skappa,        -2.048);
        m_session->LoadParameter ("Kappa",          m_Kappa,         0.0);
        m_session->LoadParameter ("mu0",            m_mu0,           1.0);
        m_session->LoadParameter ("FL",             m_FacL,          0.0);
        m_session->LoadParameter ("FH",             m_FacH,          0.0);
        m_session->LoadParameter ("hFactor",        m_hFactor,       1.0);
        m_session->LoadParameter ("epsMax",         m_eps_max,       1.0);
        m_session->LoadParameter ("C1",             m_C1,            3.0);
        m_session->LoadParameter ("C2",             m_C2,            5.0);
        m_session->LoadSolverInfo("ShockCaptureType",
                                  m_shockCaptureType,    "Off");
        m_session->LoadParameter ("thermalConductivity",
                                  m_thermalConductivity, 0.0257);
        
        m_EqTypeStr = m_session->GetSolverInfo("EQTYPE");

        m_Cp      = m_gamma / (m_gamma - 1.0) * m_gasConstant;
        m_Prandtl = m_Cp * m_mu / m_thermalConductivity;
        
        const int nPoints = m_fields[0]->GetTotPoints();
        m_un = Array<OneD, Array<OneD, NekDouble> > (m_fields.num_elements());
        for (int i = 0; i < m_fields.num_elements(); ++i)
        {
            m_un[i] = Array<OneD, NekDouble> (nPoints, 0.0);
            Vmath::Vcopy(nPoints, m_fields[i]->GetPhys(), 1, m_un[i], 1);
        }
        // Type of advection class to be used
        switch(m_projectionType)
        {
                // Continuous field
            case MultiRegions::eGalerkin:
            {
                ASSERTL0(false, "Continuous field not supported.");
                break;
            }
                // Discontinuous field
            case MultiRegions::eDiscontinuous:
            {
                string advName, diffName, riemName;
                
                // Setting up advection and diffusion operators
                m_session->LoadSolverInfo("AdvectionType", advName, "WeakDG");
                m_session->LoadSolverInfo("DiffusionType", diffName, "LDGNS");
                
                m_advection = SolverUtils::GetAdvectionFactory()
                .CreateInstance(advName, advName);
                m_diffusion = SolverUtils::GetDiffusionFactory()
                .CreateInstance(diffName, diffName);
                
                if (m_specHP_dealiasing)
                {
                    m_advection->SetFluxVector(&CompressibleFlowSystem::
                                               GetFluxVectorDeAlias, this);
                    m_diffusion->SetFluxVectorNS(
                                                 &CompressibleFlowSystem::GetViscousFluxVectorDeAlias,
                                                 this);
                }
                else
                {
                    m_advection->SetFluxVector  (&CompressibleFlowSystem::
                                                 GetFluxVector, this);
                    m_diffusion->SetFluxVectorNS(&CompressibleFlowSystem::
                                                 GetViscousFluxVector, this);
                }
                
                if (m_shockCaptureType=="Smooth" && m_EqTypeStr=="EulerADCFE")
                {
                    m_advection->SetFluxVector(&CompressibleFlowSystem::
                                               GetFluxVector, this);
                    
                    m_diffusion->SetArtificialDiffusionVector(
                                                              &CompressibleFlowSystem::GetSmoothArtificialViscosity, this);
                }
                if (m_shockCaptureType=="NonSmooth" && m_EqTypeStr=="EulerADCFE")
                {
                    m_advection->SetFluxVector(&CompressibleFlowSystem::
                                               GetFluxVector, this);
                    
                    m_diffusion->SetArtificialDiffusionVector(
                                                              &CompressibleFlowSystem::GetArtificialDynamicViscosity, this);
                }

                // Setting up Riemann solver for advection operator
                m_session->LoadSolverInfo("UpwindType", riemName, "Average");
                m_riemannSolver = SolverUtils::GetRiemannSolverFactory()
                .CreateInstance(riemName);
                
                // Setting up upwind solver for diffusion operator
                m_riemannSolverLDG = SolverUtils::GetRiemannSolverFactory()
                .CreateInstance("UpwindLDG");
                
                m_riemannSolver->SetParam ("adjointSwitch",
                                           &CompressibleFlowSystem::GetAdjointSwitch,   this);
                
                // Setting up parameters for advection operator Riemann solver
                m_riemannSolver->SetParam (
                                           "gamma",  &CompressibleFlowSystem::GetGamma,   this);
                m_riemannSolver->SetAuxVec(
                                           "vecLocs", &CompressibleFlowSystem::GetVecLocs, this);
                m_riemannSolver->SetVector(
                                           "N",      &CompressibleFlowSystem::GetNormals, this);
                
                // Setting up parameters for diffusion operator Riemann solver
                m_riemannSolverLDG->SetParam (
                                              "gamma",  &CompressibleFlowSystem::GetGamma,   this);
                m_riemannSolverLDG->SetVector(
                                              "vecLocs", &CompressibleFlowSystem::GetVecLocs, this);
                m_riemannSolverLDG->SetVector(
                                              "N",      &CompressibleFlowSystem::GetNormals, this);
                
                if (m_adjointSwitch != 0.0)
                {
                    EquationSystem::InitialisePrimalSolution();
                    
                    m_advection->SetDirectSolution(&CompressibleFlowSystem::
                                                   GetDirectSolution, this);
                    
                    m_riemannSolver->SetFwdBwdDirectSolution(
                        &CompressibleFlowSystem::GetFwdBwdDirectSolution, this);
                    
                    
                    m_riemannSolver->SetFwdBwdDIFFDirectSolution(
                    &CompressibleFlowSystem::GetFwdBwdDIFFDirectSolution, this);
                    
                    //===============Set the appropriate fluxes=================
                    //                      Advection
                    //m_advection->SetFluxVector(&CompressibleFlowSystem::
                    //         GetAdjointConvFluxVectorFromPrimitiveVar, this);
                    
                    m_advection->SetFluxVector(
                            &CompressibleFlowSystem::GetAdjointFluxVector, this);
                    
                    m_advection->SetAddFluxVector(
    &CompressibleFlowSystem::GetAdjointAddConvFluxVectorFromPrimitiveVar, this);
                    
                    m_advection->SetJacTransposeDivVector(
            &CompressibleFlowSystem::GetDerivJacVectorFromPrimitiveVar, this);
                    
                    //m_advection->SetJacTransposeDivVector(&CompressibleFlowSystem::GetJacTransposeDivVector, this);
                    
                    m_advection->SetAddJacTransposeDivVector(&CompressibleFlowSystem::GetDerivAddJacVectorFromPrimitiveVar, this);
                    
                    //                      Diffusion
                    m_diffusion->SetFluxVectorNS(&CompressibleFlowSystem::
                        GetAdjointViscousFluxVectorFromPrimitiveVar, this);
                }
                
                
                // Concluding initialisation of advection / diffusion operators
                m_advection->SetRiemannSolver   (m_riemannSolver);
                m_diffusion->SetRiemannSolver   (m_riemannSolverLDG);
                m_advection->InitObject         (m_session, m_fields);
                m_diffusion->InitObject         (m_session, m_fields);
                break;

            }
            default:
            {
                ASSERTL0(false, "Unsupported projection type.");
                break;
            }
        }
    }
    
    /**
     * @brief Destructor for CompressibleFlowSystem class.
     */
    CompressibleFlowSystem::~CompressibleFlowSystem()
    {
        
    }
    
    /**
     * @brief Print out a summary with some relevant information.
     */
    void CompressibleFlowSystem::v_GenerateSummary(SolverUtils::SummaryList& s)
    {
        UnsteadySystem::v_GenerateSummary(s);
    }
    
    /**
     * @brief Wall boundary conditions for compressible flow problems.
     */
    void CompressibleFlowSystem::WallBC(
                            int                                   bcRegion,
                            int                                   cnt,
                            Array<OneD, Array<OneD, NekDouble> > &physarray)
    {
        int i;
        int nTracePts = GetTraceTotPoints();
        int nVariables = physarray.num_elements();
        
        const Array<OneD, const int> &traceBndMap
            = m_fields[0]->GetTraceBndMap();
        
        // Get physical values of the forward trace
        Array<OneD, Array<OneD, NekDouble> > Fwd(nVariables);
        Array<OneD, Array<OneD, NekDouble> > Fwdnew(nVariables);
        for (i = 0; i < nVariables; ++i)
        {
            Fwd[i] = Array<OneD, NekDouble>(nTracePts);
            Fwdnew[i] = Array<OneD, NekDouble>(nTracePts);
            m_fields[i]->ExtractTracePhys(physarray[i], Fwd[i]);
            m_fields[i]->ExtractTracePhys(physarray[i], Fwdnew[i]);
        }
        
        // Adjust the physical values of the trace to take
        // user defined boundaries into account
        int e, id1, id2, nBCEdgePts, eMax;
        
        eMax = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExpSize();
        
        for (e = 0; e < eMax; ++e)
        {
            nBCEdgePts = m_fields[0]->GetBndCondExpansions()[bcRegion]->
            GetExp(e)->GetTotPoints();
            id1 = m_fields[0]->GetBndCondExpansions()[bcRegion]->
            GetPhys_Offset(e);
            id2 = m_fields[0]->GetTrace()->GetPhys_Offset(traceBndMap[cnt+e]);

            // Boundary condition for epsilon term.
            if (nVariables == m_spacedim+3)
            {
                NekDouble factor  = 1.0;
                NekDouble factor2 = 1.0;
                
                Array<OneD, NekDouble > tmp2(nBCEdgePts, 0.0);
                Vmath::Smul(nBCEdgePts,
                            factor,
                            &Fwd[nVariables-1][id2], 1,
                            &tmp2[0], 1);
                
                Vmath::Vsub(nBCEdgePts,
                            &Fwd[nVariables-1][id2], 1,
                            &tmp2[0], 1,
                            &Fwd[nVariables-1][id2], 1);
                
                Vmath::Smul(nBCEdgePts,
                            factor2,
                            &Fwd[nVariables-1][id2], 1,
                            &Fwd[nVariables-1][id2], 1);
                
            }
            
            // For 2D/3D, define: v* = v - 2(v.n)n
            Array<OneD, NekDouble> tmp(nBCEdgePts, 0.0);
            Array<OneD, NekDouble> tmpfwd(nBCEdgePts, 0.0);
            Array<OneD, NekDouble> tmpfwdnew(nBCEdgePts, 0.0);
            // Calculate (v.n)
            for (i = 0; i < m_spacedim; ++i)
            {
                Vmath::Vvtvp(nBCEdgePts,
                             &Fwd[1+i][id2], 1,
                             &m_traceNormals[i][id2], 1,
                             &tmp[0], 1,
                             &tmp[0], 1);
            }
            
            // Calculate 2.0(v.n)
            Vmath::Smul(nBCEdgePts, -2.0, &tmp[0], 1, &tmp[0], 1);
            
            // Calculate v* = v - 2.0(v.n)n
            for (i = 0; i < m_spacedim; ++i)
            {
                Vmath::Vvtvp(nBCEdgePts,
                             &tmp[0], 1,
                             &m_traceNormals[i][id2], 1,
                             &Fwdnew[1+i][id2], 1,
                             &Fwdnew[1+i][id2], 1);
            }
            
            // Copy boundary adjusted values into the boundary expansion
            for (i = 0; i < nVariables; ++i)
            {
                Vmath::Vcopy(nBCEdgePts, &Fwdnew[i][id2], 1,
                             &(m_fields[i]->GetBndCondExpansions()[bcRegion]->
                               UpdatePhys())[id1], 1);
            }
        }
    }
    
    /**
     * @brief Wall boundary conditions for compressible flow problems.
     */
    void CompressibleFlowSystem::AdjointWallBC(
                int                                   bcRegion,
                int                                   cnt,
                Array<OneD, Array<OneD, NekDouble> > &physarray)
    {
        int i;
        int nTracePts = GetTraceTotPoints();
        int nVariables = physarray.num_elements();
        //int nq = physarray[0].num_elements();
        
        const Array<OneD, const int> &traceBndMap
        = m_fields[0]->GetTraceBndMap();
        
        // Get physical values of the forward trace
        Array<OneD, Array<OneD, NekDouble> > Fwd(nVariables);
        Array<OneD, Array<OneD, NekDouble> > Fwdnew(nVariables);
        
        Array<OneD, Array<OneD, NekDouble> > BwdDir(nVariables);
        Array<OneD, Array<OneD, NekDouble> > FwdDir(nVariables);
        
        for (i = 0; i < nVariables; ++i)
        {
            Fwd[i] = Array<OneD, NekDouble>(nTracePts);
            Fwdnew[i] = Array<OneD, NekDouble>(nTracePts);
            m_fields[i]->ExtractTracePhys(physarray[i], Fwd[i]);
            m_fields[i]->ExtractTracePhys(physarray[i], Fwdnew[i]);
        }
        
        /*
         NekDouble Cinf = 0.5 * m_rhoInfPrimal
         * (m_uInfPrimal*m_uInfPrimal+m_vInfPrimal*m_vInfPrimal) * m_Lref;
         */
        NekDouble norm_fac = 1.0;
        
        // Adjust the physical values of the trace to take
        // user defined boundaries into account
        
        int e, id1, id2, nBCEdgePts, eMax;
        
        eMax = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExpSize();
        
        for (e = 0; e < eMax; ++e)
        {
            nBCEdgePts = m_fields[0]->GetBndCondExpansions()[bcRegion]->
            GetExp(e)->GetTotPoints();
            id1 = m_fields[0]->GetBndCondExpansions()[bcRegion]->
            GetPhys_Offset(e);
            id2 = m_fields[0]->GetTrace()->GetPhys_Offset(traceBndMap[cnt+e]);
            
            // set z_rho = 0 at wall
            Vmath::Neg(nBCEdgePts,&Fwdnew[0][id2], 1);
            
            // set z_rhou = z_rhov = 0 at wall
            Vmath::Neg(nBCEdgePts, &Fwdnew[1][id2], 1);
            Vmath::Neg(nBCEdgePts, &Fwdnew[2][id2], 1);
            
            // set z_rhoE = 0 at wall
            Vmath::Neg(nBCEdgePts,&Fwdnew[3][id2], 1);
            
            // In case of lift or drag as cost function, add forcing
            if (m_Target == "Drag")
            {
                NekDouble Dx = norm_fac * cos (m_alphaInfDir);
                NekDouble Dy = norm_fac * sin (m_alphaInfDir);
                
                Array<OneD, Array<OneD, NekDouble> > DragDir(m_spacedim);
                DragDir[0] = Array<OneD, NekDouble> (nBCEdgePts, Dx);
                DragDir[1] = Array<OneD, NekDouble> (nBCEdgePts, Dy);
                
                Vmath::Vadd(nBCEdgePts, &Fwdnew[1][id2], 1, &DragDir[0][0], 1, &Fwdnew[1][id2], 1);
                Vmath::Vadd(nBCEdgePts, &Fwdnew[2][id2], 1, &DragDir[1][0], 1, &Fwdnew[2][id2], 1);
            }
            
            if (m_Target == "Lift")
            {
                NekDouble Lx = -norm_fac * sin (m_alphaInfDir);
                NekDouble Ly =  norm_fac * cos (m_alphaInfDir);
                
                Array<OneD, Array<OneD, NekDouble> > LiftDir(m_spacedim);
                LiftDir[0] = Array<OneD, NekDouble> (nBCEdgePts, Lx);
                LiftDir[1] = Array<OneD, NekDouble> (nBCEdgePts, Ly);
                
                Vmath::Vadd(nBCEdgePts, &Fwdnew[1][id2], 1, &LiftDir[0][0], 1, &Fwdnew[1][id2], 1);
                Vmath::Vadd(nBCEdgePts, &Fwdnew[2][id2], 1, &LiftDir[1][0], 1, &Fwdnew[2][id2], 1);
            }
            // Copy boundary adjusted values into the boundary expansion
            
            for (i = 0; i < nVariables; ++i)
            {
                Vmath::Vcopy(nBCEdgePts, &Fwdnew[i][id2], 1,
                             &(m_fields[i]->GetBndCondExpansions()[bcRegion]->
                               UpdatePhys())[id1], 1);
            }
        }
    }
    
    /**
     * @brief Wall boundary conditions for viscous compressible flow problems.
     */
    void CompressibleFlowSystem::WallViscousBC(
                            int                                   bcRegion,
                            int                                   cnt,
                            Array<OneD, Array<OneD, NekDouble> > &physarray)
    {
        int i;
        int nTracePts = GetTraceTotPoints();
        int nVariables = physarray.num_elements();
        
        const Array<OneD, const int> &traceBndMap
            = m_fields[0]->GetTraceBndMap();
        
        // Get physical values of the forward trace
        Array<OneD, Array<OneD, NekDouble> > Fwd(nVariables);
        for (i = 0; i < nVariables; ++i)
        {
            Fwd[i] = Array<OneD, NekDouble>(nTracePts);
            m_fields[i]->ExtractTracePhys(physarray[i], Fwd[i]);
        }

        // Take into account that for PDE based shock capturing, eps = 0 at the
        // wall. Adjust the physical values of the trace to take user defined
        // boundaries into account
        
        int e, id1, id2, nBCEdgePts, eMax;
        
        eMax = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExpSize();
        
        for (e = 0; e < eMax; ++e)
        {
            nBCEdgePts = m_fields[0]->GetBndCondExpansions()[bcRegion]->
            GetExp(e)->GetTotPoints();
            id1  = m_fields[0]->GetBndCondExpansions()[bcRegion]->
            GetPhys_Offset(e);
            id2 = m_fields[0]->GetTrace()->GetPhys_Offset(traceBndMap[cnt+e]);

            if (nVariables == m_spacedim+3)
            {
                NekDouble factor  = 0.0;
                NekDouble factor2 = 1.0;
                
                Array<OneD, NekDouble > tmp2(nBCEdgePts, 0.0);
                Vmath::Smul(nBCEdgePts,
                            factor,
                            &Fwd[nVariables-1][id2], 1,
                            &tmp2[0], 1);
                
                Vmath::Vsub(nBCEdgePts,
                            &Fwd[nVariables-1][id2], 1,
                            &tmp2[0], 1,
                            &Fwd[nVariables-1][id2], 1);
                
                Vmath::Smul(nBCEdgePts,
                            factor2,
                            &Fwd[nVariables-1][id2], 1,
                            &Fwd[nVariables-1][id2], 1);
            }
            
            for (i = 0; i < m_spacedim; i++)
            {
                Vmath::Neg(nBCEdgePts, &Fwd[i+1][id2], 1);
            }
            
            // Copy boundary adjusted values into the boundary expansion
            for (i = 0; i < nVariables; ++i)
            {
                Vmath::Vcopy(nBCEdgePts, &Fwd[i][id2], 1,
                             &(m_fields[i]->GetBndCondExpansions()[bcRegion]->
                               UpdatePhys())[id1], 1);
            }
        }
    }
    
    /**
     * @brief Symmetry boundary conditions for compressible flow problems.
     */
    void CompressibleFlowSystem::SymmetryBC(
                int                                      bcRegion,
                int                                      cnt,
                Array<OneD, Array<OneD, NekDouble> >    &physarray)
    {
        int i;
        int nTracePts = GetTraceTotPoints();
        int nVariables = physarray.num_elements();
        
        const Array<OneD, const int> &traceBndMap
            = m_fields[0]->GetTraceBndMap();
        
        // Get physical values of the forward trace (from exp to phys)
        Array<OneD, Array<OneD, NekDouble> > Fwd(nVariables);
        for (i = 0; i < nVariables; ++i)
        {
            Fwd[i] = Array<OneD, NekDouble>(nTracePts);
            m_fields[i]->ExtractTracePhys(physarray[i], Fwd[i]);
        }
        
        // Take into account that for PDE based shock capturing, eps = 0 at the
        // wall.
        int e, id1, id2, nBCEdgePts, eMax;
        
        eMax = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExpSize();
        
        for (e = 0; e < eMax; ++e)
        {
            nBCEdgePts = m_fields[0]->GetBndCondExpansions()[bcRegion]->
            GetExp(e)->GetTotPoints();
            id1  = m_fields[0]->GetBndCondExpansions()[bcRegion]->
            GetPhys_Offset(e);
            id2 = m_fields[0]->GetTrace()->GetPhys_Offset(traceBndMap[cnt+e]);
            
            if (nVariables == m_spacedim+3)
            {
                NekDouble factor  = 0.0;
                NekDouble factor2 = 1.0;
                
                Array<OneD, NekDouble > tmp2(nBCEdgePts, 0.0);
                Vmath::Smul(nBCEdgePts,
                            factor,
                            &Fwd[nVariables-1][id2], 1,
                            &tmp2[0], 1);
                
                Vmath::Vsub(nBCEdgePts,
                            &Fwd[nVariables-1][id2], 1,
                            &tmp2[0], 1,
                            &Fwd[nVariables-1][id2], 1);
                
                Vmath::Smul(nBCEdgePts,
                            factor2,
                            &Fwd[nVariables-1][id2], 1,
                            &Fwd[nVariables-1][id2], 1);
            }
            
            // For 2D/3D, define: v* = v - 2(v.n)n
            Array<OneD, NekDouble> tmp(nBCEdgePts, 0.0);
            
            // Calculate (v.n)
            for (i = 0; i < m_spacedim; ++i)
            {
                Vmath::Vvtvp(nBCEdgePts,
                             &Fwd[1+i][id2], 1,
                             &m_traceNormals[i][id2], 1,
                             &tmp[0], 1,
                             &tmp[0], 1);
            }
            
            // Calculate 2.0(v.n)
            Vmath::Smul(nBCEdgePts, -2.0, &tmp[0], 1, &tmp[0], 1);
            
            // Calculate v* = v - 2.0(v.n)n
            for (i = 0; i < m_spacedim; ++i)
            {
                Vmath::Vvtvp(nBCEdgePts,
                             &tmp[0], 1,
                             &m_traceNormals[i][id2], 1,
                             &Fwd[1+i][id2], 1,
                             &Fwd[1+i][id2], 1);
            }
            
            // Copy boundary adjusted values into the boundary expansion
            for (i = 0; i < nVariables; ++i)
            {
                Vmath::Vcopy(nBCEdgePts, &Fwd[i][id2], 1,
                             &(m_fields[i]->GetBndCondExpansions()[bcRegion]->
                               UpdatePhys())[id1], 1);
            }
        }
    }
    
    /**
     * @brief Outflow characteristic boundary conditions for compressible
     * flow problems.
     */
    void CompressibleFlowSystem::RiemannInvariantBC(
                            int                                   bcRegion,
                            int                                   cnt,
                            Array<OneD, Array<OneD, NekDouble> > &physarray)
    {
        int i, j;
        int nTracePts = GetTraceTotPoints();
        int nVariables = physarray.num_elements();
        int nDimensions = m_spacedim;
        
        const Array<OneD, const int> &traceBndMap
            = m_fields[0]->GetTraceBndMap();
        
        NekDouble gamma            = m_gamma;
        NekDouble gammaInv         = 1.0 / gamma;
        NekDouble gammaMinusOne    = gamma - 1.0;
        NekDouble gammaMinusOneInv = 1.0 / gammaMinusOne;
        
        Array<OneD, NekDouble> tmp1 (nTracePts, 0.0);
        Array<OneD, NekDouble> tmp2 (nTracePts, 0.0);
        Array<OneD, NekDouble> VnInf(nTracePts, 0.0);
        Array<OneD, NekDouble> velInf(nDimensions, 0.0);
        
        // Computing the normal velocity for characteristics coming
        // from outside the computational domain
        velInf[0] = m_uInf;
        Vmath::Smul(nTracePts, m_uInf, m_traceNormals[0], 1, VnInf, 1);
        if (nDimensions == 2 || nDimensions == 3)
        {
            velInf[1] = m_vInf;
            Vmath::Smul(nTracePts, m_vInf, m_traceNormals[1], 1, tmp1, 1);
            Vmath::Vadd(nTracePts, VnInf, 1, tmp1, 1, VnInf, 1);
        }
        if (nDimensions == 3)
        {
            velInf[2] = m_wInf;
            Vmath::Smul(nTracePts, m_wInf, m_traceNormals[2], 1, tmp2, 1);
            Vmath::Vadd(nTracePts, VnInf, 1, tmp2, 1, VnInf, 1);
        }
        
        // Get physical values of the forward trace
        Array<OneD, Array<OneD, NekDouble> > Fwd(nVariables);
        for (i = 0; i < nVariables; ++i)
        {
            Fwd[i] = Array<OneD, NekDouble>(nTracePts);
            m_fields[i]->ExtractTracePhys(physarray[i], Fwd[i]);
        }
        
        // Computing the normal velocity for characteristics coming
        // from inside the computational domain
        Array<OneD, NekDouble > Vn (nTracePts, 0.0);
        Array<OneD, NekDouble > Vel(nTracePts, 0.0);
        for (i = 0; i < nDimensions; ++i)
        {
            Vmath::Vdiv(nTracePts, Fwd[i+1], 1, Fwd[0], 1, Vel, 1);
            Vmath::Vvtvp(nTracePts, m_traceNormals[i], 1, Vel, 1, Vn, 1, Vn, 1);
        }
        
        // Computing the absolute value of the velocity in order to compute the
        // Mach number to decide whether supersonic or subsonic
        Array<OneD, NekDouble > absVel(nTracePts, 0.0);
        for (i = 0; i < nDimensions; ++i)
        {
            Vmath::Vdiv(nTracePts, Fwd[i+1], 1, Fwd[0], 1, tmp1, 1);
            Vmath::Vmul(nTracePts, tmp1, 1, tmp1, 1, tmp1, 1);
            Vmath::Vadd(nTracePts, tmp1, 1, absVel, 1, absVel, 1);
        }
        Vmath::Vsqrt(nTracePts, absVel, 1, absVel, 1);
        
        // Get speed of sound
        Array<OneD, NekDouble > soundSpeed(nTracePts);
        Array<OneD, NekDouble > pressure  (nTracePts);
        
        for (i = 0; i < nTracePts; i++)
        {
            if (m_spacedim == 1)
            {
                pressure[i] = (gammaMinusOne) * (Fwd[2][i] -
                                                 0.5 * (Fwd[1][i] * Fwd[1][i] / Fwd[0][i]));
            }
            else if (m_spacedim == 2)
            {
                pressure[i] = (gammaMinusOne) * (Fwd[3][i] -
                                                 0.5 * (Fwd[1][i] * Fwd[1][i] / Fwd[0][i] +
                                                        Fwd[2][i] * Fwd[2][i] / Fwd[0][i]));
            }
            else
            {
                pressure[i] = (gammaMinusOne) * (Fwd[4][i] -
                                                 0.5 * (Fwd[1][i] * Fwd[1][i] / Fwd[0][i] +
                                                        Fwd[2][i] * Fwd[2][i] / Fwd[0][i] +
                                                        Fwd[3][i] * Fwd[3][i] / Fwd[0][i]));
            }
            
            soundSpeed[i] = sqrt(gamma * pressure[i] / Fwd[0][i]);
        }
        
        // Get Mach
        Array<OneD, NekDouble > Mach(nTracePts, 0.0);
        Vmath::Vdiv(nTracePts, Vn, 1, soundSpeed, 1, Mach, 1);
        Vmath::Vabs(nTracePts, Mach, 1, Mach, 1);
        
        // Auxiliary variables
        int eMax;
        int e, id1, id2, nBCEdgePts, pnt;
        NekDouble cPlus, rPlus, cMinus, rMinus, VDBC, VNBC;
        Array<OneD, NekDouble> velBC(nDimensions, 0.0);
        Array<OneD, NekDouble> rhoVelBC(nDimensions, 0.0);
        NekDouble rhoBC, EBC, cBC, sBC, pBC;
        
        eMax = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExpSize();
        
        // Loop on bcRegions
        for (e = 0; e < eMax; ++e)
        {
            nBCEdgePts = m_fields[0]->GetBndCondExpansions()[bcRegion]->
            GetExp(e)->GetNumPoints(0);
            
            id1 = m_fields[0]->GetBndCondExpansions()[bcRegion]->
            GetPhys_Offset(e);
            id2 = m_fields[0]->GetTrace()->GetPhys_Offset(traceBndMap[cnt+e]);
            
            // Loop on the points of the bcRegion
            for (i = 0; i < nBCEdgePts; i++)
            {
                pnt = id2+i;
                
                // Impose inflow Riemann invariant
                if (Vn[pnt] <= 0.0)
                {
                    // Subsonic flows
                    if (Mach[pnt] < 1.00)
                    {
                        // + Characteristic from inside
                        cPlus = sqrt(gamma * pressure[pnt] / Fwd[0][pnt]);
                        rPlus = Vn[pnt] + 2.0 * cPlus * gammaMinusOneInv;
                        
                        // - Characteristic from boundary
                        cMinus = sqrt(gamma * m_pInf / m_rhoInf);
                        rMinus = VnInf[pnt] - 2.0 * cMinus * gammaMinusOneInv;
                    }
                    else
                    {
                        // + Characteristic from inside
                        cPlus = sqrt(gamma * m_pInf / m_rhoInf);
                        rPlus = VnInf[pnt] + 2.0 * cPlus * gammaMinusOneInv;
                        
                        // + Characteristic from inside
                        cMinus = sqrt(gamma * m_pInf / m_rhoInf);
                        rMinus = VnInf[pnt] - 2.0 * cPlus * gammaMinusOneInv;
                    }
                    
                    // Riemann boundary variables
                    VNBC = 0.5 * (rPlus + rMinus);
                    cBC = 0.25 * gammaMinusOne * (rPlus - rMinus);
                    VDBC = VNBC - VnInf[pnt];
                    
                    // Thermodynamic boundary variables
                    sBC = m_pInf / (pow(m_rhoInf, gamma));
                    rhoBC = pow((cBC * cBC) / (gamma * sBC), gammaMinusOneInv);
                    pBC = rhoBC * cBC * cBC * gammaInv;
                    
                    // Kinetic energy initialiasation
                    NekDouble EkBC = 0.0;
                    
                    // Boundary velocities
                    for ( j = 0; j < nDimensions; ++j)
                    {
                        velBC[j] = velInf[j] + VDBC * m_traceNormals[j][pnt];
                        rhoVelBC[j] = rhoBC * velBC[j];
                        EkBC += 0.5 * rhoBC * velBC[j]*velBC[j];
                    }
                    
                    // Boundary energy
                    EBC = pBC * gammaMinusOneInv + EkBC;
                    
                    // Imposing Riemann Invariant boundary conditions
                    (m_fields[0]->GetBndCondExpansions()[bcRegion]->
                     UpdatePhys())[id1+i] = rhoBC;
                    for (j = 0; j < nDimensions; ++j)
                    {
                        (m_fields[j+1]->GetBndCondExpansions()[bcRegion]->
                         UpdatePhys())[id1+i] = rhoVelBC[j];
                    }
                    (m_fields[nDimensions+1]->GetBndCondExpansions()[bcRegion]->
                     UpdatePhys())[id1+i] = EBC;
                    
                }
                else // Impose outflow Riemann invariant
                {
                    // Subsonic flows
                    if (Mach[pnt] < 1.00)
                    {
                        // + Characteristic from inside
                        cPlus = sqrt(gamma * pressure[pnt] / Fwd[0][pnt]);
                        rPlus = Vn[pnt] + 2.0 * cPlus * gammaMinusOneInv;
                        
                        // - Characteristic from boundary
                        cMinus = sqrt(gamma * m_pInf / m_rhoInf);
                        rMinus = VnInf[pnt] - 2.0 * cMinus * gammaMinusOneInv;
                    }
                    else
                    {
                        // + Characteristic from inside
                        cPlus = sqrt(gamma * pressure[pnt] / Fwd[0][pnt]);
                        rPlus = Vn[pnt] + 2.0 * cPlus * gammaMinusOneInv;
                        
                        // + Characteristic from inside
                        cMinus = sqrt(gamma * pressure[pnt] / Fwd[0][pnt]);
                        rMinus = Vn[pnt] - 2.0 * cPlus * gammaMinusOneInv;
                    }
                    
                    // Riemann boundary variables
                    VNBC = 0.5 * (rPlus + rMinus);
                    cBC = 0.25 * gammaMinusOne * (rPlus - rMinus);
                    VDBC = VNBC - Vn[pnt];
                    
                    // Thermodynamic boundary variables
                    sBC = pressure[pnt] / (pow(Fwd[0][pnt], gamma));
                    rhoBC = pow((cBC * cBC) / (gamma * sBC), gammaMinusOneInv);
                    pBC = rhoBC * cBC * cBC * gammaInv;
                    
                    // Kinetic energy initialiasation
                    NekDouble EkBC = 0.0;
                    
                    // Boundary velocities
                    for ( j = 0; j < nDimensions; ++j)
                    {
                        velBC[j] = Fwd[j+1][pnt] / Fwd[0][pnt] +
                        VDBC * m_traceNormals[j][pnt];
                        rhoVelBC[j] = rhoBC * velBC[j];
                        EkBC += 0.5 * rhoBC * velBC[j]*velBC[j];
                    }
                    
                    // Boundary energy
                    EBC = pBC * gammaMinusOneInv + EkBC;
                    
                    // Imposing Riemann Invariant boundary conditions
                    (m_fields[0]->GetBndCondExpansions()[bcRegion]->
                     UpdatePhys())[id1+i] = rhoBC;
                    for (j = 0; j < nDimensions; ++j)
                    {
                        (m_fields[j+1]->GetBndCondExpansions()[bcRegion]->
                         UpdatePhys())[id1+i] = rhoVelBC[j];
                    }
                    (m_fields[nDimensions+1]->GetBndCondExpansions()[bcRegion]->
                     UpdatePhys())[id1+i] = EBC;
                }
            }
        }
    }
    
    /**
     * @brief Extrapolation of order 0 for all the variables such that,
     * at the boundaries, a trivial Riemann problem is solved.
     */
    void CompressibleFlowSystem::ExtrapOrder0BC(
                                int                                   bcRegion,
                                int                                   cnt,
                                Array<OneD, Array<OneD, NekDouble> > &physarray)
    {
        int i, j;
        int e, pnt;
        int id1, id2, nBCEdgePts;
        int nTracePts = GetTraceTotPoints();
        int nVariables = physarray.num_elements();
        int nDimensions = m_spacedim;
        
        const Array<OneD, const int> &traceBndMap
            = m_fields[0]->GetTraceBndMap();
    
        // Get physical values of the forward trace
        Array<OneD, Array<OneD, NekDouble> > Fwd(nVariables);
        for (i = 0; i < nVariables; ++i)
        {
            Fwd[i] = Array<OneD, NekDouble>(nTracePts);
            m_fields[i]->ExtractTracePhys(physarray[i], Fwd[i]);
        }
        
        int eMax;
        
        eMax = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExpSize();
        
        // Loop on bcRegions
        for (e = 0; e < eMax; ++e)
        {
            nBCEdgePts = m_fields[0]->GetBndCondExpansions()[bcRegion]->
            GetExp(e)->GetNumPoints(0);
            id1 = m_fields[0]->GetBndCondExpansions()[bcRegion]->
            GetPhys_Offset(e) ;
            id2 = m_fields[0]->GetTrace()->GetPhys_Offset(traceBndMap[cnt+e]);
            
            // Loop on points of bcRegion 'e'
            for (i = 0; i < nBCEdgePts; i++)
            {
                pnt = id2+i;
                
                // Setting up bcs for density
                (m_fields[0]->GetBndCondExpansions()[bcRegion]->
                 UpdatePhys())[id1+i] = Fwd[0][pnt];
                
                // Setting up bcs for velocities
                for (j = 1; j <=nDimensions; ++j)
                {
                    (m_fields[j]->GetBndCondExpansions()[bcRegion]->
                     UpdatePhys())[id1+i] = Fwd[j][pnt];
                }
                
                // Setting up bcs for energy
                (m_fields[nVariables-1]->GetBndCondExpansions()[bcRegion]->
                 UpdatePhys())[id1+i] = Fwd[nVariables-1][pnt];
            }
        }
    }
    
    /**
     * @brief Return the flux vector for the compressible Euler equations.
     *
     * @param i           Component of the flux vector to calculate.
     * @param physfield   Fields.
     * @param flux        Resulting flux.
     */
    void CompressibleFlowSystem::GetFluxVector(
            const Array<OneD, Array<OneD, NekDouble> >               &physfield,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &flux)
    {
        int i, j;
        int nq = physfield[0].num_elements();
        int nVariables = m_fields.num_elements();
        
        Array<OneD, NekDouble> pressure(nq);
        Array<OneD, Array<OneD, NekDouble> > velocity(m_spacedim);
        
        // Flux vector for the rho equation
        for (i = 0; i < m_spacedim; ++i)
        {
            velocity[i] = Array<OneD, NekDouble>(nq);
            Vmath::Vcopy(nq, physfield[i+1], 1, flux[0][i], 1);
        }
        
        GetVelocityVector(physfield, velocity);
        GetPressure      (physfield, velocity, pressure);
        
        // Flux vector for the velocity fields
        for (i = 0; i < m_spacedim; ++i)
        {
            for (j = 0; j < m_spacedim; ++j)
            {
                Vmath::Vmul(nq, velocity[j], 1, physfield[i+1], 1,
                            flux[i+1][j], 1);
            }
            
            // Add pressure to appropriate field
            Vmath::Vadd(nq, flux[i+1][i], 1, pressure, 1, flux[i+1][i], 1);
        }
        
        // Flux vector for energy.
        Vmath::Vadd(nq, physfield[m_spacedim+1], 1, pressure, 1,
                    pressure, 1);
        
        for (j = 0; j < m_spacedim; ++j)
        {
            Vmath::Vmul(nq, velocity[j], 1, pressure, 1,
                        flux[m_spacedim+1][j], 1);
        }
        
        // For the smooth viscosity model
        if (nVariables == m_spacedim+3)
        {
            // Add a zero row for the advective fluxes
            for (j = 0; j < m_spacedim; ++j)
            {
                Vmath::Zero(nq, flux[m_spacedim+2][j], 1);
            }
        }
    }
    
    void CompressibleFlowSystem::GetJacTransposeDivVector(
          const Array<OneD, Array<OneD, NekDouble> > &inarray,
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &outarray)
    {
        int i, p, o;
        int nq = inarray[0].num_elements();
        int nvar = m_fields.num_elements();
        
        Array<OneD, NekDouble> pressure(nq,0.0);
        Array<OneD, NekDouble> H(nq,0.0);
        Array<OneD, NekDouble> Htmp(nq,0.0);
        Array<OneD, NekDouble> U2(nq,0.0);
        Array<OneD, Array<OneD, NekDouble> > vel(m_spacedim);
        Array<OneD, Array<OneD, NekDouble> > direct_fields(nvar);
        
        Array<OneD, NekDouble> dH_dx(nq, 0.0);
        Array<OneD, NekDouble> dH_dy(nq, 0.0);
        Array<OneD, NekDouble> dH_dz(nq, 0.0);
        
        Array<OneD, NekDouble> du_dx(nq, 0.0);
        Array<OneD, NekDouble> du_dy(nq, 0.0);
        Array<OneD, NekDouble> du_dz(nq, 0.0);
        
        Array<OneD, NekDouble> dv_dx(nq, 0.0);
        Array<OneD, NekDouble> dv_dy(nq, 0.0);
        Array<OneD, NekDouble> dv_dz(nq, 0.0);
        
        Array<OneD, NekDouble> dw_dx(nq, 0.0);
        Array<OneD, NekDouble> dw_dy(nq, 0.0);
        Array<OneD, NekDouble> dw_dz(nq, 0.0);
        
        Array<OneD, NekDouble> dU2_dx(nq, 0.0);
        Array<OneD, NekDouble> dU2_dy(nq, 0.0);
        Array<OneD, NekDouble> dU2_dz(nq, 0.0);
        
        Array<OneD, NekDouble> duv_dx(nq, 0.0);
        Array<OneD, NekDouble> duv_dy(nq, 0.0);
        Array<OneD, NekDouble> duv_dz(nq, 0.0);
        
        Array<OneD, NekDouble> duU2_dx(nq, 0.0);
        Array<OneD, NekDouble> dvU2_dy(nq, 0.0);
        Array<OneD, NekDouble> dwU2_dz(nq, 0.0);
        
        Array<OneD, NekDouble> duH_dx(nq, 0.0);
        Array<OneD, NekDouble> dvH_dy(nq, 0.0);
        Array<OneD, NekDouble> dwH_dz(nq, 0.0);
        
        Array<OneD, NekDouble> du2_dx(nq, 0.0);
        Array<OneD, NekDouble> du2_dy(nq, 0.0);
        Array<OneD, NekDouble> du2_dz(nq, 0.0);
        
        Array<OneD, NekDouble> dv2_dx(nq, 0.0);
        Array<OneD, NekDouble> dv2_dy(nq, 0.0);
        Array<OneD, NekDouble> dv2_dz(nq, 0.0);
        
        Array<OneD, NekDouble> dw2_dx(nq, 0.0);
        Array<OneD, NekDouble> dw2_dy(nq, 0.0);
        Array<OneD, NekDouble> dw2_dz(nq, 0.0);
        
        NekDouble GammaMinOne     = m_gamma - 1;
        NekDouble HalfGammaMinOne = 0.5*(m_gamma - 1);
        NekDouble ThreeMinGam     = 3.0 - m_gamma;
        
        for (i = 0; i < nvar; ++i)
        {
            direct_fields[i] = Array<OneD, NekDouble>(nq, 0.0);
        }
        
        GetDirectSolution(direct_fields);
        
        for (i = 0; i < m_spacedim; ++i)
        {
            vel[i] = Array<OneD, NekDouble>(nq, 0.0);
            Vmath::Vdiv(nq, direct_fields[i+1], 1, direct_fields[0], 1, vel[i], 1);
            Vmath::Vvtvp(nq, &vel[i][0], 1, &vel[i][0], 1, &U2[0], 1, &U2[0], 1);
        }
        
        GetPressure(direct_fields, vel, pressure);
        
        // determine H = E + p/rho;
        Vmath::Vdiv(nq, direct_fields[nvar-1], 1, direct_fields[0], 1, H, 1);
        
        Vmath::Vdiv(nq, pressure, 1, direct_fields[0], 1, Htmp, 1);
        
        Vmath::Vadd(nq, H, 1, Htmp, 1, H, 1);
        
        if (m_spacedim == 1)
        {
            ASSERTL0(false, "1D adjoint solver is not yet implemented");
        }
        if (m_spacedim == 2)
        {
            m_primal[0]->PhysDeriv(vel[0], du_dx, du_dy);
            m_primal[0]->PhysDeriv(vel[1], dv_dx, dv_dy);
            m_primal[0]->PhysDeriv(H, dH_dx, dH_dy);
            
            // create du^2/dx = 2u*du/dx
            Vmath::Vmul(nq, &vel[0][0], 1, &du_dx[0], 1, &du2_dx[0], 1);
            Vmath::Smul(nq, 2.0, &du2_dx[0], 1, &du2_dx[0], 1);
            
            // create dv^2/dx = 2v*dv/dx
            Vmath::Vmul(nq, &vel[1][0], 1, &dv_dx[0], 1, &dv2_dx[0], 1);
            Vmath::Smul(nq, 2.0, &dv2_dx[0], 1, &dv2_dx[0], 1);
            
            // create du^2/dy = 2u*du/dy
            Vmath::Vmul(nq, &vel[0][0], 1, &du_dy[0], 1, &du2_dy[0], 1);
            Vmath::Smul(nq, 2.0, &du2_dy[0], 1, &du2_dy[0], 1);
            
            // create dv^2/dy = 2v*dv/dy
            Vmath::Vmul(nq, &vel[1][0], 1, &dv_dy[0], 1, &dv2_dy[0], 1);
            Vmath::Smul(nq, 2.0, &dv2_dy[0], 1, &dv2_dy[0], 1);
            
            // create dU^2/dx = du^2/dx + dv^2/dx = 2udu/dx+2vdv/dx
            Vmath::Vadd(nq, &du2_dx[0], 1, &dv2_dx[0], 1, &dU2_dx[0], 1);
            
            // create dU^2/dy = du^2/dy + dv^2/dy
            Vmath::Vadd(nq, &du2_dy[0], 1, &dv2_dy[0], 1, &dU2_dy[0], 1);
            
            // d(uU^2)/dx = d(u^3 + u*v^2)/dx = du^3/dx + d(uv^2)/dx = du^3/dx + d(u*v*v)/dx
            //  = du^3/dx + 2*u*v*dv/dx + v^2*du/dx
            
            // d(uU^2)/dx = 2*u*v*dv/dx + v^2du/dx + 3u^2du/dx
            
            Array<OneD, NekDouble> tmp_dir(nq, 0.0);
            Array<OneD, NekDouble> tmp2_dir(nq, 0.0);
            
            Vmath::Vmul(nq, &vel[0][0], 1, &vel[1][0], 1, &duU2_dx[0], 1);
            Vmath::Smul(nq, 2.0, &duU2_dx[0], 1, &duU2_dx[0], 1);
            Vmath::Vmul(nq, &duU2_dx[0], 1, &dv_dx[0], 1, &duU2_dx[0], 1);
            
            Vmath::Vmul(nq, &vel[0][0], 1, &vel[0][0], 1, &tmp_dir[0], 1);
            Vmath::Smul(nq, 3.0, &tmp_dir[0], 1, &tmp_dir[0], 1);
            Vmath::Vmul(nq, &tmp_dir[0], 1, &du_dx[0], 1, &tmp_dir[0], 1);
            
            Vmath::Vmul(nq, &vel[1][0], 1, &vel[1][0], 1, &tmp2_dir[0], 1);
            Vmath::Vmul(nq, &tmp2_dir[0], 1, &du_dx[0], 1, &tmp2_dir[0], 1);
            
            Vmath::Vadd(nq, &tmp_dir[0], 1, &duU2_dx[0], 1, &duU2_dx[0], 1);
            Vmath::Vadd(nq, &tmp2_dir[0], 1, &duU2_dx[0], 1, &duU2_dx[0], 1);
            
            Vmath::Zero(nq, &tmp_dir[0], 1);
            Vmath::Zero(nq, &tmp2_dir[0], 1);
            
            // d(vU^2)/dy = d(vu^2 + v^3)/dy = dvu^2/dy + d(v^3)/dy
            // = vdu^2/dy+u^2dv/dy + dv^3/dy
            
            // d(vU^2)/dy = 2*u*v*du/dy + u^2dv/dy + 3v^2dv/dy
            
            Vmath::Vmul(nq, &vel[0][0], 1, &vel[1][0], 1, &dvU2_dy[0], 1);
            Vmath::Smul(nq, 2.0, &dvU2_dy[0], 1, &dvU2_dy[0], 1);
            Vmath::Vmul(nq, &dvU2_dy[0], 1, &du_dy[0], 1, &dvU2_dy[0], 1);
            
            Vmath::Vmul(nq, &vel[1][0], 1, &vel[1][0], 1, &tmp_dir[0], 1);
            Vmath::Smul(nq, 3.0, &tmp_dir[0], 1, &tmp_dir[0], 1);
            Vmath::Vmul(nq, &tmp_dir[0], 1, &dv_dy[0], 1, &tmp_dir[0], 1);
            
            Vmath::Vmul(nq, &vel[0][0], 1, &vel[0][0], 1, &tmp2_dir[0], 1);
            Vmath::Vmul(nq, &tmp2_dir[0], 1, &dv_dy[0], 1, &tmp2_dir[0], 1);
            
            Vmath::Vadd(nq, &tmp_dir[0], 1, &dvU2_dy[0], 1, &dvU2_dy[0], 1);
            Vmath::Vadd(nq, &tmp2_dir[0], 1, &dvU2_dy[0], 1, &dvU2_dy[0], 1);
            
            Vmath::Zero(nq, &tmp_dir[0], 1);
            Vmath::Zero(nq, &tmp2_dir[0], 1);
            // create d(vU^2)/dy = vdU^2/dy + U^2dv/dy
            
            // create d(uH)/dx = udH/dx + Hdu/dx
            Vmath::Vmul(nq, &vel[0][0], 1, &dH_dx[0], 1, &duH_dx[0], 1);
            Vmath::Vvtvp(nq, &H[0], 1, &du_dx[0], 1, &duH_dx[0], 1, &duH_dx[0], 1);
            
            // create d(vH)/dy = vdH/dy + Hdv/dy
            Vmath::Vmul(nq, &vel[1][0], 1, &dH_dy[0], 1, &dvH_dy[0], 1);
            Vmath::Vvtvp(nq, &H[0], 1, &dv_dy[0], 1, &dvH_dy[0], 1, &dvH_dy[0], 1);
            
            // create d(uv)/dx = udv/dx + vdu/dx
            Vmath::Vmul(nq, &vel[0][0], 1, &dv_dx[0], 1, &duv_dx[0], 1);
            Vmath::Vvtvp(nq, &vel[1][0], 1, &du_dx[0], 1, &duv_dx[0], 1, &duv_dx[0], 1);
            
            // create d(uv)/dy = udv/dy + vdu/dy
            Vmath::Vmul(nq, &vel[0][0], 1, &dv_dy[0], 1, &duv_dy[0], 1);
            Vmath::Vvtvp(nq, &vel[1][0], 1, &du_dy[0], 1, &duv_dy[0], 1, &duv_dy[0], 1);
            
            // =================================================================
            //                   Create JacTransVectorX
            // =================================================================
            // ===========================Eq1===================================
            
            Array<OneD, NekDouble> tmp(nq,  0.0);
            Array<OneD, NekDouble> tmp1(nq, 0.0);
            Array<OneD, NekDouble> tmp2(nq, 0.0);
            Array<OneD, NekDouble> tmp3(nq, 0.0);
            
            // tmp = du^2_dx
            Vmath::Vcopy(nq, &du2_dx[0], 1, &tmp[0], 1);
            
            // tmp1 = (gamma-1)/2 * dU^2_dx
            Vmath::Smul(nq, HalfGammaMinOne, &dU2_dx[0], 1, &tmp1[0], 1);
            
            // tmp = (gamma-1)/2 * dU^2_dx - du^2_dx
            Vmath::Vsub(nq, &tmp1[0], 1, &tmp[0], 1, &tmp[0], 1);
            
            // tmp = z_rhou*((gamma-1)/2 * dU^2_dx - du^2_dx)
            Vmath::Vmul(nq, &inarray[1][0], 1, &tmp[0], 1, &tmp[0], 1);
            
            //------------------------------------------------------------------
            
            // tmp2 = u*v;
            Vmath::Vcopy(nq, &duv_dx[0], 1, &tmp2[0], 1);
            
            // tmp2 = -u*v;
            Vmath::Neg(nq, &tmp2[0], 1);
            
            // tmp2 = -u*v*z_rhov;
            Vmath::Vmul(nq, &inarray[2][0], 1, &tmp2[0], 1, &tmp2[0], 1);
            
            //------------------------------------------------------------------
            // tmp3 = (gamma-1)/2*duU2_dx;
            Vmath::Smul(nq, HalfGammaMinOne, &duU2_dx[0], 1, &tmp3[0], 1);
            
            // tmp3 = (gamma-1)/2*duU2_dx-duH_dx;
            Vmath::Vsub(nq, &tmp3[0], 1, &duH_dx[0], 1, &tmp3[0], 1);
            
            // tmp3 = z_rhoE*((gamma-1)/2*duU2_dx-duH_dx);
            Vmath::Vmul(nq, &inarray[3][0], 1, &tmp3[0], 1,  &tmp3[0], 1);
            //------------------------------------------------------------------
            
            Vmath::Vadd(nq, &tmp[0], 1, &tmp2[0], 1, &outarray[0][0][0], 1);
            
            Vmath::Vadd(nq, &outarray[0][0][0], 1, &tmp3[0], 1, &outarray[0][0][0], 1);
            
            // ===========================Eq2===================================
            
            Vmath::Zero(nq, &tmp[0], 1);
            Vmath::Zero(nq, &tmp1[0], 1);
            Vmath::Zero(nq, &tmp2[0], 1);
            Vmath::Zero(nq, &tmp3[0], 1);
            
            // tmp = zrho;
            //Vmath::Vcopy(nq, &inarray[0][0], 1, &tmp[0], 1);
            
            // tmp1 = (3-gamma) * du_dx
            Vmath::Smul(nq,ThreeMinGam, &du_dx[0], 1, &tmp1[0], 1);
            
            // tmp1 = (3-gamma) * du_dx * z_rhou
            Vmath::Vmul(nq, &tmp1[0], 1, &inarray[1][0], 1, &tmp1[0], 1);
            
            // tmp2 = dv_dx
            Vmath::Vcopy(nq, &dv_dx[0], 1, &tmp2[0], 1);
            
            // tmp2 = dv_dx * z_rhov
            Vmath::Vmul(nq, &tmp2[0], 1, &inarray[2][0], 1, &tmp2[0], 1);
            
            // tmp3 = (1-gamma) * du^2_dx
            Vmath::Smul(nq, -GammaMinOne, &du2_dx[0], 1, &tmp3[0], 1);
            
            // tmp3 = (1-gamma) * du^2_dx+dH_dx
            Vmath::Vadd(nq, &tmp3[0], 1, &dH_dx[0], 1, &tmp3[0], 1);
            
            // tmp3 = ((1-gamma) * du^2_dx+dH_dx)*z_rhoE
            Vmath::Vmul(nq, &inarray[3][0], 1, &tmp3[0], 1, &tmp3[0], 1);
            
            Vmath::Vcopy(nq, &tmp1[0], 1, &outarray[1][0][0], 1);
            
            Vmath::Vadd(nq, &outarray[1][0][0], 1, &tmp2[0], 1, &outarray[1][0][0], 1);
            
            Vmath::Vadd(nq, &outarray[1][0][0], 1, &tmp3[0], 1, &outarray[1][0][0], 1);
            
            // ===========================Eq3===================================
            
            Vmath::Zero(nq, &tmp[0], 1);
            Vmath::Zero(nq, &tmp1[0], 1);
            Vmath::Zero(nq, &tmp2[0], 1);
            Vmath::Zero(nq, &tmp3[0], 1);
            
            // tmp = -(gamma-1)*dv_dx
            Vmath::Smul(nq, -GammaMinOne, &dv_dx[0], 1, &tmp[0], 1);
            
            // tmp = -(gamma-1)*dv_dx * z_rhou
            Vmath::Vmul(nq, &tmp[0], 1, &inarray[1][0], 1, &tmp[0], 1);
            
            // tmp1 = du_dx * z_rhov
            Vmath::Vmul(nq, &du_dx[0], 1, &inarray[2][0], 1, &tmp1[0], 1);
            
            // tmp2 =  duv_dx
            Vmath::Vcopy(nq, &duv_dx[0], 1, &tmp2[0], 1);
            
            // tmp2 = -(gamma-1)* duv_dx
            Vmath::Smul(nq, -GammaMinOne, &tmp2[0], 1, &tmp2[0], 1);
            
            // tmp2 = -(gamma-1)* duv_dx * z_rhoE
            Vmath::Vmul(nq, &inarray[3][0], 1, &tmp2[0], 1, &tmp2[0], 1);
            
            Vmath::Vadd(nq, &tmp[0], 1, &tmp1[0], 1, &outarray[2][0][0], 1);
            Vmath::Vadd(nq, &outarray[2][0][0], 1, &tmp2[0], 1, &outarray[2][0][0], 1);
            
            // ===========================Eq4===================================
            
            Vmath::Zero(nq, &tmp[0], 1);
            Vmath::Zero(nq, &tmp1[0], 1);
            Vmath::Zero(nq, &tmp2[0], 1);
            Vmath::Zero(nq, &tmp3[0], 1);
            
            Vmath::Smul(nq, m_gamma, &du_dx[0], 1, &tmp1[0], 1);
            Vmath::Vmul(nq, &tmp1[0], 1, &inarray[3][0], 1, &outarray[3][0][0], 1);
            
            // =================================================================
            //                   Create JacTransVectorY
            // =================================================================
            // ===========================Eq1===================================
            
            Vmath::Zero(nq, &tmp[0], 1);
            Vmath::Zero(nq, &tmp1[0], 1);
            Vmath::Zero(nq, &tmp2[0], 1);
            Vmath::Zero(nq, &tmp3[0], 1);
            
            // tmp = dv^2_dy
            Vmath::Vcopy(nq, &dv2_dy[0], 1, &tmp[0], 1);
            
            // tmp1 = (gamma-1)/2 * dU^2_dy
            Vmath::Smul(nq, HalfGammaMinOne, &dU2_dy[0], 1, &tmp1[0], 1);
            
            // tmp = (gamma-1)/2 * dU^2_dy - dv^2_dy
            Vmath::Vsub(nq, &tmp1[0], 1, &tmp[0], 1, &tmp[0], 1);
            
            // tmp = z_rhov*((gamma-1)/2 * dU^2_dy - du^2_dy)
            Vmath::Vmul(nq, &inarray[2][0], 1, &tmp[0], 1, &tmp[0], 1);
            
            //--------------------------------------------------------------
            
            // tmp2 = duv_dy;
            Vmath::Vcopy(nq, &duv_dy[0], 1, &tmp2[0], 1);
            
            // tmp2 = -duv_dy;
            Vmath::Neg(nq, &tmp2[0], 1);
            
            // tmp2 = -u*v*z_rhou;
            Vmath::Vmul(nq, &inarray[1][0], 1, &tmp2[0], 1, &tmp2[0], 1);
            
            //--------------------------------------------------------------
            // tmp3 = (gamma-1)/2*dvU2_dy;
            Vmath::Smul(nq, HalfGammaMinOne, &dvU2_dy[0], 1, &tmp3[0], 1);
            
            // tmp3 = (gamma-1)/2*dvU2_dy-dvH_dy;
            Vmath::Vsub(nq, &tmp3[0], 1, &dvH_dy[0], 1, &tmp3[0], 1);
            
            // tmp3 = z_rhoE*((gamma-1)/2*dvU2_dy-dvH_dy);
            Vmath::Vmul(nq, &inarray[3][0], 1, &tmp3[0], 1,  &tmp3[0], 1);
            //--------------------------------------------------------------
            
            Vmath::Vadd(nq,&tmp[0], 1, &tmp2[0], 1, &outarray[0][1][0], 1);
            
            Vmath::Vadd(nq, &outarray[0][1][0], 1, &tmp3[0], 1, &outarray[0][1][0], 1);
            
            // ===========================Eq2===================================
            
            Vmath::Zero(nq, &tmp[0], 1);
            Vmath::Zero(nq, &tmp1[0], 1);
            Vmath::Zero(nq, &tmp2[0], 1);
            Vmath::Zero(nq, &tmp3[0], 1);
            
            // tmp = -(gamma-1)*du_dy
            Vmath::Smul(nq, -GammaMinOne, &du_dy[0], 1, &tmp[0], 1);
            
            // tmp = -(gamma-1)*du_dy * z_rhov
            Vmath::Vmul(nq, &tmp[0], 1, &inarray[2][0], 1, &tmp[0], 1);
            
            // tmp1 = dv_dy * z_rhou
            Vmath::Vmul(nq, &dv_dy[0], 1, &inarray[1][0], 1, &tmp1[0], 1);
            
            // tmp2 =  duv_dy
            Vmath::Vcopy(nq, &duv_dy[0], 1, &tmp2[0], 1);
            
            // tmp2 = -(gamma-1)* duv_dy
            Vmath::Smul(nq, -GammaMinOne, &tmp2[0], 1, &tmp2[0], 1);
            
            // tmp2 = -(gamma-1)* duv_dy * z_rhoE
            Vmath::Vmul(nq, &inarray[3][0], 1, &tmp2[0], 1, &tmp2[0], 1);
            
            Vmath::Vadd(nq, &tmp[0], 1, &tmp1[0], 1, &outarray[1][1][0], 1);
            Vmath::Vadd(nq, &outarray[1][1][0], 1, &tmp2[0], 1, &outarray[1][1][0], 1);
            
            // ===========================Eq3===============================
            
            Vmath::Zero(nq, &tmp[0], 1);
            Vmath::Zero(nq, &tmp1[0], 1);
            Vmath::Zero(nq, &tmp2[0], 1);
            Vmath::Zero(nq, &tmp3[0], 1);
            
            // tmp = zrho
            // Vmath::Vcopy(nq, &inarray[0][0], 1, &tmp[0], 1);
            
            // tmp1 = (3-gamma) * dv_dy
            Vmath::Smul(nq,ThreeMinGam, &dv_dy[0], 1, &tmp1[0], 1);
            
            // tmp1 = (3-gamma) * dv_dy * z_rhov
            Vmath::Vmul(nq, &tmp1[0], 1, &inarray[2][0], 1, &tmp1[0], 1);
            
            // tmp2 = du_dy
            Vmath::Vcopy(nq, &du_dy[0], 1, &tmp2[0], 1);
            
            // tmp2 = du_dy * z_rhou
            Vmath::Vmul(nq, &tmp2[0], 1, &inarray[1][0], 1, &tmp2[0], 1);
            
            // tmp3 = dv^2_dy
            Vmath::Vcopy(nq, &dv2_dy[0], 1, &tmp3[0], 1);
            
            // tmp3 = (1-gamma) * dv^2_dy
            Vmath::Smul(nq, -GammaMinOne, &tmp3[0], 1, &tmp3[0], 1);
            
            // tmp3 = (gamma-1) * du^2_dy+dH_dy
            Vmath::Vadd(nq, &tmp3[0], 1, &dH_dy[0], 1, &tmp3[0], 1);
            
            // tmp3 = ((gamma-1) * du^2_dy+dH_dy)*z_rhoE
            Vmath::Vmul(nq, &inarray[3][0], 1, &tmp3[0], 1, &tmp3[0], 1);
            
            Vmath::Vcopy(nq, &tmp1[0], 1, &outarray[2][1][0], 1);
            
            Vmath::Vadd(nq, &outarray[2][1][0], 1, &tmp2[0], 1, &outarray[2][1][0], 1);
            
            Vmath::Vadd(nq, &outarray[2][1][0], 1, &tmp3[0], 1, &outarray[2][1][0], 1);
            // ===========================Eq4===============================
            
            Vmath::Zero(nq, &tmp[0], 1);
            Vmath::Zero(nq, &tmp1[0], 1);
            Vmath::Zero(nq, &tmp2[0], 1);
            Vmath::Zero(nq, &tmp3[0], 1);
            
            Vmath::Smul(nq, m_gamma, &dv_dy[0], 1, &tmp1[0], 1);
            Vmath::Vmul(nq, &tmp1[0], 1, &inarray[3][0], 1, &outarray[3][1][0], 1);
            
            Vmath::Zero(nq, &tmp[0], 1);
            Vmath::Zero(nq, &tmp1[0], 1);
            Vmath::Zero(nq, &tmp2[0], 1);
            Vmath::Zero(nq, &tmp3[0], 1);
        }
        
        if (m_spacedim == 3)
        {
            ASSERTL0(false, "3D adjoint solver is not yet implemented");
        }
    }
    
    
    void CompressibleFlowSystem::GetAdjointFluxVector(
              const Array<OneD, Array<OneD, NekDouble> > &inarray,
                    Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &outarray)
    {
        
        int i;
        int nq = inarray[0].num_elements();
        int nvar = m_fields.num_elements();
        
        Array<OneD, NekDouble> pressure(nq,0.0);
        Array<OneD, NekDouble> H(nq,0.0);
        Array<OneD, NekDouble> Htmp(nq,0.0);
        Array<OneD, NekDouble> velsq(nq,0.0);
        Array<OneD, Array<OneD, NekDouble> > vel(m_spacedim);
        Array<OneD, Array<OneD, NekDouble> > direct_fields(nvar);
        
        for (i = 0; i < nvar; ++i)
        {
            direct_fields[i] = Array<OneD, NekDouble>(nq, 0.0);
        }
        
        GetDirectSolution(direct_fields);
        
        // load the physical values to obtain the adjoint vector \phi
        
        // Flux vector for the rho equation
        for (i = 0; i < m_spacedim; ++i)
        {
            vel[i] = Array<OneD, NekDouble>(nq);
        }
        
        GetVelocityVector(direct_fields, vel);
        GetPressure      (direct_fields, vel, pressure);
        
        Array<OneD, NekDouble> ones(nq, 1.0);
        
        // determine u^2+v^2+w^2;
        for (i = 0; i < m_spacedim; ++i)
        {
            Vmath::Vvtvp(nq, vel[i], 1, vel[i], 1, velsq, 1, velsq, 1);
        }
        
        // determine H = E + p/rho;
        Vmath::Vdiv(nq, direct_fields[nvar-1], 1, direct_fields[0], 1, H, 1);
        
        Vmath::Vdiv(nq, pressure, 1, direct_fields[0], 1, Htmp, 1);
        
        Vmath::Vadd(nq, H, 1, Htmp, 1, H, 1);
        
        // determine constants that are used regularly
        NekDouble GammaMinOne     = m_gamma - 1;
        NekDouble HalfGammaMinOne = 0.5*(m_gamma - 1);
        NekDouble ThreeMinGam     = 3.0 - m_gamma;
        
        if (m_spacedim == 1)
        {
            ASSERTL0(false, "1D adjoint solver is not yet implemented");
        }
        if (m_spacedim == 2)
        {
            Array<OneD, NekDouble> tmp(nq, 0.0);
            Array<OneD, NekDouble> tmp1(nq, 0.0);
            Array<OneD, NekDouble> tmp2(nq, 0.0);
            Array<OneD, NekDouble> tmp3(nq, 0.0);
            
            // tmp = u^2
            Vmath::Vmul(nq, &vel[0][0], 1, &vel[0][0], 1, &tmp[0], 1);
            
            // tmp1 = (gamma-1)/2 * velsq
            Vmath::Smul(nq, HalfGammaMinOne, &velsq[0], 1, &tmp1[0], 1);
            
            // ====================================================
            
            // tmp = (gamma-1)/2 * velsq - u^2
            Vmath::Vsub(nq, &tmp1[0], 1, &tmp[0], 1, &tmp[0], 1);
            
            // tmp = z_rhou*((gamma/2 - 1/2)*(v1^2 + v2^2) - v1^2)
            Vmath::Vmul(nq, &inarray[1][0], 1, &tmp[0], 1, &tmp[0], 1);
            
            // tmp2 = u*v;
            Vmath::Vmul(nq, &vel[0][0], 1, &vel[1][0], 1, &tmp2[0], 1);
            
            // tmp2 = -u*v;
            Vmath::Neg(nq, &tmp2[0], 1);
            
            // tmp2 = -u*v*z_rhov;
            Vmath::Vmul(nq, &inarray[2][0], 1, &tmp2[0], 1, &tmp2[0], 1);
            
            // tmp3 = (gamma-1)/2*velsq;
            Vmath::Smul(nq, HalfGammaMinOne, &velsq[0], 1, &tmp3[0], 1);
            
            // tmp3 = ((gamma-1)/2*velsq-H);
            Vmath::Vsub(nq, &tmp3[0], 1, &H[0], 1, &tmp3[0], 1);
            
            // tmp3 = u((gamma-1)/2*velsq-H);
            Vmath::Vmul(nq, &vel[0][0], 1, &tmp3[0], 1, &tmp3[0], 1);
            
            // tmp3 = z_rhoE*u((gamma-1)/2*velsq-H);
            Vmath::Vmul(nq, &inarray[3][0], 1, &tmp3[0], 1, &tmp3[0], 1);
            
            Vmath::Vadd(nq, &tmp[0], 1, &tmp2[0], 1, &outarray[0][0][0], 1);
            
            Vmath::Vadd(nq, &outarray[0][0][0], 1, &tmp3[0], 1, &outarray[0][0][0], 1);
            
            // ====================================================
            
            Vmath::Zero(nq, &tmp[0], 1);
            Vmath::Zero(nq, &tmp1[0], 1);
            Vmath::Zero(nq, &tmp2[0], 1);
            Vmath::Zero(nq, &tmp3[0], 1);
            
            // tmp = z_rho
            Vmath::Vcopy(nq, &inarray[0][0], 1, &tmp[0], 1);
            
            // tmp1 = (3-gamma) * u
            Vmath::Smul(nq, ThreeMinGam, &vel[0][0], 1, &tmp1[0], 1);
            
            // tmp1 = (3-gamma) * u * z_rhou
            Vmath::Vmul(nq, &tmp1[0], 1, &inarray[1][0], 1, &tmp1[0], 1);
            
            // tmp2 = v
            Vmath::Vcopy(nq, &vel[1][0], 1, &tmp2[0], 1);
            
            // tmp2 = v * z_rhov
            Vmath::Vmul(nq, &tmp2[0], 1, &inarray[2][0], 1, &tmp2[0], 1);
            
            // tmp3 = u^2
            Vmath::Vmul(nq, &vel[0][0], 1, &vel[0][0], 1, &tmp3[0], 1);
            
            // tmp3 = (1-gamma) * u^2
            Vmath::Smul(nq, -GammaMinOne, &tmp3[0], 1, &tmp3[0], 1);
            
            // tmp3 = (gamma-1) * u^2+H
            Vmath::Vadd(nq, &tmp3[0], 1, &H[0], 1, &tmp3[0], 1);
            
            // tmp3 = ((gamma-1) * u^2+H)*z_rhoE
            Vmath::Vmul(nq, &inarray[3][0], 1, &tmp3[0], 1, &tmp3[0], 1);
            
            
            Vmath::Vadd(nq, &tmp[0], 1, &tmp1[0], 1, &outarray[1][0][0], 1);
            Vmath::Vadd(nq, &outarray[1][0][0], 1, &tmp2[0], 1, &outarray[1][0][0], 1);
            Vmath::Vadd(nq, &outarray[1][0][0], 1, &tmp3[0], 1, &outarray[1][0][0], 1);
            
            // ====================================================
            
            Vmath::Zero(nq, &tmp[0], 1);
            Vmath::Zero(nq, &tmp1[0], 1);
            Vmath::Zero(nq, &tmp2[0], 1);
            Vmath::Zero(nq, &tmp3[0], 1);
            
            // tmp = -(gamma-1)*v
            Vmath::Smul(nq, -GammaMinOne, &vel[1][0], 1, &tmp[0], 1);
            
            // tmp = -(gamma-1)*v * z_rhou
            Vmath::Vmul(nq, &tmp[0], 1, &inarray[1][0], 1, &tmp[0], 1);
            
            // tmp1 = u * z_rhov
            Vmath::Vmul(nq, &vel[0][0], 1, &inarray[2][0], 1, &tmp1[0], 1);
            
            // tmp2 =  u * v
            Vmath::Vmul(nq, &vel[0][0], 1, &vel[1][0], 1, &tmp2[0], 1);
            
            // tmp2 = -(gamma-1) * u * v
            Vmath::Smul(nq, -GammaMinOne, &tmp2[0], 1, &tmp2[0], 1);
            
            // tmp2 = -(gamma-1) * u * v * z_rhoE
            Vmath::Vmul(nq, &inarray[3][0], 1, &tmp2[0], 1, &tmp2[0], 1);
            
            Vmath::Vadd(nq, &tmp[0], 1, &tmp1[0], 1, &outarray[2][0][0], 1);
            Vmath::Vadd(nq, &outarray[2][0][0], 1, &tmp2[0], 1, &outarray[2][0][0], 1);
            
            Vmath::Zero(nq, &tmp[0], 1);
            Vmath::Zero(nq, &tmp1[0], 1);
            Vmath::Zero(nq, &tmp2[0], 1);
            Vmath::Zero(nq, &tmp3[0], 1);
            
            Vmath::Smul(nq, GammaMinOne, &inarray[1][0], 1, &tmp[0], 1);
            
            Vmath::Smul(nq, m_gamma, &vel[0][0], 1, &tmp1[0], 1);
            Vmath::Vmul(nq, &tmp1[0], 1, &inarray[3][0], 1, &tmp1[0], 1);
            
            Vmath::Vadd(nq, &tmp[0], 1, &tmp1[0], 1, &outarray[3][0][0], 1);
            
            // =================================================================
            // =================================================================
            
            Vmath::Zero(nq, &tmp[0], 1);
            Vmath::Zero(nq, &tmp1[0], 1);
            Vmath::Zero(nq, &tmp2[0], 1);
            Vmath::Zero(nq, &tmp3[0], 1);
            
            // tmp = v^2
            Vmath::Vmul(nq, &vel[1][0], 1, &vel[1][0], 1, &tmp[0], 1);
            
            // tmp1 = (gamma-1)/2 * velsq
            Vmath::Smul(nq, HalfGammaMinOne, &velsq[0], 1, &tmp1[0], 1);
            
            
            // tmp = (gamma-1)/2 * velsq - v^2
            Vmath::Vsub(nq, &tmp1[0], 1, &tmp[0], 1, &tmp[0], 1);
            
            // tmp = z_rhov*((gamma/2 - 1/2)*(v1^2 + v2^2) - v1^2)
            Vmath::Vmul(nq, &inarray[2][0], 1, &tmp[0], 1, &tmp[0], 1);
            
            
            // ====================================================
            
            // tmp2 = u*v;
            Vmath::Vmul(nq, &vel[0][0], 1, &vel[1][0], 1, &tmp2[0], 1);
            
            // tmp2 = -u*v;
            Vmath::Neg(nq, &tmp2[0], 1);
            
            // tmp2 = -u*v*z_rhou;
            Vmath::Vmul(nq, &inarray[1][0], 1, &tmp2[0], 1, &tmp2[0], 1);
            
            // ====================================================
            
            // tmp3 = (gamma-1)/2*velsq;
            Vmath::Smul(nq, HalfGammaMinOne, &velsq[0], 1, &tmp3[0], 1);
            
            // tmp3 = ((gamma-1)/2*velsq-H);
            Vmath::Vsub(nq, &tmp3[0], 1, &H[0], 1, &tmp3[0], 1);
            
            // tmp3 = v((gamma-1)/2*velsq-H);
            Vmath::Vmul(nq,
                        &vel[1][0], 1,
                        &tmp3[0], 1,
                        &tmp3[0], 1);
            
            // tmp3 = z_rhoE*u((gamma-1)/2*velsq-H);
            Vmath::Vmul(nq, &inarray[3][0], 1, &tmp3[0], 1, &tmp3[0], 1);
            
            Vmath::Vadd(nq, &tmp[0], 1, &tmp2[0], 1, &outarray[0][1][0], 1);
            
            Vmath::Vadd(nq, &outarray[0][1][0], 1, &tmp3[0], 1, &outarray[0][1][0], 1);
            
            // ====================================================
            
            Vmath::Zero(nq, &tmp[0], 1);
            Vmath::Zero(nq, &tmp1[0], 1);
            Vmath::Zero(nq, &tmp2[0], 1);
            Vmath::Zero(nq, &tmp3[0], 1);
            
            // tmp = z_rho
            Vmath::Vcopy(nq, &inarray[0][0], 1, &tmp[0], 1);
            
            // tmp1 = -(gamma-3) * v
            Vmath::Smul(nq, ThreeMinGam, &vel[1][0], 1, &tmp1[0], 1);
            // tmp1 = -(gamma-3) * v * z_rhov
            Vmath::Vmul(nq, &tmp1[0], 1, &inarray[2][0], 1, &tmp1[0], 1);
            
            
            
            // tmp2 = u
            Vmath::Vcopy(nq, &vel[0][0], 1, &tmp2[0], 1);
            
            // tmp2 = u * z_rhou
            Vmath::Vmul(nq, &tmp2[0], 1, &inarray[1][0], 1, &tmp2[0], 1);
            
            
            
            // tmp3 = v^2
            Vmath::Vmul(nq, &vel[1][0], 1, &vel[1][0], 1, &tmp3[0], 1);
            
            // tmp3 = (1-gamma) * v^2
            Vmath::Smul(nq, -GammaMinOne, &tmp3[0], 1, &tmp3[0], 1);
            
            // tmp3 = (gamma-1) * v^2+H
            Vmath::Vadd(nq, &tmp3[0], 1, &H[0], 1, &tmp3[0], 1);
            
            // tmp3 = ((gamma-1) * v^2+H)*z_rhoE
            Vmath::Vmul(nq, &inarray[3][0], 1, &tmp3[0], 1, &tmp3[0], 1);
            
            
            
            Vmath::Vadd(nq, &tmp[0], 1, &tmp1[0], 1, &outarray[2][1][0], 1);
            
            Vmath::Vadd(nq, &outarray[2][1][0], 1, &tmp2[0], 1, &outarray[2][1][0], 1);
            
            Vmath::Vadd(nq, &outarray[2][1][0], 1, &tmp3[0], 1, &outarray[2][1][0], 1);
            
            // ====================================================
            
            Vmath::Zero(nq, &tmp[0], 1);
            Vmath::Zero(nq, &tmp1[0], 1);
            Vmath::Zero(nq, &tmp2[0], 1);
            Vmath::Zero(nq, &tmp3[0], 1);
            
            // tmp = -(gamma-1) * u
            Vmath::Smul(nq, -GammaMinOne, &vel[0][0], 1, &tmp[0], 1);
            
            // tmp = -(gamma-1) * u * z_rhov
            Vmath::Vmul(nq, &tmp[0], 1, &inarray[2][0], 1, &tmp[0], 1);
            
            // tmp1 = v * z_rhou
            Vmath::Vmul(nq, &vel[1][0], 1, &inarray[1][0], 1, &tmp1[0], 1);
            
            // tmp2 =  u * v
            Vmath::Vmul(nq, &vel[0][0], 1, &vel[1][0], 1, &tmp2[0], 1);
            
            // tmp2 = -(gamma-1)* u * v
            Vmath::Smul(nq, -GammaMinOne, &tmp2[0], 1, &tmp2[0], 1);
            
            // tmp2 = -(gamma-1)* u * v * z_rhoE
            Vmath::Vmul(nq, &inarray[3][0], 1, &tmp2[0], 1, &tmp2[0], 1);
            
            Vmath::Vadd(nq, &tmp[0], 1, &tmp1[0], 1, &outarray[1][1][0], 1);
            
            Vmath::Vadd(nq, &outarray[1][1][0], 1, &tmp2[0], 1, &outarray[1][1][0], 1);
            
            // ====================================================
            
            Vmath::Zero(nq, &tmp[0], 1);
            Vmath::Zero(nq, &tmp1[0], 1);
            Vmath::Zero(nq, &tmp2[0], 1);
            Vmath::Zero(nq, &tmp3[0], 1);
            
            Vmath::Smul(nq, GammaMinOne, &inarray[2][0], 1, &tmp[0], 1);
            
            Vmath::Smul(nq, m_gamma, &vel[1][0], 1, &tmp1[0], 1);
            Vmath::Vmul(nq, &tmp1[0], 1, &inarray[3][0], 1, &tmp1[0], 1);
            
            Vmath::Vadd(nq, &tmp[0], 1, &tmp1[0], 1, &outarray[3][1][0], 1);
            
        }
        if (m_spacedim == 3)
        {
     	    Array<OneD, NekDouble> tmp(nq, 0.0);
            Array<OneD, NekDouble> tmp1(nq, 0.0);
            Array<OneD, NekDouble> tmp2(nq, 0.0);
            Array<OneD, NekDouble> tmp3(nq, 0.0);
            Array<OneD, NekDouble> tmp4(nq, 0.0);
            
            //==================================================================
            //============================ Ax^t*z ==============================
            //==================================================================
            
            // tmp = u^2
            Vmath::Vmul(nq, &vel[0][0], 1, &vel[0][0], 1, &tmp[0], 1);
            
            // tmp1 = (gamma-1)/2 * velsq
            Vmath::Smul(nq, HalfGammaMinOne, &velsq[0], 1, &tmp1[0], 1);
            
            // ====================================================
            // tmp = (gamma-1)/2 * velsq - u^2
            Vmath::Vsub(nq, &tmp1[0], 1, &tmp[0], 1, &tmp[0], 1);
            
            // tmp = z_rhou*((gamma/2 - 1/2)*(v1^2 + v2^2) - v1^2)
            Vmath::Vmul(nq, &inarray[1][0], 1, &tmp[0], 1, &tmp[0], 1);
            
            // tmp2 = u*v;
            Vmath::Vmul(nq, &vel[0][0], 1, &vel[1][0], 1, &tmp2[0], 1);
            
            // tmp2 = -u*v;
            Vmath::Neg(nq, &tmp2[0], 1);
            
            // tmp2 = -u*v*z_rhov;
            Vmath::Vmul(nq, &inarray[2][0], 1, &tmp2[0], 1, &tmp2[0], 1);
            
            // tmp3 = u*v;
            Vmath::Vmul(nq, &vel[0][0], 1, &vel[2][0], 1, &tmp3[0], 1);
            
            // tmp3 = -u*w;
            Vmath::Neg(nq, &tmp3[0], 1);
            
            // tmp3 = -u*w*z_rhow;
            Vmath::Vmul(nq, &inarray[3][0], 1, &tmp3[0], 1, &tmp3[0], 1);
            
            // tmp3 = (gamma-1)/2*velsq;
            Vmath::Smul(nq, HalfGammaMinOne, &velsq[0], 1, &tmp4[0], 1);
            
            // tmp3 = ((gamma-1)/2*velsq-H);
            Vmath::Vsub(nq, &tmp4[0], 1, &H[0], 1, &tmp4[0], 1);
            
            // tmp3 = u((gamma-1)/2*velsq-H);
            Vmath::Vmul(nq, &vel[0][0], 1, &tmp4[0], 1, &tmp4[0], 1);
            
            // tmp3 = z_rhoE*u((gamma-1)/2*velsq-H);
            Vmath::Vmul(nq, &inarray[4][0], 1, &tmp4[0], 1, &tmp4[0], 1);
            
            Vmath::Vadd(nq, &tmp[0], 1, &tmp2[0], 1, &outarray[0][0][0], 1);
            Vmath::Vadd(nq, &outarray[0][0][0], 1, &tmp3[0], 1, &outarray[0][0][0], 1);
            Vmath::Vadd(nq, &outarray[0][0][0], 1, &tmp4[0], 1, &outarray[0][0][0], 1);
            
            // ====================================================
            
            Vmath::Zero(nq, &tmp[0], 1);
            Vmath::Zero(nq, &tmp1[0], 1);
            Vmath::Zero(nq, &tmp2[0], 1);
            Vmath::Zero(nq, &tmp3[0], 1);
            Vmath::Zero(nq, &tmp4[0], 1);
            
            // tmp = z_rho
            Vmath::Vcopy(nq, &inarray[0][0], 1, &tmp[0], 1);
            
            // tmp1 = (3-gamma) * u
            Vmath::Smul(nq, ThreeMinGam, &vel[0][0], 1, &tmp1[0], 1);
            
            // tmp1 = (3-gamma) * u * z_rhou
            Vmath::Vmul(nq, &tmp1[0], 1, &inarray[1][0], 1, &tmp1[0], 1);
            
            // tmp2 = v
            Vmath::Vcopy(nq, &vel[1][0], 1, &tmp2[0], 1);
            
            // tmp2 = v * z_rhov
            Vmath::Vmul(nq, &tmp2[0], 1, &inarray[2][0], 1, &tmp2[0], 1);
            
            // tmp3 = w
            Vmath::Vcopy(nq, &vel[2][0], 1, &tmp3[0], 1);
            
            // tmp3 = w * z_rhow
            Vmath::Vmul(nq, &tmp3[0], 1, &inarray[3][0], 1, &tmp3[0], 1);
            
            // tmp4 = u^2
            Vmath::Vmul(nq, &vel[0][0], 1, &vel[0][0], 1, &tmp4[0], 1);
            
            // tmp4 = (1-gamma) * u^2
            Vmath::Smul(nq, -GammaMinOne, &tmp4[0], 1, &tmp4[0], 1);
            
            // tmp4 = (gamma-1) * u^2+H
            Vmath::Vadd(nq, &tmp4[0], 1, &H[0], 1, &tmp4[0], 1);
            
            // tmp3 = ((gamma-1) * u^2+H)*z_rhoE
            Vmath::Vmul(nq, &inarray[4][0], 1, &tmp4[0], 1, &tmp4[0], 1);
            
            
            Vmath::Vadd(nq, &tmp[0], 1, &tmp1[0], 1, &outarray[1][0][0], 1);
            Vmath::Vadd(nq, &outarray[1][0][0], 1, &tmp2[0], 1, &outarray[1][0][0], 1);
            Vmath::Vadd(nq, &outarray[1][0][0], 1, &tmp3[0], 1, &outarray[1][0][0], 1);
            Vmath::Vadd(nq, &outarray[1][0][0], 1, &tmp4[0], 1, &outarray[1][0][0], 1);
            
            // ====================================================
            
            Vmath::Zero(nq, &tmp[0], 1);
            Vmath::Zero(nq, &tmp1[0], 1);
            Vmath::Zero(nq, &tmp2[0], 1);
            Vmath::Zero(nq, &tmp3[0], 1);
            Vmath::Zero(nq, &tmp4[0], 1);
            
            // tmp = -(gamma-1)*v
            Vmath::Smul(nq, -GammaMinOne, &vel[1][0], 1, &tmp[0], 1);
            
            // tmp = -(gamma-1)*v * z_rhou
            Vmath::Vmul(nq, &tmp[0], 1, &inarray[1][0], 1, &tmp[0], 1);
            
            // tmp1 = u * z_rhov
            Vmath::Vmul(nq, &vel[0][0], 1, &inarray[2][0], 1, &tmp1[0], 1);
            
            // tmp2 =  u * v
            Vmath::Vmul(nq, &vel[0][0], 1, &vel[1][0], 1, &tmp2[0], 1);
            
            // tmp2 = -(gamma-1) * u * v
            Vmath::Smul(nq, -GammaMinOne, &tmp2[0], 1, &tmp2[0], 1);
            
            // tmp2 = -(gamma-1) * u * v * z_rhoE
            Vmath::Vmul(nq, &inarray[3][0], 1, &tmp2[0], 1, &tmp2[0], 1);
            
            Vmath::Vadd(nq, &tmp[0], 1, &tmp1[0], 1, &outarray[2][0][0], 1);
            Vmath::Vadd(nq, &outarray[2][0][0], 1, &tmp2[0], 1, &outarray[2][0][0], 1);
            
            // ====================================================
            
            Vmath::Zero(nq, &tmp[0], 1);
            Vmath::Zero(nq, &tmp1[0], 1);
            Vmath::Zero(nq, &tmp2[0], 1);
            Vmath::Zero(nq, &tmp3[0], 1);
            Vmath::Zero(nq, &tmp4[0], 1);
            
            // tmp = -(gamma-1)*w
            Vmath::Smul(nq, -GammaMinOne, &vel[2][0], 1, &tmp[0], 1);
            
            // tmp = -(gamma-1)*w * z_rhou
            Vmath::Vmul(nq, &tmp[0], 1, &inarray[1][0], 1, &tmp[0], 1);
            
            // tmp1 = u * z_rhow
            Vmath::Vmul(nq, &vel[0][0], 1, &inarray[3][0], 1, &tmp1[0], 1);
            
            // tmp2 =  u * w
            Vmath::Vmul(nq, &vel[0][0], 1, &vel[2][0], 1, &tmp2[0], 1);
            
            // tmp2 = -(gamma-1) * u * v
            Vmath::Smul(nq, -GammaMinOne, &tmp2[0], 1, &tmp2[0], 1);
            
            // tmp2 = -(gamma-1) * u * v * z_rhoE
            Vmath::Vmul(nq, &inarray[4][0], 1, &tmp2[0], 1, &tmp2[0], 1);
            
            Vmath::Vadd(nq, &tmp[0], 1, &tmp1[0], 1, &outarray[3][0][0], 1);
            Vmath::Vadd(nq, &outarray[3][0][0], 1, &tmp2[0], 1, &outarray[3][0][0], 1);
            
            //======================================================================
            Vmath::Zero(nq, &tmp[0], 1);
            Vmath::Zero(nq, &tmp1[0], 1);
            Vmath::Zero(nq, &tmp2[0], 1);
            Vmath::Zero(nq, &tmp3[0], 1);
            
            Vmath::Smul(nq, GammaMinOne, &inarray[1][0], 1, &tmp[0], 1);
            
            Vmath::Smul(nq, m_gamma, &vel[0][0], 1, &tmp1[0], 1);
            Vmath::Vmul(nq, &tmp1[0], 1, &inarray[4][0], 1, &tmp1[0], 1);
            
            Vmath::Vadd(nq, &tmp[0], 1, &tmp1[0], 1, &outarray[4][0][0], 1);
            
            
            //==================================================================
            //============================ Ay^t*z ==============================
            //==================================================================
            
            Vmath::Zero(nq, &tmp[0], 1);
            Vmath::Zero(nq, &tmp1[0], 1);
            Vmath::Zero(nq, &tmp2[0], 1);
            Vmath::Zero(nq, &tmp3[0], 1);
            Vmath::Zero(nq, &tmp4[0], 1);
            
            // tmp = v^2
            Vmath::Vmul(nq, &vel[1][0], 1, &vel[1][0], 1, &tmp[0], 1);
            
            // tmp1 = (gamma-1)/2 * velsq
            Vmath::Smul(nq, HalfGammaMinOne, &velsq[0], 1, &tmp1[0], 1);
            
            
            // tmp = (gamma-1)/2 * velsq - v^2
            Vmath::Vsub(nq, &tmp1[0], 1, &tmp[0], 1, &tmp[0], 1);
            
            // tmp = z_rhov*((gamma/2 - 1/2)*(v1^2 + v2^2) - v1^2)
            Vmath::Vmul(nq, &inarray[2][0], 1, &tmp[0], 1, &tmp[0], 1);
            
            
            // ====================================================
            
            // tmp2 = u*v;
            Vmath::Vmul(nq, &vel[0][0], 1, &vel[1][0], 1, &tmp2[0], 1);
            
            // tmp2 = -u*v;
            Vmath::Neg(nq, &tmp2[0], 1);
            
            // tmp2 = -u*v*z_rhou;
            Vmath::Vmul(nq, &inarray[1][0], 1, &tmp2[0], 1, &tmp2[0], 1);
            
            // ====================================================
            
            // tmp3 = u*w;
            Vmath::Vmul(nq, &vel[0][0], 1, &vel[2][0], 1, &tmp3[0], 1);
            
            // tmp3 = -u*w;
            Vmath::Neg(nq, &tmp3[0], 1);
            
            // tmp2 = -u*w*z_rhou;
            Vmath::Vmul(nq, &inarray[1][0], 1, &tmp3[0], 1, &tmp3[0], 1);
            
            // ====================================================
            
            // tmp3 = (gamma-1)/2*velsq;
            Vmath::Smul(nq, HalfGammaMinOne, &velsq[0], 1, &tmp4[0], 1);
            
            // tmp3 = ((gamma-1)/2*velsq-H);
            Vmath::Vsub(nq, &tmp4[0], 1, &H[0], 1, &tmp4[0], 1);
            
            // tmp3 = v((gamma-1)/2*velsq-H);
            Vmath::Vmul(nq,
                        &vel[1][0], 1,
                        &tmp4[0], 1,
                        &tmp4[0], 1);
            
            // tmp3 = z_rhoE*u((gamma-1)/2*velsq-H);
            Vmath::Vmul(nq, &inarray[3][0], 1, &tmp4[0], 1, &tmp4[0], 1);
            
            Vmath::Vadd(nq, &tmp[0], 1, &tmp2[0], 1, &outarray[0][1][0], 1);
            Vmath::Vadd(nq, &outarray[0][1][0], 1, &tmp3[0], 1, &outarray[0][1][0], 1);
            Vmath::Vadd(nq, &outarray[0][1][0], 1, &tmp4[0], 1, &outarray[0][1][0], 1);
            
            // ====================================================
            
            Vmath::Zero(nq, &tmp[0], 1);
            Vmath::Zero(nq, &tmp1[0], 1);
            Vmath::Zero(nq, &tmp2[0], 1);
            Vmath::Zero(nq, &tmp3[0], 1);
            Vmath::Zero(nq, &tmp4[0], 1);
            
            // tmp = z_rho
            Vmath::Vcopy(nq, &inarray[0][0], 1, &tmp[0], 1);
            
            // tmp1 = -(gamma-3) * v
            Vmath::Smul(nq, ThreeMinGam, &vel[1][0], 1, &tmp1[0], 1);
            // tmp1 = -(gamma-3) * v * z_rhov
            Vmath::Vmul(nq, &tmp1[0], 1, &inarray[2][0], 1, &tmp1[0], 1);
            
            //====================================================
            // tmp2 = u
            Vmath::Vcopy(nq, &vel[0][0], 1, &tmp2[0], 1);
            
            // tmp2 = u * z_rhou
            Vmath::Vmul(nq, &tmp2[0], 1, &inarray[1][0], 1, &tmp2[0], 1);
            
            //====================================================
            // tmp3 = u
            Vmath::Vcopy(nq, &vel[2][0], 1, &tmp3[0], 1);
            
            // tmp3 = w * z_rhow
            Vmath::Vmul(nq, &tmp3[0], 1, &inarray[3][0], 1, &tmp3[0], 1);
            
            //====================================================
            // tmp4 = v^2
            Vmath::Vmul(nq, &vel[1][0], 1, &vel[1][0], 1, &tmp4[0], 1);
            
            // tmp4 = (1-gamma) * v^2
            Vmath::Smul(nq, -GammaMinOne, &tmp4[0], 1, &tmp4[0], 1);
            
            // tmp4 = (gamma-1) * v^2+H
            Vmath::Vadd(nq, &tmp4[0], 1, &H[0], 1, &tmp4[0], 1);
            
            // tmp4 = ((gamma-1) * v^2+H)*z_rhoE
            Vmath::Vmul(nq, &inarray[4][0], 1, &tmp4[0], 1, &tmp4[0], 1);
            
            Vmath::Vadd(nq, &tmp[0], 1, &tmp1[0], 1, &outarray[2][1][0], 1);
            Vmath::Vadd(nq, &outarray[2][1][0], 1, &tmp2[0], 1, &outarray[2][1][0], 1);
            Vmath::Vadd(nq, &outarray[2][1][0], 1, &tmp3[0], 1, &outarray[2][1][0], 1);
            Vmath::Vadd(nq, &outarray[2][1][0], 1, &tmp4[0], 1, &outarray[2][1][0], 1);
            // ====================================================
            
            Vmath::Zero(nq, &tmp[0], 1);
            Vmath::Zero(nq, &tmp1[0], 1);
            Vmath::Zero(nq, &tmp2[0], 1);
            Vmath::Zero(nq, &tmp3[0], 1);
            
            // tmp = -(gamma-1) * u
            Vmath::Smul(nq, -GammaMinOne, &vel[0][0], 1, &tmp[0], 1);
            
            // tmp = -(gamma-1) * u * z_rhov
            Vmath::Vmul(nq, &tmp[0], 1, &inarray[2][0], 1, &tmp[0], 1);
            
            // tmp1 = v * z_rhou
            Vmath::Vmul(nq, &vel[1][0], 1, &inarray[1][0], 1, &tmp1[0], 1);
            
            // tmp2 =  u * v
            Vmath::Vmul(nq, &vel[0][0], 1, &vel[1][0], 1, &tmp2[0], 1);
            
            // tmp2 = -(gamma-1)* u * v
            Vmath::Smul(nq, -GammaMinOne, &tmp2[0], 1, &tmp2[0], 1);
            
            // tmp2 = -(gamma-1)* u * v * z_rhoE
            Vmath::Vmul(nq, &inarray[3][0], 1, &tmp2[0], 1, &tmp2[0], 1);
            
            Vmath::Vadd(nq, &tmp[0], 1, &tmp1[0], 1, &outarray[1][1][0], 1);
            
            Vmath::Vadd(nq, &outarray[1][1][0], 1, &tmp2[0], 1, &outarray[1][1][0], 1);
            
            // ====================================================
            // ====================================================
            
            Vmath::Zero(nq, &tmp[0], 1);
            Vmath::Zero(nq, &tmp1[0], 1);
            Vmath::Zero(nq, &tmp2[0], 1);
            Vmath::Zero(nq, &tmp3[0], 1);
            
            // tmp = -(gamma-1) * u
            Vmath::Smul(nq, -GammaMinOne, &vel[0][0], 1, &tmp[0], 1);
            
            // tmp = -(gamma-1) * u * z_rhov
            Vmath::Vmul(nq, &tmp[0], 1, &inarray[2][0], 1, &tmp[0], 1);
            
            // tmp1 = v * z_rhou
            Vmath::Vmul(nq, &vel[1][0], 1, &inarray[1][0], 1, &tmp1[0], 1);
            
            // tmp2 =  u * v
            Vmath::Vmul(nq, &vel[0][0], 1, &vel[1][0], 1, &tmp2[0], 1);
            
            // tmp2 = -(gamma-1)* u * v
            Vmath::Smul(nq, -GammaMinOne, &tmp2[0], 1, &tmp2[0], 1);
            
            // tmp2 = -(gamma-1)* u * v * z_rhoE
            Vmath::Vmul(nq, &inarray[4][0], 1, &tmp2[0], 1, &tmp2[0], 1);
            
            Vmath::Vadd(nq, &tmp[0], 1, &tmp1[0], 1, &outarray[1][1][0], 1);
            
            Vmath::Vadd(nq, &outarray[1][1][0], 1, &tmp2[0], 1, &outarray[1][1][0], 1);
            
            // ====================================================
            Vmath::Zero(nq, &tmp[0], 1);
            Vmath::Zero(nq, &tmp1[0], 1);
            Vmath::Zero(nq, &tmp2[0], 1);
            Vmath::Zero(nq, &tmp3[0], 1);
            
            // tmp = -(gamma-1) * w
            Vmath::Smul(nq, -GammaMinOne, &vel[2][0], 1, &tmp[0], 1);
            
            // tmp = -(gamma-1) * w * z_rhov
            Vmath::Vmul(nq, &tmp[0], 1, &inarray[2][0], 1, &tmp[0], 1);
            
            // tmp1 = v * zrhow
            Vmath::Vmul(nq, &vel[1][0], 1, &inarray[3][0], 1, &tmp1[0], 1);
            
            // tmp2 =  v * w
            Vmath::Vmul(nq, &vel[1][0], 1, &vel[2][0], 1, &tmp2[0], 1);
            
            // tmp2 = -(gamma-1) * v * w
            Vmath::Smul(nq, -GammaMinOne, &tmp2[0], 1, &tmp2[0], 1);
            
            // tmp2 = -(gamma-1) * v * w * z_rhoE
            Vmath::Vmul(nq, &inarray[4][0], 1, &tmp2[0], 1, &tmp2[0], 1);
            
            Vmath::Vadd(nq, &tmp[0], 1, &tmp1[0], 1, &outarray[3][1][0], 1);
            Vmath::Vadd(nq, &outarray[3][1][0], 1, &tmp2[0], 1, &outarray[3][1][0], 1);
            
            //===================================================
            Vmath::Zero(nq, &tmp[0], 1);
            Vmath::Zero(nq, &tmp1[0], 1);
            Vmath::Zero(nq, &tmp2[0], 1);
            Vmath::Zero(nq, &tmp3[0], 1);
            
            Vmath::Smul(nq, GammaMinOne, &inarray[2][0], 1, &tmp[0], 1);
            
            Vmath::Smul(nq, m_gamma, &vel[1][0], 1, &tmp1[0], 1);
            Vmath::Vmul(nq, &tmp1[0], 1, &inarray[4][0], 1, &tmp1[0], 1);
            
            Vmath::Vadd(nq, &tmp[0], 1, &tmp1[0], 1, &outarray[4][1][0], 1);
            
            //==================================================
            //================= Az^t*z =========================
            //==================================================
            
            Vmath::Zero(nq, &tmp[0], 1);
            Vmath::Zero(nq, &tmp1[0], 1);
            Vmath::Zero(nq, &tmp2[0], 1);
            Vmath::Zero(nq, &tmp3[0], 1);
            Vmath::Zero(nq, &tmp4[0], 1);
            
            // tmp = v^2
            Vmath::Vmul(nq, &vel[2][0], 1, &vel[2][0], 1, &tmp[0], 1);
            
            // tmp1 = (gamma-1)/2 * velsq
            Vmath::Smul(nq, HalfGammaMinOne, &velsq[0], 1, &tmp1[0], 1);
            
            // tmp = (gamma-1)/2 * velsq - v^2
            Vmath::Vsub(nq, &tmp1[0], 1, &tmp[0], 1, &tmp[0], 1);
            
            // tmp = z_rhov*((gamma/2 - 1/2)*(v1^2 + v2^2) - v1^2)
            Vmath::Vmul(nq, &inarray[3][0], 1, &tmp[0], 1, &tmp[0], 1);
            
            // ====================================================
            
            // tmp2 = u*w;
            Vmath::Vmul(nq, &vel[0][0], 1, &vel[2][0], 1, &tmp2[0], 1);
            
            // tmp2 = -u*w;
            Vmath::Neg(nq, &tmp2[0], 1);
            
            // tmp2 = -u*w*z_rhou;
            Vmath::Vmul(nq, &inarray[1][0], 1, &tmp2[0], 1, &tmp2[0], 1);
            
            // ====================================================
            
            // tmp3 = v*w;
            Vmath::Vmul(nq, &vel[1][0], 1, &vel[2][0], 1, &tmp3[0], 1);
            
            // tmp3 = -v*w;
            Vmath::Neg(nq, &tmp3[0], 1);
            
            // tmp2 = -v*w*z_rhov;
            Vmath::Vmul(nq, &inarray[2][0], 1, &tmp3[0], 1, &tmp3[0], 1);
            
            // ====================================================
            
            // tmp3 = (gamma-1)/2*velsq;
            Vmath::Smul(nq, HalfGammaMinOne, &velsq[0], 1, &tmp4[0], 1);
            
            // tmp3 = ((gamma-1)/2*velsq-H);
            Vmath::Vsub(nq, &tmp4[0], 1, &H[0], 1, &tmp4[0], 1);
            
            // tmp3 = v((gamma-1)/2*velsq-H);
            Vmath::Vmul(nq,
                        &vel[2][0], 1,
                        &tmp4[0], 1,
                        &tmp4[0], 1);
            
            // tmp3 = z_rhoE*u((gamma-1)/2*velsq-H);
            Vmath::Vmul(nq, &inarray[3][0], 1, &tmp4[0], 1, &tmp4[0], 1);
            
            Vmath::Vadd(nq, &tmp[0], 1, &tmp2[0], 1, &outarray[0][2][0], 1);
            Vmath::Vadd(nq, &outarray[0][2][0], 1, &tmp3[0], 1, &outarray[0][2][0], 1);
            Vmath::Vadd(nq, &outarray[0][2][0], 1, &tmp4[0], 1, &outarray[0][2][0], 1);
            
            // ====================================================
            
            Vmath::Zero(nq, &tmp[0], 1);
            Vmath::Zero(nq, &tmp1[0], 1);
            Vmath::Zero(nq, &tmp2[0], 1);
            Vmath::Zero(nq, &tmp3[0], 1);
            Vmath::Zero(nq, &tmp4[0], 1);
            
            // tmp1 = -(gamma-3) * w
            Vmath::Smul(nq, ThreeMinGam, &vel[2][0], 1, &tmp1[0], 1);
            // tmp1 = -(gamma-3) * w * z_rhov
            Vmath::Vmul(nq, &tmp1[0], 1, &inarray[3][0], 1, &tmp1[0], 1);
            
            //====================================================
            // tmp2 = u
            Vmath::Vcopy(nq, &vel[0][0], 1, &tmp2[0], 1);
            
            // tmp2 = u * z_rhou
            Vmath::Vmul(nq, &tmp2[0], 1, &inarray[1][0], 1, &tmp2[0], 1);
            
            //====================================================
            // tmp3 = v
            Vmath::Vcopy(nq, &vel[1][0], 1, &tmp3[0], 1);
            
            // tmp3 = v * z_rhov
            Vmath::Vmul(nq, &tmp3[0], 1, &inarray[3][0], 1, &tmp3[0], 1);
            
            //====================================================
            // tmp4 = w^2
            Vmath::Vmul(nq, &vel[2][0], 1, &vel[2][0], 1, &tmp4[0], 1);
            
            // tmp4 = (1-gamma) * w^2
            Vmath::Smul(nq, -GammaMinOne, &tmp4[0], 1, &tmp4[0], 1);
            
            // tmp4 = (gamma-1) * w^2+H
            Vmath::Vadd(nq, &tmp4[0], 1, &H[0], 1, &tmp4[0], 1);
            
            // tmp4 = ((gamma-1) * w^2+H)*z_rhoE
            Vmath::Vmul(nq, &inarray[4][0], 1, &tmp4[0], 1, &tmp4[0], 1);
            
            Vmath::Vadd(nq, &tmp[0], 1, &tmp1[0], 1, &outarray[2][1][0], 1);
            Vmath::Vadd(nq, &outarray[3][2][0], 1, &tmp2[0], 1, &outarray[3][2][0], 1);
            Vmath::Vadd(nq, &outarray[3][2][0], 1, &tmp3[0], 1, &outarray[3][2][0], 1);
            Vmath::Vadd(nq, &outarray[3][2][0], 1, &tmp4[0], 1, &outarray[3][2][0], 1);
            // ====================================================
            
            Vmath::Zero(nq, &tmp[0], 1);
            Vmath::Zero(nq, &tmp1[0], 1);
            Vmath::Zero(nq, &tmp2[0], 1);
            Vmath::Zero(nq, &tmp3[0], 1);
            
            // tmp = -(gamma-1) * u
            Vmath::Smul(nq, -GammaMinOne, &vel[0][0], 1, &tmp[0], 1);
            
            // tmp = -(gamma-1) * u * z_rhow
            Vmath::Vmul(nq, &tmp[0], 1, &inarray[3][0], 1, &tmp[0], 1);
            
            // tmp1 = w * z_rhou
            Vmath::Vmul(nq, &vel[2][0], 1, &inarray[1][0], 1, &tmp1[0], 1);
            
            // tmp2 =  u * w
            Vmath::Vmul(nq, &vel[0][0], 1, &vel[2][0], 1, &tmp2[0], 1);
            
            // tmp2 = -(gamma-1)* u * w
            Vmath::Smul(nq, -GammaMinOne, &tmp2[0], 1, &tmp2[0], 1);
            
            // tmp2 = -(gamma-1)* u * w * z_rhoE
            Vmath::Vmul(nq, &inarray[3][0], 1, &tmp2[0], 1, &tmp2[0], 1);
            Vmath::Vadd(nq, &tmp[0], 1, &tmp1[0], 1, &outarray[1][2][0], 1);
            Vmath::Vadd(nq, &outarray[1][2][0], 1, &tmp2[0], 1, &outarray[1][2][0], 1);
            
            // ====================================================
            // ====================================================
            
            Vmath::Zero(nq, &tmp[0], 1);
            Vmath::Zero(nq, &tmp1[0], 1);
            Vmath::Zero(nq, &tmp2[0], 1);
            Vmath::Zero(nq, &tmp3[0], 1);
            
            // tmp = -(gamma-1) * v
            Vmath::Smul(nq, -GammaMinOne, &vel[1][0], 1, &tmp[0], 1);
            
            // tmp = -(gamma-1) * v * z_rhow
            Vmath::Vmul(nq, &tmp[0], 1, &inarray[3][0], 1, &tmp[0], 1);
            
            // tmp1 = w * z_rhov
            Vmath::Vmul(nq, &vel[2][0], 1, &inarray[2][0], 1, &tmp1[0], 1);
            
            // tmp2 =  v * w
            Vmath::Vmul(nq, &vel[1][0], 1, &vel[2][0], 1, &tmp2[0], 1);
            
            // tmp2 = -(gamma-1)* v * w
            Vmath::Smul(nq, -GammaMinOne, &tmp2[0], 1, &tmp2[0], 1);
            
            // tmp2 = -(gamma-1)* v * w * z_rhoE
            Vmath::Vmul(nq, &inarray[4][0], 1, &tmp2[0], 1, &tmp2[0], 1);
            
            Vmath::Vadd(nq, &tmp[0], 1, &tmp1[0], 1, &outarray[2][2][0], 1);
            Vmath::Vadd(nq, &outarray[2][2][0], 1, &tmp2[0], 1, &outarray[2][2][0], 1);
            
            //===================================================
            Vmath::Zero(nq, &tmp[0], 1);
            Vmath::Zero(nq, &tmp1[0], 1);
            Vmath::Zero(nq, &tmp2[0], 1);
            Vmath::Zero(nq, &tmp3[0], 1);
            
            Vmath::Smul(nq, GammaMinOne, &inarray[3][0], 1, &tmp[0], 1);
            Vmath::Smul(nq, m_gamma, &vel[1][0], 1, &tmp1[0], 1);
            Vmath::Vmul(nq, &tmp1[0], 1, &inarray[4][0], 1, &tmp1[0], 1);
            
            
            Vmath::Vadd(nq, &tmp[0], 1, &tmp1[0], 1, &outarray[4][2][0], 1);
        }
    }
    
    void CompressibleFlowSystem::GetDirectSolution(
                            Array<OneD, Array<OneD, NekDouble> > &directSol)
    {
        int nConvectiveFields = directSol.num_elements();
        int npts = GetTotPoints();
        int nCoeffs = m_primal[0]->GetCoeffs().num_elements();
        int npts_dir = directSol[0].num_elements();
        
        Array<OneD, Array<OneD, NekDouble> > directSol2(nConvectiveFields);
        
        for (int i = 0; i < nConvectiveFields; ++i)
        {
            m_primal[i]->BwdTrans(m_primal[i]->GetCoeffs(), directSol[i]);
        }
    }
    
    void CompressibleFlowSystem::GetFwdBwdDirectSolution(
                    Array<OneD, Array<OneD, NekDouble> > &FwdDir,
                    Array<OneD, Array<OneD, NekDouble> > &BwdDir)
    {
        int nConvectiveFields = FwdDir.num_elements();
        Array<OneD, Array<OneD, NekDouble> > directSol(nConvectiveFields);
        
        Array<OneD, Array<OneD, NekDouble> > absBwdDir(nConvectiveFields);
        Array<OneD, Array<OneD, NekDouble> > absFwdDir(nConvectiveFields);
        for (int i = 0; i < nConvectiveFields; ++i)
        {
            
            directSol[i] = Array<OneD, NekDouble>(m_primal[0]->GetPhys().num_elements());
            
            m_primal[i]->BwdTrans(m_primal[i]->GetCoeffs(), directSol[i]);
            
            m_primal[i]->GetFwdBwdTracePhys(directSol[i],
                                            FwdDir[i],
                                            BwdDir[i]);
        }
        
        // correct for boundary conditions // HACK so check for better implementation
        //int bcRegion = 0;
        int eMax     = 0;
        
        // =============================HACK!!!!================================
        int cnt = 0;
        const Array<OneD, const int> &traceBndMap =
        m_fields[0]->GetTraceBndMap();
        
        int nBCregions = m_fields[0]->GetBndConditions().num_elements();
        
        for (int k = 0; k < nBCregions; ++k)
        {
            if (m_fields[0]->GetBndConditions()[k]->GetUserDefined() ==
                SpatialDomains::eAdjointWall ||
                m_fields[0]->GetBndConditions()[k]->GetUserDefined() ==
                SpatialDomains::eWallViscous)
            {
                eMax = m_fields[0]->GetBndCondExpansions()[k]->GetExpSize();
                
                for (int e = 0; e < eMax; ++e)
                {
                    int nBCEdgePts = m_primal[0]->GetBndCondExpansions()[k]->
                    GetExp(e)->GetTotPoints();
                    int id1 = m_primal[0]->GetBndCondExpansions()[k]->
                    GetPhys_Offset(e);
                    int id2 = m_primal[0]->GetTrace()->GetPhys_Offset(traceBndMap[cnt+e]);
                    
                    Array<OneD, NekDouble> tmp(nBCEdgePts, 0.0);
                    Array<OneD, NekDouble> tmp1(nBCEdgePts, 0.0);
                    Array<OneD, NekDouble> tmpfwdnew(nBCEdgePts, 0.0);
                    Array<OneD, NekDouble> tmpfwd(nBCEdgePts, 0.0);
                    
                    //==========================================================
                    // For Euler
                    /*// Calculate (v.n) and (phi.n)
                     for (int i = 0; i < m_spacedim; ++i)
                     {
                     Vmath::Vvtvp(nBCEdgePts,
                     &FwdDir[1+i][id2], 1,
                     &m_traceNormals[i][id2], 1,
                     &tmp[0], 1,
                     &tmp[0], 1);
                     }
                     
                     // Calculate 2.0(v.n)
                     Vmath::Smul(nBCEdgePts, -2.0, &tmp[0], 1, &tmp1[0], 1);
                     
                     // Calculate v* = v - 2.0(v.n)n + 1/Cinf(phi.n)n
                     for (int i = 0; i < m_spacedim; ++i)
                     {
                     Vmath::Vvtvp(nBCEdgePts,
                     &tmp1[0], 1,
                     &m_traceNormals[i][id2], 1,
                     &BwdDir[1+i][id2], 1,
                     &BwdDir[1+i][id2], 1);
                     }*/
                    
                    //==========================================================
                    // For Navier Stokes
                    Vmath::Vcopy(nBCEdgePts,
                                 &FwdDir[0][id2], 1,
                                 &BwdDir[0][id2], 1);
                   	
                    Vmath::Vcopy(nBCEdgePts,
                                 &FwdDir[1][id2], 1,
                                 &BwdDir[1][id2], 1);
                    
                    Vmath::Vcopy(nBCEdgePts,
                                 &FwdDir[2][id2], 1,
                                 &BwdDir[2][id2], 1);
                    
                    for (int i = 0; i < m_spacedim; i++)
                    {
                        Vmath::Neg(nBCEdgePts, &BwdDir[i+1][id2], 1);
                    }
                    
                    // T = Twall at the wall hence,
                    // P = rho * R * T
                    // => T = P/(rho * R)
                    // => T = ((gamma - 1) * (rhoE - rho * q^2)/(rho * R)
                    // => rhoE = (T * (rho * R) / (gamma - 1) + rho * q^2
                    
                    NekDouble TwallFac =
                    (m_Twall * m_gasConstant) / (m_gamma - 1.0);
                    
                    Array<OneD, NekDouble> TwallVec (nBCEdgePts, TwallFac);
                    Array<OneD, NekDouble> qterm    (nBCEdgePts,     0.0);
                    
                    for (int i = 0; i < m_spacedim; ++i)
                    {
                        Vmath::Vvtvp(nBCEdgePts,
                                     &BwdDir[i+1][id2], 1,
                                     &BwdDir[i+1][id2], 1,
                                     &qterm[0],         1,
                                     &qterm[0],         1);
                    }
                    
                    Vmath::Vdiv(nBCEdgePts,
                                &qterm[0], 1,
                                &BwdDir[0][id2], 1,
                                &BwdDir[m_spacedim+1][id2], 1);
                    
                    Vmath::Vadd(nBCEdgePts,
                                &TwallVec[0], 1,
                                &BwdDir[0][id2], 1,
                                &TwallVec[0], 1);
                    
                    Vmath::Vadd(nBCEdgePts,
                                &TwallVec[0], 1,
                                &BwdDir[m_spacedim+1][id2], 1,
                                &BwdDir[m_spacedim+1][id2], 1);
                    
                    //==========================================================
                    
                    for (int i = 0; i < nConvectiveFields; ++i)
                    {
                        Vmath::Vcopy(nBCEdgePts, &BwdDir[i][id2], 1,
                                     &(m_primal[i]->GetBndCondExpansions()[k]->
                                       UpdatePhys())[id1], 1);
                    }
                }
                
                cnt += m_fields[0]->GetBndCondExpansions()[k]->GetExpSize();
            
            }
            
            if (m_fields[0]->GetBndConditions()[k]->GetUserDefined() !=
                SpatialDomains::eAdjointWall &&
                m_fields[0]->GetBndConditions()[k]->GetUserDefined() !=
                SpatialDomains::eWallViscous)
            {
                eMax = m_fields[0]->GetBndCondExpansions()[k]->GetExpSize();
                for (int e = 0; e < eMax; ++e)
                {
                    int nBCEdgePts = m_primal[0]->GetBndCondExpansions()[k]->
                    GetExp(e)->GetTotPoints();
                    int id1 = m_primal[0]->GetBndCondExpansions()[k]->
                    GetPhys_Offset(e);
                    int id2 = m_primal[0]->GetTrace()->GetPhys_Offset(traceBndMap[cnt+e]);
                    
                    
                    for (int i = 0; i < nConvectiveFields; ++i)
                    {
                        Vmath::Vcopy(nBCEdgePts,
                                     &FwdDir[i][id2], 1,
                                     &BwdDir[i][id2], 1);
                        
                        Vmath::Vcopy(nBCEdgePts, &FwdDir[i][id2], 1,
                                     &(m_primal[i]->GetBndCondExpansions()[k]->
                                       UpdatePhys())[id1], 1);
                    }
                }
                
                cnt += m_fields[0]->GetBndCondExpansions()[k]->GetExpSize();
                
            }
        }
    }
    
    void CompressibleFlowSystem::GetFwdBwdDIFFDirectSolution(
                Array<OneD, Array<OneD, NekDouble> > &FwdDir,
                Array<OneD, Array<OneD, NekDouble> > &BwdDir,
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &FwdDirDIFF,
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &BwdDirDIFF)
    {
        // FwdDirDIFF[m_spacedim][nFields][nTracePoints]
        // BwdDirDIFF[m_spacedim][nFields][nTracePoints]
        // rho u v w p;
        int nConvectiveFields = FwdDir.num_elements();
        //int nTracePoints       = FwdDir[0].num_elements();
        int i;
        
        Array<OneD, Array<OneD, NekDouble> > directSol(nConvectiveFields);
        
        Array<OneD, Array<OneD, NekDouble> > absBwdDir(nConvectiveFields);
        Array<OneD, Array<OneD, NekDouble> > absFwdDir(nConvectiveFields);
        
        for (int i = 0; i < nConvectiveFields; ++i)
        {
            
            directSol[i] = Array<OneD, NekDouble>(m_primal[0]->GetPhys().num_elements());
            
            m_primal[i]->BwdTrans(m_primal[i]->GetCoeffs(), directSol[i]);
            
            
        }
        int nq = directSol[0].num_elements();
        
        Array<OneD, NekDouble> Pressure(nq, 0.0);
        
        Array<OneD, NekDouble> drho_dX(nq, 0.0);
        Array<OneD, NekDouble> drho_dY(nq, 0.0);
        
        Array<OneD, NekDouble> dP_dX(nq, 0.0);
        Array<OneD, NekDouble> dP_dY(nq, 0.0);
        
        Array<OneD, Array<OneD, NekDouble> > vel(m_spacedim);
        
        m_primal[0]->PhysDeriv(directSol[0], drho_dX,
                               drho_dY);
        
        m_primal[0]->GetFwdBwdTracePhys(drho_dX,
                                        FwdDirDIFF[0][0],
                                        BwdDirDIFF[0][0]);
        
        m_primal[0]->GetFwdBwdTracePhys(drho_dY,
                                        FwdDirDIFF[1][0],
                                        BwdDirDIFF[1][0]);
        
        
        for (i = 0; i < m_spacedim; ++i)
        {
            Array<OneD, NekDouble> tmpx(nq, 0.0);
            Array<OneD, NekDouble> tmpy(nq, 0.0);
            
            vel[i] = Array<OneD, NekDouble> (nq, 0.0);
            
            Vmath::Vdiv(nq, directSol[i+1], 1, directSol[0], 1, vel[i], 1);
            
            m_primal[0]->PhysDeriv(vel[i],
                                   tmpx,
                                   tmpy);
            
            m_primal[0]->GetFwdBwdTracePhys(tmpx,
                                            FwdDirDIFF[0][i+1],
                                            BwdDirDIFF[0][i+1]);
            
            m_primal[0]->GetFwdBwdTracePhys(tmpy,
                                            FwdDirDIFF[1][i+1],
                                            BwdDirDIFF[1][i+1]);
        }
        
        GetPressure(directSol, vel, Pressure);
        
        m_primal[0]->PhysDeriv(Pressure,
                               dP_dX,
                               dP_dY);
        
        m_primal[0]->GetFwdBwdTracePhys(dP_dX,
                                        FwdDirDIFF[0][m_spacedim+1],
                                        BwdDirDIFF[0][m_spacedim+1]);
        
        m_primal[0]->GetFwdBwdTracePhys(dP_dY,
                                        FwdDirDIFF[1][m_spacedim+1],
                                        BwdDirDIFF[1][m_spacedim+1]);
    }
    
    /**
     * @brief Return the flux vector for the compressible Euler equations
     * by using the de-aliasing technique.
     *
     * @param i           Component of the flux vector to calculate.
     * @param physfield   Fields.
     * @param flux        Resulting flux.
     */
    void CompressibleFlowSystem::GetFluxVectorDeAlias(
            const Array<OneD, Array<OneD, NekDouble> >               &physfield,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &flux)
    {
        int i, j;
        int nq = physfield[0].num_elements();
        int nVariables = m_fields.num_elements();
        
        // Factor to rescale 1d points in dealiasing
        NekDouble OneDptscale = 2;
        nq = m_fields[0]->Get1DScaledTotPoints(OneDptscale);
        
        Array<OneD, NekDouble> pressure(nq);
        Array<OneD, Array<OneD, NekDouble> > velocity(m_spacedim);
        
        Array<OneD, Array<OneD, NekDouble> > physfield_interp(nVariables);
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > > flux_interp(
                                                                       nVariables);
        
        for (i = 0; i < nVariables; ++ i)
        {
            physfield_interp[i] = Array<OneD, NekDouble>(nq);
            flux_interp[i] = Array<OneD, Array<OneD, NekDouble> >(m_spacedim);
            m_fields[0]->PhysInterp1DScaled(
                                            OneDptscale, physfield[i], physfield_interp[i]);
            
            for (j = 0; j < m_spacedim; ++j)
            {
                flux_interp[i][j] = Array<OneD, NekDouble>(nq);
            }
        }
        
        // Flux vector for the rho equation
        for (i = 0; i < m_spacedim; ++i)
        {
            velocity[i] = Array<OneD, NekDouble>(nq);
            
            // Galerkin project solution back to original space
            m_fields[0]->PhysGalerkinProjection1DScaled(
                                                        OneDptscale, physfield_interp[i+1], flux[0][i]);
        }
        
        GetVelocityVector(physfield_interp, velocity);
        GetPressure      (physfield_interp, velocity, pressure);
        
        // Evaluation of flux vector for the velocity fields
        for (i = 0; i < m_spacedim; ++i)
        {
            for (j = 0; j < m_spacedim; ++j)
            {
                Vmath::Vmul(nq, velocity[j], 1, physfield_interp[i+1], 1,
                            flux_interp[i+1][j], 1);
            }
            
            // Add pressure to appropriate field
            Vmath::Vadd(nq, flux_interp[i+1][i], 1, pressure,1,
                        flux_interp[i+1][i], 1);
        }
        
        // Galerkin project solution back to origianl space
        for (i = 0; i < m_spacedim; ++i)
        {
            for (j = 0; j < m_spacedim; ++j)
            {
                m_fields[0]->PhysGalerkinProjection1DScaled(
                                                            OneDptscale, flux_interp[i+1][j], flux[i+1][j]);
            }
        }
        
        // Evaluation of flux vector for energy
        Vmath::Vadd(nq, physfield_interp[m_spacedim+1], 1, pressure, 1,
                    pressure, 1);
        
        for (j = 0; j < m_spacedim; ++j)
        {
            Vmath::Vmul(nq, velocity[j], 1, pressure, 1,
                        flux_interp[m_spacedim+1][j], 1);
            
            // Galerkin project solution back to origianl space
            m_fields[0]->PhysGalerkinProjection1DScaled(
                                                        OneDptscale,
                                                        flux_interp[m_spacedim+1][j],
                                                        flux[m_spacedim+1][j]);
        }
    }
    
    /**
     * @brief Return the flux vector for the LDG diffusion problem.
     * \todo Complete the viscous flux vector
     */
    void CompressibleFlowSystem::GetViscousFluxVector(
        const Array<OneD, Array<OneD, NekDouble> >               &physfield,
              Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &derivativesO1,
              Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &viscousTensor)
    {
        int j, k;
        int nVariables = m_fields.num_elements();
        int nPts       = physfield[0].num_elements();
        
        // Stokes hypotesis
        const NekDouble lambda = -2.0/3.0;
        
        // Auxiliary variables
        Array<OneD, NekDouble > mu                 (nPts, 0.0);
        Array<OneD, NekDouble > thermalConductivity(nPts, 0.0);
        Array<OneD, NekDouble > mu2                (nPts, 0.0);
        Array<OneD, NekDouble > divVel             (nPts, 0.0);

        // Variable viscosity through the Sutherland's law
        if (m_ViscosityType == "Variable")
        {
            GetDynamicViscosity(physfield[nVariables-2], mu);
            NekDouble tRa = m_Cp / m_Prandtl;
            Vmath::Smul(nPts, tRa, &mu[0], 1, &thermalConductivity[0], 1);
        }
        else
        {
            Vmath::Fill(nPts, m_mu, &mu[0], 1);
            Vmath::Fill(nPts, m_thermalConductivity,
                        &thermalConductivity[0], 1);
        }
        // Computing diagonal terms of viscous stress tensor
        Array<OneD, Array<OneD, NekDouble> > tmp(m_spacedim);
        Array<OneD, Array<OneD, NekDouble> > Sgg(m_spacedim);
        
        // mu2 = 2 * mu
        Vmath::Smul(nPts, 2.0, &mu[0], 1, &mu2[0], 1);
        
        // Velocity divergence
        for (j = 0; j < m_spacedim; ++j)
        {
            Vmath::Vadd(nPts, &divVel[0], 1, &derivativesO1[j][j][0], 1,
                        &divVel[0], 1);
        }
        
        // Velocity divergence scaled by lambda * mu
        Vmath::Smul(nPts, lambda, &divVel[0], 1, &divVel[0], 1);
        Vmath::Vmul(nPts, &mu[0], 1, &divVel[0], 1, &divVel[0], 1);
        
        // Diagonal terms of viscous stress tensor (Sxx, Syy, Szz)
        // Sjj = 2 * mu * du_j/dx_j - (2 / 3) * mu * sum_j(du_j/dx_j)
        for (j = 0; j < m_spacedim; ++j)
        {
            tmp[j] = Array<OneD, NekDouble>(nPts, 0.0);
            Sgg[j] = Array<OneD, NekDouble>(nPts, 0.0);
            
            Vmath::Vmul(nPts, &mu2[0], 1, &derivativesO1[j][j][0], 1,
                        &tmp[j][0], 1);
            
            Vmath::Vadd(nPts, &tmp[j][0], 1, &divVel[0], 1, &Sgg[j][0], 1);
        }
        
        // Extra diagonal terms of viscous stress tensor (Sxy, Sxz, Syz)
        // Note: they exist for 2D and 3D problems only
        Array<OneD, NekDouble > Sxy(nPts, 0.0);
        Array<OneD, NekDouble > Sxz(nPts, 0.0);
        Array<OneD, NekDouble > Syz(nPts, 0.0);
        
        if (m_spacedim == 2)
        {
            // Sxy = (du/dy + dv/dx)
            Vmath::Vadd(nPts, &derivativesO1[0][1][0], 1,
                        &derivativesO1[1][0][0], 1, &Sxy[0], 1);
            
            // Sxy = mu * (du/dy + dv/dx)
            Vmath::Vmul(nPts, &mu[0], 1, &Sxy[0], 1, &Sxy[0], 1);
        }
        else if (m_spacedim == 3)
        {
            // Sxy = (du/dy + dv/dx)
            Vmath::Vadd(nPts, &derivativesO1[0][1][0], 1,
                        &derivativesO1[1][0][0], 1, &Sxy[0], 1);
            
            // Sxz = (du/dz + dw/dx)
            Vmath::Vadd(nPts, &derivativesO1[0][2][0], 1,
                        &derivativesO1[2][0][0], 1, &Sxz[0], 1);
            
            // Syz = (dv/dz + dw/dy)
            Vmath::Vadd(nPts, &derivativesO1[1][2][0], 1,
                        &derivativesO1[2][1][0], 1, &Syz[0], 1);
            
            // Sxy = mu * (du/dy + dv/dx)
            Vmath::Vmul(nPts, &mu[0], 1, &Sxy[0], 1, &Sxy[0], 1);
            
            // Sxz = mu * (du/dy + dv/dx)
            Vmath::Vmul(nPts, &mu[0], 1, &Sxz[0], 1, &Sxz[0], 1);
            
            // Syz = mu * (du/dy + dv/dx)
            Vmath::Vmul(nPts, &mu[0], 1, &Syz[0], 1, &Syz[0], 1);
        }
        
        // Energy-related terms
        Array<OneD, NekDouble > STx(nPts, 0.0);
        Array<OneD, NekDouble > STy(nPts, 0.0);
        Array<OneD, NekDouble > STz(nPts, 0.0);
        
        // Building the viscous flux vector
        
        // Viscous flux vector for the rho equation
        for (k = 0; k < m_spacedim; ++k)
        {
            Vmath::Zero(nPts, viscousTensor[k][0], 1);
        }
        
        if (m_spacedim == 1)
        {
            Array<OneD, NekDouble > tmp1(nPts, 0.0);
            
            // u * Sxx
            Vmath::Vmul(nPts, &physfield[0][0], 1, &Sgg[0][0], 1, &STx[0], 1);
            
            // k * dT/dx
            Vmath::Vmul(nPts, &thermalConductivity[0], 1,
                        &derivativesO1[0][1][0], 1, &tmp1[0], 1);
            // STx = u * Sxx + (K / mu) * dT/dx
            Vmath::Vadd(nPts, &STx[0], 1, &tmp1[0], 1, &STx[0], 1);
        }
        else if (m_spacedim == 2)
        {
            Array<OneD, NekDouble > tmp1(nPts, 0.0);
            Array<OneD, NekDouble > tmp2(nPts, 0.0);
            
            // Computation of STx
            
            // u * Sxx
            Vmath::Vmul(nPts, &physfield[0][0], 1, &Sgg[0][0], 1, &STx[0], 1);
            
            // v * Sxy
            Vmath::Vmul(nPts, &physfield[1][0], 1, &Sxy[0], 1, &tmp1[0], 1);
            
            // k * dT/dx

            Vmath::Vmul(nPts, &thermalConductivity[0], 1,
                        &derivativesO1[0][2][0], 1, &tmp2[0], 1);

            // STx = u * Sxx + v * Sxy + K * dT/dx
            Vmath::Vadd(nPts, &STx[0], 1, &tmp1[0], 1, &STx[0], 1);
            Vmath::Vadd(nPts, &STx[0], 1, &tmp2[0], 1, &STx[0], 1);
            
            // Computation of STy
            
            // Re-initialise temporary arrays
            Vmath::Zero(nPts, &tmp1[0], 1);
            Vmath::Zero(nPts, &tmp2[0], 1);
            
            // v * Syy
            Vmath::Vmul(nPts, &physfield[1][0], 1, &Sgg[1][0], 1, &STy[0], 1);
            
            // u * Sxy
            Vmath::Vmul(nPts, &physfield[0][0], 1, &Sxy[0], 1, &tmp1[0], 1);
            
            // k * dT/dy

            Vmath::Vmul(nPts, &thermalConductivity[0], 1,
                        &derivativesO1[1][2][0], 1, &tmp2[0], 1);

            // STy = v * Syy + u * Sxy + K * dT/dy
            Vmath::Vadd(nPts, &STy[0], 1, &tmp1[0], 1, &STy[0], 1);
            Vmath::Vadd(nPts, &STy[0], 1, &tmp2[0], 1, &STy[0], 1);
        }
        else if (m_spacedim == 3)
        {
            Array<OneD, NekDouble > tmp1(nPts, 0.0);
            Array<OneD, NekDouble > tmp2(nPts, 0.0);
            Array<OneD, NekDouble > tmp3(nPts, 0.0);
            
            // Computation of STx
            
            // u * Sxx
            Vmath::Vmul(nPts, &physfield[0][0], 1, &Sgg[0][0], 1, &STx[0], 1);
            
            // v * Sxy
            Vmath::Vmul(nPts, &physfield[1][0], 1, &Sxy[0], 1, &tmp1[0], 1);
            
            // v * Sxz
            Vmath::Vmul(nPts, &physfield[2][0], 1, &Sxz[0], 1, &tmp2[0], 1);
            
            // k * dT/dx
            Vmath::Vmul(nPts, &thermalConductivity[0], 1,
                        &derivativesO1[0][3][0], 1, &tmp3[0], 1);
            
            // STx = u * Sxx + v * Sxy + w * Sxz + (K / mu) * dT/dx
            Vmath::Vadd(nPts, &STx[0], 1, &tmp1[0], 1, &STx[0], 1);
            Vmath::Vadd(nPts, &STx[0], 1, &tmp2[0], 1, &STx[0], 1);
            Vmath::Vadd(nPts, &STx[0], 1, &tmp3[0], 1, &STx[0], 1);
            
            // Computation of STy
            
            // Re-initialise temporary arrays
            Vmath::Zero(nPts, &tmp1[0], 1);
            Vmath::Zero(nPts, &tmp2[0], 1);
            Vmath::Zero(nPts, &tmp3[0], 1);
            
            // v * Syy
            Vmath::Vmul(nPts, &physfield[1][0], 1, &Sgg[1][0], 1, &STy[0], 1);
            
            // u * Sxy
            Vmath::Vmul(nPts, &physfield[0][0], 1, &Sxy[0], 1, &tmp1[0], 1);
            
            // w * Syz
            Vmath::Vmul(nPts, &physfield[2][0], 1, &Syz[0], 1, &tmp2[0], 1);
            
            // k * dT/dy
            Vmath::Vmul(nPts, &thermalConductivity[0], 1,
                        &derivativesO1[1][3][0], 1, &tmp3[0], 1);

            // STy = v * Syy + u * Sxy + w * Syz + K * dT/dy
            Vmath::Vadd(nPts, &STy[0], 1, &tmp1[0], 1, &STy[0], 1);
            Vmath::Vadd(nPts, &STy[0], 1, &tmp2[0], 1, &STy[0], 1);
            Vmath::Vadd(nPts, &STy[0], 1, &tmp3[0], 1, &STy[0], 1);
            
            // Computation of STz
            
            // Re-initialise temporary arrays
            Vmath::Zero(nPts, &tmp1[0], 1);
            Vmath::Zero(nPts, &tmp2[0], 1);
            Vmath::Zero(nPts, &tmp3[0], 1);
            
            // w * Szz
            Vmath::Vmul(nPts, &physfield[2][0], 1, &Sgg[2][0], 1, &STz[0], 1);
            
            // u * Sxz
            Vmath::Vmul(nPts, &physfield[0][0], 1, &Sxz[0], 1, &tmp1[0], 1);
            
            // v * Syz
            Vmath::Vmul(nPts, &physfield[1][0], 1, &Syz[0], 1, &tmp2[0], 1);
            
            // k * dT/dz
            Vmath::Vmul(nPts, &thermalConductivity[0], 1,
                        &derivativesO1[2][3][0], 1, &tmp3[0], 1);

            // STz = w * Szz + u * Sxz + v * Syz + K * dT/dz
            Vmath::Vadd(nPts, &STz[0], 1, &tmp1[0], 1, &STz[0], 1);
            Vmath::Vadd(nPts, &STz[0], 1, &tmp2[0], 1, &STz[0], 1);
            Vmath::Vadd(nPts, &STz[0], 1, &tmp3[0], 1, &STz[0], 1);
        }
        
        switch (m_spacedim)
        {
            case 1:
            {
                // f_11v = f_rho = 0
                Vmath::Zero(nPts, &viscousTensor[0][0][0], 1);
                
                // f_21v = f_rhou
                Vmath::Vcopy(nPts, &Sgg[0][0], 1, &viscousTensor[0][1][0], 1);
                
                // f_31v = f_E
                Vmath::Vcopy(nPts, &STx[0], 1, &viscousTensor[0][2][0], 1);
                break;
            }
            case 2:
            {
                // f_11v = f_rho1 = 0
                Vmath::Zero(nPts, &viscousTensor[0][0][0], 1);
                // f_12v = f_rho2 = 0
                Vmath::Zero(nPts, &viscousTensor[1][0][0], 1);
                
                // f_21v = f_rhou1
                Vmath::Vcopy(nPts, &Sgg[0][0], 1, &viscousTensor[0][1][0], 1);
                // f_22v = f_rhou2
                Vmath::Vcopy(nPts, &Sxy[0],    1, &viscousTensor[1][1][0], 1);
                
                // f_31v = f_rhov1
                Vmath::Vcopy(nPts, &Sxy[0],    1, &viscousTensor[0][2][0], 1);
                // f_32v = f_rhov2
                Vmath::Vcopy(nPts, &Sgg[1][0], 1, &viscousTensor[1][2][0], 1);
                
                // f_41v = f_E1
                Vmath::Vcopy(nPts, &STx[0], 1, &viscousTensor[0][3][0], 1);
                // f_42v = f_E2
                Vmath::Vcopy(nPts, &STy[0], 1, &viscousTensor[1][3][0], 1);
                break;
            }
            case 3:
            {
                // f_11v = f_rho1 = 0
                Vmath::Zero(nPts, &viscousTensor[0][0][0], 1);
                // f_12v = f_rho2 = 0
                Vmath::Zero(nPts, &viscousTensor[1][0][0], 1);
                // f_13v = f_rho3 = 0
                Vmath::Zero(nPts, &viscousTensor[2][0][0], 1);
                
                // f_21v = f_rhou1
                Vmath::Vcopy(nPts, &Sgg[0][0], 1, &viscousTensor[0][1][0], 1);
                // f_22v = f_rhou2
                Vmath::Vcopy(nPts, &Sxy[0],    1, &viscousTensor[1][1][0], 1);
                // f_23v = f_rhou3
                Vmath::Vcopy(nPts, &Sxz[0],    1, &viscousTensor[2][1][0], 1);
                
                // f_31v = f_rhov1
                Vmath::Vcopy(nPts, &Sxy[0],    1, &viscousTensor[0][2][0], 1);
                // f_32v = f_rhov2
                Vmath::Vcopy(nPts, &Sgg[1][0], 1, &viscousTensor[1][2][0], 1);
                // f_33v = f_rhov3
                Vmath::Vcopy(nPts, &Syz[0],    1, &viscousTensor[2][2][0], 1);
                
                // f_31v = f_rhow1
                Vmath::Vcopy(nPts, &Sxz[0],    1, &viscousTensor[0][3][0], 1);
                // f_32v = f_rhow2
                Vmath::Vcopy(nPts, &Syz[0],    1, &viscousTensor[1][3][0], 1);
                // f_33v = f_rhow3
                Vmath::Vcopy(nPts, &Sgg[2][0], 1, &viscousTensor[2][3][0], 1);
                
                // f_41v = f_E1
                Vmath::Vcopy(nPts, &STx[0], 1, &viscousTensor[0][4][0], 1);
                // f_42v = f_E2
                Vmath::Vcopy(nPts, &STy[0], 1, &viscousTensor[1][4][0], 1);
                // f_43v = f_E3
                Vmath::Vcopy(nPts, &STz[0], 1, &viscousTensor[2][4][0], 1);
                break;
            }
            default:
            {
                ASSERTL0(false, "Illegal expansion dimension");
            }
        }
    }
    
    /**
     * @brief Return the flux vector for the LDG diffusion problem.
     * \todo Complete the viscous flux vector
     */
    void CompressibleFlowSystem::GetViscousFluxVectorDeAlias(
                                                             const Array<OneD, Array<OneD, NekDouble> >               &physfield,
                                                             Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &derivativesO1,
                                                             Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &viscousTensor)
    {
#if 0
        int i, j, k;
        int nVariables = m_fields.num_elements();
        int nPts       = physfield[0].num_elements();
        
        int variables_phys = physfield.num_elements();
        
        // Factor to rescale 1d points in dealiasing.
        NekDouble OneDptscale = 2;
        
        // Get number of points to dealias a cubic non-linearity
        nPts = m_fields[0]->Get1DScaledTotPoints(OneDptscale);
        
        int nVariables_aux = derivativesO1[0].num_elements();
        
        Array<OneD, Array<OneD, NekDouble> > physfield_interp(variables_phys);
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > > derivativesO1_interp(
                                                                                m_spacedim);
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > >viscousTensor_interp(
                                                                               m_spacedim);
        
        for (i = 0; i < m_spacedim; ++ i)
        {
            viscousTensor_interp[i] = Array<OneD, Array<OneD, NekDouble> >(
                                                                           nVariables);
            for (j = 0; j < nVariables; ++j)
            {
                viscousTensor_interp[i][j] = Array<OneD, NekDouble>(nPts);
            }
        }
        
        // Stokes hypotesis
        NekDouble lambda = -2.0/3.0;
        
        // Auxiliary variables
        Array<OneD, NekDouble > mu         (nPts, 0.0);
        Array<OneD, NekDouble > mu2        (nPts, 0.0);
        Array<OneD, NekDouble > divVel     (nPts, 0.0);
        Array<OneD, NekDouble > pressure   (nPts, 0.0);
        Array<OneD, NekDouble > temperature(nPts, 0.0);
        
        for (i = 0; i < nVariables; ++i)
        {
            m_fields[0]->PhysInterp1DScaled(
                                            OneDptscale, physfield[i], fields_interp[i]);
        }
        
        for (i = 0; i < variables_phys; ++i)
        {
            physfield_interp[i] = Array<OneD, NekDouble> (nPts);
            
            // Interpolation to higher space
            m_fields[0]->PhysInterp1DScaled(OneDptscale, physfield[i],
                                            physfield_interp[i]);
        }
        
        for (i = 0; i < m_spacedim; ++i)
        {
            derivativesO1_interp[i] = Array<OneD, Array<OneD, NekDouble> >(
                                                                           nVariables_aux);
            for (j = 0; j < nVariables_aux; ++j)
            {
                derivativesO1_interp[i][j] = Array<OneD, NekDouble>(nPts);
                m_fields[0]->PhysInterp1DScaled(
                                                OneDptscale, derivativesO1[i][j],
                                                derivativesO1_interp[i][j]);
            }
        }
        
        // Thermodynamic related quantities
        GetPressure(fields_interp, pressure);
        GetTemperature(fields_interp, pressure, temperature);
        
        // Variable viscosity through the Sutherland's law
        if (m_ViscosityType == "Variable")
        {
            GetDynamicViscosity(fields_interp[variables_phys-1], mu);
        }
        else
        {
            Vmath::Sadd(nPts, m_mu, &mu[0], 1, &mu[0], 1);
        }
        
        // Computing diagonal terms of viscous stress tensor
        Array<OneD, Array<OneD, NekDouble> > tmp(m_spacedim);
        Array<OneD, Array<OneD, NekDouble> > Sgg(m_spacedim);
        
        // mu2 = 2 * mu
        Vmath::Smul(nPts, 2.0, &mu[0], 1, &mu2[0], 1);
        
        // Velocity divergence
        for (j = 0; j < m_spacedim; ++j)
        {
            Vmath::Vadd(nPts, &divVel[0], 1, &derivativesO1_interp[j][j][0], 1,
                        &divVel[0], 1);
        }
        
        // Velocity divergence scaled by lambda * mu
        Vmath::Smul(nPts, lambda, &divVel[0], 1, &divVel[0], 1);
        Vmath::Vmul(nPts, &mu[0], 1, &divVel[0], 1, &divVel[0], 1);
        
        // Digonal terms of viscous stress tensor (Sxx, Syy, Szz)
        // Sjj = 2 * mu * du_j/dx_j - (2 / 3) * mu * sum_j(du_j/dx_j)
        for (j = 0; j < m_spacedim; ++j)
        {
            tmp[j] = Array<OneD, NekDouble>(nPts, 0.0);
            Sgg[j] = Array<OneD, NekDouble>(nPts, 0.0);
            
            Vmath::Vmul(nPts, &mu2[0], 1, &derivativesO1_interp[j][j][0], 1,
                        &tmp[j][0], 1);
            
            Vmath::Vadd(nPts, &tmp[j][0], 1, &divVel[0], 1, &Sgg[j][0], 1);
        }
        
        // Extra diagonal terms of viscous stress tensor (Sxy, Sxz, Syz)
        // Note: they exist for 2D and 3D problems only
        Array<OneD, NekDouble > Sxy(nPts, 0.0);
        Array<OneD, NekDouble > Sxz(nPts, 0.0);
        Array<OneD, NekDouble > Syz(nPts, 0.0);
        
        if (m_spacedim == 2)
        {
            // Sxy = (du/dy + dv/dx)
            Vmath::Vadd(nPts, &derivativesO1_interp[0][1][0], 1,
                        &derivativesO1_interp[1][0][0], 1, &Sxy[0], 1);
            
            // Sxy = mu * (du/dy + dv/dx)
            Vmath::Vmul(nPts, &mu[0], 1, &Sxy[0], 1, &Sxy[0], 1);
        }
        else if (m_spacedim == 3)
        {
            // Sxy = (du/dy + dv/dx)
            Vmath::Vadd(nPts, &derivativesO1_interp[0][1][0], 1,
                        &derivativesO1_interp[1][0][0], 1, &Sxy[0], 1);
            
            // Sxz = (du/dz + dw/dx)
            Vmath::Vadd(nPts, &derivativesO1_interp[0][2][0], 1,
                        &derivativesO1_interp[2][0][0], 1, &Sxz[0], 1);
            
            // Syz = (dv/dz + dw/dy)
            Vmath::Vadd(nPts, &derivativesO1_interp[1][2][0], 1,
                        &derivativesO1_interp[2][1][0], 1, &Syz[0], 1);
            
            // Sxy = mu * (du/dy + dv/dx)
            Vmath::Vmul(nPts, &mu[0], 1, &Sxy[0], 1, &Sxy[0], 1);
            
            // Sxz = mu * (du/dy + dv/dx)
            Vmath::Vmul(nPts, &mu[0], 1, &Sxz[0], 1, &Sxz[0], 1);
            
            // Syz = mu * (du/dy + dv/dx)
            Vmath::Vmul(nPts, &mu[0], 1, &Syz[0], 1, &Syz[0], 1);
        }
        
        // Energy-related terms
        Array<OneD, NekDouble > STx(nPts, 0.0);
        Array<OneD, NekDouble > STy(nPts, 0.0);
        Array<OneD, NekDouble > STz(nPts, 0.0);
        // Building the viscous flux vector
        if (i == 0)
        {
            // Viscous flux vector for the rho equation
            for (k = 0; k < m_spacedim; ++k)
            {
                Vmath::Zero(nPts, viscousTensor_interp[k][i], 1);
            }
        }
        
        if (m_spacedim == 1)
        {
            Array<OneD, NekDouble > tmp1(nPts, 0.0);
            
            // u * Sxx
            Vmath::Vmul(nPts, &physfield_interp[0][0], 1,
                        &Sgg[0][0], 1, &STx[0], 1);
            
            // k * dT/dx
            Vmath::Smul(nPts, m_thermalConductivity,
                        &derivativesO1_interp[0][1][0], 1, &tmp1[0], 1);
            
            // STx = u * Sxx + (K / mu) * dT/dx
            Vmath::Vadd(nPts, &STx[0], 1, &tmp1[0], 1, &STx[0], 1);
        }
        else if (m_spacedim == 2)
        {
            Array<OneD, NekDouble > tmp1(nPts, 0.0);
            Array<OneD, NekDouble > tmp2(nPts, 0.0);
            
            // Computation of STx
            
            // u * Sxx
            Vmath::Vmul(nPts, &physfield_interp[0][0], 1,
                        &Sgg[0][0], 1, &STx[0], 1);
            
            // v * Sxy
            Vmath::Vmul(nPts, &physfield_interp[1][0], 1,
                        &Sxy[0], 1, &tmp1[0], 1);
            
            // k * dT/dx
            Vmath::Smul(nPts, m_thermalConductivity,
                        &derivativesO1_interp[0][2][0], 1, &tmp2[0], 1);
            
            // STx = u * Sxx + v * Sxy + K * dT/dx
            Vmath::Vadd(nPts, &STx[0], 1, &tmp1[0], 1, &STx[0], 1);
            Vmath::Vadd(nPts, &STx[0], 1, &tmp2[0], 1, &STx[0], 1);
            
            // Computation of STy
            
            // Re-initialise temporary arrays
            Vmath::Zero(nPts, &tmp1[0], 1);
            Vmath::Zero(nPts, &tmp2[0], 1);
            
            // v * Syy
            Vmath::Vmul(nPts, &physfield_interp[1][0], 1,
                        &Sgg[1][0], 1, &STy[0], 1);
            
            // u * Sxy
            Vmath::Vmul(nPts, &physfield_interp[0][0], 1,
                        &Sxy[0], 1, &tmp1[0], 1);
            
            // k * dT/dy
            Vmath::Smul(nPts, m_thermalConductivity,
                        &derivativesO1_interp[1][2][0], 1, &tmp2[0], 1);
            
            // STy = v * Syy + u * Sxy + K * dT/dy
            Vmath::Vadd(nPts, &STy[0], 1, &tmp1[0], 1, &STy[0], 1);
            Vmath::Vadd(nPts, &STy[0], 1, &tmp2[0], 1, &STy[0], 1);
        }
        else if (m_spacedim == 3)
        {
            Array<OneD, NekDouble > tmp1(nPts, 0.0);
            Array<OneD, NekDouble > tmp2(nPts, 0.0);
            Array<OneD, NekDouble > tmp3(nPts, 0.0);
            
            // Computation of STx
            
            // u * Sxx
            Vmath::Vmul(nPts, &physfield_interp[0][0], 1,
                        &Sgg[0][0], 1, &STx[0], 1);
            
            // v * Sxy
            Vmath::Vmul(nPts, &physfield_interp[1][0], 1,
                        &Sxy[0], 1, &tmp1[0], 1);
            
            // v * Sxy
            Vmath::Vmul(nPts, &physfield_interp[2][0], 1,
                        &Sxz[0], 1, &tmp2[0], 1);
            
            // k * dT/dx
            Vmath::Smul(nPts, m_thermalConductivity,
                        &derivativesO1_interp[0][3][0], 1, &tmp3[0], 1);
            
            // STx = u * Sxx + v * Sxy + w * Sxz + (K / mu) * dT/dx
            Vmath::Vadd(nPts, &STx[0], 1, &tmp1[0], 1, &STx[0], 1);
            Vmath::Vadd(nPts, &STx[0], 1, &tmp2[0], 1, &STx[0], 1);
            Vmath::Vadd(nPts, &STx[0], 1, &tmp3[0], 1, &STx[0], 1);
            
            // Computation of STy
            
            // Re-initialise temporary arrays
            Vmath::Zero(nPts, &tmp1[0], 1);
            Vmath::Zero(nPts, &tmp2[0], 1);
            Vmath::Zero(nPts, &tmp3[0], 1);
            
            // v * Syy
            Vmath::Vmul(nPts, &physfield_interp[1][0], 1,
                        &Sgg[1][0], 1, &STy[0], 1);
            
            // u * Sxy
            Vmath::Vmul(nPts, &physfield_interp[0][0], 1,
                        &Sxy[0], 1, &tmp1[0], 1);
            
            // w * Syz
            Vmath::Vmul(nPts, &physfield_interp[2][0], 1,
                        &Syz[0], 1, &tmp2[0], 1);
            
            // k * dT/dy
            Vmath::Smul(nPts, m_thermalConductivity,
                        &derivativesO1_interp[1][3][0], 1, &tmp3[0], 1);
            
            // STy = v * Syy + u * Sxy + w * Syz + K * dT/dy
            Vmath::Vadd(nPts, &STy[0], 1, &tmp1[0], 1, &STy[0], 1);
            Vmath::Vadd(nPts, &STy[0], 1, &tmp2[0], 1, &STy[0], 1);
            Vmath::Vadd(nPts, &STy[0], 1, &tmp3[0], 1, &STy[0], 1);
            
            // Computation of STz
            
            // Re-initialise temporary arrays
            Vmath::Zero(nPts, &tmp1[0], 1);
            Vmath::Zero(nPts, &tmp2[0], 1);
            Vmath::Zero(nPts, &tmp3[0], 1);
            
            // w * Szz
            Vmath::Vmul(nPts, &physfield_interp[2][0], 1,
                        &Sgg[2][0], 1, &STz[0], 1);
            
            // u * Sxz
            Vmath::Vmul(nPts, &physfield_interp[0][0], 1,
                        &Sxz[0], 1, &tmp1[0], 1);
            
            // v * Syz
            Vmath::Vmul(nPts, &physfield_interp[1][0], 1,
                        &Syz[0], 1, &tmp2[0], 1);
            
            // k * dT/dz
            Vmath::Smul(nPts, m_thermalConductivity,
                        &derivativesO1_interp[2][3][0], 1, &tmp3[0], 1);
            
            // STz = w * Szz + u * Sxz + v * Syz + K * dT/dz
            Vmath::Vadd(nPts, &STz[0], 1, &tmp1[0], 1, &STz[0], 1);
            Vmath::Vadd(nPts, &STz[0], 1, &tmp2[0], 1, &STz[0], 1);
            Vmath::Vadd(nPts, &STz[0], 1, &tmp3[0], 1, &STz[0], 1);
        }
        
        switch (m_spacedim)
        {
            case 1:
            {
                
                int nq = physfield[0].num_elements();
                // f_11v = f_rho = 0
                Vmath::Zero(nq, &viscousTensor_interp[0][0][0], 1);
                
                // f_21v = f_rhou
                Vmath::Vcopy(nPts, &Sgg[0][0], 1, &viscousTensor_interp[0][1][0], 1);
                
                // f_31v = f_E
                Vmath::Vcopy(nPts, &STx[0], 1, &viscousTensor_interp[0][2][0], 1);
                break;
            }
            case 2:
            {
                int nq = physfield[0].num_elements();
                // f_11v = f_rho1 = 0
                Vmath::Zero(nq, &viscousTensor_interp[0][0][0], 1);
                // f_12v = f_rho2 = 0
                Vmath::Zero(nq, &viscousTensor_interp[1][0][0], 1);
                
                // f_21v = f_rhou1
                Vmath::Vcopy(nPts, &Sgg[0][0], 1, &viscousTensor_interp[0][1][0], 1);
                // f_22v = f_rhou2
                Vmath::Vcopy(nPts, &Sxy[0],    1, &viscousTensor_interp[1][1][0], 1);
                
                // f_31v = f_rhov1
                Vmath::Vcopy(nPts, &Sxy[0],    1, &viscousTensor_interp[0][2][0], 1);
                // f_32v = f_rhov2
                Vmath::Vcopy(nPts, &Sgg[1][0], 1, &viscousTensor_interp[1][2][0], 1);
                
                // f_41v = f_E1
                Vmath::Vcopy(nPts, &STx[0], 1, &viscousTensor_interp[0][3][0], 1);
                // f_42v = f_E2
                Vmath::Vcopy(nPts, &STy[0], 1, &viscousTensor_interp[1][3][0], 1);
                break;
            }
            case 3:
            {
                int nq = physfield[0].num_elements();
                // f_11v = f_rho1 = 0
                Vmath::Zero(nq, &viscousTensor_interp[0][0][0], 1);
                // f_12v = f_rho2 = 0
                Vmath::Zero(nq, &viscousTensor_interp[1][0][0], 1);
                // f_13v = f_rho3 = 0
                Vmath::Zero(nq, &viscousTensor_interp[2][0][0], 1);
                
                // f_21v = f_rhou1
                Vmath::Vcopy(nPts, &Sgg[0][0], 1, &viscousTensor_interp[0][1][0], 1);
                // f_22v = f_rhou2
                Vmath::Vcopy(nPts, &Sxy[0],    1, &viscousTensor_interp[1][1][0], 1);
                // f_23v = f_rhou3
                Vmath::Vcopy(nPts, &Sxz[0],    1, &viscousTensor_interp[2][1][0], 1);
                
                // f_31v = f_rhov1
                Vmath::Vcopy(nPts, &Sxy[0],    1, &viscousTensor_interp[0][2][0], 1);
                // f_32v = f_rhov2
                Vmath::Vcopy(nPts, &Sgg[1][0], 1, &viscousTensor_interp[1][2][0], 1);
                // f_33v = f_rhov3
                Vmath::Vcopy(nPts, &Syz[0],    1, &viscousTensor_interp[2][2][0], 1);
                
                // f_31v = f_rhow1
                Vmath::Vcopy(nPts, &Sxz[0],    1, &viscousTensor_interp[0][3][0], 1);
                // f_32v = f_rhow2
                Vmath::Vcopy(nPts, &Syz[0],    1, &viscousTensor_interp[1][3][0], 1);
                // f_33v = f_rhow3
                Vmath::Vcopy(nPts, &Sgg[2][0], 1, &viscousTensor_interp[2][3][0], 1);
                
                // f_41v = f_E1
                Vmath::Vcopy(nPts, &STx[0], 1, &viscousTensor_interp[0][4][0], 1);
                // f_42v = f_E2
                Vmath::Vcopy(nPts, &STy[0], 1, &viscousTensor_interp[1][4][0], 1);
                // f_43v = f_E3
                Vmath::Vcopy(nPts, &STz[0], 1, &viscousTensor_interp[2][4][0], 1);
                break;
            }
            default:
            {
                ASSERTL0(false, "Illegal expansion dimension");
            }
        }
        
        for (i = 0; i < m_spacedim; ++i)
        {
            for (j = 1; j < nVariables; ++j)
            {
                // Galerkin project solution back to origianl space
                m_fields[0]->PhysGalerkinProjection1DScaled(
                                                            OneDptscale,
                                                            viscousTensor_interp[i][j],
                                                            viscousTensor[i][j]);
            }
        }
#endif
    }
    
    /**
     * @brief Calculate the pressure field \f$ p =
     * (\gamma-1)(E-\frac{1}{2}\rho\| \mathbf{v} \|^2) \f$ assuming an ideal
     * gas law.
     *
     * @param physfield  Input momentum.
     * @param pressure   Computed pressure field.
     */
    void CompressibleFlowSystem::GetPressure(
                                             const Array<OneD, const Array<OneD, NekDouble> > &physfield,
                                             Array<OneD,                   NekDouble>   &pressure)
    {
        int       nBCEdgePts  = physfield[0].num_elements();
        NekDouble alpha = -0.5;
        
        // Calculate ||rho v||^2
        Vmath::Vmul(nBCEdgePts, physfield[1], 1, physfield[1], 1, pressure, 1);
        for (int i = 1; i < m_spacedim; ++i)
        {
            Vmath::Vvtvp(nBCEdgePts, physfield[1+i], 1, physfield[1+i], 1,
                         pressure,       1, pressure,       1);
        }
        // Divide by rho to get rho*||v||^2
        Vmath::Vdiv (nBCEdgePts, pressure, 1, physfield[0], 1, pressure, 1);
        // pressure <- E - 0.5*pressure
        Vmath::Svtvp(nBCEdgePts, alpha,
                     pressure, 1, physfield[m_spacedim+1], 1, pressure, 1);
        // Multiply by (gamma-1)
        Vmath::Smul (nBCEdgePts, m_gamma-1, pressure, 1, pressure, 1);
    }
    
    /**
     * @brief Compute the enthalpy term \f$ H = E + p/rho \$.
     */
    void CompressibleFlowSystem::GetEnthalpy(
            const Array<OneD, const Array<OneD, NekDouble> > &physfield,
                  Array<OneD,                   NekDouble>   &pressure,
                  Array<OneD,                   NekDouble>   &enthalpy)
    {
        int npts  = m_fields[0]->GetTotPoints();
        Array<OneD, NekDouble> tmp(npts, 0.0);

        // Calculate E = rhoE/rho
        Vmath::Vdiv(npts, physfield[m_spacedim+1], 1, physfield[0], 1, tmp, 1);
        // Calculate p/rho
        Vmath::Vdiv(npts, pressure, 1, physfield[0], 1, enthalpy, 1);
        // Calculate H = E + p/rho
        Vmath::Vadd(npts, tmp, 1, enthalpy, 1, enthalpy, 1);
    }

    /**
     * @brief Calculate the pressure field \f$ p =
     * (\gamma-1)(E-\frac{1}{2}\rho\| \mathbf{v} \|^2) \f$ assuming an ideal
     * gas law.
     *
     * This is a slightly optimised way to calculate the pressure field which
     * avoids division by the density field if the velocity field has already
     * been calculated.
     *
     * @param physfield  Input momentum.
     * @param velocity   Velocity vector.
     * @param pressure   Computed pressure field.
     */
    void CompressibleFlowSystem::GetPressure(
                                             const Array<OneD, const Array<OneD, NekDouble> > &physfield,
                                             const Array<OneD, const Array<OneD, NekDouble> > &velocity,
                                             Array<OneD,                   NekDouble>   &pressure)
    {
        int nBCEdgePts = physfield[0].num_elements();
        NekDouble alpha = -0.5;
        
        // Calculate ||\rho v||^2.
        Vmath::Vmul (nBCEdgePts, velocity[0], 1, physfield[1], 1, pressure, 1);
        for (int i = 1; i < m_spacedim; ++i)
        {
            Vmath::Vvtvp(nBCEdgePts, velocity[i], 1, physfield[1+i], 1,
                         pressure,    1, pressure,       1);
        }
        // pressure <- E - 0.5*pressure
        Vmath::Svtvp(nBCEdgePts,     alpha,
                     pressure, 1, physfield[m_spacedim+1], 1, pressure, 1);
        // Multiply by (gamma-1).
        Vmath::Smul (nBCEdgePts, m_gamma-1, pressure, 1, pressure, 1);
    }
    
    /**
     * @brief Compute the velocity field \f$ \mathbf{v} \f$ given the momentum
     * \f$ \rho\mathbf{v} \f$.
     *
     * @param physfield  Momentum field.
     * @param velocity   Velocity field.
     */
    void CompressibleFlowSystem::GetVelocityVector(
                                                   const Array<OneD, Array<OneD, NekDouble> > &physfield,
                                                   Array<OneD, Array<OneD, NekDouble> > &velocity)
    {
        const int nBCEdgePts = physfield[0].num_elements();
        
        for (int i = 0; i < m_spacedim; ++i)
        {
            Vmath::Vdiv(nBCEdgePts, physfield[1+i], 1, physfield[0], 1,
                        velocity[i],    1);
        }
    }
    
    /**
     * @brief Compute the temperature \f$ T = p/\rho R \f$.
     *
     * @param physfield    Input physical field.
     * @param pressure     The pressure field corresponding to physfield.
     * @param temperature  The resulting temperature \f$ T \f$.
     */
    void CompressibleFlowSystem::GetTemperature(
                                                const Array<OneD, const Array<OneD, NekDouble> > &physfield,
                                                Array<OneD,                         NekDouble  > &pressure,
                                                Array<OneD,                         NekDouble  > &temperature)
    {
        const int nq = physfield[0].num_elements();
        
        Vmath::Vdiv(nq, pressure, 1, physfield[0], 1, temperature, 1);
        Vmath::Smul(nq, 1.0/m_gasConstant, temperature, 1, temperature, 1);
    }
    
    /**
     * @brief Compute the sound speed \f$ c = sqrt(\gamma p/\rho) \f$.
     *
     * @param physfield    Input physical field.
     * @param pressure     The pressure field corresponding to physfield.
     * @param soundspeed   The resulting sound speed \f$ c \f$.
     */
    void CompressibleFlowSystem::GetSoundSpeed(
                                               const Array<OneD, Array<OneD, NekDouble> > &physfield,
                                               Array<OneD,             NekDouble  > &pressure,
                                               Array<OneD,             NekDouble  > &soundspeed)
    {
        const int nq = m_fields[0]->GetTotPoints();
        Vmath::Vdiv (nq, pressure, 1, physfield[0], 1, soundspeed, 1);
        Vmath::Smul (nq, m_gamma, soundspeed, 1, soundspeed, 1);
        Vmath::Vsqrt(nq, soundspeed, 1, soundspeed, 1);
    }
    
    /**
     * @brief Compute the mach number \f$ M = \| \mathbf{v} \|^2 / c \f$.
     *
     * @param physfield    Input physical field.
     * @param soundfield   The speed of sound corresponding to physfield.
     * @param mach         The resulting mach number \f$ M \f$.
     */
    void CompressibleFlowSystem::GetMach(
            Array<OneD, Array<OneD, NekDouble> > &physfield,
            Array<OneD,             NekDouble  > &soundspeed,
            Array<OneD,             NekDouble  > &mach)
    {

        const int nq = m_fields[0]->GetTotPoints();
        
        Vmath::Vmul(nq, physfield[1], 1, physfield[1], 1, mach, 1);
        
        for (int i = 1; i < m_spacedim; ++i)
        {
            Vmath::Vvtvp(nq, physfield[1+i], 1, physfield[1+i], 1,
                             mach,           1, mach,           1);
        }
        
        Vmath::Vdiv(nq, mach, 1, physfield[0], 1, mach, 1);
        Vmath::Vdiv(nq, mach, 1, physfield[0], 1, mach, 1);
        Vmath::Vsqrt(nq, mach, 1, mach, 1);
        
        Vmath::Vdiv(nq, mach, 1, soundspeed,   1, mach, 1);
    }
    
    /**
     * @brief Compute the dynamic viscosity using the Sutherland's law
     * \f$ \mu = \mu_star * (T / T_star)^3/2 * (T_star + 110) / (T + 110) \f$,
     * where:   \mu_star = 1.7894 * 10^-5 Kg / (m * s)
     *          T_star   = 288.15 K
     *
     * @param physfield    Input physical field.
     * @param mu           The resulting dynamic viscosity.
     */
    void CompressibleFlowSystem::GetDynamicViscosity(
                                                     const Array<OneD, const NekDouble> &temperature,
                                                     Array<OneD,       NekDouble> &mu)
    {
        const int nPts    = temperature.num_elements();
        NekDouble mu_star = m_mu;
        NekDouble T_star  = m_pInf / (m_rhoInf * m_gasConstant);
        NekDouble ratio;
        
        for (int i = 0; i < nPts; ++i)
        {
            ratio = temperature[i] / T_star;
            mu[i] = mu_star * ratio * sqrt(ratio) *
            (T_star + 110.0) / (temperature[i] + 110.0);
        }
    }
    
    /**
     * @brief Calcualte entropy.
     */
    void CompressibleFlowSystem::GetEntropy(
        const Array<OneD, const Array<OneD, NekDouble> > &physfield,
        const Array<OneD, const NekDouble>               &pressure,
        const Array<OneD, const NekDouble>               &temperature,
              Array<OneD,       NekDouble>               &entropy)
    {
        NekDouble entropy_sum = 0.0;
        const int npts = m_fields[0]->GetTotPoints();
        const NekDouble temp_inf = m_pInf/(m_rhoInf*m_gasConstant);;
        Array<OneD, NekDouble> L2entropy(npts, 0.0);

        for (int i = 0; i < npts; ++i)
        {
            entropy[i] = m_gamma/(m_gamma-1.0)*m_gasConstant*log(temperature[i]/temp_inf) -
                m_gasConstant*log(pressure[i]/m_pInf);
        }
        
        Vmath::Vmul(npts,entropy,1,entropy,1,L2entropy,1);
        
        entropy_sum = Vmath::Vsum(npts, L2entropy, 1);
        
        entropy_sum = sqrt(entropy_sum);
        
        std::ofstream m_file( "L2entropy.txt", std::ios_base::app);
        
        m_file << setprecision(16) << scientific << entropy_sum << endl;
        //m_file << Vmath::Vmax(entropy.num_elements(),entropy,1) << endl;
        
        m_file.close();
    }

    /**
     * @brief Calculate the maximum timestep subject to CFL restrictions.
     */
    NekDouble CompressibleFlowSystem::v_GetTimeStep(
                const Array<OneD, const Array<OneD, NekDouble> > &inarray)
    {
        int n;
        int nElements = m_fields[0]->GetExpSize();
        const Array<OneD, int> ExpOrder = GetNumExpModesPerExp();
        
        Array<OneD, NekDouble> tstep      (nElements, 0.0);
        Array<OneD, NekDouble> stdVelocity(nElements);
        
        // Get standard velocity to compute the time-step limit
        GetStdVelocity(inarray, stdVelocity);
        
        // Factors to compute the time-step limit
        NekDouble minLength;
        NekDouble alpha   = MaxTimeStepEstimator();
        NekDouble cLambda = 0.2; // Spencer book-317
        
        // Loop over elements to compute the time-step limit for each element
        for(n = 0; n < nElements; ++n)
        {
            int npoints = m_fields[0]->GetExp(n)->GetTotPoints();
            Array<OneD, NekDouble> one2D(npoints, 1.0);
            NekDouble Area = m_fields[0]->GetExp(n)->Integral(one2D);

            if (m_fields[0]->GetExp(n)->as<LocalRegions::TriExp>())
            {
                minLength = 2.0 * sqrt(Area);
            }
            else if (m_fields[0]->GetExp(n)->as<LocalRegions::QuadExp>())
            {
                minLength = sqrt(Area);
            }
            else if (m_fields[0]->GetExp(n)->as<LocalRegions::HexExp>())
            {
                minLength = sqrt(Area);
            }
            
            tstep[n] = m_cflSafetyFactor * alpha * minLength
            / (stdVelocity[n] * cLambda
               * (ExpOrder[n] - 1) * (ExpOrder[n] - 1));
        }
        
        // Get the minimum time-step limit and return the time-step
        NekDouble TimeStep = Vmath::Vmin(nElements, tstep, 1);
        m_comm->AllReduce(TimeStep, LibUtilities::ReduceMin);
        return TimeStep;
    }
    
    /**
     * @brief Compute the advection velocity in the standard space
     * for each element of the expansion.
     *
     * @param inarray    Momentum field.
     * @param stdV       Standard velocity field.
     */
    void CompressibleFlowSystem::GetStdVelocity(
                    const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                          Array<OneD,                   NekDouble>   &stdV)
    {
        int nTotQuadPoints = GetTotPoints();
        int n_element      = m_fields[0]->GetExpSize();
        int nBCEdgePts           = 0;
        
        // Getting the velocity vector on the 2D normal space
        Array<OneD, Array<OneD, NekDouble> > velocity   (m_spacedim);
        Array<OneD, Array<OneD, NekDouble> > stdVelocity(m_spacedim);
        Array<OneD, NekDouble>               pressure   (nTotQuadPoints);
        Array<OneD, NekDouble>               soundspeed (nTotQuadPoints);
        LibUtilities::PointsKeyVector        ptsKeys;
        
        // Zero output array
        Vmath::Zero(stdV.num_elements(), stdV, 1);
        
        for (int i = 0; i < m_spacedim; ++i)
        {
            velocity   [i] = Array<OneD, NekDouble>(nTotQuadPoints);
            stdVelocity[i] = Array<OneD, NekDouble>(nTotQuadPoints, 0.0);
        }

        GetVelocityVector(inarray, velocity);
        GetPressure      (inarray, velocity, pressure);
        GetSoundSpeed    (inarray, pressure, soundspeed);
        
        for(int el = 0; el < n_element; ++el)
        {
            ptsKeys = m_fields[0]->GetExp(el)->GetPointsKeys();
            
            // Possible bug: not multiply by jacobian??
            const SpatialDomains::GeomFactorsSharedPtr metricInfo =
            m_fields[0]->GetExp(el)->GetGeom()->GetMetricInfo();
            const Array<TwoD, const NekDouble> &gmat =
            m_fields[0]->GetExp(el)->GetGeom()->GetMetricInfo()
            ->GetDerivFactors(ptsKeys);
            
            int nq = m_fields[0]->GetExp(el)->GetTotPoints();
            
            if(metricInfo->GetGtype() == SpatialDomains::eDeformed)
            {
                // d xi/ dx = gmat = 1/J * d x/d xi
                for (int i = 0; i < m_spacedim; ++i)
                {
                    Vmath::Vmul(nq, gmat[i], 1, velocity[0], 1,
                                stdVelocity[i], 1);
                    for (int j = 1; j < m_spacedim; ++j)
                    {
                        Vmath::Vvtvp(nq, gmat[m_spacedim*j+i], 1, velocity[j],
                                     1, stdVelocity[i], 1, stdVelocity[i], 1);
                    }
                }
            }
            else
            {
                for (int i = 0; i < m_spacedim; ++i)
                {
                    Vmath::Smul(nq, gmat[i][0], velocity[0], 1,
                                stdVelocity[i], 1);
                    for (int j = 1; j < m_spacedim; ++j)
                    {
                        Vmath::Svtvp(nq, gmat[m_spacedim*j+i][0], velocity[j],
                                     1, stdVelocity[i], 1, stdVelocity[i], 1);
                    }
                }
            }
            
            for (int i = 0; i < nq; ++i)
            {
                NekDouble pntVelocity = 0.0;
                for (int j = 0; j < m_spacedim; ++j)
                {
                    pntVelocity += stdVelocity[j][i]*stdVelocity[j][i];
                }
                pntVelocity = sqrt(pntVelocity) + soundspeed[nBCEdgePts];
                if (pntVelocity > stdV[el])
                {
                    stdV[el] = pntVelocity;
                }
                nBCEdgePts++;
            }
        }
    }
    
    /**
     * @brief Set the denominator to compute the time step when a cfl
     * control is employed. This function is no longer used but is still
     * here for being utilised in the future.
     *
     * @param n   Order of expansion element by element.
     */
    NekDouble CompressibleFlowSystem::GetStabilityLimit(int n)
    {
        ASSERTL0(n <= 20, "Illegal modes dimension for CFL calculation "
                 "(P has to be less then 20)");
        
        NekDouble CFLDG[21] = {  2.0000,   6.0000,  11.8424,  19.1569,
            27.8419,  37.8247,  49.0518,  61.4815,
            75.0797,  89.8181, 105.6700, 122.6200,
            140.6400, 159.7300, 179.8500, 201.0100,
            223.1800, 246.3600, 270.5300, 295.6900,
            321.8300}; //CFLDG 1D [0-20]
        NekDouble CFL = 0.0;
        
        if (m_projectionType == MultiRegions::eDiscontinuous)
        {
            CFL = CFLDG[n];
        }
        else
        {
            ASSERTL0(false, "Continuous Galerkin stability coefficients "
                     "not introduced yet.");
        }
        
        return CFL;
    }
    
    /**
     * @brief Compute the vector of denominators to compute the time step
     * when a cfl control is employed. This function is no longer used but
     * is still here for being utilised in the future.
     *
     * @param ExpOrder   Order of expansion element by element.
     */
    Array<OneD, NekDouble> CompressibleFlowSystem::GetStabilityLimitVector(
                                                                           const Array<OneD,int> &ExpOrder)
    {
        int i;
        Array<OneD,NekDouble> returnval(m_fields[0]->GetExpSize(), 0.0);
        for (i =0; i<m_fields[0]->GetExpSize(); i++)
        {
            returnval[i] = GetStabilityLimit(ExpOrder[i]);
        }
        return returnval;
    }
    // Determined matrix to change from conservative (U) to primitive (V) variables
    // where U = [rho, rhou, rhov, rhow, rhoE]^t  and V = [rho, u, v, w, p]^t
    
    void CompressibleFlowSystem::GetConservToPrimVariableMat(
                  const Array<OneD, Array<OneD, NekDouble> > &inarray,
                        Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &dUdV)
    {
        int i;
        int nq = inarray[0].num_elements();
        int nvar = m_fields.num_elements();
        int nDimensions = m_fields.num_elements();
        
        NekDouble GamMinOneScalar           = m_gamma - 1.0;
        NekDouble OneMinGamScalar           = 1.0 - m_gamma;
        NekDouble GamMinOneOverTwoScalar    = 0.5 * (m_gamma - 1.0);
        NekDouble OneOverGamMinOneScalar    = 1.0 / (m_gamma - 1.0);
        Array<OneD, NekDouble> OneOverGamMinOne(nq, OneOverGamMinOneScalar);
        Array<OneD, NekDouble> GamMinOneVec(nq, GamMinOneScalar);
        
        Array<OneD, NekDouble> pressure(nq, 0.0);
        Array<OneD, NekDouble> ones(nq,  1.0);
        Array<OneD, NekDouble> zeros(nq, 0.0);
        Array<OneD, NekDouble> velsq(nq,0.0);
        Array<OneD, NekDouble> HalfVelsq(nq, 0.0);
        
        Array<OneD, Array<OneD, NekDouble> > vel(m_spacedim);
        Array<OneD, Array<OneD, NekDouble> > direct_fields(nvar);
        
        for (i = 0; i < nvar; ++i)
        {
            direct_fields[i] = Array<OneD, NekDouble>(nq, 0.0);
        }
        
        GetDirectSolution(direct_fields);
        
        for (i = 0; i < m_spacedim; ++i)
        {
            vel[i] = Array<OneD, NekDouble>(nq);
        }
        
        GetVelocityVector(direct_fields, vel);
        GetPressure      (direct_fields, vel, pressure);
        
        // determine u^2+v^2+w^2;
        for (i = 0; i < m_spacedim; ++i)
        {
            Vmath::Vvtvp(nq, vel[i], 1, vel[i], 1, velsq, 1, velsq, 1);
        }
        
        Vmath::Smul(nq, 0.5, &velsq[0], 1, &HalfVelsq[0], 1);
        
        if (m_spacedim == 2)
        {
            //
            /* =================================================================
             M =
             [      1,      0,            0,            0]
             [      u,    rho,            0,            0]
             [      v,      0,          rho,            0]
             [  U^2/2,  rho*u,        rho*v,  1/(gamma-1)]
             */
            // =================================================================
            // row one
            Vmath::Vcopy(nq, &ones[0],              1, &dUdV[0][0][0], 1);
            Vmath::Vcopy(nq, &zeros[0],             1, &dUdV[0][1][0], 1);
            Vmath::Vcopy(nq, &zeros[0],             1, &dUdV[0][2][0], 1);
            Vmath::Vcopy(nq, &zeros[0],             1, &dUdV[0][3][0], 1);
            // row two
            Vmath::Vcopy(nq, &vel[0][0],            1, &dUdV[1][0][0], 1);
            Vmath::Vcopy(nq, &inarray[0][0],        1, &dUdV[1][1][0], 1);
            Vmath::Vcopy(nq, &zeros[0],             1, &dUdV[1][2][0], 1);
            Vmath::Vcopy(nq, &zeros[0],             1, &dUdV[1][3][0], 1);
            // row three
            Vmath::Vcopy(nq, &vel[1][0],            1, &dUdV[2][0][0], 1);
            Vmath::Vcopy(nq, &zeros[0],             1, &dUdV[2][1][0], 1);
            Vmath::Vcopy(nq, &inarray[0][0],        1, &dUdV[2][2][0], 1);
            Vmath::Vcopy(nq, &zeros[0],             1, &dUdV[2][3][0], 1);
            // row four
            Vmath::Vcopy(nq, &HalfVelsq[0],         1, &dUdV[3][0][0], 1);
            Vmath::Vcopy(nq, &inarray[1][0],        1, &dUdV[3][1][0], 1);
            Vmath::Vcopy(nq, &inarray[2][0],        1, &dUdV[3][2][0], 1);
            Vmath::Vcopy(nq, &OneOverGamMinOne[0],  1, &dUdV[3][3][0], 1);
        }
        
        // Needs to be checked still
        if (m_spacedim == 3)
        {
            /* =================================================================
             M =
             [      1,      0,            0,              0,              0]
             [      u,    rho,            0,              0,              0]
             [      v,      0,          rho,              0,              0]
             [      W,      0,            0,            rho,              0]
             [  U^2/2,  rho*u,        rho*v,          rho*w,    1/(gamma-1)]
             */
            // =================================================================
            // row one
            Vmath::Vcopy(nq, &ones[0],              1, &dUdV[0][0][0], 1);
            Vmath::Vcopy(nq, &zeros[0],             1, &dUdV[0][1][0], 1);
            Vmath::Vcopy(nq, &zeros[0],             1, &dUdV[0][2][0], 1);
            Vmath::Vcopy(nq, &zeros[0],             1, &dUdV[0][3][0], 1);
            Vmath::Vcopy(nq, &zeros[0],             1, &dUdV[0][4][0], 1);
            // row two
            Vmath::Vcopy(nq, &vel[0][0],            1, &dUdV[1][0][0], 1);
            Vmath::Vcopy(nq, &inarray[0][0],        1, &dUdV[1][1][0], 1);
            Vmath::Vcopy(nq, &zeros[0],             1, &dUdV[1][2][0], 1);
            Vmath::Vcopy(nq, &zeros[0],             1, &dUdV[1][3][0], 1);
            Vmath::Vcopy(nq, &zeros[0],             1, &dUdV[1][4][0], 1);
            // row three
            Vmath::Vcopy(nq, &vel[1][0],            1, &dUdV[2][0][0], 1);
            Vmath::Vcopy(nq, &zeros[0],             1, &dUdV[2][1][0], 1);
            Vmath::Vcopy(nq, &inarray[0][0],        1, &dUdV[2][2][0], 1);
            Vmath::Vcopy(nq, &zeros[0],             1, &dUdV[2][3][0], 1);
            Vmath::Vcopy(nq, &zeros[0],             1, &dUdV[2][4][0], 1);
            // row four
            Vmath::Vcopy(nq, &vel[2][0],            1, &dUdV[3][0][0], 1);
            Vmath::Vcopy(nq, &zeros[0],             1, &dUdV[3][1][0], 1);
            Vmath::Vcopy(nq, &zeros[0],             1, &dUdV[3][2][0], 1);
            Vmath::Vcopy(nq, &inarray[0][0],        1, &dUdV[3][3][0], 1);
            Vmath::Vcopy(nq, &zeros[0],             1, &dUdV[3][4][0], 1);
            // row five
            Vmath::Vcopy(nq, &HalfVelsq[0],         1, &dUdV[4][0][0], 1);
            Vmath::Vcopy(nq, &inarray[1][0],        1, &dUdV[4][1][0], 1);
            Vmath::Vcopy(nq, &inarray[2][0],        1, &dUdV[4][2][0], 1);
            Vmath::Vcopy(nq, &inarray[3][0],        1, &dUdV[4][3][0], 1);
            Vmath::Vcopy(nq, &OneOverGamMinOne[0],  1, &dUdV[4][4][0], 1);
        }
    }
    
    void CompressibleFlowSystem::GetConservToPrimVariableInvMat(
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &dUdVInv)
    {
        int i;
        int nq = dUdVInv[0][0].num_elements();
        int nvar = m_fields.num_elements();
        int nDimensions = m_fields.num_elements();
        
        NekDouble GamMinOneScalar           = m_gamma - 1.0;
        NekDouble OneMinGamScalar           = 1.0 - m_gamma;
        NekDouble GamMinOneOverTwoScalar    = 0.5 * (m_gamma - 1.0);
        NekDouble OneOverGamMinOneScalar    = 1.0 / (m_gamma - 1.0);
        
        Array<OneD, NekDouble> ones(nq,  1.0);
        
        Array<OneD, NekDouble> pressure(nq, 0.0);
        Array<OneD, NekDouble> zeros(nq, 0.0);
        Array<OneD, NekDouble> oneOverRho(nq, 0.0);
        Array<OneD, NekDouble> GamMinOneOverTwoVelsq(nq, 0.0);
        Array<OneD, NekDouble> HalfVelsq(nq, 0.0);
        Array<OneD, NekDouble> velsq(nq, 0.0);
        
        Array<OneD, NekDouble> GamMinOneVec(nq, GamMinOneScalar);
        Array<OneD, NekDouble> OneMinGamVec(nq, OneMinGamScalar);
        Array<OneD, NekDouble> OneOverGamMinOne(nq, OneOverGamMinOneScalar);
        Array<OneD, NekDouble> GamMinOneOverTwoVec(nq, GamMinOneOverTwoScalar);
        
        Array<OneD, Array<OneD, NekDouble> > vel(m_spacedim);
        Array<OneD, Array<OneD, NekDouble> > OneMinGamVel(m_spacedim);
        Array<OneD, Array<OneD, NekDouble> > velOverRho(m_spacedim);
        Array<OneD, Array<OneD, NekDouble> > direct_fields(nvar);
        
        for (i = 0; i < nvar; ++i)
        {
            direct_fields[i] = Array<OneD, NekDouble>(nq, 0.0);
        }
        
        GetDirectSolution(direct_fields);
        
        for (i = 0; i < m_spacedim; ++i)
        {
            vel[i]          = Array<OneD, NekDouble>(nq);
            velOverRho[i]   = Array<OneD, NekDouble>(nq);
            OneMinGamVel[i] = Array<OneD, NekDouble>(nq);
        }
        
        GetVelocityVector(direct_fields, vel);
        GetPressure      (direct_fields, vel, pressure);
        
        // determine u^2+v^2+w^2;
        for (i = 0; i < m_spacedim; ++i)
        {
            Vmath::Vvtvp(nq, vel[i], 1, vel[i], 1, velsq, 1, velsq, 1);
            Vmath::Vdiv(nq, vel[i], 1, direct_fields[0], 1, velOverRho[i], 1);
            Vmath::Neg(nq, &velOverRho[i][0], 1);
            Vmath::Vmul(nq, &OneMinGamVec[0], 1,
                        &vel[i][0], 1,
                        &OneMinGamVel[i][0], 1);
        }
        
        Vmath::Smul(nq, 0.5, &velsq[0], 1, &HalfVelsq[0], 1);
        
        Vmath::Vmul(nq, &GamMinOneOverTwoVec[0], 1,
                    &velsq[0], 1,
                    &GamMinOneOverTwoVelsq[0], 1);
        
        Vmath::Vdiv(nq, &ones[0], 1, &direct_fields[0][0], 1, &oneOverRho[0], 1);
        
        if (m_spacedim == 2)
        {
            /* =================================================================
             M^-1 =
             [             1,             0,            0,            0]
             [        -u/rho,         1/rho,            0,            0]
             [        -v/rho,             0,        1/rho,            0]
             [(gamma-1)U^2/2,   (1-gamma)*u,  (1-gamma)*v,    (gamma-1)]
             */
            // =================================================================
            // row one
            Vmath::Vcopy(nq, &ones[0],                  1, &dUdVInv[0][0][0], 1);
            Vmath::Vcopy(nq, &zeros[0],                 1, &dUdVInv[1][0][0], 1);
            Vmath::Vcopy(nq, &zeros[0],                 1, &dUdVInv[2][0][0], 1);
            Vmath::Vcopy(nq, &zeros[0],                 1, &dUdVInv[3][0][0], 1);
            // row two
            Vmath::Vcopy(nq, &velOverRho[0][0],         1, &dUdVInv[0][1][0], 1);
            Vmath::Vcopy(nq, &oneOverRho[0],            1, &dUdVInv[1][1][0], 1);
            Vmath::Vcopy(nq, &zeros[0],                 1, &dUdVInv[2][1][0], 1);
            Vmath::Vcopy(nq, &zeros[0],                 1, &dUdVInv[3][1][0], 1);
            // row three
            Vmath::Vcopy(nq, &velOverRho[1][0],         1, &dUdVInv[0][2][0], 1);
            Vmath::Vcopy(nq, &zeros[0],                 1, &dUdVInv[1][2][0], 1);
            Vmath::Vcopy(nq, &oneOverRho[0],            1, &dUdVInv[2][2][0], 1);
            Vmath::Vcopy(nq, &zeros[0],                 1, &dUdVInv[3][2][0], 1);
            // row four
            Vmath::Vcopy(nq, &GamMinOneOverTwoVelsq[0], 1, &dUdVInv[0][3][0], 1);
            Vmath::Vcopy(nq, &OneMinGamVel[0][0],       1, &dUdVInv[1][3][0], 1);
            Vmath::Vcopy(nq, &OneMinGamVel[1][0],       1, &dUdVInv[2][3][0], 1);
            Vmath::Vcopy(nq, &GamMinOneVec[0],          1, &dUdVInv[3][3][0], 1);
            
        }
        
        // Needs to be checked still
        if (m_spacedim == 3)
        {
            /* =================================================================
             M^-1 =
             [             1,           0,          0,          0,          0]
             [        -u/rho,       1/rho,          0,          0,          0]
             [        -v/rho,           0,      1/rho,          0,          0]
             [        -w/rho,           0,          0,      1/rho,          0]
             [(gamma-1)U^2/2, (1-gamma)*u,(1-gamma)*v,(gamma-1)*w,  (gamma-1)]
             */
            // =================================================================
            // row one
            Vmath::Vcopy(nq, &ones[0],                  1, &dUdVInv[0][0][0], 1);
            Vmath::Vcopy(nq, &zeros[0],                 1, &dUdVInv[1][1][0], 1);
            Vmath::Vcopy(nq, &zeros[0],                 1, &dUdVInv[2][2][0], 1);
            Vmath::Vcopy(nq, &zeros[0],                 1, &dUdVInv[3][3][0], 1);
            Vmath::Vcopy(nq, &zeros[0],                 1, &dUdVInv[4][4][0], 1);
            // row two
            Vmath::Vcopy(nq, &velOverRho[0][0],         1, &dUdVInv[1][0][0], 1);
            Vmath::Vcopy(nq, &oneOverRho[0],            1, &dUdVInv[1][1][0], 1);
            Vmath::Vcopy(nq, &zeros[0],                 1, &dUdVInv[1][2][0], 1);
            Vmath::Vcopy(nq, &zeros[0],                 1, &dUdVInv[1][3][0], 1);
            Vmath::Vcopy(nq, &zeros[0],                 1, &dUdVInv[1][4][0], 1);
            // row three
            Vmath::Vcopy(nq, &velOverRho[1][0],         1, &dUdVInv[2][0][0], 1);
            Vmath::Vcopy(nq, &zeros[0],                 1, &dUdVInv[2][1][0], 1);
            Vmath::Vcopy(nq, &oneOverRho[0],            1, &dUdVInv[2][2][0], 1);
            Vmath::Vcopy(nq, &zeros[0],                 1, &dUdVInv[2][3][0], 1);
            Vmath::Vcopy(nq, &zeros[0],                 1, &dUdVInv[2][4][0], 1);
            // row four
            Vmath::Vcopy(nq, &velOverRho[2][0],         1, &dUdVInv[3][0][0], 1);
            Vmath::Vcopy(nq, &zeros[0],                 1, &dUdVInv[3][1][0], 1);
            Vmath::Vcopy(nq, &zeros[0],                 1, &dUdVInv[3][2][0], 1);
            Vmath::Vcopy(nq, &oneOverRho[0],            1, &dUdVInv[3][3][0], 1);
            Vmath::Vcopy(nq, &zeros[0],                 1, &dUdVInv[3][4][0], 1);
            // row five
            Vmath::Vcopy(nq, &GamMinOneOverTwoVelsq[0], 1, &dUdVInv[4][0][0], 1);
            Vmath::Vcopy(nq, &OneMinGamVel[0][0],       1, &dUdVInv[4][1][0], 1);
            Vmath::Vcopy(nq, &OneMinGamVel[1][0],       1, &dUdVInv[4][2][0], 1);
            Vmath::Vcopy(nq, &OneMinGamVel[2][0],       1, &dUdVInv[4][3][0], 1);
            Vmath::Vcopy(nq, &GamMinOneVec[0],          1, &dUdVInv[4][4][0], 1);
        }
    }
    
    
    void CompressibleFlowSystem::GetConservToPrimVariableInvMatDiv(
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &dUdVInv,
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &dUdVInvdX,
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &dUdVInvdY)
    {
        int i, j, k;
        int nq = dUdVInvdX[0][0].num_elements();
        int nvar = m_fields.num_elements();
        int nDimensions = m_fields.num_elements();
        
        
        NekDouble GamMinOneScalar           = m_gamma - 1.0;
        NekDouble OneMinGamScalar           = 1.0 - m_gamma;
        NekDouble GamMinOneOverTwoScalar    = 0.5 * (m_gamma - 1.0);
        NekDouble OneOverGamMinOneScalar    = 1.0 / (m_gamma - 1.0);
        
        Array<OneD, NekDouble> ones(nq,  1.0);
        Array<OneD, NekDouble> tmp(nq, 0.0);
        Array<OneD, NekDouble> tmp1(nq, 0.0);
        Array<OneD, NekDouble> pressure(nq, 0.0);
        Array<OneD, NekDouble> zeros(nq, 0.0);
        Array<OneD, NekDouble> oneOverRho(nq, 0.0);
        Array<OneD, NekDouble> GamMinOneOverTwoVelsqdX(nq, 0.0);
        Array<OneD, NekDouble> GamMinOneOverTwoVelsqdY(nq, 0.0);
        Array<OneD, NekDouble> HalfVelsq(nq, 0.0);
        Array<OneD, NekDouble> velsq(nq, 0.0);
        
        Array<OneD, NekDouble> GamMinOneVec(nq, GamMinOneScalar);
        Array<OneD, NekDouble> OneMinGamVec(nq, OneMinGamScalar);
        Array<OneD, NekDouble> OneOverGamMinOne(nq, OneOverGamMinOneScalar);
        Array<OneD, NekDouble> GamMinOneOverTwoVec(nq, GamMinOneOverTwoScalar);
        
        Array<OneD, Array<OneD, NekDouble> > vel(m_spacedim);
        
        Array<OneD, Array<OneD, NekDouble> > dui_dx(m_spacedim);
        Array<OneD, Array<OneD, NekDouble> > dui_dy(m_spacedim);
        
        Array<OneD, Array<OneD, NekDouble> > du2i_dx(m_spacedim);
        Array<OneD, Array<OneD, NekDouble> > du2i_dy(m_spacedim);
        
        Array<OneD, Array<OneD, NekDouble> > duiOverRho_dx(m_spacedim);
        Array<OneD, Array<OneD, NekDouble> > duiOverRho_dy(m_spacedim);
        Array<OneD, Array<OneD, NekDouble> > d1overRho_dxi(m_spacedim);
        
        Array<OneD, Array<OneD, NekDouble> > dU2_dxi(m_spacedim);
        Array<OneD, Array<OneD, NekDouble> > drho_dxi(m_spacedim);
        
        Array<OneD, Array<OneD, NekDouble> > OneMinGamVeldX(m_spacedim);
        Array<OneD, Array<OneD, NekDouble> > OneMinGamVeldY(m_spacedim);
        
        Array<OneD, Array<OneD, NekDouble> > velOverRho(m_spacedim);
        Array<OneD, Array<OneD, NekDouble> > dirFlds(nvar);
        
        for (i = 0; i < nvar; ++i)
        {
            dirFlds[i] = Array<OneD, NekDouble>(nq, 0.0);
        }
        
        GetDirectSolution(dirFlds);
        
        for (i = 0; i < m_spacedim; ++i)
        {
            vel[i]            = Array<OneD, NekDouble>(nq, 0.0);
            velOverRho[i]     = Array<OneD, NekDouble>(nq, 0.0);
            OneMinGamVeldX[i] = Array<OneD, NekDouble>(nq, 0.0);
            OneMinGamVeldY[i] = Array<OneD, NekDouble>(nq, 0.0);
            dui_dx[i]         = Array<OneD, NekDouble>(nq, 0.0);
            dui_dy[i]         = Array<OneD, NekDouble>(nq, 0.0);
            
            duiOverRho_dx[i]  = Array<OneD, NekDouble>(nq, 0.0);
            duiOverRho_dy[i]  = Array<OneD, NekDouble>(nq, 0.0);
            
            dU2_dxi[i]        = Array<OneD, NekDouble>(nq, 0.0);
            drho_dxi[i]       = Array<OneD, NekDouble>(nq, 0.0);
            d1overRho_dxi[i]  = Array<OneD, NekDouble>(nq, 0.0);
            
            du2i_dx[i]        = Array<OneD, NekDouble>(nq, 0.0);
            du2i_dy[i]        = Array<OneD, NekDouble>(nq, 0.0);
        }
        
        GetVelocityVector(dirFlds, vel);
        GetPressure      (dirFlds, vel, pressure);
        
        // determine u^2+v^2+w^2;
        for (i = 0; i < m_spacedim; ++i)
        {
            Vmath::Vvtvp(nq, vel[i], 1, vel[i], 1, velsq, 1, velsq, 1);
            Vmath::Vdiv(nq, vel[i], 1, dirFlds[0], 1, velOverRho[i], 1);
            Vmath::Neg(nq, &velOverRho[i][0], 1);
        }
        
        Vmath::Smul(nq, 0.5, &velsq[0], 1, &HalfVelsq[0], 1);
        
        Vmath::Vdiv(nq, &ones[0], 1, &dirFlds[0][0], 1, &oneOverRho[0], 1);
        
        if (m_spacedim == 2)
        {
            
            /* =================================================================
             if
             M^-1 =
             [             1,             0,            0,            0]
             [        -u/rho,         1/rho,            0,            0]
             [        -v/rho,             0,        1/rho,            0]
             [(gamma-1)U^2/2,   (1-gamma)*u,  (1-gamma)*v,    (gamma-1)]
             
             then
             
             dUdVInvdX = dM^-1/dx =
             [                    0,                 0,              0,  0]
             [         d(-u/rho)/dx,       d(1/rho)/dx,              0,  0]
             [         d(-v/rho)/dx,                 0,    d(1/rho)/dx,  0]
             [(gamma-1)/2*d(U^2)/dx,   (1-gamma)*du/dx,(1-gamma)*dv/dx,  0]
             
             
             dUdVInvdY = dM^-1/dy =
             [                    0,                 0,              0,  0]
             [         d(-u/rho)/dy,       d(1/rho)/dy,              0,  0]
             [         d(-v/rho)/dy,                 0,    d(1/rho)/dy,  0]
             [(gamma-1)/2*d(U^2)/dy,   (1-gamma)*du/dy,(1-gamma)*dv/dy,  0]
             */
            // =================================================================
            /*
             Array<OneD, NekDouble> tmp(nq,0.0);
             */
            for (i = 0; i < nvar; ++i)
            {
                for (j = 0; j < nvar; ++j)
                {
                    m_primal[0]->PhysDeriv(dUdVInv[i][j],
                                           dUdVInvdX[i][j],
                                           dUdVInvdY[i][j]);
                }
            }
        }
        
        // Needs to be checked still
        if (m_spacedim == 3)
        {
            ASSERTL0(false, "3D adjoint solver is not yet implemented");
        }
    }
    
    
    void CompressibleFlowSystem::GetJacobianConvFluxPrimitiveVar(
        Array<OneD, Array<OneD, Array<OneD, Array<OneD, NekDouble> > > > &Jac)
    {
        int i;
        
        int nq = Jac[0][0][0].num_elements();
        int nvar = m_fields.num_elements();
        NekDouble GamOverGamMinOneScalar = m_gamma / (m_gamma - 1.0);
        
        Array<OneD, NekDouble> pressure(nq, 0.0);
        Array<OneD, NekDouble> G_over_GminOne(nq,GamOverGamMinOneScalar);
        Array<OneD, NekDouble> H(nq, 0.0);
        Array<OneD, NekDouble> tmp(nq, 0.0);
        Array<OneD, NekDouble> tmp1(nq, 0.0);
        Array<OneD, NekDouble> tmp2(nq, 0.0);
        Array<OneD, NekDouble> tmp3(nq, 0.0);
        Array<OneD, NekDouble> zeros(nq, 0.0);
        Array<OneD, NekDouble> ones(nq, 1.0);
        Array<OneD, NekDouble> Htmp(nq, 0.0);
        Array<OneD, NekDouble> velsq(nq, 0.0);
        Array<OneD, Array<OneD, NekDouble> > vel(m_spacedim);
        Array<OneD, Array<OneD, NekDouble> > dirFlds(nvar);
        
        Array<OneD, Array<OneD, NekDouble> > tmpx(nvar);
        Array<OneD, Array<OneD, NekDouble> > tmpy(nvar);
        
        for (i = 0; i < nvar; ++i)
        {
            dirFlds[i] = Array<OneD, NekDouble>(nq, 0.0);
        }
        
        GetDirectSolution(dirFlds);
        
        for (i = 0; i < m_spacedim; ++i)
        {
            vel[i]          = Array<OneD, NekDouble>(nq);
        }
        
        GetVelocityVector(dirFlds, vel);
        GetPressure      (dirFlds, vel, pressure);
        
        for (i = 0; i < m_spacedim; ++i)
        {
            Vmath::Vvtvp(nq, vel[i], 1, vel[i], 1, velsq, 1, velsq, 1);
        }
        
        // determine H = E + p/rho;
        Vmath::Vdiv(nq, dirFlds[nvar-1], 1, dirFlds[0], 1, H, 1);
        
        Vmath::Vdiv(nq, pressure, 1, dirFlds[0], 1, Htmp, 1);
        
        Vmath::Vadd(nq, H, 1, Htmp, 1, H, 1);
        
        if (m_spacedim == 1)
        {
            ASSERTL0(false, "1D adjoint solver is not yet implemented");
        }
        if (m_spacedim == 2)
        {
            // construction of the Ax Matrix
            /* =================================================================
             Jac[0][nvar][nvar][nq] = Ax^t =
             [  u,      u^2,          u*v,            (u*(u^2 + v^2))/2]
             [  0,    rho*u,      2*rho*v,            rho*u^2 + rho * H]
             [rho,        0,        rho*v,                      rho*u*v]
             [  0,        1,            0,            (gam*v)/(gam - 1)]
             */
            // =================================================================
            // row one
            Vmath::Vcopy(nq, &vel[0][0], 1, &Jac[0][0][0][0], 1);
            // u^2
            Vmath::Vmul(nq, &vel[0][0], 1, &vel[0][0], 1, &Jac[0][0][1][0], 1);
            // u * v
            Vmath::Vmul(nq, &vel[0][0], 1, &vel[1][0], 1, &Jac[0][0][2][0], 1);
            // 1/2 * u * (u^2 + v^2)
            Vmath::Vmul(nq, &vel[0][0], 1, &velsq[0], 1, &Jac[0][0][3][0], 1);
            Vmath::Smul(nq, 0.5, &Jac[0][0][3][0], 1, &Jac[0][0][3][0], 1);
            
            // row two
            // rho
            Vmath::Vcopy(nq, &dirFlds[0][0], 1, &Jac[0][1][0][0], 1);
            // 2 * rho * u
            Vmath::Vmul(nq, &dirFlds[0][0], 1, &vel[0][0], 1, &Jac[0][1][1][0], 1);
            Vmath::Smul(nq, 2.0, &Jac[0][1][1][0], 1, &Jac[0][1][1][0], 1);
            // rho * v
            Vmath::Vmul(nq, &dirFlds[0][0], 1, &vel[1][0], 1, &Jac[0][1][2][0], 1);
            // rho * (H + u^2)
            Vmath::Zero(nq, &tmp[0], 1);
            Vmath::Vmul(nq, &vel[0][0], 1, &vel[0][0], 1, &tmp[0], 1);
            Vmath::Vadd(nq, &tmp[0], 1, &H[0], 1, &tmp[0], 1);
            Vmath::Vmul(nq, &dirFlds[0][0], 1, &tmp[0], 1, &Jac[0][1][3][0], 1);
            
            // row three
            // 0
            Vmath::Vcopy(nq, &zeros[0], 1, &Jac[0][2][0][0], 1);
            // 0
            Vmath::Vcopy(nq, &zeros[0], 1, &Jac[0][2][1][0], 1);
            // rho * u
            Vmath::Vmul(nq, &dirFlds[0][0], 1, &vel[0][0], 1, &Jac[0][2][2][0], 1);
            // rho * u * v
            Vmath::Vmul(nq, &dirFlds[0][0], 1, &vel[0][0], 1, &Jac[0][2][3][0], 1);
            Vmath::Vmul(nq, &Jac[0][2][3][0], 1, &vel[1][0], 1, &Jac[0][2][3][0], 1);
            
            // row four
            // 0
            Vmath::Vcopy(nq, &zeros[0], 1, &Jac[0][3][0][0], 1);
            // 1
            Vmath::Vcopy(nq, &ones[0], 1, &Jac[0][3][1][0], 1);
            // 0
            Vmath::Vcopy(nq, &zeros[0], 1, &Jac[0][3][2][0], 1);
            // u * gamma / (gamma - 1.0)
            Vmath::Vmul(nq, &vel[0][0], 1, &G_over_GminOne[0], 1, &Jac[0][3][3][0], 1);
            
            // construction of the Ay Matrix
            /* =================================================================
             Jac[1][nvar][nvar][nq] = Ay^t =
             [  v,      u*v,          v^2,            (v*(u^2 + v^2))/2]
             [  0,    rho*v,            0,                      rho*u*v]
             [rho,    rho*u,      2*rho*v,            rho*v^2 + rho * H]
             [  0,        0,            1,            (gam*v)/(gam - 1)]
             */
            // =================================================================
            
            Vmath::Zero(nq, &tmp[0], 1);
            // row one
            // v
            Vmath::Vcopy(nq, &vel[1][0], 1, &Jac[1][0][0][0], 1);
            // u * v
            Vmath::Vmul(nq, &vel[1][0], 1, &vel[0][0], 1, &Jac[1][0][1][0], 1);
            // v^2
            Vmath::Vmul(nq, &vel[1][0], 1, &vel[1][0], 1, &Jac[1][0][2][0], 1);
            // 1/2 * v * (u^2 + v^2)
            Vmath::Vmul(nq, &vel[1][0], 1, &velsq[0], 1, &Jac[1][0][3][0], 1);
            Vmath::Smul(nq, 0.5, &Jac[1][0][3][0], 1, &Jac[1][0][3][0], 1);
            
            // row two
            // 0
            Vmath::Vcopy(nq, &zeros[0], 1, &Jac[1][1][0][0], 1);
            // rho * v
            Vmath::Vmul(nq, &dirFlds[0][0], 1, &vel[1][0], 1, &Jac[1][1][1][0], 1);
            // 0
            Vmath::Vcopy(nq, &zeros[0], 1, &Jac[1][1][2][0], 1);
            // rho * u * v
            Vmath::Vmul(nq, &dirFlds[0][0], 1, &vel[0][0], 1, &Jac[1][1][3][0], 1);
            Vmath::Vmul(nq, &Jac[1][1][3][0], 1, &vel[1][0], 1, &Jac[1][1][3][0], 1);
            
            // row three
            // rho
            Vmath::Vcopy(nq, &dirFlds[0][0], 1, &Jac[1][2][0][0], 1);
            // rho * u
            Vmath::Vmul(nq, &dirFlds[0][0], 1, &vel[0][0], 1, &Jac[1][2][1][0], 1);
            // 2 * rho * v
            Vmath::Vmul(nq, &dirFlds[0][0], 1, &vel[1][0], 1, &Jac[1][2][2][0], 1);
            Vmath::Smul(nq, 2.0, &Jac[1][2][2][0], 1, &Jac[1][2][2][0], 1);
            // rho * (H + v^2)
            Vmath::Zero(nq, &tmp[0], 1);
            Vmath::Vmul(nq, &vel[1][0], 1, &vel[1][0], 1, &tmp[0], 1);
            Vmath::Vadd(nq, &tmp[0], 1, &H[0], 1, &tmp[0], 1);
            Vmath::Vmul(nq, &dirFlds[0][0], 1, &tmp[0], 1, &Jac[1][2][3][0], 1);
            
            // row four
            // 0
            Vmath::Vcopy(nq, &zeros[0], 1, &Jac[1][3][0][0], 1);
            // 0
            Vmath::Vcopy(nq, &zeros[0], 1, &Jac[1][3][1][0], 1);
            // 1
            Vmath::Vcopy(nq, &ones[0], 1, &Jac[1][3][2][0], 1);
            // v * gamma / (gamma - 1.0)
            Vmath::Vmul(nq, &vel[1][0], 1, &G_over_GminOne[0], 1, &Jac[1][3][3][0], 1);
        }
        if (m_spacedim == 3)
        {
            ASSERTL0(false, "3D adjoint solver is not yet implemented");
        }
    }
    
    void CompressibleFlowSystem::GetAdjointConvFluxVectorFromPrimitiveVar(
                                                                          const Array<OneD, Array<OneD, NekDouble> > &inarray,
                                                                          Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &outarray)
    {
        int i,j,k;
        
        int nq = inarray[0].num_elements();
        int nvar = m_fields.num_elements();
        
        Array<OneD, NekDouble> tmp(nq, 0.0);
        Array<OneD, NekDouble> tmp1(nq, 0.0);
        Array<OneD, Array<OneD, NekDouble> > tmpx(nvar);
        Array<OneD, Array<OneD, NekDouble> > tmpy(nvar);
        
        
        Array<OneD, Array<OneD, Array<OneD,
        Array<OneD, NekDouble > > > > Jac(m_spacedim);
        Array<OneD, Array<OneD, Array<OneD,
        Array<OneD, NekDouble > > > >JacCons(m_spacedim);
        
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > > dVdU(nvar);
        
        // span the required arrays
        for (i = 0; i < m_spacedim; ++i)
        {
            Jac[i]     = Array<OneD, Array<OneD,
            Array<OneD, NekDouble > > >(nvar);
            JacCons[i] = Array<OneD, Array<OneD,
            Array<OneD, NekDouble > > >(nvar);
            
            for (j = 0; j < nvar; ++j)
            {
                Jac[i][j] = Array<OneD, Array<OneD, NekDouble> >(nvar);
                JacCons[i][j] = Array<OneD, Array<OneD, NekDouble> >(nvar);
                
                for (k = 0; k < nvar; ++k)
                {
                    Jac[i][j][k] = Array<OneD, NekDouble>(nq, 0.0);
                    JacCons[i][j][k] = Array<OneD, NekDouble>(nq, 0.0);
                }
            }
        }
        
        for (i = 0; i < nvar; ++i)
        {
            dVdU[i] = Array<OneD, Array<OneD, NekDouble > >(nvar);
            
            for (j = 0; j < nvar; ++j)
            {
                dVdU[i][j] = Array<OneD, NekDouble>(nq, 0.0);
            }
        }
        
        GetJacobianConvFluxPrimitiveVar(Jac);
        
        GetConservToPrimVariableInvMat(dVdU);
        
        // Perform tensor product dVdU*dFdV
        for (i = 0; i < nvar; ++i)
        {
            for (j = 0; j < nvar; ++j)
            {
                Vmath::Zero(nq, &tmp[0], 1);
                Vmath::Zero(nq, &tmp1[0], 1);
                
                for (k = 0; k < nvar; ++k)
                {
                    Vmath::Vvtvp(nq,
                                 &dVdU[i][k][0], 1,
                                 &Jac[0][k][j][0], 1,
                                 &tmp[0], 1,
                                 &tmp[0], 1);
                    
                    Vmath::Vcopy(nq,
                                 &tmp[0], 1,
                                 &JacCons[0][i][j][0], 1);
                    
                    Vmath::Vvtvp(nq,
                                 &dVdU[i][k][0], 1,
                                 &Jac[1][k][j][0], 1,
                                 &tmp1[0], 1,
                                 &tmp1[0], 1);
                    
                    Vmath::Vcopy(nq,
                                 &tmp1[0], 1,
                                 &JacCons[1][i][j][0], 1);
                }
            }
        }
        
        for (i = 0; i < nvar; ++i)
        {
            tmpx[i] = Array<OneD, NekDouble>(nq, 0.0);
            tmpy[i] = Array<OneD, NekDouble>(nq, 0.0);
            
            for (j = 0; j < nvar; ++j)
            {
                if (m_spacedim == 2)
                {
                    Vmath::Vvtvp(nq,
                                 &JacCons[0][i][j][0], 1,
                                 &inarray[j][0], 1,
                                 &tmpx[i][0], 1,
                                 &tmpx[i][0], 1);
                    
                    Vmath::Vvtvp(nq,
                                 &JacCons[1][i][j][0], 1,
                                 &inarray[j][0], 1,
                                 &tmpy[i][0], 1,
                                 &tmpy[i][0], 1);
                    
                }
                
                if (m_spacedim == 3)
                {
                    ASSERTL0(false, "3D adjoint solver is ot implemented yet")
                }
                
            }
            
            Vmath::Vcopy(nq, &tmpx[i][0], 1, &outarray[i][0][0], 1);
            Vmath::Vcopy(nq, &tmpy[i][0], 1, &outarray[i][1][0], 1);
        }
    }
    
    
    void CompressibleFlowSystem::GetDerivJacVectorFromPrimitiveVar(
                const Array<OneD, Array<OneD, NekDouble> > &inarray,
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &outarray)
    {
        int i,j,k;
        
        int nq = inarray[0].num_elements();
        int nvar = m_fields.num_elements();
        
        Array<OneD, NekDouble> tmp(nq, 0.0);
        Array<OneD, NekDouble> tmp1(nq, 0.0);
        Array<OneD, NekDouble> tmp2(nq, 0.0);
        Array<OneD, NekDouble> tmp3(nq, 0.0);
        Array<OneD, NekDouble> tmp4(nq, 0.0);
        Array<OneD, NekDouble> tmp5(nq, 0.0);
        
        Array<OneD, Array<OneD, NekDouble> > tmpx(nvar);
        Array<OneD, Array<OneD, NekDouble> > tmpy(nvar);
        
        Array<OneD, Array<OneD, Array<OneD,
        Array<OneD, NekDouble > > > > Jac(m_spacedim);
        Array<OneD, Array<OneD, Array<OneD,
        Array<OneD, NekDouble > > > > JacCons(m_spacedim);
        Array<OneD, Array<OneD, Array<OneD,
        Array<OneD, NekDouble > > > > JacDiv(m_spacedim);
        Array<OneD, Array<OneD, Array<OneD,
        Array<OneD, NekDouble > > > >JacConsD(m_spacedim);
        Array<OneD, Array<OneD, Array<OneD,
        Array<OneD, NekDouble > > > >JacConsDiv(m_spacedim);
        
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > > dVdU(nvar);
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > > dVdUdX(nvar);
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > > dVdUdY(nvar);
        
        // span the required arrays
        for (i = 0; i < m_spacedim; ++i)
        {
            Jac[i]     = Array<OneD, Array<OneD,
            Array<OneD, NekDouble > > >(nvar);
            JacDiv[i]     = Array<OneD, Array<OneD,
            Array<OneD, NekDouble > > >(nvar);
            JacConsDiv[i] = Array<OneD, Array<OneD,
            Array<OneD, NekDouble > > >(nvar);
            
            for (j = 0; j < nvar; ++j)
            {
                Jac[i][j] = Array<OneD, Array<OneD, NekDouble> >(nvar);
                JacDiv[i][j] = Array<OneD, Array<OneD, NekDouble> >(nvar);
                JacConsDiv[i][j] = Array<OneD, Array<OneD, NekDouble> >(nvar);
                
                for (k = 0; k < nvar; ++k)
                {
                    Jac[i][j][k] = Array<OneD, NekDouble>(nq, 0.0);
                    JacDiv[i][j][k] = Array<OneD, NekDouble>(nq, 0.0);
                    JacConsDiv[i][j][k] = Array<OneD, NekDouble>(nq, 0.0);
                }
            }
        }
        
        for (i = 0; i < nvar; ++i)
        {
            dVdU[i] = Array<OneD, Array<OneD, NekDouble > >(nvar);
            dVdUdX[i] = Array<OneD, Array<OneD, NekDouble > >(nvar);
            dVdUdY[i] = Array<OneD, Array<OneD, NekDouble > >(nvar);
            
            for (j = 0; j < nvar; ++j)
            {
                dVdU[i][j] = Array<OneD, NekDouble>(nq, 0.0);
                dVdUdX[i][j] = Array<OneD, NekDouble>(nq, 0.0);
                dVdUdY[i][j] = Array<OneD, NekDouble>(nq, 0.0);
            }
        }
        
        GetJacobianConvFluxPrimitiveVar(Jac);
        
        GetDerivJacobian(Jac, JacDiv);
        
        GetConservToPrimVariableInvMat(dVdU);
        
        GetConservToPrimVariableInvMatDiv(dVdU, dVdUdX, dVdUdY);
        
        for (i = 0; i < nvar; ++i)
        {
            for (j = 0; j < nvar; ++j)
            {
                Vmath::Zero(nq, &tmp[0], 1);
                Vmath::Zero(nq, &tmp1[0], 1);
                Vmath::Zero(nq, &tmp2[0], 1);
                Vmath::Zero(nq, &tmp3[0], 1);
                Vmath::Zero(nq, &tmp4[0], 1);
                Vmath::Zero(nq, &tmp5[0], 1);
                
                for (k = 0; k < nvar; ++k)
                {
                    //
                    Vmath::Vvtvp(nq,
                                 &dVdUdX[i][k][0], 1,
                                 &Jac[0][k][j][0], 1,
                                 &tmp[0], 1,
                                 &tmp[0], 1);
                    
                    Vmath::Vvtvp(nq,
                                 &dVdUdY[i][k][0], 1,
                                 &Jac[1][k][j][0], 1,
                                 &tmp1[0], 1,
                                 &tmp1[0], 1);
                    
                    //=============
                    
                    Vmath::Vvtvp(nq,
                                 &dVdU[i][k][0], 1,
                                 &JacDiv[0][k][j][0], 1,
                                 &tmp2[0], 1,
                                 &tmp2[0], 1);
                    
                    Vmath::Vvtvp(nq,
                                 &dVdU[i][k][0], 1,
                                 &JacDiv[1][k][j][0], 1,
                                 &tmp3[0], 1,
                                 &tmp3[0], 1);
                    
                    Vmath::Vadd(nq,
                                &tmp[0], 1,
                                &tmp2[0], 1,
                                &JacConsDiv[0][i][j][0], 1);
                    
                    Vmath::Vadd(nq,
                                &tmp1[0], 1,
                                &tmp3[0], 1,
                                &JacConsDiv[1][i][j][0], 1);
                }
            }
        }
        
        for (i = 0; i < nvar; ++i)
        {
            tmpx[i] = Array<OneD, NekDouble>(nq, 0.0);
            tmpy[i] = Array<OneD, NekDouble>(nq, 0.0);
            
            for (j = 0; j < nvar; ++j)
            {
                Vmath::Vvtvp(nq,
                             &JacConsDiv[0][i][j][0], 1,
                             &inarray[j][0], 1,
                             &tmpx[i][0], 1,
                             &tmpx[i][0], 1);
                
                Vmath::Vvtvp(nq,
                             &JacConsDiv[1][i][j][0], 1,
                             &inarray[j][0], 1,
                             &tmpy[i][0], 1,
                             &tmpy[i][0], 1);
                
            }
            
            Vmath::Vcopy(nq, &tmpx[i][0], 1, &outarray[i][0][0], 1);
            Vmath::Vcopy(nq, &tmpy[i][0], 1, &outarray[i][1][0], 1);
        }
    }
    
    void CompressibleFlowSystem::GetAdjointViscousFluxVectorFromPrimitiveVar(
          const Array<OneD, Array<OneD, NekDouble> > &inarray,
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &derivatives,
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &outarray)
    {
        int i,j,k;
        
        int nq = inarray[0].num_elements();
        int nvar = m_fields.num_elements();
        
        Array<OneD, NekDouble> tmp(nq, 0.0);
        Array<OneD, NekDouble> tmp1(nq, 0.0);
        Array<OneD, NekDouble> tmp2(nq, 0.0);
        Array<OneD, NekDouble> tmp3(nq, 0.0);
        
        Array<OneD, Array<OneD, NekDouble> > tmpx1(nvar);
        Array<OneD, Array<OneD, NekDouble> > tmpy1(nvar);
        Array<OneD, Array<OneD, NekDouble> > tmpx2(nvar);
        Array<OneD, Array<OneD, NekDouble> > tmpy2(nvar);
        
        Array<OneD, Array<OneD, NekDouble> > tmpxtot(nvar);
        Array<OneD, Array<OneD, NekDouble> > tmpytot(nvar);
        
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > > dVdU(nvar);
        
        for (i = 0; i < nvar; ++i)
        {
            dVdU[i] = Array<OneD, Array<OneD, NekDouble > >(nvar);
            
            for (j = 0; j < nvar; ++j)
            {
                dVdU[i][j] = Array<OneD, NekDouble>(nq, 0.0);
            }
        }
        
        Array<OneD, Array<OneD, NekDouble> > dirFlds(nvar);
        for (i = 0; i < nvar; ++i)
        {
            dirFlds[i] = Array<OneD, NekDouble>(nq, 0.0);
        }
        
        GetDirectSolution(dirFlds);
        
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > > Dxx(nvar);
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > > Dyx(nvar);
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > > Dxy(nvar);
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > > Dyy(nvar);
        
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > > DxxCons(nvar);
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > > DyxCons(nvar);
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > > DxyCons(nvar);
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > > DyyCons(nvar);
        
        for (i = 0; i < nvar; ++i)
        {
            Dxx[i] = Array<OneD, Array<OneD, NekDouble> >(nvar);
            Dyx[i] = Array<OneD, Array<OneD, NekDouble> >(nvar);
            Dxy[i] = Array<OneD, Array<OneD, NekDouble> >(nvar);
            Dyy[i] = Array<OneD, Array<OneD, NekDouble> >(nvar);
            
            DxxCons[i] = Array<OneD, Array<OneD, NekDouble> >(nvar);
            DyxCons[i] = Array<OneD, Array<OneD, NekDouble> >(nvar);
            DxyCons[i] = Array<OneD, Array<OneD, NekDouble> >(nvar);
            DyyCons[i] = Array<OneD, Array<OneD, NekDouble> >(nvar);
            
            for (j = 0; j < nvar; ++j)
            {
                Dxx[i][j] = Array<OneD, NekDouble>(nq, 0.0);
                Dyx[i][j] = Array<OneD, NekDouble>(nq, 0.0);
                Dxy[i][j] = Array<OneD, NekDouble>(nq, 0.0);
                Dyy[i][j] = Array<OneD, NekDouble>(nq, 0.0);
                
                DxxCons[i][j] = Array<OneD, NekDouble>(nq, 0.0);
                DyxCons[i][j] = Array<OneD, NekDouble>(nq, 0.0);
                DxyCons[i][j] = Array<OneD, NekDouble>(nq, 0.0);
                DyyCons[i][j] = Array<OneD, NekDouble>(nq, 0.0);
            }
        }
        
        GetJacobianViscousFluxFromDiffPrimitiveVar(dirFlds,
                                                   Dxx, Dxy,
                                                   Dyx, Dyy);
        
        GetConservToPrimVariableInvMat(dVdU);
        
        // Perform tensor product dVdU*dFdV
        /*for (i = 0; i < nvar; ++i)
         {
         for (j = 0; j < nvar; ++j)
         {
         Vmath::Zero(nq, &tmp[0], 1);
         Vmath::Zero(nq, &tmp1[0], 1);
         Vmath::Zero(nq, &tmp2[0], 1);
         Vmath::Zero(nq, &tmp3[0], 1);
         
         for (k = 0; k < nvar; ++k)
         {
         Vmath::Vvtvp(nq,
         &dVdU[i][k][0], 1,
         &Dxx[k][j][0], 1,
         &tmp[0], 1,
         &tmp[0], 1);
         
         Vmath::Vcopy(nq,
         &tmp[0], 1,
         &DxxCons[i][j][0], 1);
         
         Vmath::Vvtvp(nq,
         &dVdU[i][k][0], 1,
         &Dyx[k][j][0], 1,
         &tmp1[0], 1,
         &tmp1[0], 1);
         
         Vmath::Vcopy(nq,
         &tmp1[0], 1,
         &DyxCons[i][j][0], 1);
         
         Vmath::Vvtvp(nq,
         &dVdU[i][k][0], 1,
         &Dxy[k][j][0], 1,
         &tmp2[0], 1,
         &tmp2[0], 1);
         
         Vmath::Vcopy(nq,
         &tmp2[0], 1,
         &DxyCons[i][j][0], 1);
         
         Vmath::Vvtvp(nq,
         &dVdU[i][k][0], 1,
         &Dyy[k][j][0], 1,
         &tmp3[0], 1,
         &tmp3[0], 1);
         
         Vmath::Vcopy(nq,
         &tmp3[0], 1,
         &DyyCons[i][j][0], 1);
         }
         }
         }
         
         cout << "====================================================" << endl;
         
         cout << Vmath::Vmax(nq, &DxxCons[0][0][0], 1) << " "
         << Vmath::Vmax(nq, &DxxCons[0][1][0], 1) << " "
         << Vmath::Vmax(nq, &DxxCons[0][2][0], 1) << " "
         << Vmath::Vmax(nq, &DxxCons[0][3][0], 1) << endl;
         cout << Vmath::Vmax(nq, &DxxCons[1][0][0], 1) << " "
         << Vmath::Vmax(nq, &DxxCons[1][1][0], 1) << " "
         << Vmath::Vmax(nq, &DxxCons[1][2][0], 1) << " "
         << Vmath::Vmax(nq, &DxxCons[1][3][0], 1) << endl;
         cout << Vmath::Vmax(nq, &DxxCons[2][0][0], 1) << " "
         << Vmath::Vmax(nq, &DxxCons[2][1][0], 1) << " "
         << Vmath::Vmax(nq, &DxxCons[2][2][0], 1) << " "
         << Vmath::Vmax(nq, &DxxCons[2][3][0], 1) << endl;
         cout << Vmath::Vmax(nq, &DxxCons[3][0][0], 1) << " "
         << Vmath::Vmax(nq, &DxxCons[3][1][0], 1) <<  " "
         << Vmath::Vmax(nq, &DxxCons[3][2][0], 1) << " "
         << Vmath::Vmax(nq, &DxxCons[3][3][0], 1) << endl;
         
         cout << "====================================================" << endl;
         */
        if (m_spacedim == 2)
        {
            for (i = 0; i < nvar; ++i)
            {
                tmpx1[i] = Array<OneD, NekDouble>(nq, 0.0);
                tmpy1[i] = Array<OneD, NekDouble>(nq, 0.0);
                
                tmpx2[i] = Array<OneD, NekDouble>(nq, 0.0);
                tmpy2[i] = Array<OneD, NekDouble>(nq, 0.0);
                
                for (j = 0; j < nvar; ++j)
                {
                    Vmath::Vvtvp(nq,
                                 &Dxx[i][j][0], 1,
                                 &derivatives[0][j][0], 1,
                                 &tmpx1[i][0], 1,
                                 &tmpx1[i][0], 1);
                    
                    Vmath::Vvtvp(nq,
                                 &Dyx[i][j][0], 1,
                                 &derivatives[1][j][0], 1,
                                 &tmpx2[i][0], 1,
                                 &tmpx2[i][0], 1);
                    
                    Vmath::Vvtvp(nq,
                                 &Dxy[i][j][0], 1,
                                 &derivatives[0][j][0], 1,
                                 &tmpy1[i][0], 1,
                                 &tmpy1[i][0], 1);
                    
                    Vmath::Vvtvp(nq,
                                 &Dyy[i][j][0], 1,
                                 &derivatives[1][j][0], 1,
                                 &tmpy2[i][0], 1,
                                 &tmpy2[i][0], 1);
                }
            }
            
            Array<OneD, Array<OneD, NekDouble> > tmpnewx(nvar);
            Array<OneD, Array<OneD, NekDouble> > tmpnewy(nvar);
            
            for (k = 0; k < nvar; ++k)
            {
                tmpnewx[k] = Array<OneD, NekDouble>(nq, 0.0);
                tmpnewy[k] = Array<OneD, NekDouble>(nq, 0.0);
                Vmath::Vadd(nq,
                            &tmpx1[k][0], 1,
                            &tmpx2[k][0], 1,
                            &tmpnewx[k][0], 1);
                
                Vmath::Vadd(nq,
                            &tmpy1[k][0], 1,
                            &tmpy2[k][0], 1,
                            &tmpnewy[k][0], 1);
            }
            
            
            
            for (i = 0; i < nvar; ++i)
            {
                
                
                for (j = 0; j < nvar; ++j)
                {
                    Vmath::Vvtvp(nq,
                                 &dVdU[i][j][0], 1,
                                 &tmpnewx[j][0], 1,
                                 &outarray[0][i][0], 1,
                                 &outarray[0][i][0], 1);
                    
                    Vmath::Vvtvp(nq,
                                 &dVdU[i][j][0], 1,
                                 &tmpnewy[j][0], 1,
                                 &outarray[1][i][0], 1,
                                 &outarray[1][i][0], 1);
                }
            }
        }
        if (m_spacedim == 3)
        {
            ASSERTL0(false, "3D adjoint solver is not yet implemented");
        }
    }
    
    
    void CompressibleFlowSystem::GetJacobianViscousFluxFromDiffPrimitiveVar(
                        Array<OneD, Array<OneD, NekDouble> > inarray,
                        Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &Dxx,
                        Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &Dxy,
                        Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &Dyx,
                        Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &Dyy)
    {
        // Jac[m_spacedim][m_spacedim][nvar][nvar][nq];
        int i, j, k;
        int nq = inarray[0].num_elements();
        int nvar = m_fields.num_elements();
        
        Array<OneD, NekDouble> tmp(nq, 0.0);
        Array<OneD, NekDouble> zeros(nq, 0.0);
        Array<OneD, NekDouble> ones(nq, 1.0);
        Array<OneD, NekDouble> muVec(nq, m_mu);
        Array<OneD, NekDouble> Prandtl(nq, 0.0);
        Array<OneD, NekDouble> Pressure(nq, 0.0);
        Array<OneD, NekDouble> kPrho2Vec(nq, 0.0);
        
        Array<OneD, Array<OneD, NekDouble> > vel(m_spacedim);
        Array<OneD, Array<OneD, NekDouble> > MuVel(m_spacedim);
        Array<OneD, Array<OneD, NekDouble> > tmpx(nvar);
        Array<OneD, Array<OneD, NekDouble> > tmpy(nvar);
        
        for (i = 0; i < m_spacedim; ++i)
        {
            vel[i]          = Array<OneD, NekDouble>(nq, 0.0);
            MuVel[i]          = Array<OneD, NekDouble>(nq, 0.0);
        }
        
        GetVelocityVector(inarray, vel);
        
        for (i = 0; i < m_spacedim; ++i)
        {
            Vmath::Smul(nq, m_mu, &vel[i][0], 1, &MuVel[i][0], 1);
        }
        
        GetPressure(inarray, vel, Pressure);
        
        NekDouble FourOThreeMuSclr = 4.0/3.0*m_mu;
        NekDouble TwoOThreeMuSclr  = 2.0/3.0*m_mu;
        NekDouble GamMu_O_GamMinOnePrSclr = m_gamma*m_mu/((m_gamma-1)*m_Prandtl);
        
        Array<OneD, NekDouble> FourOThreeMuVector(nq, FourOThreeMuSclr);
        Array<OneD, NekDouble> TwoOThreeMuVector(nq, TwoOThreeMuSclr);
        Array<OneD, NekDouble> GamMu_O_GamMinOnePrVector(nq, GamMu_O_GamMinOnePrSclr);
        
        Vmath::Vmul(nq,
                    &GamMu_O_GamMinOnePrVector[0], 1,
                    &Pressure[0], 1,
                    &kPrho2Vec[0], 1);
        
        Vmath::Vmul(nq, &inarray[0][0], 1, &inarray[0][0], 1, &tmp[0], 1);
        Vmath::Vdiv(nq, &kPrho2Vec[0], 1, &tmp[0], 1, &kPrho2Vec[0], 1);
        
        Vmath::Zero(nq, &tmp[0], 1);
        
        
        if (m_spacedim == 2)
        {
            /* =============================================================
             Dxx =
             
             [ 0,        0,  0,    -gamma*mu/(Pr*(gamma-1))*p/rho^2]
             [ 0,   4/3*mu,  0,                            4/3*u*mu]
             [ 0,        0, mu,                                v*mu]
             [ 0,        0,  0,         gamma/(Pr*(gamma-1))*mu/rho]
             
             =============================================================*/
            // row 1;
            Vmath::Zero(nq, &Dxx[0][0][0], 1);
            Vmath::Zero(nq, &Dxx[0][1][0], 1);
            Vmath::Zero(nq, &Dxx[0][2][0], 1);
            Vmath::Vcopy(nq, &kPrho2Vec[0], 1, &Dxx[0][3][0], 1);
            Vmath::Neg(nq, &Dxx[0][3][0], 1);
            
            // row 2;
            Vmath::Zero(nq, &Dxx[1][0][0], 1);
            Vmath::Vcopy(nq, &FourOThreeMuVector[0], 1, &Dxx[1][1][0], 1);
            Vmath::Zero(nq, &Dxx[1][2][0], 1);
            Vmath::Smul(nq,
                        FourOThreeMuSclr,
                        &vel[0][0], 1,
                        &Dxx[1][3][0], 1);
            
            // row 3;
            Vmath::Zero(nq, &Dxx[2][0][0], 1);
            Vmath::Zero(nq, &Dxx[2][1][0], 1);
            Vmath::Vcopy(nq, &muVec[0], 1, &Dxx[2][2][0], 1);
            Vmath::Vcopy(nq, &MuVel[1][0], 1, &Dxx[2][3][0], 1);
            
            // row 4;
            Vmath::Zero(nq, &Dxx[3][0][0], 1);
            Vmath::Zero(nq, &Dxx[3][1][0], 1);
            Vmath::Zero(nq, &Dxx[3][2][0], 1);
            Vmath::Vdiv(nq,
                        &ones[0], 1,
                        &inarray[0][0], 1,
                        &Dxx[3][3][0], 1);
            
            Array<OneD, NekDouble> rhotest(nq, 0.0);
            Vmath::Vcopy(nq, &inarray[0][0], 1, &rhotest[0], 1);
            Vmath::Vabs(nq,  &rhotest[0], 1,  &rhotest[0], 1);
            
            Vmath::Vmul(nq,
                        &GamMu_O_GamMinOnePrVector[0], 1,
                        &Dxx[3][3][0], 1,
                        &Dxx[3][3][0], 1);
            
            /* =============================================================
             Dyx =
             
             [ 0,        0,       0,            0]
             [ 0,        0, -2/3*mu,    -2/3*v*mu]
             [ 0,       mu,       0,         u*mu]
             [ 0,        0,       0,            0]
             
             =============================================================*/
            // row 1;
            Vmath::Zero(nq, &Dyx[0][0][0], 1);
            Vmath::Zero(nq, &Dyx[0][1][0], 1);
            Vmath::Zero(nq, &Dyx[0][2][0], 1);
            Vmath::Zero(nq, &Dyx[0][3][0], 1);
            
            // row 2;
            Vmath::Zero(nq, &Dyx[1][0][0], 1);
            Vmath::Zero(nq, &Dyx[1][1][0], 1);
            Vmath::Vcopy(nq, &TwoOThreeMuVector[0], 1, &Dyx[1][2][0], 1);
            Vmath::Neg(nq, &Dyx[1][2][0], 1);
            
            Vmath::Vmul(nq,
                        &TwoOThreeMuVector[0], 1,
                        &vel[1][0], 1,
                        &Dyx[1][3][0], 1);
            
            Vmath::Neg(nq, &Dyx[1][3][0], 1);
            
            // row 3;
            Vmath::Zero(nq, &Dyx[2][0][0], 1);
            Vmath::Vcopy(nq, &muVec[0], 1, &Dyx[2][1][0], 1);
            Vmath::Zero(nq, &Dyx[2][2][0], 1);
            Vmath::Vcopy(nq, &MuVel[0][0], 1, &Dyx[2][3][0], 1);
            
            // row 4;
            Vmath::Zero(nq, &Dyx[3][0][0], 1);
            Vmath::Zero(nq, &Dyx[3][1][0], 1);
            Vmath::Zero(nq, &Dyx[3][2][0], 1);
            Vmath::Zero(nq, &Dyx[3][3][0], 1);
            /* =============================================================
             Dxy =
             
             [ 0,        0,       0,            0]
             [ 0,        0,      mu,         v*mu]
             [ 0,       -2/3*mu,  0,    -2/3*u*mu]
             [ 0,        0,       0,            0]
             
             =============================================================*/
            // row 1;
            Vmath::Zero(nq, &Dxy[0][0][0], 1);
            Vmath::Zero(nq, &Dxy[0][1][0], 1);
            Vmath::Zero(nq, &Dxy[0][2][0], 1);
            Vmath::Zero(nq, &Dxy[0][3][0], 1);
            
            // row 2;
            Vmath::Zero(nq, &Dxy[1][0][0], 1);
            Vmath::Zero(nq, &Dxy[1][1][0], 1);
            Vmath::Vcopy(nq, &muVec[0], 1, &Dxy[1][2][0], 1);
            Vmath::Vcopy(nq, &MuVel[1][0], 1, &Dxy[1][3][0], 1);
            
            // row 3;
            Vmath::Zero(nq, &Dxy[2][0][0], 1);
            Vmath::Vcopy(nq, &TwoOThreeMuVector[0], 1, &Dxy[2][1][0], 1);
            Vmath::Neg(nq, &Dxy[2][1][0], 1);
            Vmath::Zero(nq, &Dxy[2][2][0], 1);
            Vmath::Vmul(nq,
                        &TwoOThreeMuVector[0], 1,
                        &vel[1][0], 1,
                        &Dxy[2][3][0], 1);
            
            Vmath::Neg(nq, &Dxy[2][3][0], 1);
            
            // row 4;
            Vmath::Zero(nq, &Dxy[3][0][0], 1);
            Vmath::Zero(nq, &Dxy[3][1][0], 1);
            Vmath::Zero(nq, &Dxy[3][2][0], 1);
            Vmath::Zero(nq, &Dxy[3][3][0], 1);
            
            /* =============================================================
             Dyy =
             
             [ 0,        0,      0,    -gamma*mu/(Pr*(gamma-1))*p/rho^2]
             [ 0,       mu,      0,                                u*mu]
             [ 0,        0, 4/3*mu,                            4/3*v*mu]
             [ 0,        0,      0,         gamma/(Pr*(gamma-1))*mu/rho]
             
             ===============================================================*/
            
            // row 1;
            Vmath::Zero(nq, &Dyy[0][0][0], 1);
            Vmath::Zero(nq, &Dyy[0][1][0], 1);
            Vmath::Zero(nq, &Dyy[0][2][0], 1);
            Vmath::Vcopy(nq, &kPrho2Vec[0], 1, &Dyy[0][3][0], 1);
            Vmath::Neg(nq, &Dyy[0][3][0], 1);
            
            // row 2;
            Vmath::Zero(nq, &Dyy[1][0][0], 1);
            Vmath::Vcopy(nq, &muVec[0], 1, &Dyy[1][1][0], 1);
            Vmath::Zero(nq, &Dyy[1][2][0], 1);
            Vmath::Vcopy(nq, &MuVel[0][0], 1, &Dyy[1][3][0], 1);
            
            // row 3;
            Vmath::Zero(nq, &Dyy[2][0][0], 1);
            Vmath::Zero(nq, &Dyy[2][1][0], 1);
            Vmath::Vcopy(nq, &FourOThreeMuVector[0], 1, &Dyy[2][2][0], 1);
            Vmath::Vmul(nq,
                        &FourOThreeMuVector[0], 1,
                        &vel[1][0], 1,
                        &Dyy[2][3][0], 1);
            
            // row 4;
            Vmath::Zero(nq, &Dyy[3][0][0], 1);
            Vmath::Zero(nq, &Dyy[3][1][0], 1);
            Vmath::Zero(nq, &Dyy[3][2][0], 1);
            Vmath::Vdiv(nq,
                        &ones[0], 1,
                        &inarray[0][0], 1,
                        &Dyy[3][3][0], 1);
            
            Vmath::Vmul(nq,
                        &GamMu_O_GamMinOnePrVector[0], 1,
                        &Dyy[3][3][0], 1,
                        &Dyy[3][3][0], 1);
        }
        
        if (m_spacedim == 3)
        {
            ASSERTL0(false, "3D adjoint solver is not yet implemented");
        }
    }
    
    
    void CompressibleFlowSystem::
    GetJacobianAddConvFluxFromPrimitiveVar(
                                           Array<OneD, Array<OneD, Array<OneD, Array<OneD, NekDouble> > > > &Jac)
    {
        int i, j;
        int nq = Jac[0][0][0].num_elements();
        int nvar = m_fields.num_elements();
        
        Array<OneD, NekDouble> zeros(nq, 0.0);
        Array<OneD, NekDouble> PoRho2(nq, 0.0);
        Array<OneD, NekDouble> OneORho(nq, 0.0);
        Array<OneD, NekDouble> ones(nq, 1.0);
        Array<OneD, NekDouble> Pressure(nq, 0.0);
        Array<OneD, NekDouble> kPrhoVec(nq, 0.0);
        
        NekDouble FourOThreeMuSclr = 4.0/3.0*m_mu;
        NekDouble TwoOThreeMuSclr  = 2.0/3.0*m_mu;
        NekDouble GamMu_O_GamMinOnePrSclr = m_gamma*m_mu/((m_gamma-1)*m_Prandtl);
        
        Array<OneD, NekDouble> FourOThreeMuVector(nq, FourOThreeMuSclr);
        Array<OneD, NekDouble> TwoOThreeMuVector(nq, TwoOThreeMuSclr);
        Array<OneD, NekDouble> GamMu_O_GamMinOnePrVector(nq,GamMu_O_GamMinOnePrSclr);
        
        Array<OneD, Array<OneD, NekDouble> > dirFlds(nvar);
        
        for (i = 0; i < nvar; ++i)
        {
            dirFlds[i] = Array<OneD, NekDouble>(nq, 0.0);
        }
        
        GetDirectSolution(dirFlds);
        
        
        if (m_ViscosityType == "Constant")
        {
            if (m_spacedim == 2)
            {
                
                Array<OneD, Array<OneD, NekDouble> > vel(m_spacedim);
                
                Array<OneD, Array<OneD, NekDouble> > dui_dx(m_spacedim);
                Array<OneD, Array<OneD, NekDouble> > dui_dy(m_spacedim);
                Array<OneD, Array<OneD, NekDouble> > drho_dxi(m_spacedim);
                Array<OneD, Array<OneD, NekDouble> > dP_dxi(m_spacedim);
                
                Array<OneD, Array<OneD, NekDouble> > dPoRho2_dxi(m_spacedim);
                Array<OneD, Array<OneD, NekDouble> > dPoRho2_dxi2(m_spacedim);
                Array<OneD, Array<OneD, NekDouble> > dOneORho_dxi(m_spacedim);
                Array<OneD, Array<OneD, NekDouble> > dOneORho_dxi2(m_spacedim);
                
                for (i = 0; i < m_spacedim; ++i)
                {
                    vel[i]           = Array<OneD, NekDouble>(nq, 0.0);
                    dui_dx[i]        = Array<OneD, NekDouble>(nq, 0.0);
                    dui_dy[i]        = Array<OneD, NekDouble>(nq, 0.0);
                    drho_dxi[i]      = Array<OneD, NekDouble>(nq, 0.0);
                    dP_dxi[i]        = Array<OneD, NekDouble>(nq, 0.0);
                    dPoRho2_dxi[i]   = Array<OneD, NekDouble>(nq, 0.0);
                    dPoRho2_dxi2[i]   = Array<OneD, NekDouble>(nq, 0.0);
                    dOneORho_dxi[i]   = Array<OneD, NekDouble>(nq, 0.0);
                    dOneORho_dxi2[i]   = Array<OneD, NekDouble>(nq, 0.0);
                }
                
                GetVelocityVector(dirFlds, vel);
                GetPressure      (dirFlds, vel, Pressure);
                
                for (i = 0; i < m_spacedim; ++i)
                {
                    m_primal[0]->PhysDeriv(vel[i], dui_dx[i], dui_dy[i]);
                    
                }
                
                m_fields[0]->PhysDeriv(dirFlds[0], drho_dxi[0], drho_dxi[1]);
                m_fields[0]->PhysDeriv(Pressure, dP_dxi[0], dP_dxi[1]);
                
                // Calculate d(P/rho^2)/dx = (ƒDer∂rho*dP/dx-2*P*drho/dx)/rho^3
                /*Array<OneD, NekDouble> tmp(nq, 0.0);
                 Vmath::Vmul(nq, &dP_dxi[0][0], 1, &dirFlds[0][0], 1, &tmp[0], 1);
                 Vmath::Vmul(nq, &drho_dxi[0][0], 1, &Pressure[0], 1, &dPoRho2_dxi[0][0], 1);
                 Vmath::Smul(nq, 2.0,  &dPoRho2_dxi[0][0], 1,  &dPoRho2_dxi[0][0], 1);
                 Vmath::Vsub(nq, &tmp[0], 1,  &dPoRho2_dxi[0][0], 1,  &dPoRho2_dxi[0][0], 1);
                 Vmath::Zero(nq, &tmp[0], 1);
                 Vmath::Vmul(nq, &dirFlds[0][0], 1, &dirFlds[0][0], 1, &tmp[0], 1);
                 Vmath::Vmul(nq, &tmp[0], 1, &dirFlds[0][0], 1, &tmp[0], 1);
                 Vmath::Vdiv(nq, &dPoRho2_dxi[0][0], 1, &tmp[0], 1, &dPoRho2_dxi[0][0], 1);
                 
                 // Calculate d(P/rho^2)/dy = (rho*dP/dy-2*P*drho/dy)/rho^3
                 Vmath::Zero(nq, &tmp[0], 1);
                 Vmath::Vmul(nq, &dP_dxi[1][0], 1, &dirFlds[0][0], 1, &tmp[0], 1);
                 Vmath::Vmul(nq, &drho_dxi[1][0], 1, &Pressure[0], 1, &dPoRho2_dxi[1][0], 1);
                 Vmath::Smul(nq, 2.0,  &dPoRho2_dxi[1][0], 1,  &dPoRho2_dxi[1][0], 1);
                 Vmath::Vsub(nq, &tmp[0], 1,  &dPoRho2_dxi[1][0], 1,  &dPoRho2_dxi[1][0], 1);
                 Vmath::Zero(nq, &tmp[0], 1);
                 Vmath::Vmul(nq, &dirFlds[0][0], 1, &dirFlds[0][0], 1, &tmp[0], 1);
                 Vmath::Vmul(nq, &tmp[0], 1, &dirFlds[0][0], 1, &tmp[0], 1);
                 Vmath::Vdiv(nq, &dPoRho2_dxi[1][0], 1, &tmp[0], 1, &dPoRho2_dxi[1][0], 1);
                 Vmath::Zero(nq, &tmp[0], 1);
                 
                 // Calculate d(1/rho)/dx
                 Vmath::Vmul(nq, &dirFlds[0][0], 1, &dirFlds[0][0], 1, &tmp[0], 1);
                 Vmath::Vdiv(nq, &drho_dxi[0][0], 1, &tmp[0], 1, &dOneORho_dxi[0][0], 1);
                 Vmath::Neg(nq, &dOneORho_dxi[0][0], 1);
                 Vmath::Zero(nq, &tmp[0], 1);
                 // Calculate d(1/rho)/dy
                 Vmath::Vmul(nq, &dirFlds[0][0], 1, &dirFlds[0][0], 1, &tmp[0], 1);
                 Vmath::Vdiv(nq, &drho_dxi[1][0], 1, &tmp[0], 1, &dOneORho_dxi[1][0], 1);
                 Vmath::Neg(nq, &dOneORho_dxi[1][0], 1);*/
                
                Vmath::Vdiv(nq, &ones[0], 1, &dirFlds[0][0], 1, &OneORho[0], 1);
                m_primal[0]->PhysDeriv(OneORho, dOneORho_dxi[0], dOneORho_dxi[1]);
                
                Vmath::Vmul(nq, &dirFlds[0][0], 1, &dirFlds[0][0], 1, &PoRho2[0], 1);
                Vmath::Vdiv(nq, &Pressure[0], 1, &PoRho2[0], 1, &PoRho2[0], 1);
                m_primal[0]->PhysDeriv(PoRho2, dPoRho2_dxi[0], dPoRho2_dxi[1]);
                /*
                 Array<OneD, NekDouble> diffx(nq, 0.0);
                 Array<OneD, NekDouble> diffy(nq, 0.0);
                 
                 Vmath::Vsub(nq, dOneORho_dxi2[0], 1, dOneORho_dxi[0], 1, diffx, 1);
                 Vmath::Vsub(nq, dOneORho_dxi2[1], 1, dOneORho_dxi[1], 1, diffy, 1);*/
                
                //cout << Vmath::Vmax(nq, diffx, 1) << "  " << Vmath::Vmax(nq, diffy, 1) << endl;
                
                /*
                 for (i = 0; i < nq; ++i)
                 {
                 if (diffx[i] > 1.0)
                 {
                 cout << i << " " << diffx[i] << " " << diffy[i] << " " << dirFlds[0][i] << endl;
                 }
                 if (diffy[i] > 1.0)
                 {
                 cout << i << " " << diffx[i] << " " << diffy[i] << " " << dirFlds[0][i] << endl;
                 }
                 }
                 */
                
                //cout << Vmath::Vmin(nq, &rhoabs[0], 1) << endl;
                
                /*cout << "MAX = " << Vmath::Vmax(nq, diffx, 1)
                 << " " << Vmath::Vmax(nq, diffy, 1) << endl;
                 cout << "MIN = " << Vmath::Vmin(nq, diffx, 1)
                 << " " << Vmath::Vmin(nq, diffy, 1) << endl;
                 cout << "AVG = " << avgx << " " << avgy << endl;*/
                
                /*cout << "Pressure/rho^2 = "
                 << Vmath::Vmax(nq, dPoRho2_dxi[0], 1) << " "
                 << Vmath::Vmax(nq, dPoRho2_dxi[1], 1) << endl;
                 
                 
                 cout << "1/rho = "
                 << Vmath::Vmax(nq, dOneORho_dxi[0], 1) << " "
                 << Vmath::Vmax(nq, dOneORho_dxi[1], 1) <<  endl;*/
                
                /* =============================================================
                 Acvx =
                 
                 [ 0,        0,  0, -gamma/((gamma-1)*Pr)*mu*d(P/rho^2)/dx)]
                 [ 0,        0,  0,                (2*mu*(2*dudx - dvdy))/3]
                 [ 0,        0,  0,                        mu*(dudy + dvdx)]
                 [ 0,        0,  0,    gamma/((gamma-1)*Pr)*mu*d(1/rho)/dx)]
                 
                 =============================================================*/
                
                Vmath::Zero(nq, &Jac[0][0][0][0], 1);
                Vmath::Zero(nq, &Jac[0][0][1][0], 1);
                Vmath::Zero(nq, &Jac[0][0][2][0], 1);
                
                Array<OneD, NekDouble> tmp(nq, 0.0);
                Array<OneD, NekDouble> tmp1(nq, 0.0);
                
                //(k*(dpdx - drhodx*p))/rho^2
                /*Vmath::Vmul(nq,
                 &inarray[0][0], 1,
                 &inarray[0][0], 1,
                 &tmp[0], 1);
                 Vmath::Vmul(nq,
                 &Pressure[0], 1,
                 &drho_dxi[0][0], 1,
                 &tmp1[0], 1);
                 
                 Vmath::Vsub(nq,
                 &dP_dxi[0][0], 1,
                 &tmp1[0], 1,
                 &tmp1[0], 1);
                 
                 Vmath::Vdiv(nq, &tmp1[0], 1, &tmp[0], 1, &Avx[0][3][0], 1);
                 Vmath::Smul(nq,
                 m_thermalConductivity,
                 &Avx[0][3][0], 1,
                 &Avx[0][3][0], 1);*/
                
                //-gamma/((gamma-1)*Pr)*mu*d(P/rho^2)/dx)
                Vmath::Vmul(nq,
                            &GamMu_O_GamMinOnePrVector[0], 1,
                            &dPoRho2_dxi[0][0], 1,
                            &Jac[0][0][3][0], 1);
                
                Vmath::Zero(nq, &Jac[0][1][0][0], 1);
                Vmath::Zero(nq, &Jac[0][1][1][0], 1);
                Vmath::Zero(nq, &Jac[0][1][2][0], 1);
                
                //(2*mu*(2*dudx - dvdy))/3
                Vmath::Zero(nq, &tmp[0], 1);
                Vmath::Zero(nq, &tmp1[0], 1);
                
                Vmath::Smul(nq, 2.0, &dui_dx[0][0], 1, &tmp[0], 1);
                Vmath::Vsub(nq, &tmp[0], 1, &dui_dy[1][0], 1, &tmp[0], 1);
                Vmath::Smul(nq,
                            TwoOThreeMuSclr,
                            &tmp[0], 1,
                            &Jac[0][1][3][0], 1);
                
                Vmath::Neg(nq, &Jac[0][1][3][0], 1);
                
                Vmath::Zero(nq, &Jac[0][2][0][0], 1);
                Vmath::Zero(nq, &Jac[0][2][1][0], 1);
                Vmath::Zero(nq, &Jac[0][2][2][0], 1);
                
                //mu*(dudy + dvdx)
                Vmath::Zero(nq, &tmp[0], 1);
                Vmath::Zero(nq, &tmp1[0], 1);
                
                Vmath::Vadd(nq,
                            &dui_dy[0][0], 1,
                            &dui_dx[1][0], 1,
                            &Jac[0][2][3][0], 1);
                Vmath::Smul(nq, m_mu, &Jac[0][2][3][0], 1, &Jac[0][2][3][0], 1);
                
                Vmath::Neg(nq, &Jac[0][2][3][0], 1);
                
                Vmath::Zero(nq, &Jac[0][3][0][0], 1);
                Vmath::Zero(nq, &Jac[0][3][1][0], 1);
                Vmath::Zero(nq, &Jac[0][3][2][0], 1);
                
                // drho/dx*k/rho
                /*Vmath::Zero(nq, &tmp[0], 1);
                 Vmath::Zero(nq, &tmp1[0], 1);
                 Vmath::Vcopy(nq, &drho_dxi[0][0], 1, &Avx[3][3][0], 1);
                 Vmath::Vdiv(nq,
                 &Avx[3][3][0], 1,
                 &inarray[0][0], 1,
                 &Avx[3][3][0], 1);
                 Vmath::Smul(nq,
                 m_thermalConductivity,
                 &Avx[3][3][0], 1,
                 &Avx[3][3][0], 1);*/
                
                // gamma/((gamma-1)*Pr)*mu*d(1/rho)/dx)
                
                Vmath::Vmul(nq,
                            &GamMu_O_GamMinOnePrVector[0], 1,
                            &dOneORho_dxi[0][0], 1,
                            &Jac[0][3][3][0], 1);
                
                Vmath::Neg(nq, &Jac[0][3][3][0], 1);
                
                /* =============================================================
                 Acvy =
                 
                 [ 0,        0,  0,         (k*(dpdy - drhody*p))/rho^2]
                 [ 0,        0,  0,                    mu*(dudy + dvdx)]
                 [ 0,        0,  0,           -(2*mu*(dudx - 2*dvdy))/3]
                 [ 0,        0,  0,                      (drhody*k)/rho]
                 
                 =============================================================*/
                
                Vmath::Zero(nq, &Jac[1][0][0][0], 1);
                Vmath::Zero(nq, &Jac[1][0][1][0], 1);
                Vmath::Zero(nq, &Jac[1][0][2][0], 1);
                
                Vmath::Zero(nq, &tmp[0], 1);
                Vmath::Zero(nq, &tmp1[0], 1);
                
                /*Vmath::Vmul(nq,
                 &inarray[0][0], 1,
                 &inarray[0][0], 1,
                 &tmp[0], 1);
                 
                 Vmath::Vmul(nq,
                 &Pressure[0], 1,
                 &drho_dxi[0][0], 1,
                 &tmp1[0], 1);
                 
                 Vmath::Vsub(nq,
                 &dP_dxi[0][0], 1,
                 &tmp1[0], 1,
                 &tmp1[0], 1);
                 
                 Vmath::Vdiv(nq, &tmp1[0], 1, &tmp[0], 1, &Avy[0][3][0], 1);
                 
                 Vmath::Smul(nq,
                 m_thermalConductivity,
                 &Avy[0][3][0], 1,
                 &Avy[0][3][0], 1);*/
                
                Vmath::Vmul(nq,
                            &GamMu_O_GamMinOnePrVector[0], 1,
                            &dPoRho2_dxi[1][0], 1,
                            &Jac[1][0][3][0], 1);
                
                Vmath::Zero(nq, &Jac[1][1][0][0], 1);
                Vmath::Zero(nq, &Jac[1][1][1][0], 1);
                Vmath::Zero(nq, &Jac[1][1][2][0], 1);
                
                //mu*(dudy + dvdx)
                Vmath::Zero(nq, &tmp[0], 1);
                Vmath::Zero(nq, &tmp1[0], 1);
                
                Vmath::Vadd(nq,
                            &dui_dy[0][0], 1,
                            &dui_dx[1][0], 1,
                            &Jac[1][1][3][0], 1);
                
                Vmath::Smul(nq, m_mu, &Jac[1][1][3][0], 1, &Jac[1][1][3][0], 1);
                
                Vmath::Neg(nq, &Jac[1][1][3][0], 1);
                
                Vmath::Zero(nq, &Jac[1][2][0][0], 1);
                Vmath::Zero(nq, &Jac[1][2][1][0], 1);
                Vmath::Zero(nq, &Jac[1][2][2][0], 1);
                
                //(2*mu*(dudx - 2*dvdy))/3
                Vmath::Zero(nq, &tmp[0], 1);
                Vmath::Zero(nq, &tmp1[0], 1);
                
                Vmath::Smul(nq, 2.0, &dui_dy[1][0], 1, &tmp[0], 1);
                Vmath::Vsub(nq, &dui_dx[0][0], 1, &tmp[0], 1, &tmp[0], 1);
                Vmath::Vmul(nq,
                            &TwoOThreeMuVector[0], 1,
                            &tmp[0], 1,
                            &Jac[1][2][3][0], 1);
                
                Vmath::Neg(nq, &Jac[1][2][3][0], 1);
                
                Vmath::Zero(nq, &Jac[1][3][0][0], 1);
                Vmath::Zero(nq, &Jac[1][3][1][0], 1);
                Vmath::Zero(nq, &Jac[1][3][2][0], 1);
                
                // drho/dy*k/rho
                /*Vmath::Zero(nq, &tmp[0], 1);
                 Vmath::Zero(nq, &tmp1[0], 1);
                 Vmath::Vcopy(nq, &drho_dxi[1][0], 1, &Avy[3][3][0], 1);
                 Vmath::Vdiv(nq,
                 &Avy[3][3][0], 1,
                 &inarray[0][0], 1,
                 &Avy[3][3][0], 1);
                 Vmath::Smul(nq,
                 m_thermalConductivity,
                 &Avy[3][3][0], 1,
                 &Avy[3][3][0], 1);*/
                
                // gamma/((gamma-1)*Pr)*mu*d(1/rho)/dy)
                
                Vmath::Vmul(nq,
                            &GamMu_O_GamMinOnePrVector[0], 1,
                            &dOneORho_dxi[1][0], 1,
                            &Jac[1][3][3][0], 1);
                
                Vmath::Neg(nq, &Jac[1][3][3][0], 1);
            }
        }
    }
    
    void CompressibleFlowSystem::GetAdjointAddConvFluxVectorFromPrimitiveVar(
                                                                             const Array<OneD, Array<OneD, NekDouble> > &inarray,
                                                                             Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &outarray)
    {
        int i,j,k;
        
        int nq = inarray[0].num_elements();
        int nvar = m_fields.num_elements();
        
        Array<OneD, NekDouble> tmp(nq, 0.0);
        Array<OneD, NekDouble> tmp1(nq, 0.0);
        Array<OneD, Array<OneD, NekDouble> > tmpx(nvar);
        Array<OneD, Array<OneD, NekDouble> > tmpy(nvar);
        
        
        Array<OneD, Array<OneD, Array<OneD,
        Array<OneD, NekDouble > > > > Jac(m_spacedim);
        Array<OneD, Array<OneD, Array<OneD,
        Array<OneD, NekDouble > > > >JacCons(m_spacedim);
        
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > > dVdU(nvar);
        
        // span the required arrays
        for (i = 0; i < m_spacedim; ++i)
        {
            Jac[i]     = Array<OneD, Array<OneD,
            Array<OneD, NekDouble > > >(nvar);
            JacCons[i] = Array<OneD, Array<OneD,
            Array<OneD, NekDouble > > >(nvar);
            
            for (j = 0; j < nvar; ++j)
            {
                Jac[i][j] = Array<OneD, Array<OneD, NekDouble> >(nvar);
                JacCons[i][j] = Array<OneD, Array<OneD, NekDouble> >(nvar);
                
                for (k = 0; k < nvar; ++k)
                {
                    Jac[i][j][k] = Array<OneD, NekDouble>(nq, 0.0);
                    JacCons[i][j][k] = Array<OneD, NekDouble>(nq, 0.0);
                }
            }
        }
        
        for (i = 0; i < nvar; ++i)
        {
            dVdU[i] = Array<OneD, Array<OneD, NekDouble > >(nvar);
            
            for (j = 0; j < nvar; ++j)
            {
                dVdU[i][j] = Array<OneD, NekDouble>(nq, 0.0);
            }
        }
        
        GetJacobianAddConvFluxFromPrimitiveVar(Jac);
        
        GetConservToPrimVariableInvMat(dVdU);
        
        // Perform tensor product dVdU*dFdV
        for (i = 0; i < nvar; ++i)
        {
            for (j = 0; j < nvar; ++j)
            {
                Vmath::Zero(nq, &tmp[0], 1);
                Vmath::Zero(nq, &tmp1[0], 1);
                
                for (k = 0; k < nvar; ++k)
                {
                    Vmath::Vvtvp(nq,
                                 &dVdU[i][k][0], 1,
                                 &Jac[0][k][j][0], 1,
                                 &tmp[0], 1,
                                 &tmp[0], 1);
                    
                    Vmath::Vcopy(nq,
                                 &tmp[0], 1,
                                 &JacCons[0][i][j][0], 1);
                    
                    Vmath::Vvtvp(nq,
                                 &dVdU[i][k][0], 1,
                                 &Jac[1][k][j][0], 1,
                                 &tmp1[0], 1,
                                 &tmp1[0], 1);
                    
                    Vmath::Vcopy(nq,
                                 &tmp1[0], 1,
                                 &JacCons[1][i][j][0], 1);
                }
            }
        }
        
        for (i = 0; i < nvar; ++i)
        {
            tmpx[i] = Array<OneD, NekDouble>(nq, 0.0);
            tmpy[i] = Array<OneD, NekDouble>(nq, 0.0);
            
            for (j = 0; j < nvar; ++j)
            {
                if (m_spacedim == 2)
                {
                    Vmath::Vvtvp(nq,
                                 &JacCons[0][i][j][0], 1,
                                 &inarray[j][0], 1,
                                 &tmpx[i][0], 1,
                                 &tmpx[i][0], 1);
                    
                    Vmath::Vvtvp(nq,
                                 &JacCons[1][i][j][0], 1,
                                 &inarray[j][0], 1,
                                 &tmpy[i][0], 1,
                                 &tmpy[i][0], 1);
                    
                }
                
                if (m_spacedim == 3)
                {
                    ASSERTL0(false, "3D adjoint solver is ot implemented yet")
                }
                
            }
            
            Vmath::Vcopy(nq, &tmpx[i][0], 1, &outarray[i][0][0], 1);
            Vmath::Vcopy(nq, &tmpy[i][0], 1, &outarray[i][1][0], 1);
        }
    }
    
    
    void CompressibleFlowSystem::GetDerivJacobian(
                                                  Array<OneD, Array<OneD, Array<OneD, Array<OneD, NekDouble> > > > &Jac,
                                                  Array<OneD, Array<OneD, Array<OneD, Array<OneD, NekDouble> > > > &JacDiv)
    {
        int i,j,k;
        
        int nq = JacDiv[0][0][0].num_elements();
        int nvar = m_fields.num_elements();
        
        
        if (m_spacedim == 2)
        {
            Array<OneD, NekDouble> tmp(nq,0.0);
            
            for (i = 0; i < nvar; ++i)
            {
                for (j = 0; j < nvar; ++j)
                {
                    m_primal[0]->PhysDeriv(Jac[0][i][j],
                                           JacDiv[0][i][j],
                                           tmp);
                    
                    m_primal[0]->PhysDeriv(Jac[1][i][j],
                                           tmp,
                                           JacDiv[1][i][j]);
                }
            }
        }
        if (m_spacedim == 3)
        {
            
            Array<OneD, NekDouble> tmp1(nq,0.0);
            Array<OneD, NekDouble> tmp2(nq,0.0);
            
            for (i = 0; i < nvar; ++i)
            {
                for (j = 0; j < nvar; ++j)
                {
                    m_primal[0]->PhysDeriv(Jac[0][i][j],
                                           JacDiv[0][i][j],
                                           tmp1,
                                           tmp2);
                    
                    m_primal[0]->PhysDeriv(Jac[1][i][j],
                                           tmp1,
                                           JacDiv[1][i][j],
                                           tmp2);
                    
                    m_primal[0]->PhysDeriv(Jac[2][i][j],
                                           tmp1,
                                           tmp2,
                                           JacDiv[2][i][j]);
                }
            }
        }
    }
    
    void CompressibleFlowSystem::GetDerivAddJacVectorFromPrimitiveVar(
                                                                      const Array<OneD, Array<OneD, NekDouble> > &inarray,
                                                                      Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &outarray)
    {
        int i,j,k;
        
        int nq = inarray[0].num_elements();
        int nvar = m_fields.num_elements();
        
        Array<OneD, NekDouble> tmp(nq, 0.0);
        Array<OneD, NekDouble> tmp1(nq, 0.0);
        Array<OneD, NekDouble> tmp2(nq, 0.0);
        Array<OneD, NekDouble> tmp3(nq, 0.0);
        Array<OneD, NekDouble> tmp4(nq, 0.0);
        Array<OneD, NekDouble> tmp5(nq, 0.0);
        
        Array<OneD, Array<OneD, NekDouble> > tmpx(nvar);
        Array<OneD, Array<OneD, NekDouble> > tmpy(nvar);
        
        Array<OneD, Array<OneD, Array<OneD,
        Array<OneD, NekDouble > > > > Jac(m_spacedim);
        Array<OneD, Array<OneD, Array<OneD,
        Array<OneD, NekDouble > > > > JacCons(m_spacedim);
        Array<OneD, Array<OneD, Array<OneD,
        Array<OneD, NekDouble > > > > JacDiv(m_spacedim);
        Array<OneD, Array<OneD, Array<OneD,
        Array<OneD, NekDouble > > > >JacConsD(m_spacedim);
        Array<OneD, Array<OneD, Array<OneD,
        Array<OneD, NekDouble > > > >JacConsDiv(m_spacedim);
        
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > > dVdU(nvar);
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > > dVdUdX(nvar);
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > > dVdUdY(nvar);
        
        // span the required arrays
        for (i = 0; i < m_spacedim; ++i)
        {
            Jac[i]     = Array<OneD, Array<OneD,
            Array<OneD, NekDouble > > >(nvar);
            JacDiv[i]     = Array<OneD, Array<OneD,
            Array<OneD, NekDouble > > >(nvar);
            JacConsDiv[i] = Array<OneD, Array<OneD,
            Array<OneD, NekDouble > > >(nvar);
            
            for (j = 0; j < nvar; ++j)
            {
                Jac[i][j] = Array<OneD, Array<OneD, NekDouble> >(nvar);
                JacDiv[i][j] = Array<OneD, Array<OneD, NekDouble> >(nvar);
                JacConsDiv[i][j] = Array<OneD, Array<OneD, NekDouble> >(nvar);
                
                for (k = 0; k < nvar; ++k)
                {
                    Jac[i][j][k] = Array<OneD, NekDouble>(nq, 0.0);
                    JacDiv[i][j][k] = Array<OneD, NekDouble>(nq, 0.0);
                    JacConsDiv[i][j][k] = Array<OneD, NekDouble>(nq, 0.0);
                }
            }
        }
        
        for (i = 0; i < nvar; ++i)
        {
            dVdU[i] = Array<OneD, Array<OneD, NekDouble > >(nvar);
            dVdUdX[i] = Array<OneD, Array<OneD, NekDouble > >(nvar);
            dVdUdY[i] = Array<OneD, Array<OneD, NekDouble > >(nvar);
            
            for (j = 0; j < nvar; ++j)
            {
                dVdU[i][j] = Array<OneD, NekDouble>(nq, 0.0);
                dVdUdX[i][j] = Array<OneD, NekDouble>(nq, 0.0);
                dVdUdY[i][j] = Array<OneD, NekDouble>(nq, 0.0);
            }
        }
        
        GetJacobianAddConvFluxFromPrimitiveVar(Jac);
        
        GetDerivJacobian(Jac, JacDiv);
        
        GetConservToPrimVariableInvMat(dVdU);
        
        GetConservToPrimVariableInvMatDiv(dVdU, dVdUdX, dVdUdY);
        
        for (i = 0; i < nvar; ++i)
        {
            for (j = 0; j < nvar; ++j)
            {
                Vmath::Zero(nq, &tmp[0], 1);
                Vmath::Zero(nq, &tmp1[0], 1);
                Vmath::Zero(nq, &tmp2[0], 1);
                Vmath::Zero(nq, &tmp3[0], 1);
                Vmath::Zero(nq, &tmp4[0], 1);
                Vmath::Zero(nq, &tmp5[0], 1);
                
                for (k = 0; k < nvar; ++k)
                {
                    //
                    Vmath::Vvtvp(nq,
                                 &dVdUdX[i][k][0], 1,
                                 &Jac[0][k][j][0], 1,
                                 &tmp[0], 1,
                                 &tmp[0], 1);
                    
                    Vmath::Vvtvp(nq,
                                 &dVdUdY[i][k][0], 1,
                                 &Jac[1][k][j][0], 1,
                                 &tmp1[0], 1,
                                 &tmp1[0], 1);
                    
                    //=============
                    
                    Vmath::Vvtvp(nq,
                                 &dVdU[i][k][0], 1,
                                 &JacDiv[0][k][j][0], 1,
                                 &tmp2[0], 1,
                                 &tmp2[0], 1);
                    
                    Vmath::Vvtvp(nq,
                                 &dVdU[i][k][0], 1,
                                 &JacDiv[1][k][j][0], 1,
                                 &tmp3[0], 1,
                                 &tmp3[0], 1);
                    
                    Vmath::Vadd(nq,
                                &tmp[0], 1,
                                &tmp2[0], 1,
                                &JacConsDiv[0][i][j][0], 1);
                    
                    Vmath::Vadd(nq,
                                &tmp1[0], 1,
                                &tmp3[0], 1,
                                &JacConsDiv[1][i][j][0], 1);
                }
            }
        }
        
        for (i = 0; i < nvar; ++i)
        {
            tmpx[i] = Array<OneD, NekDouble>(nq, 0.0);
            tmpy[i] = Array<OneD, NekDouble>(nq, 0.0);
            
            for (j = 0; j < nvar; ++j)
            {
                Vmath::Vvtvp(nq,
                             &JacConsDiv[0][i][j][0], 1,
                             &inarray[j][0], 1,
                             &tmpx[i][0], 1,
                             &tmpx[i][0], 1);
                
                Vmath::Vvtvp(nq,
                             &JacConsDiv[1][i][j][0], 1,
                             &inarray[j][0], 1,
                             &tmpy[i][0], 1,
                             &tmpy[i][0], 1);
                
            }
            
            Vmath::Vcopy(nq, &tmpx[i][0], 1, &outarray[i][0][0], 1);
            Vmath::Vcopy(nq, &tmpy[i][0], 1, &outarray[i][1][0], 1);
        }
    }

    void CompressibleFlowSystem::GetSensor(
        const Array<OneD, const Array<OneD, NekDouble> > &physarray,
              Array<OneD,                   NekDouble>   &Sensor,
              Array<OneD,                   NekDouble>   &SensorKappa)
    {
        int e, NumModesElement, nQuadPointsElement;
        int nTotQuadPoints  = GetTotPoints();
        int nElements       = m_fields[0]->GetExpSize();
        
        // Find solution (SolP) at p = P;
        // The input array (physarray) is the solution at p = P;
        
        Array<OneD,int> ExpOrderElement = GetNumExpModesPerExp();
        
        Array<OneD, NekDouble> SolP(nTotQuadPoints,0.0);
        Array<OneD, NekDouble> SolPmOne(nTotQuadPoints,0.0);
        Array<OneD, NekDouble> SolNorm(nTotQuadPoints,0.0);
        
        Vmath::Vcopy(nTotQuadPoints,physarray[0],1,SolP,1);
        
        int CoeffsCount = 0;
        
        for (e = 0; e < nElements; e++)
        {
            NumModesElement         = ExpOrderElement[e];
            
            int nQuadPointsElement = m_fields[0]->GetExp(e)->GetTotPoints();
            int nCoeffsElement = m_fields[0]->GetExp(e)->GetNcoeffs();
            int numCutOff = NumModesElement - 1;
            
            // Set-up of the Orthogonal basis for a Quadrilateral element which is
            // needed to obtain thesolution at P =  p - 1;
            
            Array<OneD, NekDouble> SolPElementPhys(nQuadPointsElement,0.0);
            Array<OneD, NekDouble> SolPElementCoeffs(nCoeffsElement,0.0);
            
            Array<OneD, NekDouble> SolPmOneElementPhys(nQuadPointsElement,0.0);
            Array<OneD, NekDouble> SolPmOneElementCoeffs(nCoeffsElement,0.0);
            
            // create vector the save the solution points per element at P = p;
            
            for (int i = 0; i < nQuadPointsElement; i++)
            {
                SolPElementPhys[i] = SolP[CoeffsCount+i];
            }

            m_fields[0]->GetExp(e)->FwdTrans(SolPElementPhys,
                                             SolPElementCoeffs);
            
            // ReduceOrderCoeffs reduces the polynomial order of the solution that
            // is represented by the coeffs given as an inarray. This is done by
            // projecting the higher order solution onto the orthogonal basis and
            // padding the higher order coefficients with zeros.
            
            m_fields[0]->GetExp(e)->ReduceOrderCoeffs(numCutOff,
                                                      SolPElementCoeffs,
                                                      SolPmOneElementCoeffs);
            
            m_fields[0]->GetExp(e)->BwdTrans(SolPmOneElementCoeffs,
                                             SolPmOneElementPhys);
            
            for (int i = 0; i < nQuadPointsElement; i++)
            {
                SolPmOne[CoeffsCount+i] = SolPmOneElementPhys[i];
            }
            
            NekDouble SolPmeanNumerator = 0.0;
            NekDouble SolPmeanDenumerator = 0.0;
            
            // Determining the norm of the numerator of the Sensor
            
            Vmath::Vsub(nQuadPointsElement,
                        SolPElementPhys, 1,
                        SolPmOneElementPhys, 1,
                        SolNorm, 1);
            
            Vmath::Vmul(nQuadPointsElement,
                        SolNorm, 1,
                        SolNorm, 1,
                        SolNorm, 1);
            
            for (int i = 0; i < nQuadPointsElement; i++)
            {
                SolPmeanNumerator   += SolNorm[i];
                SolPmeanDenumerator += SolPElementPhys[i];
            }
            
            for (int i = 0; i < nQuadPointsElement; ++i)
            {
                Sensor[CoeffsCount+i] = sqrt(SolPmeanNumerator/nQuadPointsElement);
                
                Sensor[CoeffsCount+i] = log10(Sensor[CoeffsCount+i]);
            }
            CoeffsCount += nQuadPointsElement;
        }
    
        CoeffsCount = 0.0;
        
        for (e = 0; e < nElements; e++)
        {
            NumModesElement         = ExpOrderElement[e];
            NekDouble ThetaS        = m_mu0;
            NekDouble Phi0          = m_Skappa;
            NekDouble DeltaPhi      = m_Kappa;
            nQuadPointsElement      = m_fields[0]->GetExp(e)->GetTotPoints();
            
            for (int i = 0; i < nQuadPointsElement; i++)
            {
                if (Sensor[CoeffsCount+i] <= (Phi0 - DeltaPhi))
                {
                    SensorKappa[CoeffsCount+i] = 0;
                }
                else if(Sensor[CoeffsCount+i] >= (Phi0 + DeltaPhi))
                {
                    SensorKappa[CoeffsCount+i] = ThetaS;
                }
                else if(abs(Sensor[CoeffsCount+i]-Phi0) < DeltaPhi)
                {
                    SensorKappa[CoeffsCount+i] = ThetaS/2*(1+sin(M_PI*
                                        (Sensor[CoeffsCount+i]-Phi0)
                                                    /(2*DeltaPhi)));
                }
            }
            
            CoeffsCount += nQuadPointsElement;
        }
    }
    
    void CompressibleFlowSystem::GetForcingTerm(
                const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                      Array<OneD, Array<OneD, NekDouble> > outarrayForcing)
    {
        const int nPts = m_fields[0]->GetTotPoints();
        const int nvariables = m_fields.num_elements();
        const int nElements = m_fields[0]->GetExpSize();
        
        Array<OneD,  NekDouble>  Sensor(nPts, 0.0);
        Array<OneD,  NekDouble>  SensorKappa(nPts, 0.0);
        Array <OneD, NekDouble > Lambda(nPts, 0.0);
        Array <OneD, NekDouble > Tau(nPts, 1.0);
        Array <OneD, NekDouble > soundspeed(nPts, 0.0);
        Array <OneD, NekDouble > pressure(nPts, 0.0);
        Array <OneD, NekDouble > temperature(nPts, 0.0);
        Array <OneD, NekDouble > absVelocity(nPts, 0.0);
        Array <OneD, NekDouble > hel(nPts, 0.0);
        Array <OneD, NekDouble > h_minmin(m_spacedim, 0.0);
        
        Array<OneD,int> pOrderElmt = GetNumExpModesPerExp();
        Array<OneD, NekDouble> pOrder (nPts, 0.0);
        // Thermodynamic related quantities
        GetPressure(inarray, pressure);
        GetTemperature(inarray, pressure, temperature);
        GetSoundSpeed(inarray, pressure, soundspeed);
        GetAbsoluteVelocity(inarray, absVelocity);
        GetSensor(inarray, Sensor, SensorKappa);
        
        // Determine the maximum wavespeed
        Vmath::Vadd(nPts, absVelocity, 1, soundspeed, 1, Lambda, 1);
        
        NekDouble LambdaMax = Vmath::Vmax(nPts, Lambda, 1);
        
        int PointCount = 0;
        
        for (int e = 0; e < nElements; e++)
        {
            int nQuadPointsElement = m_fields[0]->GetExp(e)->GetTotPoints();
            
            for (int n = 0; n < nQuadPointsElement; n++)
            {
                pOrder[n + PointCount] = pOrderElmt[e];
                
                Tau[n + PointCount] = 1.0/(m_C1*pOrder[n + PointCount]*LambdaMax); // order 1.0e-06
                
                outarrayForcing[nvariables-1][n + PointCount] = 1/Tau[n + PointCount]*(m_hFactor*LambdaMax/pOrder[n + PointCount]*SensorKappa[n + PointCount]-inarray[nvariables-1][n + PointCount]);
            }
            PointCount += nQuadPointsElement;
        }
    }

    void CompressibleFlowSystem::GetElementDimensions(
                          Array<OneD,       Array<OneD, NekDouble> > &outarray,
                          Array<OneD,       NekDouble > &hmin)
    {
        // So far, this function is only implemented for quads
        const int nElements = m_fields[0]->GetExpSize();
        
        SpatialDomains::HexGeomSharedPtr    ElHexGeom;
        
        NekDouble hx = 0.0;
        NekDouble hy = 0.0;
        
        for (int e = 0; e < nElements; e++)
        {
            NekDouble nedges = m_fields[0]->GetExp(e)->GetNedges();
            Array <OneD, NekDouble> L1(nedges, 0.0);
            
            for (int j = 0; j < nedges; ++j)
            {
                
                NekDouble x0 = 0.0;
                NekDouble y0 = 0.0;
                NekDouble z0 = 0.0;
                
                NekDouble x1 = 0.0;
                NekDouble y1 = 0.0;
                NekDouble z1 = 0.0;
                
                if (boost::dynamic_pointer_cast<SpatialDomains::QuadGeom>(m_fields[0]->GetExp(e)->GetGeom()))
                {
                    SpatialDomains::QuadGeomSharedPtr ElQuadGeom = boost::dynamic_pointer_cast<SpatialDomains::QuadGeom>(m_fields[0]->GetExp(e)->GetGeom());
                    
                    ElQuadGeom->GetEdge(j)->GetVertex(0)->GetCoords(x0,y0,z0);
                    ElQuadGeom->GetEdge(j)->GetVertex(1)->GetCoords(x1,y1,z1);
                    
                    L1[j] = sqrt(pow((x0-x1),2)+pow((y0-y1),2)+pow((z0-z1),2));
                }
                else
                {
                    ASSERTL0(false, "GetElementDimensions() is only implemented for quadrilateral elements");
                }
            }
            // determine the minimum length in x and y direction
            // still have to find a better estimate when dealing with unstructured meshes
            if(boost::dynamic_pointer_cast<SpatialDomains::QuadGeom>(m_fields[0]->GetExp(e)->GetGeom()))
            {
                hx = min(L1[0], L1[2]);
                hy = min(L1[1], L1[3]);
                    
                outarray[0][e] = hx;
                outarray[1][e] = hy;
            }
        }
    
        hmin[0] = Vmath::Vmin(outarray[0].num_elements(), outarray[0], 1);
        hmin[1] = Vmath::Vmin(outarray[1].num_elements(), outarray[1], 1);
    }
    
    void CompressibleFlowSystem::GetAbsoluteVelocity(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
              Array<OneD,                   NekDouble>   &Vtot)
    {
        int nTotQuadPoints = GetTotPoints();
        
        // Getting the velocity vector on the 2D normal space
        Array<OneD, Array<OneD, NekDouble> > velocity   (m_spacedim);
        
        Vmath::Zero(Vtot.num_elements(), Vtot, 1);
        
        for (int i = 0; i < m_spacedim; ++i)
        {
            velocity[i] = Array<OneD, NekDouble>(nTotQuadPoints);
        }
        
        GetVelocityVector(inarray, velocity);
        
        for (int i = 0; i < m_spacedim; ++i)
        {
            Vmath::Vvtvp(nTotQuadPoints,
                         velocity[i], 1,
                         velocity[i], 1,
                         Vtot, 1,
                         Vtot, 1);
        }
        
        Vmath::Vsqrt(Vtot.num_elements(),Vtot,1,Vtot,1);
    }
    
    void CompressibleFlowSystem::GetSmoothArtificialViscosity(
                        const Array<OneD, Array<OneD, NekDouble> > &physfield,
                              Array<OneD,             NekDouble  > &eps_bar)
    {
        int nvariables = physfield.num_elements();
        int nPts       = m_fields[0]->GetTotPoints();
        
        Array<OneD, NekDouble > pressure            (nPts, 0.0);
        Array<OneD, NekDouble > temperature         (nPts, 0.0);
        Array <OneD, NekDouble > sensor             (nPts, 0.0);
        Array <OneD, NekDouble > SensorKappa        (nPts, 0.0);
        Array <OneD, NekDouble > absVelocity        (nPts, 0.0);
        Array <OneD, NekDouble > soundspeed         (nPts, 0.0);
        Array <OneD, NekDouble > Lambda             (nPts, 0.0);
        Array <OneD, NekDouble > mu_var             (nPts, 0.0);
        Array <OneD, NekDouble > h_minmin           (m_spacedim, 0.0);
        Vmath::Zero(nPts, eps_bar, 1);
        
        // Thermodynamic related quantities
        GetPressure(physfield, pressure);
        GetTemperature(physfield, pressure, temperature);
        GetSoundSpeed(physfield, pressure, soundspeed);
        GetAbsoluteVelocity(physfield, absVelocity);
        GetSensor(physfield, sensor, SensorKappa);
        
        // Determine the maximum wavespeed
        Vmath::Vadd(nPts, absVelocity, 1, soundspeed, 1, Lambda, 1);
        
        // Determine hbar = hx_i/h
        Array<OneD,int> pOrderElmt = GetNumExpModesPerExp();
        
        NekDouble ThetaH = m_FacH;
        NekDouble ThetaL = m_FacL;
        
        NekDouble Phi0     = (ThetaH+ThetaL)/2;
        NekDouble DeltaPhi = ThetaH-Phi0;
    
        Vmath::Zero(eps_bar.num_elements(), eps_bar, 1);

	    /*Vmath::Smul(eps_bar.num_elements(),
                    m_eps_max,
                    &physfield[nvariables-1][0], 1,
                    &eps_bar[0], 1);*/

        for (int e = 0; e < eps_bar.num_elements(); e++)
        {
            if (physfield[nvariables-1][e] <= (Phi0 - DeltaPhi))
            {
                eps_bar[e] = 0;
            }
            else if(physfield[nvariables-1][e] >= (Phi0 + DeltaPhi))
            {
                eps_bar[e] = m_mu0;
            }
            else if(abs(physfield[nvariables-1][e]-Phi0) < DeltaPhi)
            {
                eps_bar[e] = m_mu0/2*(1+sin(M_PI*
                (physfield[nvariables-1][e]-Phi0)/(2*DeltaPhi)));
            }
        }	
	
    }
    
    void CompressibleFlowSystem::GetArtificialDynamicViscosity(
        const Array<OneD, Array<OneD, NekDouble> > &physfield,
              Array<OneD,             NekDouble  > &mu_var)
    {
        const int nElements  = m_fields[0]->GetExpSize();
    
        int PointCount = 0;
        int nTotQuadPoints  = GetTotPoints();
        
        Array<OneD, NekDouble> S_e        (nTotQuadPoints, 0.0);
        Array<OneD, NekDouble> se         (nTotQuadPoints, 0.0);
        Array<OneD, NekDouble> Sensor     (nTotQuadPoints, 0.0);
        Array<OneD, NekDouble> SensorKappa(nTotQuadPoints, 0.0);
        Array<OneD, NekDouble> absVelocity(nTotQuadPoints, 0.0);
        Array<OneD, NekDouble> soundspeed (nTotQuadPoints, 0.0);
        Array<OneD, NekDouble> pressure   (nTotQuadPoints, 0.0);
        
        GetAbsoluteVelocity(physfield, absVelocity);
        GetPressure        (physfield, pressure);
        GetSoundSpeed      (physfield, pressure, soundspeed);
        GetSensor          (physfield, Sensor, SensorKappa);
        
        Array<OneD, int> pOrderElmt = GetNumExpModesPerExp();
        Array<OneD, NekDouble> Lambda(nTotQuadPoints, 1.0);
        Vmath::Vadd(nTotQuadPoints, absVelocity, 1, soundspeed, 1, Lambda, 1);
        
        for (int e = 0; e < nElements; e++)
        {
            // Threshold value specified in C. Biottos thesis.  Based on a 1D
            // shock tube problem S_k = log10(1/p^4). See G.E. Barter and
            // D.L. Darmofal. Shock Capturing with PDE-based artificial
            // diffusion for DGFEM: Part 1 Formulation, Journal of Computational
            // Physics 229 (2010) 1810-1827 for further reference

            // Adjustable depending on the coarsness of the mesh. Might want to
            // move this variable into the session file

            int nQuadPointsElement = m_fields[0]->GetExp(e)->GetTotPoints();
            Array <OneD, NekDouble> one2D(nQuadPointsElement, 1.0);
            
            for (int n = 0; n < nQuadPointsElement; n++)
            {
                NekDouble mu_0 = m_mu0;
                
                if (Sensor[n+PointCount] < (m_Skappa-m_Kappa))
                {
                    mu_var[n+PointCount] = 0;
                }
                else if(Sensor[n+PointCount] >= (m_Skappa-m_Kappa)
                        && Sensor[n+PointCount] <= (m_Skappa+m_Kappa))
                {
                    mu_var[n+PointCount] = mu_0*(0.5*(1+sin(
                                           M_PI*(Sensor[n+PointCount]-m_Skappa-m_Kappa)/(2*m_Kappa))));
                }
                else if(Sensor[n+PointCount] > (m_Skappa+m_Kappa))
                {
                    mu_var[n+PointCount] = mu_0;
                }
            }
            
            PointCount += nQuadPointsElement;
        }
    }
    bool CompressibleFlowSystem::v_PostIntegrate(int step)
    {
        NekDouble maxL2 = CalcSteadyState();
        return false;
    }
    
    
    
    // Calculate if the solution reached a steady state
    NekDouble CompressibleFlowSystem::CalcSteadyState()
    {
        /*
         int nPoints = GetTotPoints();
         Array<OneD, NekDouble>        L2   (m_fields.num_elements());
         Array<OneD, NekDouble>        U2np1(m_fields.num_elements());
         static Array<OneD, NekDouble> U2n  (m_fields.num_elements());
         
         // Calculate L2 discrete summation
         for (int i = 0; i < m_fields.num_elements(); ++i)
         {
         U2np1[i]  = 0.0;
         U2np1[i] += Vmath::Dot(nPoints, m_fields[i]->GetPhys(), 1,
         m_fields[i]->GetPhys(), 1);
         m_comm->AllReduce(U2np1[i], LibUtilities::ReduceSum);
         
         L2[i]  = sqrt(fabs(U2np1[i] - U2n[i]) / fabs(U2np1[i]));
         U2n[i] = U2np1[i];
         }
         */
        
        int nPoints = GetTotPoints();
        Array<OneD, NekDouble> L2   (m_fields.num_elements());
        Array<OneD, NekDouble> numer(m_fields.num_elements());
        Array<OneD, NekDouble> denom(m_fields.num_elements());
        Array<OneD, Array<OneD, NekDouble> > unp1(m_fields.num_elements());
        Array<OneD, Array<OneD, NekDouble> > diff(m_fields.num_elements());
        Array<OneD, Array<OneD, NekDouble> > diff2(m_fields.num_elements());
        Array<OneD, Array<OneD, NekDouble> > u2np1(m_fields.num_elements());
        
        
        for (int i = 0; i < m_fields.num_elements(); ++i)
        {
            unp1[i] = Array<OneD, NekDouble>(nPoints, 0.0);
            diff[i] = Array<OneD, NekDouble>(nPoints, 0.0);
            diff2[i] = Array<OneD, NekDouble>(nPoints, 0.0);
            u2np1[i] = Array<OneD, NekDouble>(nPoints, 0.0);
            
            Vmath::Vcopy(nPoints, m_fields[i]->GetPhys(), 1, unp1[i], 1);
            Vmath::Vsub(nPoints, unp1[i], 1, m_un[i], 1, diff[i], 1);
            Vmath::Vmul(nPoints, diff[i], 1, diff[i], 1, diff2[i], 1);
            numer[i] = Vmath::Vsum(nPoints, diff2[i], 1);
            m_comm->AllReduce(numer[i], LibUtilities::ReduceSum);
            
            Vmath::Vmul(nPoints, unp1[i], 1, unp1[i], 1, u2np1[i], 1);
            denom[i] = Vmath::Vsum(nPoints, u2np1[i], 1);
            m_comm->AllReduce(denom[i], LibUtilities::ReduceSum);
            
            L2[i] = sqrt(numer[i]/denom[i]);
            
            Vmath::Vcopy(nPoints, unp1[i], 1, m_un[i], 1);
        }
        
        NekDouble maxL2 = Vmath::Vmax(m_fields.num_elements(), L2, 1);
        
        if (m_fields.num_elements() == 3)
        {
            if (m_comm->GetRank() == 0)
            {
                cout
                /*<< "L2_rho = "  << L2[0] << "    "
                << "L2_rhou = " << L2[1] << "    "
                << "L2_E = "    << L2[2] << "    "*/
                << "L2_max = "  << maxL2 << endl;
            }
        }
        else if (m_fields.num_elements() == 4)
        {
            if (m_comm->GetRank() == 0)
            {
                cout
                /*<< "L2_rho = "  << L2[0] << "    "
                << "L2_rhou = " << L2[1] << "    "
                << "L2_rhov = " << L2[2] << "    "
                << "L2_E = "    << L2[3] << "    "*/
                << "L2_max = "  << maxL2 << endl;
            }
        }
        else if (m_fields.num_elements() == 5)
        {
            if (m_comm->GetRank() == 0)
            {
                cout
                /*<< "L2_rho = "  << L2[0] << "    "
                << "L2_rhou = " << L2[1] << "    "
                << "L2_rhov = " << L2[2] << "    "
                << "L2_rhow = " << L2[3] << "    "
                << "L2_E = "    << L2[4] << "    "*/
                << "L2_max = "  << maxL2 << endl;
            }
        }
        
        std::ofstream myfile;
        myfile.open ("Convergence.txt", std::ios_base::app);
        
        myfile << maxL2 << endl;
        
        myfile.close();
        
        return maxL2;
    }

    void CompressibleFlowSystem::SetVarPOrderElmt(
                const Array<OneD, const Array<OneD, NekDouble> > &physfield,
                      Array<OneD,                    NekDouble  > &PolyOrder)
    {
        int e;
        NekDouble s_ds, s_sm, s_fl;
        
        int nElements  = m_fields[0]->GetExpSize();
        int npts       = m_fields[0]->GetTotPoints();
        
        Array<OneD, NekDouble > Sensor           (npts, 0.0);
        Array<OneD, NekDouble > SensorKappa      (npts, 0.0);
        Array<OneD, NekDouble > se               (npts,0.0);

        GetSensor(physfield, Sensor, SensorKappa);

        int nQuadPointsElement  = 0;
        int npCount             = 0;
        int MinOrder            = 2;
        int MaxOrder            = 12;
        int MinOrderShock       = 4;
        
        std::ofstream m_file( "VariablePComposites.txt", std::ios_base::app);
        for (int e = 0; e < nElements; e++)
        {
            m_file << "<C ID=\"" << e+1 << "\"> Q[" << e << "] </C>"<< endl;
        }
        m_file.close();
        
        std::ofstream m_file2( "VariablePExpansions.txt", std::ios_base::app);
        
        for (e = 0; e < nElements; e++)
        {
            nQuadPointsElement = m_fields[0]->GetExp(e)->GetTotPoints();
            
            // Define thresholds
            // Ideally, these threshold values could be given in the Session File
            
            s_ds =  -5.0;
            //s_ds = s_0*log10(PolyOrder[e]);
            s_sm = -6;
            s_fl = -7;
            
            
            for (int i = 0; i < nQuadPointsElement; i++)
            {
                se[npCount + i] = (Sensor[npCount + i]);
                
                if (se[npCount + i] > s_ds)
                {
                    if (PolyOrder[npCount + i] > MinOrderShock)
                    {
                        PolyOrder[npCount + i] = PolyOrder[npCount + i] - 1;
                    }
                    else if(PolyOrder[e] < MinOrderShock)
                    {
                        PolyOrder[npCount + i] = PolyOrder[npCount + i] + 1;
                    }
                    
                }
                else if(se[npCount + i] > s_sm && se[npCount + i] < s_ds)
                {
                    if (PolyOrder[npCount + i] < MaxOrder)
                    {
                        PolyOrder[npCount + i] = PolyOrder[npCount + i] + 2;
                    }
                }
                else if(se[npCount + i] > s_fl && se[npCount + i] < s_sm)
                {
                    PolyOrder[npCount + i] = PolyOrder[npCount + i] + 1;
                }
                else if(se[npCount + i] < s_fl)
                {
                    if (PolyOrder[npCount + i] > MinOrder)
                    {
                            PolyOrder[npCount + i] = PolyOrder[npCount + i];
                    }
                }
            }
            m_file2 << "<E COMPOSITE= \"C[" << e+1
                    << "]\" NUMMODES=\"" << PolyOrder[npCount + 1]
                    << "\" TYPE=\"MODIFIED\" FIELDS=\"rho,rhou,rhov,rhow,E\" />"
                    << endl;
            npCount += nQuadPointsElement;
        }
        
        m_file2.close();
    }

    void CompressibleFlowSystem::v_ExtraFldOutput(
        std::vector<Array<OneD, NekDouble> > &fieldcoeffs,
        std::vector<std::string>             &variables)
    {
        if (m_adjointSwitch == 0.0)
        {
            const int nPhys   = m_fields[0]->GetNpoints();
            const int nCoeffs = m_fields[0]->GetNcoeffs();
            Array<OneD, Array<OneD, NekDouble> > tmp(m_fields.num_elements());
            
            for (int i = 0; i < m_fields.num_elements(); ++i)
            {
                tmp[i] = m_fields[i]->GetPhys();
            }
            
            Array<OneD, NekDouble> pressure(nPhys), soundspeed(nPhys), mach(nPhys), sensor(nPhys), SensorKappa(nPhys), smooth(nPhys);
            
            GetPressure  (tmp, pressure);
            GetSoundSpeed(tmp, pressure, soundspeed);
            GetMach      (tmp, soundspeed, mach);
            GetSensor    (tmp, sensor, SensorKappa);
            GetSmoothArtificialViscosity    (tmp, smooth);
            
            Array<OneD, NekDouble> pFwd(nCoeffs), sFwd(nCoeffs), mFwd(nCoeffs), sensFwd(nCoeffs), smoothFwd(nCoeffs);
            
            m_fields[0]->FwdTrans(pressure,   pFwd);
            m_fields[0]->FwdTrans(soundspeed, sFwd);
            m_fields[0]->FwdTrans(mach,       mFwd);
            m_fields[0]->FwdTrans(sensor,     sensFwd);
            m_fields[0]->FwdTrans(smooth,     smoothFwd);
            
            variables.push_back  ("p");
            variables.push_back  ("a");
            variables.push_back  ("Mach");
            variables.push_back  ("Sensor");
            variables.push_back  ("SmoothVisc");
            fieldcoeffs.push_back(pFwd);
            fieldcoeffs.push_back(sFwd);
            fieldcoeffs.push_back(mFwd);
            fieldcoeffs.push_back(sensFwd);
            fieldcoeffs.push_back(smoothFwd);
        }
    }
}

