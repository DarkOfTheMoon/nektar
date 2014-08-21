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

#include <CompressibleFlowSolver/EquationSystems/CompressibleFlowSystem.h>
#include <LocalRegions/TriExp.h>
#include <LocalRegions/QuadExp.h>
#include <LocalRegions/HexExp.h>
#include <MultiRegions/ExpList.h>
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
        m_velLoc = Array<OneD, NekDouble>(m_spacedim);
        for (int i = 0; i < m_spacedim; ++i)
        {
            m_velLoc[i] = i+1;
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
        m_session->LoadSolverInfo("Target",         m_Target,        "Lift");
        m_session->LoadParameter ("mu",             m_mu,            1.78e-05);
        m_session->LoadParameter ("thermalConductivity",
                                                    m_thermalConductivity, 0.0257);
        m_session->LoadParameter ("adjointSwitch",  m_adjointSwitch, 0.0);
        m_session->LoadParameter ("rhoInfPrimal",   m_rhoInfPrimal, 0.0);
        m_session->LoadParameter ("uInfPrimal",     m_uInfPrimal, 1.0);
        m_session->LoadParameter ("vInfPrimal",     m_vInfPrimal, 1.0);
        m_session->LoadParameter ("pInfPrimal",     m_pInfPrimal, 1.0);
        m_session->LoadParameter ("alphaInfPrimal", m_alphaInfDir, 0.0);
        m_session->LoadParameter ("Skappa",         m_Skappa,       0.0);
        m_session->LoadParameter ("Kappa",          m_Kappa,         0.0);
        m_session->LoadParameter ("mu0",            m_mu0,           1.0);
        m_session->LoadParameter ("Lref",           m_Lref, 1.0);
        m_session->LoadParameter ("alpha",          m_alpha, 0.0);
        
        m_Cp      = m_gamma / (m_gamma - 1.0) * m_gasConstant;
        m_Prandtl = m_Cp * m_mu / m_thermalConductivity;

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
                // Setting up flux vector for advection operator
                if (m_specHP_dealiasing)
                {
                    m_advection->SetFluxVector(&CompressibleFlowSystem::
                                               GetFluxVectorDeAlias, this);
                }
                else
                {
                    m_advection->SetFluxVector(&CompressibleFlowSystem::
                                               GetFluxVector, this);
                }

                // Setting up flux vector for diffusion operator
                if (m_specHP_dealiasing)
                {
                    m_diffusion->SetFluxVectorNS(
                                            &CompressibleFlowSystem::
                                            GetViscousFluxVectorDeAlias, this);
                }
                else
                {
                    m_diffusion->SetFluxVectorNS(&CompressibleFlowSystem::
                                                 GetViscousFluxVector, this);
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
                m_riemannSolver->SetAuxiliary(
                    "velLoc", &CompressibleFlowSystem::GetVelLoc,  this);
                m_riemannSolver->SetVector(
                    "N",      &CompressibleFlowSystem::GetNormals, this);

                // Setting up parameters for diffusion operator Riemann solver
                m_riemannSolverLDG->SetParam (
                    "gamma",  &CompressibleFlowSystem::GetGamma,   this);
                m_riemannSolverLDG->SetAuxiliary(
                    "velLoc", &CompressibleFlowSystem::GetVelLoc,  this);
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
                
                //m_diffusion->SetArtificialDiffusionVector(&CompressibleFlowSystem::GetArtificialDynamicViscosity, this);
                
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
        int nq = physarray[0].num_elements();
        
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
        
        NekDouble Cinf = 0.5 * m_rhoInfPrimal
                       * (m_uInfPrimal*m_uInfPrimal+m_vInfPrimal*m_vInfPrimal) * m_Lref;
        
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
            
            Array<OneD, NekDouble> tmp(nBCEdgePts, 0.0);
            Array<OneD, NekDouble> tmp1(nBCEdgePts, 0.0);
            Array<OneD, NekDouble> tmp2(nBCEdgePts, 0.0);
            Array<OneD, NekDouble> tmp3(nBCEdgePts, 0.0);
            
            Array<OneD, Array<OneD, NekDouble> > DragDir(m_spacedim);
            Array<OneD, Array<OneD, NekDouble> > LiftDir(m_spacedim);
            
	    if(m_spacedim == 2)
	    {
		NekDouble Dx = norm_fac * cos (m_alphaInfDir);
		NekDouble Dy = norm_fac * sin (m_alphaInfDir);
            
		NekDouble Lx = -norm_fac * sin (m_alphaInfDir);
		NekDouble Ly =  norm_fac * cos (m_alphaInfDir);
            
		DragDir[0] = Array<OneD, NekDouble> (nBCEdgePts, Dx);
		DragDir[1] = Array<OneD, NekDouble> (nBCEdgePts, Dy);
            
		LiftDir[0] = Array<OneD, NekDouble> (nBCEdgePts, Lx);
		LiftDir[1] = Array<OneD, NekDouble> (nBCEdgePts, Ly);
	    }
	    
            if(m_spacedim == 3)
	    {
		NekDouble Dx = 0.0;
		NekDouble Dy = 0.0;
                NekDouble Dz = 0.0;
            
		NekDouble Lx = 0.0;
		NekDouble Ly = 0.0;
		NekDouble Lz = 0.0;
            
		DragDir[0] = Array<OneD, NekDouble> (nBCEdgePts, Dx);
		DragDir[1] = Array<OneD, NekDouble> (nBCEdgePts, Dy);
		DragDir[2] = Array<OneD, NekDouble> (nBCEdgePts, Dz);
            
		LiftDir[0] = Array<OneD, NekDouble> (nBCEdgePts, Lx);
		LiftDir[1] = Array<OneD, NekDouble> (nBCEdgePts, Ly);
		DragDir[2] = Array<OneD, NekDouble> (nBCEdgePts, Lz);
	    }
            
            
            for (i = 0; i < m_spacedim; ++i)
            {

                Vmath::Vvtvp(nBCEdgePts,
                             &Fwd[1+i][id2], 1,
                             &m_traceNormals[i][id2], 1,
                             &tmp[0], 1,
                             &tmp[0], 1);
                
                Vmath::Vvtvp(nBCEdgePts,
                             &DragDir[i][0], 1,
                             &m_traceNormals[i][id2], 1,
                             &tmp1[0], 1,
                             &tmp1[0], 1);
                
                
                Vmath::Vvtvp(nBCEdgePts,
                             &LiftDir[i][0], 1,
                             &m_traceNormals[i][id2], 1,
                             &tmp2[0], 1,
                             &tmp2[0], 1);
            }
            //

            Vmath::Vsub(nBCEdgePts, &tmp1[0], 1, &tmp[0], 1, &tmp3[0], 1);
            // Calculate 2.0(1/Cinf(phi.n)-v.n)
            Vmath::Smul(nBCEdgePts, 2.0, &tmp3[0], 1, &tmp3[0], 1);
            // Calculate 1/Cinf(phi.n)
            
            // Calculate v* = v - 2.0(1/Cinf(phi.n)-v.n)n
            for (i = 0; i < m_spacedim; ++i)
            {
                Vmath::Vvtvp(nBCEdgePts,
                             &tmp3[0], 1,
                             &m_traceNormals[i][id2], 1,
                             &Fwdnew[1+i][id2], 1,
                             &Fwdnew[1+i][id2], 1);
            }
	    
	    Array<OneD, NekDouble> ones(nBCEdgePts, 1.0);
            
            //Vmath::Neg(nBCEdgePts,&Fwdnew[0][id2], 1);
            Vmath::Neg(nBCEdgePts,&Fwdnew[nVariables-1][id2], 1);
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

        // Adjust the physical values of the trace to
        // take user defined boundaries into account
        int e, id1, id2, nBCEdgePts, eMax;

        eMax = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExpSize();

        for (e = 0; e < eMax; ++e)
        {
            nBCEdgePts = m_fields[0]->GetBndCondExpansions()[bcRegion]->
                GetExp(e)->GetTotPoints();
            id1  = m_fields[0]->GetBndCondExpansions()[bcRegion]->
                GetPhys_Offset(e);
            id2 = m_fields[0]->GetTrace()->GetPhys_Offset(traceBndMap[cnt+e]);

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
     * @brief Simmetry boundary conditions for compressible flow problems.
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

        int e, id1, id2, nBCEdgePts, eMax;

        eMax = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExpSize();

        for (e = 0; e < eMax; ++e)
        {
            nBCEdgePts = m_fields[0]->GetBndCondExpansions()[bcRegion]->
                GetExp(e)->GetTotPoints();
            id1  = m_fields[0]->GetBndCondExpansions()[bcRegion]->
                GetPhys_Offset(e);
            id2 = m_fields[0]->GetTrace()->GetPhys_Offset(traceBndMap[cnt+e]);

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
            Vmath::Smul(nTracePts, m_vInf, m_traceNormals[0], 1, tmp1, 1);
            Vmath::Vadd(nTracePts, VnInf, 1, tmp1, 1, VnInf, 1);
        }
        if (nDimensions == 3)
        {
            velInf[2] = m_wInf;
            Vmath::Smul(nTracePts, m_wInf, m_traceNormals[0], 1, tmp2, 1);
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

    
    void CompressibleFlowSystem::AdjointPressureOutflow(
                            int                                   bcRegion,
                            int                                   cnt,
                            Array<OneD, Array<OneD, NekDouble> > &physarray)
    {
        int i;
        int nTracePts = GetTraceTotPoints();
        int nVariables = physarray.num_elements();
        int nq = physarray[0].num_elements();
        
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
            
            Array<OneD, NekDouble> qterm    (nBCEdgePts,     0.0);
            
            for (int i = 0; i < m_spacedim; ++i)
            {
                Vmath::Vvtvp(nBCEdgePts,
                             &Fwd[i+1][id2], 1,
                             &Fwd[i+1][id2], 1,
                             &qterm[0],         1,
                             &qterm[0],         1);
            }
            //Calculate Bwd[1]^2/Bwd[0] + Bwd[2]^2/Bwd[0]
            /*Vmath::Vdiv(nBCEdgePts,
                        &qterm[0], 1,
                        &Fwd[0][id2], 1,
                        &Fwdnew[m_spacedim+1][id2], 1);*/
            
            for (i = 0; i < nVariables; ++i)
            {
                Vmath::Vcopy(nBCEdgePts, &Fwdnew[i][id2], 1,
                             &(m_fields[i]->GetBndCondExpansions()[bcRegion]->
                               UpdatePhys())[id1], 1);
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
    }
    
    
    
    void CompressibleFlowSystem::GetJacTransposeDivVector(
            const Array<OneD, Array<OneD, NekDouble> > &inarray,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &outarray)
    {
        int i, j, k, n, p, o, t;
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
        
        int i, j, k, n, p, o, t;
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
        int bcRegion = 0;
        int eMax     = 0;
        
        // =============================HACK!!!!================================
        int cnt = 0;
        const Array<OneD, const int> &traceBndMap =
                                                  m_fields[0]->GetTraceBndMap();
        
        int nBCregions = m_fields[0]->GetBndConditions().num_elements();
        
        for (int k = 0; k < nBCregions; ++k)
        {
            if (m_fields[0]->GetBndConditions()[k]->GetUserDefined() ==
                SpatialDomains::eAdjointWall)
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
            
            if (m_fields[0]->GetBndConditions()[k]->GetUserDefined() ==
                SpatialDomains::eAdjointPressureOutflow)
            {
                eMax = m_fields[0]->GetBndCondExpansions()[k]->GetExpSize();
                
                for (int e = 0; e < eMax; ++e)
                {
                    int nBCEdgePts = m_primal[0]->GetBndCondExpansions()[k]->
                    GetExp(e)->GetTotPoints();
                    int id1 = m_primal[0]->GetBndCondExpansions()[k]->
                    GetPhys_Offset(e);
                    int id2 = m_primal[0]->GetTrace()->GetPhys_Offset(traceBndMap[cnt+e]);
                    
                    // T = Pinf at outflow condition,
                    // P = (gamma - 1) (rhoE - rho * q^2)
                    // => rhoE = Pinf/ (gamma - 1) + rho * q^2;
                    
                    // at pressure outflow for 2 D:
                    
                    // rho =>  Bwd[0] = Fwd[0]
                    // rhou => Bwd[1] = Fwd[1]
                    // rhov => Bwd[2] = Fwd[2]
                    // rhoE => Bwd[3] = (Pinf + Bwd[1]^2/Bwd[0] + Bwd[2]^2/Bwd[0]) / (gamma - 1);
                    
                    Vmath::Vcopy(nBCEdgePts,
                                 &FwdDir[0][id2], 1,
                                 &BwdDir[0][id2], 1);
                    
                    Vmath::Vcopy(nBCEdgePts,
                                 &FwdDir[1][id2], 1,
                                 &BwdDir[1][id2], 1);
                    
                    Vmath::Vcopy(nBCEdgePts,
                                 &FwdDir[2][id2], 1,
                                 &BwdDir[2][id2], 1);
                    
                    NekDouble pInfInv = m_pInfPrimal / (m_gamma - 1.0);
                    
                    Array<OneD, NekDouble> pInfFac (nBCEdgePts, pInfInv);
                    Array<OneD, NekDouble> qterm    (nBCEdgePts,     0.0);
                    
                    // Calculate qterm = Bwd[1]^2 + Bwd[2]^2
                    for (int i = 0; i < m_spacedim; ++i)
                    {
                        Vmath::Vvtvp(nBCEdgePts,
                                     &BwdDir[i+1][id2], 1,
                                     &BwdDir[i+1][id2], 1,
                                     &qterm[0],         1,
                                     &qterm[0],         1);
                    }
                    //Calculate Bwd[1]^2/Bwd[0] + Bwd[2]^2/Bwd[0]
                    Vmath::Vdiv(nBCEdgePts,
                                &qterm[0], 1,
                                &BwdDir[0][id2], 1,
                                &BwdDir[m_spacedim+1][id2], 1);
                    //Calculate pInf/(gamma-1) + Bwd[1]^2/Bwd[0] + Bwd[2]^2/Bwd[0]
                    Vmath::Vadd(nBCEdgePts,
                                &pInfFac[0], 1,
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
                SpatialDomains::eAdjointPressureOutflow && m_fields[0]->GetBndConditions()[k]->GetUserDefined() !=
                SpatialDomains::eAdjointWall)
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
                    Array<OneD, NekDouble> ones(nBCEdgePts, 1.0);
                    
                    NekDouble rhouInfPrimal = m_rhoInfPrimal*m_uInfPrimal;
                    NekDouble rhovInfPrimal = m_rhoInfPrimal*m_vInfPrimal;
                    
                    NekDouble rhoEInfPrimal = m_pInfPrimal/(m_gamma - 1)
                    + 0.5 * m_rhoInfPrimal * (m_uInfPrimal * m_uInfPrimal + m_vInfPrimal * m_vInfPrimal);
                    
                    Vmath::Smul(nBCEdgePts,
                                m_rhoInfPrimal,
                                &ones[0], 1,
                                &BwdDir[0][id2], 1);
                    Vmath::Smul(nBCEdgePts,
                                rhouInfPrimal,
                                &ones[0], 1,
                                &BwdDir[1][id2], 1);
                    Vmath::Smul(nBCEdgePts,
                                rhovInfPrimal,
                                &ones[0], 1,
                                &BwdDir[2][id2], 1);
                    Vmath::Smul(nBCEdgePts,
                                rhoEInfPrimal,
                                &ones[0], 1,
                                &BwdDir[3][id2], 1);
                    
                    for (int i = 0; i < nConvectiveFields; ++i)
                    {
                        Vmath::Vcopy(nBCEdgePts,
                                     &FwdDir[i][id2], 1,
                                     &BwdDir[i][id2], 1);
                        
                    }
                    
                    for (int i = 0; i < nConvectiveFields; ++i)
                    {
                        Vmath::Vcopy(nBCEdgePts, &BwdDir[i][id2], 1,
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
        int nTracePoints       = FwdDir[0].num_elements();
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
        

        
        
        /*int nConvectiveFields = FwdDir.num_elements();
        int nTracePoints       = FwdDir[0].num_elements();
        int i;
        Array<OneD, Array<OneD, NekDouble> > velFwd(m_spacedim);
        Array<OneD, Array<OneD, NekDouble> > velBwd(m_spacedim);
        
        Array<OneD, NekDouble> PresFwd(nTracePoints, 0.0);
        Array<OneD, NekDouble> PresBwd(nTracePoints, 0.0);
        
        if (m_spacedim == 2)
        {
            for (i = 0; i < m_spacedim; ++i)
            {
                velFwd[i] = Array<OneD, NekDouble>(nTracePoints, 0.0);
                velBwd[i] = Array<OneD, NekDouble>(nTracePoints, 0.0);
                
                Vmath::Vdiv(nTracePoints,
                            FwdDir[i+1], 1,
                            FwdDir[0], 1,
                            velFwd[i], 1);
                Vmath::Vdiv(nTracePoints,
                            BwdDir[i+1], 1,
                            BwdDir[0], 1,
                            velBwd[i], 1);
                
                m_primal[0]->PhysDeriv(velFwd[i],
                                       FwdDirDIFF[0][i],
                                       FwdDirDIFF[1][i]);
                
                m_primal[0]->PhysDeriv(velBwd[i],
                                       BwdDirDIFF[0][i],
                                       BwdDirDIFF[1][i]);
            }
            
            GetPressure(FwdDir, velFwd, PresFwd);
            GetPressure(FwdDir, velFwd, PresFwd);
            
            // density gradients
            m_primal[0]->PhysDeriv(FwdDir[0], FwdDirDIFF[0][0],
                                              FwdDirDIFF[1][0]);
            
            m_primal[0]->PhysDeriv(BwdDir[0], BwdDirDIFF[0][0],
                                              BwdDirDIFF[1][0]);
            
            // pressure gradients
            m_primal[0]->PhysDeriv(PresFwd, FwdDirDIFF[0][3],
                                            FwdDirDIFF[1][3]);
            
            m_primal[0]->PhysDeriv(PresBwd, BwdDirDIFF[0][3],
                                            BwdDirDIFF[1][3]);
            
        }
        if (m_spacedim == 3)
        {
            for (i = 0; i < m_spacedim; ++i)
            {
                velFwd[i] = Array<OneD, NekDouble>(nTracePoints, 0.0);
                velBwd[i] = Array<OneD, NekDouble>(nTracePoints, 0.0);
                
                Vmath::Vdiv(nTracePoints,
                            FwdDir[i+1], 1,
                            FwdDir[0], 1,
                            velFwd[i], 1);
                Vmath::Vdiv(nTracePoints,
                            BwdDir[i+1], 1,
                            BwdDir[0], 1,
                            velBwd[i], 1);
                
                m_primal[0]->PhysDeriv(velFwd[i], FwdDirDIFF[0][i],
                                                  FwdDirDIFF[1][i],
                                                  BwdDirDIFF[2][i]);
                
                m_primal[0]->PhysDeriv(velBwd[i], BwdDirDIFF[0][i],
                                                  BwdDirDIFF[1][i],
                                                  BwdDirDIFF[2][i]);
            }
            
            GetPressure(FwdDir, velFwd, PresFwd);
            GetPressure(FwdDir, velFwd, PresFwd);
            
            // density gradients
            m_primal[0]->PhysDeriv(FwdDir[0], FwdDirDIFF[0][0],
                                              FwdDirDIFF[1][0],
                                              FwdDirDIFF[2][0]);
            
            m_primal[0]->PhysDeriv(BwdDir[0], BwdDirDIFF[0][0],
                                              BwdDirDIFF[1][0],
                                              BwdDirDIFF[2][0]);
            
            
            // pressure gradients
            m_primal[0]->PhysDeriv(PresFwd, FwdDirDIFF[0][4],
                                            FwdDirDIFF[1][4],
                                            FwdDirDIFF[2][4]);
            
            m_primal[0]->PhysDeriv(PresBwd, BwdDirDIFF[0][4],
                                            BwdDirDIFF[1][4],
                                            BwdDirDIFF[2][4]);
            
        }
        
        */
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
        Array<OneD, NekDouble > mu    (nPts, 0.0);
        Array<OneD, NekDouble > mu2   (nPts, 0.0);
        Array<OneD, NekDouble > divVel(nPts, 0.0);

        // Variable viscosity through the Sutherland's law
        if (m_ViscosityType == "Variable")
        {
            GetDynamicViscosity(physfield[nVariables-2], mu);
        }
        else
        {
            Vmath::Fill(nPts, m_mu, &mu[0], 1);
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
            Vmath::Smul(nPts, m_thermalConductivity, &derivativesO1[0][1][0], 1,
                        &tmp1[0], 1);

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
            Vmath::Smul(nPts, m_thermalConductivity, &derivativesO1[0][2][0], 1,
                        &tmp2[0], 1);

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
            Vmath::Smul(nPts, m_thermalConductivity, &derivativesO1[1][2][0], 1,
                        &tmp2[0], 1);

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
            Vmath::Smul(nPts, m_thermalConductivity, &derivativesO1[0][3][0], 1,
                        &tmp3[0], 1);

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
            Vmath::Smul(nPts, m_thermalConductivity, &derivativesO1[1][3][0], 1,
                        &tmp3[0], 1);

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
            Vmath::Smul(nPts, m_thermalConductivity, &derivativesO1[2][3][0], 1,
                        &tmp3[0], 1);

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
            GetDynamicViscosity(fields_interp, mu);
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
        const int nBCEdgePts = m_fields[0]->GetTotPoints();

        Vmath::Vmul(nBCEdgePts, physfield[1], 1, physfield[1], 1, mach, 1);

        for (int i = 1; i < m_spacedim; ++i)
        {
            Vmath::Vvtvp(nBCEdgePts, physfield[1+i], 1, physfield[1+i], 1,
                               mach,           1, mach,           1);
        }

        Vmath::Vdiv(nBCEdgePts, mach, 1, physfield[0], 1, mach, 1);
        Vmath::Vdiv(nBCEdgePts, mach, 1, physfield[0], 1, mach, 1);
        Vmath::Vdiv(nBCEdgePts, mach, 1, soundspeed,   1, mach, 1);
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

            if (boost::dynamic_pointer_cast<LocalRegions::TriExp>(
                    m_fields[0]->GetExp(n)))
            {
                minLength = 2.0 * sqrt(Area);
            }

            else if (boost::dynamic_pointer_cast<LocalRegions::QuadExp>(
                         m_fields[0]->GetExp(n)))
            {
                minLength = sqrt(Area);
            }
            else if (boost::dynamic_pointer_cast<LocalRegions::HexExp>(
                         m_fields[0]->GetExp(n)))
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
        int i, j, k, n, p, o, t;
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
        int i, j, k, n, p, o, t;
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
        int i, j, k, n, p, o, t;
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
            
            /*Vmath::Zero(nq, dUdVInvdY[0][0], 1);
            Vmath::Zero(nq, dUdVInvdY[0][1], 1);
            Vmath::Zero(nq, dUdVInvdY[0][2], 1);
            Vmath::Zero(nq, dUdVInvdY[0][3], 1);
            
            Vmath::Zero(nq, dUdVInvdY[1][0], 1);
            Vmath::Zero(nq, dUdVInvdY[1][1], 1);
            Vmath::Zero(nq, dUdVInvdY[1][2], 1);
            Vmath::Zero(nq, dUdVInvdY[1][3], 1);
            
            Vmath::Zero(nq, dUdVInvdY[2][0], 1);
            Vmath::Zero(nq, dUdVInvdY[2][1], 1);
            Vmath::Zero(nq, dUdVInvdY[2][2], 1);
            Vmath::Zero(nq, dUdVInvdY[2][3], 1);
            
            Vmath::Zero(nq, dUdVInvdY[3][0], 1);
            Vmath::Zero(nq, dUdVInvdY[3][1], 1);
            Vmath::Zero(nq, dUdVInvdY[3][2], 1);
            Vmath::Zero(nq, dUdVInvdY[3][3], 1);
            
            Vmath::Zero(nq, dUdVInvdX[0][0], 1);
            Vmath::Zero(nq, dUdVInvdX[0][1], 1);
            Vmath::Zero(nq, dUdVInvdX[0][2], 1);
            Vmath::Zero(nq, dUdVInvdX[0][3], 1);
            
            Vmath::Zero(nq, dUdVInvdX[1][0], 1);
            Vmath::Zero(nq, dUdVInvdX[1][1], 1);
            Vmath::Zero(nq, dUdVInvdX[1][2], 1);
            Vmath::Zero(nq, dUdVInvdX[1][3], 1);
            
            Vmath::Zero(nq, dUdVInvdX[2][0], 1);
            Vmath::Zero(nq, dUdVInvdX[2][1], 1);
            Vmath::Zero(nq, dUdVInvdX[2][2], 1);
            Vmath::Zero(nq, dUdVInvdX[2][3], 1);
            
            Vmath::Zero(nq, dUdVInvdX[3][0], 1);
            Vmath::Zero(nq, dUdVInvdX[3][1], 1);
            Vmath::Zero(nq, dUdVInvdX[3][2], 1);
            Vmath::Zero(nq, dUdVInvdX[3][3], 1);

            for (i = 0; i < m_spacedim; ++i)
            {
                m_primal[0]->PhysDeriv(vel[i], dui_dx[i], dui_dy[i]);
                
                // Calculating the derivatives of dui2/dxi
                
                Vmath::Vmul(nq,
                            &vel[i][0], 1,
                            &dui_dx[i][0], 1,
                            &du2i_dx[i][0], 1);
                
                Vmath::Smul(nq,
                            2.0,
                            &du2i_dx[i][0], 1,
                            &du2i_dx[i][0], 1);
                
                Vmath::Vmul(nq,
                            &vel[i][0], 1,
                            &dui_dy[i][0], 1,
                            &du2i_dy[i][0], 1);
                
                Vmath::Smul(nq,
                            2.0,
                            &du2i_dy[i][0], 1,
                            &du2i_dy[i][0], 1);
            }
            
            Vmath::Vadd(nq, &du2i_dx[0][0], 1, &du2i_dx[1][0], 1, &dU2_dxi[0][0], 1);
            Vmath::Vadd(nq, &du2i_dy[0][0], 1, &du2i_dy[1][0], 1, &dU2_dxi[1][0], 1);
            // Calculating the derivatives of dU2/dxi
            for (j = 0; j < m_spacedim; ++j)
            {
                /*Vmath::Vadd(nq,
                 &du2i_dx[j][0], 1,
                 &du2i_dx[j][0], 1,
                 &dU2_dxi[0][0], 1);
                 
                 Vmath::Vadd(nq,
                 &du2i_dy[j][0], 1,
                 &du2i_dy[j][0], 1,
                 &dU2_dxi[1][0], 1);*/
           /*
                Vmath::Vmul(nq,
                            &OneMinGamVec[0], 1,
                            &dui_dx[j][0], 1,
                            &OneMinGamVeldX[j][0], 1);
                
                Vmath::Vmul(nq,
                            &OneMinGamVec[0], 1,
                            &dui_dy[j][0], 1,
                            &OneMinGamVeldY[j][0], 1);
            }
            
            m_primal[0]->PhysDeriv(velsq, dU2_dxi[0], dU2_dxi[1]);
            
            Vmath::Vmul(nq,
                        &GamMinOneOverTwoVec[0], 1,
                        &dU2_dxi[0][0], 1,
                        &GamMinOneOverTwoVelsqdX[0], 1);
            
            Vmath::Vmul(nq,
                        &GamMinOneOverTwoVec[0], 1,
                        &dU2_dxi[1][0], 1,
                        &GamMinOneOverTwoVelsqdY[0], 1);
            
            m_primal[0]->PhysDeriv(dirFlds[0], drho_dxi[0], drho_dxi[1]);
            */
            /* =================================================================
             if
             M^-1 =
             [             1,             0,            0,            0]
             [        -u/rho,         1/rho,            0,            0]
             [        -v/rho,             0,        1/rho,            0]
             [(gamma-1)U^2/2,   (1-gamma)*u,  (1-gamma)*v,    (gamma-1)]
             
             then
             
             dM^-1/dx =
             [                    0,                 0,              0,  0]
             [         d(-u/rho)/dx,       d(1/rho)/dx,              0,  0]
             [         d(-v/rho)/dx,                 0,    d(1/rho)/dx,  0]
             [(gamma-1)/2*d(U^2)/dx,   (1-gamma)*du/dx,(1-gamma)*dv/dx,  0]
             */
            // =================================================================
            
            // calculate d(1/rho)/dx
            /*Vmath::Zero(nq, &tmp[0], 1);
            Vmath::Zero(nq, &tmp1[0], 1);
            
            Vmath::Vmul(nq, &dirFlds[0][0], 1, &dirFlds[0][0], 1, &tmp[0], 1);
            
            Vmath::Vcopy(nq, &drho_dxi[0][0], 1, &d1overRho_dxi[0][0], 1);
            
            Vmath::Vdiv(nq,
                        &d1overRho_dxi[0][0], 1,
                        &tmp[0], 1,
                        &d1overRho_dxi[0][0], 1);
            
            Vmath::Neg(nq, &d1overRho_dxi[0][0], 1);
            
            // calculate d(-u/rho)/dx
            Vmath::Zero(nq, &tmp[0], 1);
            Vmath::Zero(nq, &tmp1[0], 1);
            
            Vmath::Vmul(nq,
                        &dirFlds[0][0], 1,
                        &dirFlds[0][0], 1,
                        &tmp[0], 1);
            
            Vmath::Vmul(nq,
                        &vel[0][0], 1,
                        &drho_dxi[0][0], 1,
                        &tmp1[0], 1);
            
            Vmath::Vmul(nq,
                        &dirFlds[0][0], 1,
                        &dui_dx[0][0], 1,
                        &duiOverRho_dx[0][0], 1);
            
            Vmath::Vsub(nq,
                        &tmp1[0], 1,
                        &duiOverRho_dx[0][0], 1,
                        &duiOverRho_dx[0][0], 1);
            
            Vmath::Vdiv(nq,
                        &duiOverRho_dx[0][0], 1,
                        &tmp[0], 1,
                        &duiOverRho_dx[0][0], 1);
            
            // calculate d(-v/rho)/dx
            Vmath::Zero(nq, &tmp[0], 1);
            Vmath::Zero(nq, &tmp1[0], 1);
            
            Vmath::Vmul(nq,
                        &dirFlds[0][0], 1,
                        &dirFlds[0][0], 1,
                        &tmp[0], 1);
            
            Vmath::Vmul(nq, &vel[1][0], 1, &drho_dxi[0][0], 1, &tmp1[0], 1);
            
            Vmath::Vmul(nq,
                        &dirFlds[0][0], 1,
                        &dui_dx[1][0], 1,
                        &duiOverRho_dx[1][0], 1);
            
            Vmath::Vsub(nq,
                        &tmp1[0], 1,
                        &duiOverRho_dx[1][0], 1,
                        &duiOverRho_dx[1][0], 1);
            
            Vmath::Vdiv(nq,
                        &duiOverRho_dx[1][0], 1,
                        &tmp[0], 1,
                        &duiOverRho_dx[1][0], 1);
            
            // Set up the matrix for d(M^-1)/dx
            // row one
            Vmath::Vcopy(nq, &zeros[0],              1, &dUdVInvdX[0][0][0], 1);
            Vmath::Vcopy(nq, &zeros[0],              1, &dUdVInvdX[1][0][0], 1);
            Vmath::Vcopy(nq, &zeros[0],              1, &dUdVInvdX[2][0][0], 1);
            Vmath::Vcopy(nq, &zeros[0],              1, &dUdVInvdX[3][0][0], 1);
            // row two
            Vmath::Vcopy(nq, &duiOverRho_dx[0][0],   1, &dUdVInvdX[0][1][0], 1);
            Vmath::Vcopy(nq, &d1overRho_dxi[0][0],   1, &dUdVInvdX[1][1][0], 1);
            Vmath::Vcopy(nq, &zeros[0],              1, &dUdVInvdX[2][1][0], 1);
            Vmath::Vcopy(nq, &zeros[0],              1, &dUdVInvdX[3][1][0], 1);
            // row three
            Vmath::Vcopy(nq, &duiOverRho_dx[1][0],   1, &dUdVInvdX[0][2][0], 1);
            Vmath::Vcopy(nq, &zeros[0],              1, &dUdVInvdX[1][2][0], 1);
            Vmath::Vcopy(nq, &d1overRho_dxi[0][0],   1, &dUdVInvdX[2][2][0], 1);
            Vmath::Vcopy(nq, &zeros[0],              1, &dUdVInvdX[3][2][0], 1);
            // row four
            Vmath::Vcopy(nq, &GamMinOneOverTwoVelsqdX[0], 1, &dUdVInvdX[0][3][0], 1);
            Vmath::Vcopy(nq, &OneMinGamVeldX[0][0],  1, &dUdVInvdX[1][3][0], 1);
            Vmath::Vcopy(nq, &OneMinGamVeldX[1][0],  1, &dUdVInvdX[2][3][0], 1);
            Vmath::Vcopy(nq, &zeros[0],              1, &dUdVInvdX[3][3][0], 1);
            */
            /* =================================================================
             if
             M^-1 =
             [             1,             0,            0,            0]
             [        -u/rho,         1/rho,            0,            0]
             [        -v/rho,             0,        1/rho,            0]
             [(gamma-1)U^2/2,   (1-gamma)*u,  (1-gamma)*v,    (gamma-1)]
             
             then
             
             dM^-1/dy =
             [                    0,                 0,              0,  0]
             [         d(-u/rho)/dy,       d(1/rho)/dy,              0,  0]
             [         d(-v/rho)/dy,                 0,    d(1/rho)/dy,  0]
             [(gamma-1)/2*d(U^2)/dy,   (1-gamma)*du/dy,(1-gamma)*dv/dy,  0]
             */
            // =================================================================
            
            // calculate d(1/rho)/dy
            /*Vmath::Zero(nq, &tmp[0], 1);
            Vmath::Zero(nq, &tmp1[0], 1);
            
            Vmath::Vmul(nq, &dirFlds[0][0], 1, &dirFlds[0][0], 1, &tmp[0], 1);
            Vmath::Vcopy(nq, &drho_dxi[1][0], 1, &d1overRho_dxi[1][0], 1);
            
            Vmath::Vdiv(nq,
                        &d1overRho_dxi[1][0], 1,
                        &tmp[0], 1,
                        &d1overRho_dxi[1][0], 1);
            
            Vmath::Neg(nq, &d1overRho_dxi[1][0], 1);
            
            // calculate d(-u/rho)/dy
            Vmath::Zero(nq, &tmp[0], 1);
            Vmath::Zero(nq, &tmp1[0], 1);
            Vmath::Vmul(nq, &dirFlds[0][0], 1, &dirFlds[0][0], 1, &tmp[0], 1);
            Vmath::Vmul(nq, &vel[0][0], 1, &drho_dxi[1][0], 1, &tmp1[0], 1);
            
            Vmath::Vmul(nq,
                        &dirFlds[0][0], 1,
                        &dui_dy[0][0], 1,
                        &duiOverRho_dy[0][0], 1);
            
            Vmath::Vsub(nq,
                        &tmp1[0], 1,
                        &duiOverRho_dy[0][0], 1,
                        &duiOverRho_dy[0][0], 1);
            
            Vmath::Vdiv(nq,
                        &duiOverRho_dy[0][0], 1,
                        &tmp[0], 1,
                        &duiOverRho_dy[0][0], 1);
            
            // calculate d(-v/rho)/dy
            Vmath::Zero(nq, &tmp[0], 1);
            Vmath::Zero(nq, &tmp1[0], 1);
            Vmath::Vmul(nq,
                        &dirFlds[0][0], 1,
                        &dirFlds[0][0], 1,
                        &tmp[0], 1);
            
            Vmath::Vmul(nq, &vel[1][0], 1, &drho_dxi[1][0], 1, &tmp1[0], 1);
            
            Vmath::Vmul(nq,
                        &dirFlds[0][0], 1,
                        &dui_dy[1][0], 1,
                        &duiOverRho_dy[1][0], 1);
            
            Vmath::Vsub(nq,
                        &tmp1[0], 1,
                        &duiOverRho_dy[1][0], 1,
                        &duiOverRho_dy[1][0], 1);
            
            Vmath::Vdiv(nq,
                        &duiOverRho_dy[1][0], 1,
                        &tmp[0], 1,
                        &duiOverRho_dy[1][0], 1);
            
            // Set up the matrix for d(M^-1)/dy
            // row one
            Vmath::Vcopy(nq, &zeros[0],              1, &dUdVInvdY[0][0][0], 1);
            Vmath::Vcopy(nq, &zeros[0],              1, &dUdVInvdY[1][0][0], 1);
            Vmath::Vcopy(nq, &zeros[0],              1, &dUdVInvdY[2][0][0], 1);
            Vmath::Vcopy(nq, &zeros[0],              1, &dUdVInvdY[3][0][0], 1);
            // row two
            Vmath::Vcopy(nq, &duiOverRho_dy[0][0],   1, &dUdVInvdY[0][1][0], 1);
            Vmath::Vcopy(nq, &d1overRho_dxi[1][0],   1, &dUdVInvdY[1][1][0], 1);
            Vmath::Vcopy(nq, &zeros[0],              1, &dUdVInvdY[2][1][0], 1);
            Vmath::Vcopy(nq, &zeros[0],              1, &dUdVInvdY[3][1][0], 1);
            // row three
            Vmath::Vcopy(nq, &duiOverRho_dy[1][0],   1, &dUdVInvdY[0][2][0], 1);
            Vmath::Vcopy(nq, &zeros[0],              1, &dUdVInvdY[1][2][0], 1);
            Vmath::Vcopy(nq, &d1overRho_dxi[1][0],   1, &dUdVInvdY[2][2][0], 1);
            Vmath::Vcopy(nq, &zeros[0],              1, &dUdVInvdY[3][2][0], 1);
            // row four
            Vmath::Vcopy(nq, &GamMinOneOverTwoVelsqdY[0], 1, &dUdVInvdY[0][3][0], 1);
            Vmath::Vcopy(nq, &OneMinGamVeldY[0][0],  1, &dUdVInvdY[1][3][0], 1);
            Vmath::Vcopy(nq, &OneMinGamVeldY[1][0],  1, &dUdVInvdY[2][3][0], 1);
            Vmath::Vcopy(nq, &zeros[0],              1, &dUdVInvdY[3][3][0], 1);*/
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
        int i,j,k;
        
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
        for (i = 0; i < nvar; ++i)
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
                                 &DxxCons[i][j][0], 1,
                                 &derivatives[0][j][0], 1,
                                 &tmpx1[i][0], 1,
                                 &tmpx1[i][0], 1);
                    
                    Vmath::Vvtvp(nq,
                                 &DyxCons[i][j][0], 1,
                                 &derivatives[1][j][0], 1,
                                 &tmpx2[i][0], 1,
                                 &tmpx2[i][0], 1);
                    
                    Vmath::Vvtvp(nq,
                                 &DxyCons[i][j][0], 1,
                                 &derivatives[0][j][0], 1,
                                 &tmpy1[i][0], 1,
                                 &tmpy1[i][0], 1);
                    
                    Vmath::Vvtvp(nq,
                                 &DyyCons[i][j][0], 1,
                                 &derivatives[1][j][0], 1,
                                 &tmpy2[i][0], 1,
                                 &tmpy2[i][0], 1);
                }
                
                Vmath::Vadd(nq,
                            &tmpx1[i][0], 1,
                            &tmpx2[i][0], 1,
                            &outarray[0][i][0], 1);
                
                Vmath::Vadd(nq,
                            &tmpy1[i][0], 1,
                            &tmpy2[i][0], 1,
                            &outarray[1][i][0], 1);
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
        int i, j, k;
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
        
        int i, e, nCoeffsElement, NumModesElement, NumModesCuttOff, NumModesCuttOff_Dir1, NumModesCuttOff_Dir2, nQuadPointsElement;
        NekDouble SensorNumerator, SensorDenominator;
        
        int nVariables      = m_fields.num_elements();
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
            
            m_primal[0]->GetExp(e)->FwdTrans(SolPElementPhys,
                                             SolPElementCoeffs);
            
            // ReduceOrderCoeffs reduces the polynomial order of the solution that
            // is represented by the coeffs given as an inarray. This is done by
            // projecting the higher order solution onto the orthogonal basis and
            // padding the higher order coefficients with zeros.
            
            m_primal[0]->GetExp(e)->ReduceOrderCoeffs(numCutOff,
                                                      SolPElementCoeffs,
                                                      SolPmOneElementCoeffs);
            
            m_primal[0]->GetExp(e)->BwdTrans(SolPmOneElementCoeffs,
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
                Sensor[CoeffsCount+i] = sqrt(SolPmeanNumerator/nQuadPointsElement)
                /sqrt(SolPmeanDenumerator/nQuadPointsElement);
                
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
                            (Sensor[CoeffsCount+i]-Phi0)/(2*DeltaPhi)));
                }
            }
            
            CoeffsCount += nQuadPointsElement;
        }
     }
    
    void CompressibleFlowSystem::GetArtificialDynamicViscosity(
                const Array<OneD, Array<OneD, NekDouble> > &physfield,
                      Array<OneD,             NekDouble  > &mu_var)
    {
        const int npts       = m_fields[0]->GetTotPoints();
        const int nElements  = m_fields[0]->GetExpSize();
        const int ploc       = m_fields[0]->GetExpSize();
        int nVariables       = physfield.num_elements();
        
        int PointCount = 0;
        int nTotQuadPoints  = GetTotPoints();
        
        Array <OneD , NekDouble > S_e               (nTotQuadPoints, 0.0);
        Array <OneD , NekDouble > se                (nTotQuadPoints, 0.0);
        Array <OneD, NekDouble > Sensor             (nTotQuadPoints, 0.0);
        Array <OneD, NekDouble > SensorKappa        (nTotQuadPoints, 0.0);
        Array <OneD, NekDouble > absVelocity        (nTotQuadPoints, 0.0);
        Array <OneD, NekDouble > soundspeed         (nTotQuadPoints, 0.0);
        Array <OneD, NekDouble > pressure           (nTotQuadPoints, 0.0);
        
        Array <OneD, Array <OneD, NekDouble > > direct_fields     (nVariables);
        
        for (int i = 0; i < nVariables; ++i)
        {
            direct_fields[i] = Array<OneD, NekDouble>(nTotQuadPoints, 0.0);
        }
        
        GetDirectSolution(direct_fields);
        
        //GetAbsoluteVelocity   (physfield, absVelocity);
        GetPressure      (direct_fields, pressure);
        GetSoundSpeed    (direct_fields, pressure, soundspeed);
        GetSensor        (direct_fields, Sensor, SensorKappa);
        
        Array<OneD,int> pOrderElmt = GetNumExpModesPerExp();
        
        Array <OneD , NekDouble > Lambda(nTotQuadPoints, 1.0);
        Vmath::Vadd(nTotQuadPoints,absVelocity,1,soundspeed,1,Lambda,1);
        
        // Determining the maximum wave speed in the element
        NekDouble LambdaMax = Vmath::Vmax(nTotQuadPoints,Lambda,1);
        NekDouble StdVMax = Vmath::Vmax(nTotQuadPoints,absVelocity,1);
        
        for (int e = 0; e < nElements; e++)
        {
            // Threshold value specified in C. Biottos thesis.
            // Based on a 1D shock tube problem S_k = log10(1/p^4)
            // See  G.E. Barter and D.L. Darmofal. Shock Capturing with PDE-based
            // artificial diffusion for DGFEM: Part 1 Formulation, Journal od
            // Computational Physics 229 (2010) 1810-1827 for further reference
            
            //NekDouble S_Kappa = -4.0*log10(pOrderElmt[e]);
            //const double Kappa = 0.5;
            
            // Adjustable depending on the coarsness of the mesh. Might want to
            // move this variable into the session file
            
            int nQuadPointsElement = m_fields[0]->GetExp(e)->GetTotPoints();
            Array <OneD, NekDouble> one2D(nQuadPointsElement, 1.0);
            NekDouble Area = m_fields[0]->GetExp(e)->Integral(one2D);
            
            for (int n = 0; n < nQuadPointsElement; n++)
            {
                NekDouble mu_0 = m_mu0;
                
                if (Sensor[n+PointCount] < (m_Skappa-m_Kappa))
                {
                    mu_var[n+PointCount] = 0;
                }
                if(Sensor[n+PointCount] >= (m_Skappa-m_Kappa)
                        && Sensor[n+PointCount] <= (m_Skappa+m_Kappa))
                {
                    mu_var[n+PointCount] = mu_0*(0.5*(1+sin(
                    M_PI*(Sensor[n+PointCount]-m_Skappa-m_Kappa)/(2*m_Kappa))));
                }
                if(Sensor[n+PointCount] > (m_Skappa+m_Kappa))
                {
                    mu_var[n+PointCount] = mu_0;
                }
            }
            
            PointCount += nQuadPointsElement;
        }
    }
}
