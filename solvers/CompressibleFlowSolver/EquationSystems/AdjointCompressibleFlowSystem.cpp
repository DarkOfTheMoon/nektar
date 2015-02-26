///////////////////////////////////////////////////////////////////////////////
//
// File AdjointCompressibleFlowSystem.cpp
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
#include <CompressibleFlowSolver/EquationSystems/AdjointCompressibleFlowSystem.h>
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
    string AdjointCompressibleFlowSystem::className =
    SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
                    "AdjointCompressibleFlowSystem",
                    AdjointCompressibleFlowSystem::create,
                    "Auxiliary functions for the compressible flow system.");
    
    AdjointCompressibleFlowSystem::AdjointCompressibleFlowSystem(
                    const LibUtilities::SessionReaderSharedPtr& pSession)
    : UnsteadySystem(pSession)
    {
    }
    /**
     * @brief Initialization object for AdjointCompressibleFlowSystem class.
     */
    void AdjointCompressibleFlowSystem::v_InitObject()
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
        
        m_cnt = 0;
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
        m_session->LoadParameter ("GasConstant",    m_gasConstant,     287.058);
        m_session->LoadParameter ("Twall",          m_Twall,            300.15);
        m_session->LoadSolverInfo("ViscosityType",  m_ViscosityType,"Constant");
        m_session->LoadSolverInfo("Target",         m_Target,            "NoN");
        m_session->LoadParameter ("mu",             m_mu,             1.78e-05);
        m_session->LoadParameter ("thermalConductivity",
                                           m_thermalConductivity,       0.0257);
        m_session->LoadParameter ("rhoInfBase",   m_rhoInfBase,            0.0);
        m_session->LoadParameter ("uInfBase",     m_uInfBase,              1.0);
        m_session->LoadParameter ("vInfBase",     m_vInfBase,              1.0);
        m_session->LoadParameter ("pInfBase",     m_pInfBase,              1.0);
        m_session->LoadParameter ("alphaInfBase", m_alphaInfBase,          0.0);
        m_session->LoadParameter ("Lref",           m_Lref,                1.0);
        m_session->LoadParameter ("alpha",          m_alpha,               0.0);
        m_session->LoadParameter ("numCheck",       m_numCheck,            0.0);
        m_session->LoadParameter ("finalCheck",     m_finalCheck,          0.0);
        m_session->LoadParameter ("thermalConductivity",
                                                 m_thermalConductivity, 0.0257);
        
        m_EqTypeStr = m_session->GetSolverInfo("EQTYPE");

        m_Cp      = m_gamma / (m_gamma - 1.0) * m_gasConstant;
        m_Prandtl = m_Cp * m_mu / m_thermalConductivity;
        
        const int nPoints = m_fields[0]->GetTotPoints();
        const int nVar    = m_fields.num_elements();
        m_un = Array<OneD, Array<OneD, NekDouble> > (m_fields.num_elements());
        for (int i = 0; i < m_fields.num_elements(); ++i)
        {
            m_un[i] = Array<OneD, NekDouble> (nPoints, 0.0);
            Vmath::Vcopy(nPoints, m_fields[i]->GetPhys(), 1, m_un[i], 1);
        }
        
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
            
                // Setting up Riemann solver for advection operator
                m_session->LoadSolverInfo("UpwindType", riemName, "Average");
                
                m_riemannSolver = SolverUtils::GetRiemannSolverFactory()
                .CreateInstance(riemName);
                
                // Setting up upwind solver for diffusion operator
                m_riemannSolverLDG = SolverUtils::GetRiemannSolverFactory()
                .CreateInstance("UpwindLDG");
                
                // Setting up parameters for advection operator Riemann solver
                m_riemannSolver->SetParam (
                    "gamma",   &AdjointCompressibleFlowSystem::GetGamma,   this);
                m_riemannSolver->SetAuxVec(
                    "vecLocs", &AdjointCompressibleFlowSystem::GetVecLocs, this);
                m_riemannSolver->SetVector(
                    "N",       &AdjointCompressibleFlowSystem::GetNormals, this);
                
                // Setting up parameters for diffusion operator Riemann solver
                m_riemannSolverLDG->SetParam (
                    "gamma",   &AdjointCompressibleFlowSystem::GetGamma,   this);
                m_riemannSolverLDG->SetVector(
                    "vecLocs", &AdjointCompressibleFlowSystem::GetVecLocs, this);
                m_riemannSolverLDG->SetVector(
                    "N",       &AdjointCompressibleFlowSystem::GetNormals, this);
                
                m_riemannSolver->SetFwdBwdBaseFlow(
                        &AdjointCompressibleFlowSystem::GetFwdBwdBaseFlow, this);
                // Concluding initialisation of advection / diffusion operators
                m_advection->SetRiemannSolver   (m_riemannSolver);
                m_diffusion->SetRiemannSolver   (m_riemannSolverLDG);
                m_advection->InitObject         (m_session, m_fields);
                m_diffusion->InitObject         (m_session, m_fields);
                
                m_baseflow = Array<OneD, Array<OneD, NekDouble> > (nVar);
                for (int i = 0; i < nVar; ++i)
                {
                    m_baseflow[i] = Array<OneD, NekDouble> (nPoints, 0.0);
                }
                //
                string basefile = m_session->GetFunctionFilename("BaseFlow", 0);
                m_advection->ImportFldBase(basefile, m_fields);
                
                m_advection->GetBaseFlow(m_baseflow);
            
                cout << Vmath::Vmax(m_baseflow[0].num_elements(), m_baseflow[0], 1) << endl;
                cout << Vmath::Vmax(m_baseflow[1].num_elements(), m_baseflow[0], 1) << endl;
                cout << Vmath::Vmax(m_baseflow[2].num_elements(), m_baseflow[0], 1) << endl;
                cout << Vmath::Vmax(m_baseflow[3].num_elements(), m_baseflow[0], 1) << endl;
                
                m_dVdUdXi = Array<OneD, Array<OneD, Array<OneD,
                Array<OneD, NekDouble> > > > (m_spacedim);
                
                m_JacPrim = Array<OneD, Array<OneD, Array<OneD,
                Array<OneD, NekDouble> > > > (m_spacedim);
            
                m_JacDivPrim = Array<OneD, Array<OneD, Array<OneD,
                Array<OneD, NekDouble> > > > (m_spacedim);
                
                m_JacDivCons = Array<OneD, Array<OneD, Array<OneD,
                Array<OneD, NekDouble> > > > (m_spacedim);
                
                m_JacAddPrim = Array<OneD, Array<OneD, Array<OneD,
                Array<OneD, NekDouble> > > > (m_spacedim);
                m_JacAddDivPrim = Array<OneD, Array<OneD, Array<OneD,
                Array<OneD, NekDouble> > > > (m_spacedim);
                
                m_JacAddDivCons = Array<OneD, Array<OneD, Array<OneD,
                Array<OneD, NekDouble> > > > (m_spacedim);
                
                m_JacAddCons = Array<OneD, Array<OneD, Array<OneD,
                Array<OneD, NekDouble> > > > (m_spacedim);
                
                m_JacVisc = Array<OneD, Array<OneD, Array<OneD,
                Array<OneD, Array<OneD, NekDouble> > > > > (m_spacedim);
                
                for (int i = 0; i < m_spacedim; ++i)
                {
                    m_dVdUdXi[i] = Array<OneD, Array<OneD,
                    Array<OneD, NekDouble> > > (nVar);
                    
                    m_JacPrim[i] = Array<OneD, Array<OneD,
                    Array<OneD, NekDouble> > > (nVar);
                    
                    m_JacDivPrim[i] = Array<OneD, Array<OneD,
                    Array<OneD, NekDouble> > > (nVar);
                    
                    m_JacDivCons[i] = Array<OneD, Array<OneD,
                    Array<OneD, NekDouble> > > (nVar);
                    
                    for (int j = 0; j < nVar; ++j)
                    {
                        m_dVdUdXi[i][j] = Array<OneD, Array<OneD,
                        NekDouble > > (nVar);
                        
                        m_JacPrim[i][j] = Array<OneD, Array<OneD,
                        NekDouble > > (nVar);
                        
                        m_JacDivPrim[i][j] = Array<OneD, Array<OneD,
                        NekDouble > > (nVar);
                        
                        m_JacDivCons[i][j] = Array<OneD, Array<OneD,
                        NekDouble > > (nVar);
                        
                        for (int k = 0; k < nVar; k++)
                        {
                            m_dVdUdXi[i][j][k] = Array<OneD,
                            NekDouble> (nPoints, 0.0);
                            
                            m_JacPrim[i][j][k] = Array<OneD,
                            NekDouble> (nPoints, 0.0);
                            
                            m_JacDivPrim[i][j][k] = Array<OneD,
                            NekDouble> (nPoints, 0.0);
                            
                            m_JacDivCons[i][j][k] = Array<OneD,
                            NekDouble> (nPoints, 0.0);
                        }
                    }
                }
                
                m_dVdU = Array<OneD, Array<OneD,
                Array<OneD, NekDouble > > > (nVar);
                
                for (int i = 0; i < nVar; ++i)
                {
                    m_dVdU[i] = Array<OneD, Array<OneD,
                    NekDouble > > (nVar);
                    
                    for (int j = 0; j < nVar; j++)
                    {
                        m_dVdU[i][j] = Array<OneD,
                        NekDouble > (nPoints, 0.0);
                    }
                }
                if (m_EqTypeStr == "AdjointEulerCFE")
                {
                    // Setting up the matrices for the viscous jacobians.
                    // since dim{dFv/dU}
                    //       = [m_spacedim][nVar][nVar][nPoints]
                    
                    GetConservToPrimVariableInvMat(m_dVdU);
                    GetConservToPrimVariableInvMatDiv(m_dVdU, m_dVdUdXi);
                    // Get J^c_v = dFcdV
                    // Not working in 3D yet
                    GetJacobianConvFluxPrim(m_JacPrim);
                    // Get d(J^c_v)/dxi = d(dFcdV)/dxi
                    GetDerivJacobian(m_JacPrim, m_JacDivPrim);
                    
                    if (m_spacedim == 2)
                    {
                        // Allocate memory
                        Array<OneD, NekDouble> tmpConsX1(nPoints, 0.0);
                        Array<OneD, NekDouble> tmpConsX2(nPoints, 0.0);
                        Array<OneD, NekDouble> tmpConsY1(nPoints, 0.0);
                        Array<OneD, NekDouble> tmpConsY2(nPoints, 0.0);
        
                        for (int i = 0; i < nVar; ++i)
                        {
                            for (int j = 0; j < nVar; ++j)
                            {
                                // Allocate memory
                                Vmath::Zero(nPoints, &tmpConsX1[0], 1);
                                Vmath::Zero(nPoints, &tmpConsX2[0], 1);
                                Vmath::Zero(nPoints, &tmpConsY1[0], 1);
                                Vmath::Zero(nPoints, &tmpConsY2[0], 1);
                                
                                for (int k = 0; k < nVar; ++k)
                                {
                                    
                                    Vmath::Vvtvp(nPoints,
                                                 &m_dVdUdXi[0][i][k][0], 1,
                                                 &m_JacPrim[0][k][j][0], 1,
                                                 &tmpConsX1[0], 1,
                                                 &tmpConsX1[0], 1);
                                    
                                    Vmath::Vvtvp(nPoints,
                                                 &m_dVdUdXi[1][i][k][0], 1,
                                                 &m_JacPrim[1][k][j][0], 1,
                                                 &tmpConsX2[0], 1,
                                                 &tmpConsX2[0], 1);
                                    
                                    //=============
                                    
                                    Vmath::Vvtvp(nPoints,
                                                 &m_dVdU[i][k][0], 1,
                                                 &m_JacDivPrim[0][k][j][0], 1,
                                                 &tmpConsY1[0], 1,
                                                 &tmpConsY1[0], 1);
                                    
                                    Vmath::Vvtvp(nPoints,
                                                 &m_dVdU[i][k][0], 1,
                                                 &m_JacDivPrim[1][k][j][0], 1,
                                                 &tmpConsY2[0], 1,
                                                 &tmpConsY2[0], 1);
                                    
                                    Vmath::Vadd(nPoints,
                                                &tmpConsX1[0], 1,
                                                &tmpConsY1[0], 1,
                                                &m_JacDivCons[0][i][j][0], 1);
                                    
                                    Vmath::Vadd(nPoints,
                                                &tmpConsX2[0], 1,
                                                &tmpConsY2[0], 1,
                                                &m_JacDivCons[1][i][j][0], 1);
                                }
                            }
                        }
                    }
                    if (m_spacedim == 3)
                    {
                        
                        // Allocate memory
                        Array<OneD, NekDouble> tmpConsX1(nPoints, 0.0);
                        Array<OneD, NekDouble> tmpConsX2(nPoints, 0.0);
                        Array<OneD, NekDouble> tmpConsX3(nPoints, 0.0);
                        
                        Array<OneD, NekDouble> tmpConsY1(nPoints, 0.0);
                        Array<OneD, NekDouble> tmpConsY2(nPoints, 0.0);
                        Array<OneD, NekDouble> tmpConsY3(nPoints, 0.0);
                        
                        for (int i = 0; i < nVar; ++i)
                        {
                            for (int j = 0; j < nVar; ++j)
                            {
                                Vmath::Zero(nPoints, &tmpConsX1[0], 1);
                                Vmath::Zero(nPoints, &tmpConsX2[0], 1);
                                Vmath::Zero(nPoints, &tmpConsX3[0], 1);
                                
                                Vmath::Zero(nPoints, &tmpConsY1[0], 1);
                                Vmath::Zero(nPoints, &tmpConsY2[0], 1);
                                Vmath::Zero(nPoints, &tmpConsY3[0], 1);
                                
                                for (int k = 0; k < nVar; ++k)
                                {
                                    
                                    Vmath::Vvtvp(nPoints,
                                                 &m_dVdUdXi[0][i][k][0], 1,
                                                 &m_JacPrim[0][k][j][0], 1,
                                                 &tmpConsX1[0], 1,
                                                 &tmpConsX1[0], 1);
                                    
                                    Vmath::Vvtvp(nPoints,
                                                 &m_dVdUdXi[1][i][k][0], 1,
                                                 &m_JacPrim[1][k][j][0], 1,
                                                 &tmpConsX2[0], 1,
                                                 &tmpConsX2[0], 1);
                                    
                                    Vmath::Vvtvp(nPoints,
                                                 &m_dVdUdXi[2][i][k][0], 1,
                                                 &m_JacPrim[2][k][j][0], 1,
                                                 &tmpConsX3[0], 1,
                                                 &tmpConsX3[0], 1);
                                    
                                    //=============
                                    
                                    Vmath::Vvtvp(nPoints,
                                                 &m_dVdU[i][k][0], 1,
                                                 &m_JacDivPrim[0][k][j][0], 1,
                                                 &tmpConsY1[0], 1,
                                                 &tmpConsY1[0], 1);
                                    
                                    Vmath::Vvtvp(nPoints,
                                                 &m_dVdU[i][k][0], 1,
                                                 &m_JacDivPrim[1][k][j][0], 1,
                                                 &tmpConsY2[0], 1,
                                                 &tmpConsY2[0], 1);
                                    
                                    Vmath::Vvtvp(nPoints,
                                                 &m_dVdU[i][k][0], 1,
                                                 &m_JacDivPrim[2][k][j][0], 1,
                                                 &tmpConsY3[0], 1,
                                                 &tmpConsY3[0], 1);
                                    
                                    //=============
                                    
                                    
                                    Vmath::Vadd(nPoints,
                                                &tmpConsX1[0], 1,
                                                &tmpConsY1[0], 1,
                                                &m_JacDivCons[0][i][j][0], 1);
                                    
                                    
                                    Vmath::Vadd(nPoints,
                                                &tmpConsX2[0], 1,
                                                &tmpConsY2[0], 1,
                                                &m_JacDivCons[1][i][j][0], 1);
                                    
                                    
                                    Vmath::Vadd(nPoints,
                                                &tmpConsX3[0], 1,
                                                &tmpConsY3[0], 1,
                                                &m_JacDivCons[2][i][j][0], 1);
                                }
                            }
                        }
                    }
                    m_advection->SetAdjointFluxVector(
                    &AdjointCompressibleFlowSystem::GetAdjointFluxVector, this);
                    
                    m_advection->SetJacTransposeDivVector(
                    &AdjointCompressibleFlowSystem::GetAdjointDerivJacVector, this);
                }
                if (m_EqTypeStr == "AdjointNavierStokesCFE")
                {
                    // Setting up the matrices for the viscous jacobians.
                    // since dim{dFv/dU}
                    //       = [m_spacedim][nVar][nVar][nPoints]
                    for (int i = 0; i < m_spacedim; ++i)
                    {
                        m_JacAddPrim[i] = Array<OneD, Array<OneD,
                        Array<OneD, NekDouble> > > (nVar);
                        
                        m_JacAddDivPrim[i] = Array<OneD, Array<OneD,
                        Array<OneD, NekDouble> > > (nVar);
                        
                        m_JacAddDivCons[i] = Array<OneD, Array<OneD,
                        Array<OneD, NekDouble> > > (nVar);
                        
                        m_JacAddCons[i] = Array<OneD, Array<OneD,
                        Array<OneD, NekDouble> > > (nVar);
                        
                        for (int j = 0; j < nVar; ++j)
                        {
                            m_JacAddPrim[i][j] = Array<OneD, Array<OneD,
                            NekDouble > > (nVar);
                            
                            m_JacAddDivPrim[i][j] = Array<OneD, Array<OneD,
                            NekDouble > > (nVar);
                            
                            m_JacAddDivCons[i][j] = Array<OneD, Array<OneD,
                            NekDouble > > (nVar);
                            
                            m_JacAddCons[i][j] = Array<OneD, Array<OneD,
                            NekDouble > > (nVar);
                            
                            for (int k = 0; k < nVar; k++)
                            {
                                m_JacAddPrim[i][j][k] = Array<OneD,
                                NekDouble> (nPoints, 0.0);
                                
                                m_JacAddDivPrim[i][j][k] = Array<OneD,
                                NekDouble> (nPoints, 0.0);
                                
                                m_JacAddDivCons[i][j][k] = Array<OneD,
                                NekDouble> (nPoints, 0.0);
                                
                                m_JacAddCons[i][j][k] = Array<OneD,
                                NekDouble> (nPoints, 0.0);
                            }
                        }
                    }
                    
                    m_dVdU = Array<OneD, Array<OneD,
                    Array<OneD, NekDouble > > > (nVar);
                    
                    for (int i = 0; i < nVar; ++i)
                    {
                        m_dVdU[i] = Array<OneD, Array<OneD,
                        NekDouble > > (nVar);
                        
                        for (int j = 0; j < nVar; j++)
                        {
                            m_dVdU[i][j] = Array<OneD,
                            NekDouble > (nPoints, 0.0);
                        }
                    }
                    
                    // Setting up the matrices for the viscous jacobians.
                    // since dim{dFv/d(dUdxi)}
                    //       = [m_spacedim][m_spacedim][nVar][nVar][nPoints]
                    for (int i = 0; i < m_spacedim; ++i)
                    {
                        m_JacVisc[i] = Array<OneD, Array<OneD,
                        Array<OneD, Array<OneD, NekDouble> > > > (m_spacedim);
                        
                        for (int j = 0; j < m_spacedim; ++j)
                        {
                            m_JacVisc[i][j] = Array<OneD, Array<OneD,
                            Array<OneD, NekDouble> > > (nVar);
                            
                            for (int k = 0; k < nVar; ++k)
                            {
                                m_JacVisc[i][j][k] = Array<OneD, Array<OneD,
                                NekDouble> > (nVar);
                                
                                for (int l = 0; l < nVar; ++l)
                                {
                                    m_JacVisc[i][j][k][l] = Array<OneD, NekDouble> (nPoints, 0.0);
                                }
                            }
                        }
                    }
                    
                    // Get Change of Variable matrices
                    GetConservToPrimVariableInvMat(m_dVdU);
                    GetConservToPrimVariableInvMatDiv(m_dVdU, m_dVdUdXi);
                    
                    // Get J^c_v = dFcdV
                    // Not working in 3D yet
                    GetJacobianConvFluxPrim(m_JacPrim);
                    // Get d(J^c_v)/dxi = d(dFcdV)/dxi
                    GetDerivJacobian(m_JacPrim, m_JacDivPrim);
                    // Get J^v_v =dFvdV
                    // Not working in 3D yet
                    GetJacobianAddConvFlux(m_JacAddPrim);
                    // Get d(J^v_v)/dxi = d(dFv/dV)/dxi
                    GetDerivJacobian(m_JacAddPrim, m_JacAddDivPrim);
                    // Get the viscous jacobians
                    //GetJacobianViscousFluxPrim(m_JacVisc);
                    // Obtain the jacobians in conservative form
                    
                    if (m_spacedim == 2)
                    {
                        // Allocate memory
                        Array<OneD, NekDouble> tmpConsX1(nPoints, 0.0);
                        Array<OneD, NekDouble> tmpConsX2(nPoints, 0.0);
                        Array<OneD, NekDouble> tmpConsY1(nPoints, 0.0);
                        Array<OneD, NekDouble> tmpConsY2(nPoints, 0.0);
                        
                        Array<OneD, NekDouble> tmpAddConsX1(nPoints, 0.0);
                        Array<OneD, NekDouble> tmpAddConsX2(nPoints, 0.0);
                        Array<OneD, NekDouble> tmpAddConsY1(nPoints, 0.0);
                        Array<OneD, NekDouble> tmpAddConsY2(nPoints, 0.0);
                        
                        for (int i = 0; i < nVar; ++i)
                        {
                            for (int j = 0; j < nVar; ++j)
                            {
                                // Allocate memory
                                Vmath::Zero(nPoints, &tmpConsX1[0], 1);
                                Vmath::Zero(nPoints, &tmpConsX2[0], 1);
                                Vmath::Zero(nPoints, &tmpConsY1[0], 1);
                                Vmath::Zero(nPoints, &tmpConsY2[0], 1);
                                
                                Vmath::Zero(nPoints, &tmpAddConsX1[0], 1);
                                Vmath::Zero(nPoints, &tmpAddConsX2[0], 1);
                                Vmath::Zero(nPoints, &tmpAddConsY1[0], 1);
                                Vmath::Zero(nPoints, &tmpAddConsY2[0], 1);
                                
                                for (int k = 0; k < nVar; ++k)
                                {
                                    
                                    Vmath::Vvtvp(nPoints,
                                                 &m_dVdUdXi[0][i][k][0], 1,
                                                 &m_JacPrim[0][k][j][0], 1,
                                                 &tmpConsX1[0], 1,
                                                 &tmpConsX1[0], 1);
                                    
                                    Vmath::Vvtvp(nPoints,
                                                 &m_dVdUdXi[1][i][k][0], 1,
                                                 &m_JacPrim[1][k][j][0], 1,
                                                 &tmpConsX2[0], 1,
                                                 &tmpConsX2[0], 1);
                                    
                                    //=============
                                    
                                    Vmath::Vvtvp(nPoints,
                                                 &m_dVdU[i][k][0], 1,
                                                 &m_JacDivPrim[0][k][j][0], 1,
                                                 &tmpConsY1[0], 1,
                                                 &tmpConsY1[0], 1);
                                    
                                    Vmath::Vvtvp(nPoints,
                                                 &m_dVdU[i][k][0], 1,
                                                 &m_JacDivPrim[1][k][j][0], 1,
                                                 &tmpConsY2[0], 1,
                                                 &tmpConsY2[0], 1);
                                    
                                    Vmath::Vadd(nPoints,
                                                &tmpConsX1[0], 1,
                                                &tmpConsY1[0], 1,
                                                &m_JacDivCons[0][i][j][0], 1);
                                    
                                    Vmath::Vadd(nPoints,
                                                &tmpConsX2[0], 1,
                                                &tmpConsY2[0], 1,
                                                &m_JacDivCons[1][i][j][0], 1);
                                    
                                    //======================================
                                    
                                    //
                                    Vmath::Vvtvp(nPoints,
                                                 &m_dVdUdXi[0][i][k][0], 1,
                                                 &m_JacAddPrim[0][k][j][0], 1,
                                                 &tmpAddConsX1[0], 1,
                                                 &tmpAddConsX1[0], 1);
                                    
                                    Vmath::Vvtvp(nPoints,
                                                 &m_dVdUdXi[1][i][k][0], 1,
                                                 &m_JacAddPrim[1][k][j][0], 1,
                                                 &tmpAddConsX2[0], 1,
                                                 &tmpAddConsX2[0], 1);
                                    
                                    //=============
                                    
                                    Vmath::Vvtvp(nPoints,
                                                 &m_dVdU[i][k][0], 1,
                                                 &m_JacAddDivPrim[0][k][j][0], 1,
                                                 &tmpAddConsY1[0], 1,
                                                 &tmpAddConsY1[0], 1);
                                    
                                    Vmath::Vvtvp(nPoints,
                                                 &m_dVdU[i][k][0], 1,
                                                 &m_JacAddDivPrim[1][k][j][0], 1,
                                                 &tmpAddConsY2[0], 1,
                                                 &tmpAddConsY2[0], 1);
                                    
                                    Vmath::Vadd(nPoints,
                                                &tmpAddConsX1[0], 1,
                                                &tmpAddConsY1[0], 1,
                                                &m_JacAddDivCons[0][i][j][0], 1);
                                    
                                    Vmath::Vadd(nPoints,
                                                &tmpAddConsX2[0], 1,
                                                &tmpAddConsY2[0], 1,
                                                &m_JacAddDivCons[1][i][j][0], 1);
                                }
                            }
                        }
                        for (int i = 0; i < nVar; ++i)
                        {
                            for (int j = 0; j < nVar; ++j)
                            {
                                for (int k = 0; k < nVar; ++k)
                                {
                                    Vmath::Vvtvp(nPoints,
                                                 &m_dVdU[i][k][0], 1,
                                                 &m_JacAddPrim[0][k][j][0], 1,
                                                 &m_JacAddCons[0][i][j][0], 1,
                                                 &m_JacAddCons[0][i][j][0], 1);
                                    
                                    Vmath::Vvtvp(nPoints,
                                                 &m_dVdU[i][k][0], 1,
                                                 &m_JacAddPrim[1][k][j][0], 1,
                                                 &m_JacAddCons[1][i][j][0], 1,
                                                 &m_JacAddCons[1][i][j][0], 1);
                                }
                            }
                        }
                    }
                    if (m_spacedim == 3)
                    {
                        
                        // Allocate memory
                        Array<OneD, NekDouble> tmpConsX1(nPoints, 0.0);
                        Array<OneD, NekDouble> tmpConsX2(nPoints, 0.0);
                        Array<OneD, NekDouble> tmpConsX3(nPoints, 0.0);
                        
                        Array<OneD, NekDouble> tmpConsY1(nPoints, 0.0);
                        Array<OneD, NekDouble> tmpConsY2(nPoints, 0.0);
                        Array<OneD, NekDouble> tmpConsY3(nPoints, 0.0);
                        
                        Array<OneD, NekDouble> tmpConsZ1(nPoints, 0.0);
                        Array<OneD, NekDouble> tmpConsZ2(nPoints, 0.0);
                        Array<OneD, NekDouble> tmpConsZ3(nPoints, 0.0);
                        
                        Array<OneD, NekDouble> tmpAddConsX1(nPoints, 0.0);
                        Array<OneD, NekDouble> tmpAddConsX2(nPoints, 0.0);
                        Array<OneD, NekDouble> tmpAddConsX3(nPoints, 0.0);
                        
                        Array<OneD, NekDouble> tmpAddConsY1(nPoints, 0.0);
                        Array<OneD, NekDouble> tmpAddConsY2(nPoints, 0.0);
                        Array<OneD, NekDouble> tmpAddConsY3(nPoints, 0.0);
                        
                        Array<OneD, NekDouble> tmpAddConsZ1(nPoints, 0.0);
                        Array<OneD, NekDouble> tmpAddConsZ2(nPoints, 0.0);
                        Array<OneD, NekDouble> tmpAddConsZ3(nPoints, 0.0);
                        
                        for (int i = 0; i < nVar; ++i)
                        {
                            for (int j = 0; j < nVar; ++j)
                            {
                                
                                // Allocate memory
                                Vmath::Zero(nPoints, &tmpConsX1[0], 1);
                                Vmath::Zero(nPoints, &tmpConsX2[0], 1);
                                Vmath::Zero(nPoints, &tmpConsX3[0], 1);
                                
                                Vmath::Zero(nPoints, &tmpConsY1[0], 1);
                                Vmath::Zero(nPoints, &tmpConsY2[0], 1);
                                Vmath::Zero(nPoints, &tmpConsY3[0], 1);
                                
                                Vmath::Zero(nPoints, &tmpAddConsX1[0], 1);
                                Vmath::Zero(nPoints, &tmpAddConsX2[0], 1);
                                Vmath::Zero(nPoints, &tmpAddConsX3[0], 1);
                                
                                Vmath::Zero(nPoints, &tmpAddConsY1[0], 1);
                                Vmath::Zero(nPoints, &tmpAddConsY2[0], 1);
                                Vmath::Zero(nPoints, &tmpAddConsY3[0], 1);
                                
                                for (int k = 0; k < nVar; ++k)
                                {
                                    
                                    Vmath::Vvtvp(nPoints,
                                                 &m_dVdUdXi[0][i][k][0], 1,
                                                 &m_JacPrim[0][k][j][0], 1,
                                                 &tmpConsX1[0], 1,
                                                 &tmpConsX1[0], 1);
                                    
                                    Vmath::Vvtvp(nPoints,
                                                 &m_dVdUdXi[1][i][k][0], 1,
                                                 &m_JacPrim[1][k][j][0], 1,
                                                 &tmpConsX2[0], 1,
                                                 &tmpConsX2[0], 1);
                                    
                                    Vmath::Vvtvp(nPoints,
                                                 &m_dVdUdXi[2][i][k][0], 1,
                                                 &m_JacPrim[2][k][j][0], 1,
                                                 &tmpConsX3[0], 1,
                                                 &tmpConsX3[0], 1);
                                    
                                    //=============
                                    
                                    Vmath::Vvtvp(nPoints,
                                                 &m_dVdU[i][k][0], 1,
                                                 &m_JacDivPrim[0][k][j][0], 1,
                                                 &tmpConsY1[0], 1,
                                                 &tmpConsY1[0], 1);
                                    
                                    Vmath::Vvtvp(nPoints,
                                                 &m_dVdU[i][k][0], 1,
                                                 &m_JacDivPrim[1][k][j][0], 1,
                                                 &tmpConsY2[0], 1,
                                                 &tmpConsY2[0], 1);
                                    
                                    Vmath::Vvtvp(nPoints,
                                                 &m_dVdU[i][k][0], 1,
                                                 &m_JacDivPrim[2][k][j][0], 1,
                                                 &tmpConsY3[0], 1,
                                                 &tmpConsY3[0], 1);
                                    
                                    //=============
                                    
                                    
                                    Vmath::Vadd(nPoints,
                                                &tmpConsX1[0], 1,
                                                &tmpConsY1[0], 1,
                                                &m_JacDivCons[0][i][j][0], 1);
                                    
                                    
                                    Vmath::Vadd(nPoints,
                                                &tmpConsX2[0], 1,
                                                &tmpConsY2[0], 1,
                                                &m_JacDivCons[1][i][j][0], 1);
                                    
                                    
                                    Vmath::Vadd(nPoints,
                                                &tmpConsX3[0], 1,
                                                &tmpConsY3[0], 1,
                                                &m_JacDivCons[2][i][j][0], 1);
                                    
                                    //======================================
                                    
                                    //
                                    Vmath::Vvtvp(nPoints,
                                                 &m_dVdUdXi[0][i][k][0], 1,
                                                 &m_JacAddPrim[0][k][j][0], 1,
                                                 &tmpAddConsX1[0], 1,
                                                 &tmpAddConsX1[0], 1);
                                    
                                    Vmath::Vvtvp(nPoints,
                                                 &m_dVdUdXi[1][i][k][0], 1,
                                                 &m_JacAddPrim[1][k][j][0], 1,
                                                 &tmpAddConsX2[0], 1,
                                                 &tmpAddConsX2[0], 1);
                                    
                                    Vmath::Vvtvp(nPoints,
                                                 &m_dVdUdXi[1][i][k][0], 1,
                                                 &m_JacAddPrim[1][k][j][0], 1,
                                                 &tmpAddConsX3[0], 1,
                                                 &tmpAddConsX3[0], 1);
                                    
                                    //=============
                                    
                                    Vmath::Vvtvp(nPoints,
                                                 &m_dVdU[i][k][0], 1,
                                                 &m_JacAddDivPrim[0][k][j][0], 1,
                                                 &tmpAddConsY1[0], 1,
                                                 &tmpAddConsY1[0], 1);
                                    
                                    Vmath::Vvtvp(nPoints,
                                                 &m_dVdU[i][k][0], 1,
                                                 &m_JacAddDivPrim[1][k][j][0], 1,
                                                 &tmpAddConsY2[0], 1,
                                                 &tmpAddConsY2[0], 1);
                                    
                                    Vmath::Vvtvp(nPoints,
                                                 &m_dVdU[i][k][0], 1,
                                                 &m_JacAddDivPrim[2][k][j][0], 1,
                                                 &tmpAddConsY3[0], 1,
                                                 &tmpAddConsY3[0], 1);
                                    
                                    Vmath::Vadd(nPoints,
                                                &tmpAddConsX1[0], 1,
                                                &tmpAddConsY1[0], 1,
                                                &m_JacAddDivCons[0][i][j][0], 1);
                                    
                                    Vmath::Vadd(nPoints,
                                                &tmpAddConsX2[0], 1,
                                                &tmpAddConsY2[0], 1,
                                                &m_JacAddDivCons[1][i][j][0], 1);
                                    
                                    
                                    Vmath::Vadd(nPoints,
                                                &tmpAddConsX3[0], 1,
                                                &tmpAddConsY3[0], 1,
                                                &m_JacAddDivCons[2][i][j][0], 1);
                                }
                            }
                        }
                        for (int i = 0; i < nVar; ++i)
                        {
                            for (int j = 0; j < nVar; ++j)
                            {
                                for (int k = 0; k < nVar; ++k)
                                {
                                    Vmath::Vvtvp(nPoints,
                                                 &m_dVdU[i][k][0], 1,
                                                 &m_JacAddPrim[0][k][j][0], 1,
                                                 &m_JacAddCons[0][i][j][0], 1,
                                                 &m_JacAddCons[0][i][j][0], 1);
                                    
                                    Vmath::Vvtvp(nPoints,
                                                 &m_dVdU[i][k][0], 1,
                                                 &m_JacAddPrim[1][k][j][0], 1,
                                                 &m_JacAddCons[1][i][j][0], 1,
                                                 &m_JacAddCons[1][i][j][0], 1);
                                    
                                    Vmath::Vvtvp(nPoints,
                                                 &m_dVdU[i][k][0], 1,
                                                 &m_JacAddPrim[2][k][j][0], 1,
                                                 &m_JacAddCons[2][i][j][0], 1,
                                                 &m_JacAddCons[2][i][j][0], 1);
                                }
                            }
                        }
                    }
                    
                    m_advection->SetAdjointFluxVector(
                    &AdjointCompressibleFlowSystem::GetAdjointFluxVector, this);
                    
                    m_advection->SetJacTransposeDivVector(
                    &AdjointCompressibleFlowSystem::GetAdjointDerivJacVector, this);
                    
                    m_advection->SetAddFluxVector(
                    &AdjointCompressibleFlowSystem::GetAdjointAddConvFluxVector, this);
                    
                    m_advection->SetAddJacTransposeDivVector(
                    &AdjointCompressibleFlowSystem::GetAdjointDerivAddJacVector, this);
                    
                    m_diffusion->SetFluxVectorNS(&AdjointCompressibleFlowSystem::GetAdjointViscousFluxVector, this);
                }
                
                break;
            }
            default:
            {
                ASSERTL0(false, "Unsupported projection type.");
                break;
            }
        }
        // Type of advection class to be used
    }
    
    /**
     * @brief Destructor for AdjointCompressibleFlowSystem class.
     */
    AdjointCompressibleFlowSystem::~AdjointCompressibleFlowSystem()
    {
        
    }
    
    void AdjointCompressibleFlowSystem::v_DoInitialise(void)
    {
        const int nPoints = m_fields[0]->GetTotPoints();
        const int nVar    = m_fields.num_elements();
        //
        if (m_numCheck != 0)
        {
            string basefile = m_session->GetFunctionFilename("BaseFlow", 0);
            std::string chkout = std::to_string(m_finalCheck-m_cnt);
            std::string basename = basefile+"_"+chkout+".chk";
            
            cout << "Imprt Initialise" << endl;
            m_advection->ImportFldBase(basename, m_fields);
            
            if (m_session->GetComm()->GetRank() == 0)
            {
                if (m_session->GetComm()->GetRank() == 0)
                {
                    cout << endl;
                    cout << "=================================================="<< endl;
                    cout << "=================== Base flow ====================" << endl;
                    cout << "=================================================="<< endl;
                    cout << endl;
                }
                for (int i = 0; i < m_fields.num_elements(); ++i)
                {
                    std::string varName = m_session->GetVariable(i);
                    cout << "  - Field " << varName << ": "
                    << "from file " << basefile+"_"+chkout+".chk"
                    << endl;
                }
            }
            if (m_cnt != 0)
            {
                CheckForRestart(m_time);
            }
            
            SetBoundaryConditions(m_time);
            SetInitialConditions(m_time);
            
            m_cnt++;
        }
        else
        {
            UnsteadySystem::v_DoInitialise();
        }
    }
    
    void AdjointCompressibleFlowSystem::CheckForRestart(NekDouble &time)
    {
        if (m_session->DefinesFunction("InitialConditions"))
        {
            for (int i = 0; i < m_fields.num_elements(); ++i)
            {
                LibUtilities::FunctionType vType;
                
                vType = m_session->GetFunctionType(
                                "InitialConditions", m_session->GetVariable(i));
                
                if (vType == LibUtilities::eFunctionTypeFile)
                {
                    std::string filename
                    = m_session->GetFunctionFilename(
                                "InitialConditions", m_session->GetVariable(i));
                    
                    std::string chkout = std::to_string(m_cnt);
                    
                    if (m_cnt != 0)
                    {
                        filename = filename+"_"+chkout+".chk";
                    }
                    
                    fs::path pfilename(filename);
                    
                    // redefine path for parallel file which is in directory
                    if(fs::is_directory(pfilename))
                    {
                        fs::path metafile("Info.xml");
                        fs::path fullpath = pfilename / metafile;
                        filename = LibUtilities::PortablePath(fullpath);
                    }
                    m_fld->ImportFieldMetaData(filename, m_fieldMetaDataMap);
                    
                    // check to see if time defined
                    if (m_fieldMetaDataMap !=
                        LibUtilities::NullFieldMetaDataMap)
                    {
                        LibUtilities::FieldMetaDataMap::iterator iter;
                        
                        iter = m_fieldMetaDataMap.find("Time");
                        if (iter != m_fieldMetaDataMap.end())
                        {
                            time = boost::lexical_cast<NekDouble>(iter->second);
                        }
                    }
                    
                    break;
                }
            }
        }
    }
    
    void AdjointCompressibleFlowSystem::SymmetryBC(
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
    
    NekDouble AdjointCompressibleFlowSystem::v_GetTimeStep(
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
    
    void AdjointCompressibleFlowSystem::GetAdjointFluxVector(
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
        
        // load the physical values to obtain the adjoint vector \phi
        
        // Flux vector for the rho equation
        for (i = 0; i < m_spacedim; ++i)
        {
            vel[i] = Array<OneD, NekDouble>(nq);
        }
        
        GetVelocityVector(m_baseflow, vel);
        GetPressure      (m_baseflow, vel, pressure);
        
        Array<OneD, NekDouble> ones(nq, 1.0);
        
        // determine u^2+v^2+w^2;
        for (i = 0; i < m_spacedim; ++i)
        {
            Vmath::Vvtvp(nq, vel[i], 1, vel[i], 1, velsq, 1, velsq, 1);
        }
        
        // determine H = E + p/rho;
        Vmath::Vdiv(nq, m_baseflow[nvar-1], 1, m_baseflow[0], 1, H, 1);
        
        Vmath::Vdiv(nq, pressure, 1, m_baseflow[0], 1, Htmp, 1);
        
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
    
    void AdjointCompressibleFlowSystem::GetAdjointDerivJacVector(
          const Array<OneD, Array<OneD, NekDouble> > &inarray,
                Array<OneD, Array<OneD, Array<OneD, NekDouble > > > &outarray)
    {
        int i,j;
        
        int nq = inarray[0].num_elements();
        int nvar = m_fields.num_elements();
        
        if (m_spacedim == 2)
        {
            Array<OneD, Array<OneD, NekDouble> > tmpx(nvar);
            Array<OneD, Array<OneD, NekDouble> > tmpy(nvar);
            
            for (i = 0; i < nvar; ++i)
            {
                tmpx[i] = Array<OneD, NekDouble>(nq, 0.0);
                tmpy[i] = Array<OneD, NekDouble>(nq, 0.0);
                
                for (j = 0; j < nvar; ++j)
                {
                    Vmath::Vvtvp(nq,
                                 &m_JacDivCons[0][i][j][0], 1,
                                 &inarray[j][0], 1,
                                 &tmpx[i][0], 1,
                                 &tmpx[i][0], 1);
                    
                    Vmath::Vvtvp(nq,
                                 &m_JacDivCons[1][i][j][0], 1,
                                 &inarray[j][0], 1,
                                 &tmpy[i][0], 1,
                                 &tmpy[i][0], 1);
                    
                }
                
                Vmath::Vcopy(nq, &tmpx[i][0], 1, &outarray[i][0][0], 1);
                Vmath::Vcopy(nq, &tmpy[i][0], 1, &outarray[i][1][0], 1);
            }
        }
        if (m_spacedim == 3)
        {
            Array<OneD, Array<OneD, NekDouble> > tmpx(nvar);
            Array<OneD, Array<OneD, NekDouble> > tmpy(nvar);
            Array<OneD, Array<OneD, NekDouble> > tmpz(nvar);
            
            for (i = 0; i < nvar; ++i)
            {
                tmpx[i] = Array<OneD, NekDouble>(nq, 0.0);
                tmpy[i] = Array<OneD, NekDouble>(nq, 0.0);
                tmpz[i] = Array<OneD, NekDouble>(nq, 0.0);
                
                for (j = 0; j < nvar; ++j)
                {
                    Vmath::Vvtvp(nq,
                                 &m_JacDivCons[0][i][j][0], 1,
                                 &inarray[j][0], 1,
                                 &tmpx[i][0], 1,
                                 &tmpx[i][0], 1);
                    
                    Vmath::Vvtvp(nq,
                                 &m_JacDivCons[1][i][j][0], 1,
                                 &inarray[j][0], 1,
                                 &tmpy[i][0], 1,
                                 &tmpy[i][0], 1);
                    
                    Vmath::Vvtvp(nq,
                                 &m_JacDivCons[2][i][j][0], 1,
                                 &inarray[j][0], 1,
                                 &tmpz[i][0], 1,
                                 &tmpz[i][0], 1);
                }
                
                Vmath::Vcopy(nq, &tmpx[i][0], 1, &outarray[i][0][0], 1);
                Vmath::Vcopy(nq, &tmpy[i][0], 1, &outarray[i][1][0], 1);
                Vmath::Vcopy(nq, &tmpz[i][0], 1, &outarray[i][2][0], 1);
            }
        }
    }
    
    void AdjointCompressibleFlowSystem::GetAdjointAddConvFluxVector(
              const Array<OneD, Array<OneD, NekDouble> > &inarray,
                    Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &outarray)
    {
        int i,j;
        
        int nq = inarray[0].num_elements();
        int nvar = m_fields.num_elements();
        
        if (m_spacedim == 2)
        {
            Array<OneD, Array<OneD, NekDouble> > tmpx(nvar);
            Array<OneD, Array<OneD, NekDouble> > tmpy(nvar);
            
            for (i = 0; i < nvar; ++i)
            {
                tmpx[i] = Array<OneD, NekDouble>(nq, 0.0);
                tmpy[i] = Array<OneD, NekDouble>(nq, 0.0);
                
                for (j = 0; j < nvar; ++j)
                {
                    
                    Vmath::Vvtvp(nq,
                                 &m_JacAddCons[0][i][j][0], 1,
                                 &inarray[j][0], 1,
                                 &tmpx[i][0], 1,
                                 &tmpx[i][0], 1);
                    
                    Vmath::Vvtvp(nq,
                                 &m_JacAddCons[1][i][j][0], 1,
                                 &inarray[j][0], 1,
                                 &tmpy[i][0], 1,
                                 &tmpy[i][0], 1);
                }
                
                Vmath::Vcopy(nq, &tmpx[i][0], 1, &outarray[i][0][0], 1);
                Vmath::Vcopy(nq, &tmpy[i][0], 1, &outarray[i][1][0], 1);
                
            }
        }
        if (m_spacedim == 3)
        {
            Array<OneD, Array<OneD, NekDouble> > tmpx(nvar);
            Array<OneD, Array<OneD, NekDouble> > tmpy(nvar);
            Array<OneD, Array<OneD, NekDouble> > tmpz(nvar);
            
            for (i = 0; i < nvar; ++i)
            {
                tmpx[i] = Array<OneD, NekDouble>(nq, 0.0);
                tmpy[i] = Array<OneD, NekDouble>(nq, 0.0);
                tmpz[i] = Array<OneD, NekDouble>(nq, 0.0);
                
                for (j = 0; j < nvar; ++j)
                {
                    Vmath::Vvtvp(nq,
                                 &m_JacAddCons[0][i][j][0], 1,
                                 &inarray[j][0], 1,
                                 &tmpx[i][0], 1,
                                 &tmpx[i][0], 1);
                    
                    Vmath::Vvtvp(nq,
                                 &m_JacAddCons[1][i][j][0], 1,
                                 &inarray[j][0], 1,
                                 &tmpy[i][0], 1,
                                 &tmpy[i][0], 1);
                    
                    Vmath::Vvtvp(nq,
                                 &m_JacAddCons[2][i][j][0], 1,
                                 &inarray[j][0], 1,
                                 &tmpz[i][0], 1,
                                 &tmpz[i][0], 1);
                }
                
                Vmath::Vcopy(nq, &tmpx[i][0], 1, &outarray[i][0][0], 1);
                Vmath::Vcopy(nq, &tmpy[i][0], 1, &outarray[i][1][0], 1);
                Vmath::Vcopy(nq, &tmpz[i][0], 1, &outarray[i][2][0], 1);
            }
        }
    }
    
    void AdjointCompressibleFlowSystem::GetAdjointDerivAddJacVector(
            const Array<OneD, Array<OneD, NekDouble> > &inarray,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &outarray)
    {
        int i,j;
        
        int nq = inarray[0].num_elements();
        int nvar = m_fields.num_elements();
        
        if (m_spacedim == 2)
        {
            Array<OneD, Array<OneD, NekDouble> > tmpx(nvar);
            Array<OneD, Array<OneD, NekDouble> > tmpy(nvar);
            
            for (i = 0; i < nvar; ++i)
            {
                tmpx[i] = Array<OneD, NekDouble>(nq, 0.0);
                tmpy[i] = Array<OneD, NekDouble>(nq, 0.0);
                
                for (j = 0; j < nvar; ++j)
                {
                    Vmath::Vvtvp(nq,
                                 &m_JacAddDivCons[0][i][j][0], 1,
                                 &inarray[j][0], 1,
                                 &tmpx[i][0], 1,
                                 &tmpx[i][0], 1);
                    
                    Vmath::Vvtvp(nq,
                                 &m_JacAddDivCons[1][i][j][0], 1,
                                 &inarray[j][0], 1,
                                 &tmpy[i][0], 1,
                                 &tmpy[i][0], 1);
                    
                }
                
                Vmath::Vcopy(nq, &tmpx[i][0], 1, &outarray[i][0][0], 1);
                Vmath::Vcopy(nq, &tmpy[i][0], 1, &outarray[i][1][0], 1);
            }
        }
        if (m_spacedim == 3)
        {
            Array<OneD, Array<OneD, NekDouble> > tmpx(nvar);
            Array<OneD, Array<OneD, NekDouble> > tmpy(nvar);
            Array<OneD, Array<OneD, NekDouble> > tmpz(nvar);
            
            for (i = 0; i < nvar; ++i)
            {
                tmpx[i] = Array<OneD, NekDouble>(nq, 0.0);
                tmpy[i] = Array<OneD, NekDouble>(nq, 0.0);
                tmpz[i] = Array<OneD, NekDouble>(nq, 0.0);
                
                for (j = 0; j < nvar; ++j)
                {
                    Vmath::Vvtvp(nq,
                                 &m_JacAddDivCons[0][i][j][0], 1,
                                 &inarray[j][0], 1,
                                 &tmpx[i][0], 1,
                                 &tmpx[i][0], 1);
                    
                    Vmath::Vvtvp(nq,
                                 &m_JacAddDivCons[1][i][j][0], 1,
                                 &inarray[j][0], 1,
                                 &tmpy[i][0], 1,
                                 &tmpy[i][0], 1);
                    
                    Vmath::Vvtvp(nq,
                                 &m_JacAddDivCons[2][i][j][0], 1,
                                 &inarray[j][0], 1,
                                 &tmpz[i][0], 1,
                                 &tmpz[i][0], 1);
                }
                
                Vmath::Vcopy(nq, &tmpx[i][0], 1, &outarray[i][0][0], 1);
                Vmath::Vcopy(nq, &tmpy[i][0], 1, &outarray[i][1][0], 1);
                Vmath::Vcopy(nq, &tmpz[i][0], 1, &outarray[i][2][0], 1);
            }
        }
    }
    
    void AdjointCompressibleFlowSystem::GetAdjointViscousFluxVector(
      const Array<OneD, Array<OneD, NekDouble> > &inarray,
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &z_derivatives,
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &outarray)
    {
        
        // Jac[m_spacedim][m_spacedim][nvar][nvar][nq];
        int i,j,k;
        int nq = inarray[0].num_elements();
        int nvar = m_fields.num_elements();
        
        if (m_spacedim == 2)
        {
            Array<OneD, Array<OneD, NekDouble> > tmpx1(nvar);
            Array<OneD, Array<OneD, NekDouble> > tmpx2(nvar);
            
            Array<OneD, Array<OneD, NekDouble> > tmpy1(nvar);
            Array<OneD, Array<OneD, NekDouble> > tmpy2(nvar);
            
            for (i = 0; i < nvar; ++i)
            {
                tmpx1[i] = Array<OneD, NekDouble>(nq, 0.0);
                tmpy1[i] = Array<OneD, NekDouble>(nq, 0.0);
                
                tmpx2[i] = Array<OneD, NekDouble>(nq, 0.0);
                tmpy2[i] = Array<OneD, NekDouble>(nq, 0.0);
                
                for (j = 0; j < nvar; ++j)
                {
                    Vmath::Vvtvp(nq,
                                 &m_JacVisc[0][0][i][j][0], 1,
                                 &z_derivatives[0][j][0], 1,
                                 &tmpx1[i][0], 1,
                                 &tmpx1[i][0], 1);
                    
                    Vmath::Vvtvp(nq,
                                 &m_JacVisc[0][1][i][j][0], 1,
                                 &z_derivatives[1][j][0], 1,
                                 &tmpx2[i][0], 1,
                                 &tmpx2[i][0], 1);
                    
                    Vmath::Vvtvp(nq,
                                 &m_JacVisc[1][0][i][j][0], 1,
                                 &z_derivatives[0][j][0], 1,
                                 &tmpy1[i][0], 1,
                                 &tmpy1[i][0], 1);
                    
                    Vmath::Vvtvp(nq,
                                 &m_JacVisc[1][1][i][j][0], 1,
                                 &z_derivatives[1][j][0], 1,
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
            
            // Correct?
            for (i = 0; i < nvar; ++i)
            {
                for (j = 0; j < nvar; ++j)
                {
                    Vmath::Vvtvp(nq,
                                 &m_dVdU[i][j][0], 1,
                                 &tmpnewx[j][0], 1,
                                 &outarray[0][i][0], 1,
                                 &outarray[0][i][0], 1);
                    
                    Vmath::Vvtvp(nq,
                                 &m_dVdU[i][j][0], 1,
                                 &tmpnewy[j][0], 1,
                                 &outarray[1][i][0], 1,
                                 &outarray[1][i][0], 1);
                }
            }
        }
        if (m_spacedim == 3)
        {
            ASSERTL0(false, "3D LAMINAR adjoint solver is not yet implemented");
        }
    }
    
    void AdjointCompressibleFlowSystem::GetJacobianViscousFluxPrim(
                Array<OneD, Array<OneD, Array<OneD, Array<OneD,
                Array<OneD, NekDouble> > > > >        JacVisc)
    {
        
        // Jac[m_spacedim][m_spacedim][nvar][nvar][nq];
        int i;
        int nq = m_baseflow[0].num_elements();
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
        
        GetVelocityVector(m_baseflow, vel);
        
        for (i = 0; i < m_spacedim; ++i)
        {
            Vmath::Smul(nq, m_mu, &vel[i][0], 1, &MuVel[i][0], 1);
        }
        
        GetPressure(m_baseflow, vel, Pressure);
        
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
        
        Vmath::Vmul(nq, &m_baseflow[0][0], 1, &m_baseflow[0][0], 1, &tmp[0], 1);
        Vmath::Vdiv(nq, &kPrho2Vec[0], 1, &tmp[0], 1, &kPrho2Vec[0], 1);
        
        Vmath::Zero(nq, &tmp[0], 1);
        
        if (m_spacedim == 2)
        {
            /* =============================================================
             JacVisc[0][0] =
             
             [ 0,        0,  0,    -gamma*mu/(Pr*(gamma-1))*p/rho^2]
             [ 0,   4/3*mu,  0,                            4/3*u*mu]
             [ 0,        0, mu,                                v*mu]
             [ 0,        0,  0,         gamma/(Pr*(gamma-1))*mu/rho]
             
             =============================================================*/
            // row 1;
            Vmath::Zero(nq, &JacVisc[0][0][0][0][0], 1);
            Vmath::Zero(nq, &JacVisc[0][0][0][1][0], 1);
            Vmath::Zero(nq, &JacVisc[0][0][0][2][0], 1);
            
            Vmath::Vcopy(nq,
                         &kPrho2Vec[0], 1,
                         &JacVisc[0][0][0][3][0], 1);
            
            Vmath::Neg(nq, &JacVisc[0][0][0][3][0], 1);
            
            // row 2;
            Vmath::Zero(nq, &JacVisc[0][0][1][0][0], 1);
            Vmath::Vcopy(nq,
                         &FourOThreeMuVector[0], 1,
                         &JacVisc[0][0][1][1][0], 1);
            
            Vmath::Zero(nq, &JacVisc[0][0][1][2][0], 1);
            Vmath::Smul(nq,
                        FourOThreeMuSclr,
                        &vel[0][0], 1,
                        &JacVisc[0][0][1][3][0], 1);
            
            // row 3;
            Vmath::Zero(nq, &JacVisc[0][0][2][0][0], 1);
            Vmath::Zero(nq, &JacVisc[0][0][2][1][0], 1);
            Vmath::Vcopy(nq, &muVec[0], 1, &JacVisc[0][0][2][2][0], 1);
            Vmath::Vcopy(nq, &MuVel[1][0], 1, &JacVisc[0][0][2][3][0], 1);
            
            // row 4;
            Vmath::Zero(nq, &JacVisc[0][0][3][0][0], 1);
            Vmath::Zero(nq, &JacVisc[0][0][3][1][0], 1);
            Vmath::Zero(nq, &JacVisc[0][0][3][2][0], 1);
            
            Vmath::Vdiv(nq,
                        &ones[0], 1,
                        &m_baseflow[0][0], 1,
                        &JacVisc[0][0][3][3][0], 1);
            
            Vmath::Vmul(nq,
                        &GamMu_O_GamMinOnePrVector[0], 1,
                        &JacVisc[0][0][3][3][0], 1,
                        &JacVisc[0][0][3][3][0], 1);
            
            /* =============================================================
             JacVisc[0][1] =
             
             [ 0,        0,       0,            0]
             [ 0,        0, -2/3*mu,    -2/3*v*mu]
             [ 0,       mu,       0,         u*mu]
             [ 0,        0,       0,            0]
             
             =============================================================*/
            // row 1;
            Vmath::Zero(nq, &JacVisc[0][1][0][0][0], 1);
            Vmath::Zero(nq, &JacVisc[0][1][0][1][0], 1);
            Vmath::Zero(nq, &JacVisc[0][1][0][2][0], 1);
            Vmath::Zero(nq, &JacVisc[0][1][0][3][0], 1);
            
            // row 2;
            Vmath::Zero(nq, &JacVisc[0][1][1][0][0], 1);
            Vmath::Zero(nq, &JacVisc[0][1][1][1][0], 1);
            Vmath::Vcopy(nq, &TwoOThreeMuVector[0], 1, &JacVisc[0][1][1][2][0], 1);
            Vmath::Neg(nq, &JacVisc[0][1][1][2][0], 1);
            
            Vmath::Vmul(nq,
                        &TwoOThreeMuVector[0], 1,
                        &vel[1][0], 1,
                        &JacVisc[0][1][1][3][0], 1);
            
            Vmath::Neg(nq, &JacVisc[0][1][1][3][0], 1);
            
            // row 3;
            Vmath::Zero(nq, &JacVisc[0][1][2][0][0], 1);
            Vmath::Vcopy(nq, &muVec[0], 1, &JacVisc[0][1][2][1][0], 1);
            Vmath::Zero(nq, &JacVisc[0][1][2][2][0], 1);
            Vmath::Vcopy(nq, &MuVel[0][0], 1, &JacVisc[0][1][2][3][0], 1);
            
            // row 4;
            Vmath::Zero(nq, &JacVisc[0][1][3][0][0], 1);
            Vmath::Zero(nq, &JacVisc[0][1][3][1][0], 1);
            Vmath::Zero(nq, &JacVisc[0][1][3][2][0], 1);
            Vmath::Zero(nq, &JacVisc[0][1][3][3][0], 1);
            /* =============================================================
             JacVisc[1][0] =
             
             [ 0,        0,       0,            0]
             [ 0,        0,      mu,         v*mu]
             [ 0,       -2/3*mu,  0,    -2/3*u*mu]
             [ 0,        0,       0,            0]
             
             =============================================================*/
            // row 1;
            Vmath::Zero(nq, &JacVisc[1][0][0][0][0], 1);
            Vmath::Zero(nq, &JacVisc[1][0][0][1][0], 1);
            Vmath::Zero(nq, &JacVisc[1][0][0][2][0], 1);
            Vmath::Zero(nq, &JacVisc[1][0][0][3][0], 1);
            
            // row 2;
            Vmath::Zero(nq, &JacVisc[1][0][1][0][0], 1);
            Vmath::Zero(nq, &JacVisc[1][0][1][1][0], 1);
            Vmath::Vcopy(nq, &muVec[0], 1, &JacVisc[1][0][1][2][0], 1);
            Vmath::Vcopy(nq, &MuVel[1][0], 1, &JacVisc[1][0][1][3][0], 1);
            
            // row 3;
            Vmath::Zero(nq, &JacVisc[1][0][2][0][0], 1);
            Vmath::Vcopy(nq, &TwoOThreeMuVector[0], 1, &JacVisc[1][0][2][1][0], 1);
            Vmath::Neg(nq, &JacVisc[1][0][2][1][0], 1);
            Vmath::Zero(nq, &JacVisc[1][0][2][2][0], 1);
            Vmath::Vmul(nq,
                        &TwoOThreeMuVector[0], 1,
                        &vel[1][0], 1,
                        &JacVisc[1][0][2][3][0], 1);
            
            Vmath::Neg(nq, &JacVisc[1][0][2][3][0], 1);
            
            // row 4;
            Vmath::Zero(nq, &JacVisc[1][0][3][0][0], 1);
            Vmath::Zero(nq, &JacVisc[1][0][3][1][0], 1);
            Vmath::Zero(nq, &JacVisc[1][0][3][2][0], 1);
            Vmath::Zero(nq, &JacVisc[1][0][3][3][0], 1);
            
            /* =============================================================
             JacVisc[1][1] =
             
             [ 0,        0,      0,    -gamma*mu/(Pr*(gamma-1))*p/rho^2]
             [ 0,       mu,      0,                                u*mu]
             [ 0,        0, 4/3*mu,                            4/3*v*mu]
             [ 0,        0,      0,         gamma/(Pr*(gamma-1))*mu/rho]
             
             ===============================================================*/
            
            // row 1;
            Vmath::Zero(nq, &JacVisc[1][1][0][0][0], 1);
            Vmath::Zero(nq, &JacVisc[1][1][0][1][0], 1);
            Vmath::Zero(nq, &JacVisc[1][1][0][2][0], 1);
            Vmath::Vcopy(nq, &kPrho2Vec[0], 1, &JacVisc[1][1][0][3][0], 1);
            Vmath::Neg(nq, &JacVisc[1][1][0][3][0], 1);
            
            // row 2;
            Vmath::Zero(nq, &JacVisc[1][1][1][0][0], 1);
            Vmath::Vcopy(nq, &muVec[0], 1, &JacVisc[1][1][1][1][0], 1);
            Vmath::Zero(nq, &JacVisc[1][1][1][2][0], 1);
            Vmath::Vcopy(nq, &MuVel[0][0], 1, &JacVisc[1][1][1][3][0], 1);
            
            // row 3;
            Vmath::Zero(nq, &JacVisc[1][1][2][0][0], 1);
            Vmath::Zero(nq, &JacVisc[1][1][2][1][0], 1);
            Vmath::Vcopy(nq,
                         &FourOThreeMuVector[0], 1,
                         &JacVisc[1][1][2][2][0], 1);
            
            Vmath::Vmul(nq,
                        &FourOThreeMuVector[0], 1,
                        &vel[1][0], 1,
                        &JacVisc[1][1][2][3][0], 1);
            
            // row 4;
            Vmath::Zero(nq, &JacVisc[1][1][3][0][0], 1);
            Vmath::Zero(nq, &JacVisc[1][1][3][1][0], 1);
            Vmath::Zero(nq, &JacVisc[1][1][3][2][0], 1);
            
            Vmath::Vdiv(nq,
                        &ones[0], 1,
                        &m_baseflow[0][0], 1,
                        &JacVisc[1][1][3][3][0], 1);
            
            Vmath::Vmul(nq,
                        &GamMu_O_GamMinOnePrVector[0], 1,
                        &JacVisc[1][1][3][3][0], 1,
                        &JacVisc[1][1][3][3][0], 1);
        }
        
        if (m_spacedim == 3)
        {
            /* =============================================================
             JacVisc[0][0] =
             
             [ 0,        0,  0,   0,     -gamma*mu/(Pr*(gamma-1))*p/rho^2]
             [ 0,   4/3*mu,  0,   0,                             4/3*u*mu]
             [ 0,        0, mu,   0,                                 v*mu]
             [ 0,        0,  0,   mu,                                w*mu]
             [ 0,        0,  0,   0,          gamma/(Pr*(gamma-1))*mu/rho]
             
             =============================================================*/
            // row 1;
            Vmath::Zero(nq, &JacVisc[0][0][0][0][0], 1);
            Vmath::Zero(nq, &JacVisc[0][0][0][1][0], 1);
            Vmath::Zero(nq, &JacVisc[0][0][0][2][0], 1);
            Vmath::Zero(nq, &JacVisc[0][0][0][3][0], 1);
            
            Vmath::Vcopy(nq,
                         &kPrho2Vec[0], 1,
                         &JacVisc[0][0][0][4][0], 1);
            
            Vmath::Neg(nq, &JacVisc[0][0][0][4][0], 1);
            
            // row 2;
            Vmath::Zero(nq, &JacVisc[0][0][1][0][0], 1);
            Vmath::Vcopy(nq,
                         &FourOThreeMuVector[0], 1,
                         &JacVisc[0][0][1][1][0], 1);
            
            Vmath::Zero(nq, &JacVisc[0][0][1][2][0], 1);
            Vmath::Zero(nq, &JacVisc[0][0][1][3][0], 1);
            Vmath::Smul(nq,
                        FourOThreeMuSclr,
                        &vel[0][0], 1,
                        &JacVisc[0][0][1][4][0], 1);
            
            // row 3;
            Vmath::Zero(nq, &JacVisc[0][0][2][0][0], 1);
            Vmath::Zero(nq, &JacVisc[0][0][2][1][0], 1);
            Vmath::Vcopy(nq, &muVec[0], 1, &JacVisc[0][0][2][2][0], 1);
            Vmath::Zero(nq, &JacVisc[0][0][2][3][0], 1);
            Vmath::Vcopy(nq, &MuVel[1][0], 1, &JacVisc[0][0][2][4][0], 1);
            
            // row 4;
            Vmath::Zero(nq, &JacVisc[0][0][2][0][0], 1);
            Vmath::Zero(nq, &JacVisc[0][0][2][1][0], 1);
            Vmath::Zero(nq, &JacVisc[0][0][2][2][0], 1);
            Vmath::Vcopy(nq, &muVec[0], 1, &JacVisc[0][0][2][3][0], 1);
            Vmath::Vcopy(nq, &MuVel[2][0], 1, &JacVisc[0][0][2][4][0], 1);
            
            // row 5;
            Vmath::Zero(nq, &JacVisc[0][0][3][0][0], 1);
            Vmath::Zero(nq, &JacVisc[0][0][3][1][0], 1);
            Vmath::Zero(nq, &JacVisc[0][0][3][2][0], 1);
            Vmath::Zero(nq, &JacVisc[0][0][3][3][0], 1);
            
            Vmath::Vdiv(nq,
                        &ones[0], 1,
                        &m_baseflow[0][0], 1,
                        &JacVisc[0][0][3][4][0], 1);
            
            Vmath::Vmul(nq,
                        &GamMu_O_GamMinOnePrVector[0], 1,
                        &JacVisc[0][0][3][4][0], 1,
                        &JacVisc[0][0][3][4][0], 1);
            
            /* =============================================================
             JacVisc[0][1] =
             
             [ 0,        0,       0,    0,          0]
             [ 0,        0, -2/3*mu,    0,  -2/3*v*mu]
             [ 0,       mu,       0,    0,       u*mu]
             [ 0,        0,       0,    0,          0]
             [ 0,        0,       0,    0,          0]
             
             =============================================================*/
            // row 1;
            Vmath::Zero(nq, &JacVisc[0][1][0][0][0], 1);
            Vmath::Zero(nq, &JacVisc[0][1][0][1][0], 1);
            Vmath::Zero(nq, &JacVisc[0][1][0][2][0], 1);
            Vmath::Zero(nq, &JacVisc[0][1][0][3][0], 1);
            Vmath::Zero(nq, &JacVisc[0][1][0][4][0], 1);
            
            // row 2;
            Vmath::Zero(nq, &JacVisc[0][1][1][0][0], 1);
            Vmath::Zero(nq, &JacVisc[0][1][1][1][0], 1);
            Vmath::Vcopy(nq, &TwoOThreeMuVector[0], 1, &JacVisc[0][1][1][2][0], 1);
            Vmath::Neg(nq, &JacVisc[0][1][1][2][0], 1);
            Vmath::Zero(nq, &JacVisc[0][1][1][3][0], 1);
            Vmath::Vmul(nq,
                        &TwoOThreeMuVector[0], 1,
                        &vel[1][0], 1,
                        &JacVisc[0][1][1][4][0], 1);
            
            Vmath::Neg(nq, &JacVisc[0][1][1][4][0], 1);
            
            // row 3;
            Vmath::Zero(nq, &JacVisc[0][1][2][0][0], 1);
            Vmath::Vcopy(nq, &muVec[0], 1, &JacVisc[0][1][2][1][0], 1);
            Vmath::Zero(nq, &JacVisc[0][1][2][2][0], 1);
            Vmath::Zero(nq, &JacVisc[0][1][2][3][0], 1);
            Vmath::Vcopy(nq, &MuVel[0][0], 1, &JacVisc[0][1][2][4][0], 1);
            
            // row 4;
            Vmath::Zero(nq, &JacVisc[0][1][3][0][0], 1);
            Vmath::Zero(nq, &JacVisc[0][1][3][1][0], 1);
            Vmath::Zero(nq, &JacVisc[0][1][3][2][0], 1);
            Vmath::Zero(nq, &JacVisc[0][1][3][3][0], 1);
            Vmath::Zero(nq, &JacVisc[0][1][3][4][0], 1);
            
            
            /* =============================================================
             JacVisc[0][2] =
             
             [ 0,        0,       0,    0,             0]
             [ 0,        0,       0, -2/3*mu,  -2/3*w*mu]
             [ 0,        0,       0,    0,             0]
             [ 0,       mu,       0,    0,          u*mu]
             [ 0,        0,       0,    0,             0]
             
             =============================================================*/
            // row 1;
            Vmath::Zero(nq, &JacVisc[0][2][0][0][0], 1);
            Vmath::Zero(nq, &JacVisc[0][2][0][1][0], 1);
            Vmath::Zero(nq, &JacVisc[0][2][0][2][0], 1);
            Vmath::Zero(nq, &JacVisc[0][2][0][3][0], 1);
            Vmath::Zero(nq, &JacVisc[0][2][0][4][0], 1);
            
            // row 2;
            Vmath::Zero(nq, &JacVisc[0][2][1][0][0], 1);
            Vmath::Zero(nq, &JacVisc[0][2][1][1][0], 1);
            Vmath::Zero(nq, &JacVisc[0][2][1][2][0], 1);
            Vmath::Vcopy(nq,
                         &TwoOThreeMuVector[0], 1,
                         &JacVisc[0][2][1][3][0], 1);
            
            Vmath::Neg(nq, &JacVisc[0][2][1][3][0], 1);
            
            Vmath::Vmul(nq,
                        &TwoOThreeMuVector[0], 1,
                        &vel[1][0], 1,
                        &JacVisc[0][1][2][4][0], 1);
            
            Vmath::Neg(nq, &JacVisc[0][2][1][4][0], 1);
            
            // row 3;
            Vmath::Zero(nq, &JacVisc[0][2][2][0][0], 1);
            Vmath::Zero(nq, &JacVisc[0][2][2][1][0], 1);
            Vmath::Zero(nq, &JacVisc[0][2][2][2][0], 1);
            Vmath::Zero(nq, &JacVisc[0][2][2][3][0], 1);
            Vmath::Zero(nq, &JacVisc[0][2][2][4][0], 1);
            
            // row 4;
            Vmath::Zero(nq, &JacVisc[0][2][2][0][0], 1);
            Vmath::Vcopy(nq, &muVec[0], 1, &JacVisc[0][2][2][1][0], 1);
            Vmath::Zero(nq, &JacVisc[0][2][2][2][0], 1);
            Vmath::Zero(nq, &JacVisc[0][2][2][3][0], 1);
            Vmath::Vcopy(nq, &MuVel[0][0], 1, &JacVisc[0][2][2][4][0], 1);
            
            // row 5;
            Vmath::Zero(nq, &JacVisc[0][2][3][0][0], 1);
            Vmath::Zero(nq, &JacVisc[0][2][3][1][0], 1);
            Vmath::Zero(nq, &JacVisc[0][2][3][2][0], 1);
            Vmath::Zero(nq, &JacVisc[0][2][3][3][0], 1);
            Vmath::Zero(nq, &JacVisc[0][2][3][4][0], 1);
            /* =============================================================
             Avxy = JacVisc[1][0] =
             
             [ 0,        0,       0,    0,                    0]
             [ 0,        0,      mu,    0,                 v*mu]
             [ 0,       -2/3*mu,  0,    0,            -2/3*u*mu]
             [ 0,        0,       0,    0,                    0]
             [ 0,        0,       0,    0,                    0]
             
             =============================================================*/
            // row 1;
            Vmath::Zero(nq, &JacVisc[1][0][0][0][0], 1);
            Vmath::Zero(nq, &JacVisc[1][0][0][1][0], 1);
            Vmath::Zero(nq, &JacVisc[1][0][0][2][0], 1);
            Vmath::Zero(nq, &JacVisc[1][0][0][3][0], 1);
            Vmath::Zero(nq, &JacVisc[1][0][0][4][0], 1);
            
            // row 2;
            Vmath::Zero(nq, &JacVisc[1][0][1][0][0], 1);
            Vmath::Zero(nq, &JacVisc[1][0][1][1][0], 1);
            Vmath::Vcopy(nq, &muVec[0], 1, &JacVisc[1][0][1][2][0], 1);
            Vmath::Zero(nq, &JacVisc[1][0][1][3][0], 1);
            Vmath::Vcopy(nq, &MuVel[1][0], 1, &JacVisc[1][0][1][4][0], 1);
            
            // row 3;
            Vmath::Zero(nq, &JacVisc[1][0][2][0][0], 1);
            Vmath::Vcopy(nq, &TwoOThreeMuVector[0], 1, &JacVisc[1][0][2][1][0], 1);
            Vmath::Neg(nq, &JacVisc[1][0][2][1][0], 1);
            Vmath::Zero(nq, &JacVisc[1][0][2][2][0], 1);
            Vmath::Zero(nq, &JacVisc[1][0][2][3][0], 1);
            Vmath::Vmul(nq,
                        &TwoOThreeMuVector[0], 1,
                        &vel[1][0], 1,
                        &JacVisc[1][0][2][4][0], 1);
            
            Vmath::Neg(nq, &JacVisc[1][0][2][4][0], 1);
            
            // row 4;
            Vmath::Zero(nq, &JacVisc[1][0][3][0][0], 1);
            Vmath::Zero(nq, &JacVisc[1][0][3][1][0], 1);
            Vmath::Zero(nq, &JacVisc[1][0][3][2][0], 1);
            Vmath::Zero(nq, &JacVisc[1][0][3][3][0], 1);
            Vmath::Zero(nq, &JacVisc[1][0][3][4][0], 1);
            
            // row 5;
            Vmath::Zero(nq, &JacVisc[1][0][3][0][0], 1);
            Vmath::Zero(nq, &JacVisc[1][0][3][1][0], 1);
            Vmath::Zero(nq, &JacVisc[1][0][3][2][0], 1);
            Vmath::Zero(nq, &JacVisc[1][0][3][3][0], 1);
            Vmath::Zero(nq, &JacVisc[1][0][3][4][0], 1);
            
            /* =============================================================
             Avyy = JacVisc[1][1] =
             
             [ 0,        0,      0,    0,     -gamma*mu/(Pr*(gamma-1))*p/rho^2]
             [ 0,       mu,      0,    0,                                 u*mu]
             [ 0,        0, 4/3*mu,    0,                             4/3*v*mu]
             [ 0,        0,      0,   mu,                                 w*mu]
             [ 0,        0,      0,    0,          gamma/(Pr*(gamma-1))*mu/rho]
             
             ===============================================================*/
            
            // row 1;
            Vmath::Zero(nq, &JacVisc[1][1][0][0][0], 1);
            Vmath::Zero(nq, &JacVisc[1][1][0][1][0], 1);
            Vmath::Zero(nq, &JacVisc[1][1][0][2][0], 1);
            Vmath::Zero(nq, &JacVisc[1][1][0][3][0], 1);
            Vmath::Vcopy(nq, &kPrho2Vec[0], 1, &JacVisc[1][1][0][4][0], 1);
            Vmath::Neg(nq, &JacVisc[1][1][0][4][0], 1);
            
            // row 2;
            Vmath::Zero(nq, &JacVisc[1][1][1][0][0], 1);
            Vmath::Vcopy(nq, &muVec[0], 1, &JacVisc[1][1][1][1][0], 1);
            Vmath::Zero(nq, &JacVisc[1][1][1][2][0], 1);
            Vmath::Zero(nq, &JacVisc[1][1][1][3][0], 1);
            Vmath::Vcopy(nq, &MuVel[0][0], 1, &JacVisc[1][1][1][4][0], 1);
            
            // row 3;
            Vmath::Zero(nq, &JacVisc[1][1][2][0][0], 1);
            Vmath::Zero(nq, &JacVisc[1][1][2][1][0], 1);
            Vmath::Vcopy(nq,
                         &FourOThreeMuVector[0], 1,
                         &JacVisc[1][1][2][2][0], 1);
            
            Vmath::Zero(nq, &JacVisc[1][1][2][3][0], 1);
            
            Vmath::Vmul(nq,
                        &FourOThreeMuVector[0], 1,
                        &vel[1][0], 1,
                        &JacVisc[1][1][2][4][0], 1);
            
            // row 4;
            Vmath::Zero(nq, &JacVisc[1][1][3][0][0], 1);
            Vmath::Zero(nq, &JacVisc[1][1][3][1][0], 1);
            Vmath::Zero(nq, &JacVisc[1][1][3][2][0], 1);
            Vmath::Vcopy(nq, &muVec[0], 1, &JacVisc[1][1][3][3][0], 1);
            Vmath::Vcopy(nq, &MuVel[2][0], 1, &JacVisc[1][1][3][4][0], 1);
            
            // row 5;
            Vmath::Zero(nq, &JacVisc[1][1][3][0][0], 1);
            Vmath::Zero(nq, &JacVisc[1][1][3][1][0], 1);
            Vmath::Zero(nq, &JacVisc[1][1][3][2][0], 1);
            Vmath::Zero(nq, &JacVisc[1][1][3][3][0], 1);
            
            Vmath::Vdiv(nq,
                        &ones[0], 1,
                        &m_baseflow[0][0], 1,
                        &JacVisc[1][1][3][4][0], 1);
            
            Vmath::Vmul(nq,
                        &GamMu_O_GamMinOnePrVector[0], 1,
                        &JacVisc[1][1][3][4][0], 1,
                        &JacVisc[1][1][3][4][0], 1);
            
            
            /* =============================================================
             Avzy = JacVisc[1][2] =
             [ 0,        0,       0,          0,                    0]
             [ 0,        0,       0,          0,                    0]
             [ 0,        0,       0,    -2/3*mu,            -2/3*w*mu]
             [ 0,        0,      mu,          0,                 v*mu]
             [ 0,        0,       0,          0,                    0]
             =============================================================*/
            
            // row 1;
            Vmath::Zero(nq, &JacVisc[1][2][0][0][0], 1);
            Vmath::Zero(nq, &JacVisc[1][2][0][1][0], 1);
            Vmath::Zero(nq, &JacVisc[1][2][0][2][0], 1);
            Vmath::Zero(nq, &JacVisc[1][2][0][3][0], 1);
            Vmath::Zero(nq, &JacVisc[1][2][0][4][0], 1);
            
            // row 2;
            Vmath::Zero(nq, &JacVisc[1][2][1][0][0], 1);
            Vmath::Zero(nq, &JacVisc[1][2][1][1][0], 1);
            Vmath::Zero(nq, &JacVisc[1][2][1][2][0], 1);
            Vmath::Zero(nq, &JacVisc[1][2][1][3][0], 1);
            Vmath::Zero(nq, &JacVisc[1][2][1][4][0], 1);
            
            // row 3;
            Vmath::Zero(nq, &JacVisc[1][2][2][0][0], 1);
            Vmath::Zero(nq, &JacVisc[1][2][2][1][0], 1);
            Vmath::Zero(nq, &JacVisc[1][2][2][2][0], 1);
            Vmath::Vcopy(nq,
                         &TwoOThreeMuVector[0], 1,
                         &JacVisc[1][2][2][3][0], 1);
            Vmath::Neg(nq, &JacVisc[1][2][2][3][0], 1);
            
            Vmath::Vmul(nq,
                        &TwoOThreeMuVector[0], 1,
                        &vel[2][0], 1,
                        &JacVisc[1][2][2][4][0], 1);
            
            Vmath::Neg(nq, &JacVisc[1][2][2][4][0], 1);
            
            // row 4;
            Vmath::Zero(nq, &JacVisc[1][2][3][0][0], 1);
            Vmath::Zero(nq, &JacVisc[1][2][3][1][0], 1);
            Vmath::Vcopy(nq, &muVec[0], 1, &JacVisc[1][2][3][2][0], 1);
            Vmath::Zero(nq, &JacVisc[1][2][2][3][0], 1);
            Vmath::Vcopy(nq, &MuVel[1][0], 1, &JacVisc[1][2][3][4][0], 1);
            
            // row 5;
            Vmath::Zero(nq, &JacVisc[1][2][4][0][0], 1);
            Vmath::Zero(nq, &JacVisc[1][2][4][1][0], 1);
            Vmath::Zero(nq, &JacVisc[1][2][4][2][0], 1);
            Vmath::Zero(nq, &JacVisc[1][2][4][3][0], 1);
            Vmath::Zero(nq, &JacVisc[1][2][4][4][0], 1);
            
            
            /* =============================================================
             Avxz = JacVisc[2][0] =
             [ 0,        0,       0,          0,                        0]
             [ 0,        0,       0,         mu,                     mu*w]
             [ 0,        0,       0,          0,                        0]
             [ 0,  -2/3*mu,       0,          0,                -2/3*u*mu]
             [ 0,        0,       0,          0,                        0]
             =============================================================*/
            
            // row 1;
            Vmath::Zero(nq, &JacVisc[2][0][0][0][0], 1);
            Vmath::Zero(nq, &JacVisc[2][0][0][1][0], 1);
            Vmath::Zero(nq, &JacVisc[2][0][0][2][0], 1);
            Vmath::Zero(nq, &JacVisc[2][0][0][3][0], 1);
            Vmath::Zero(nq, &JacVisc[2][0][0][4][0], 1);
            
            // row 2;
            Vmath::Zero(nq, &JacVisc[2][0][1][0][0], 1);
            Vmath::Zero(nq, &JacVisc[2][0][1][1][0], 1);
            Vmath::Zero(nq, &JacVisc[2][0][1][2][0], 1);
            Vmath::Vcopy(nq, &muVec[0], 1, &JacVisc[2][0][1][3][0], 1);
            Vmath::Vcopy(nq, &MuVel[2][0], 1, &JacVisc[2][0][1][4][0], 1);
            
            // row 3;
            Vmath::Zero(nq, &JacVisc[2][0][2][0][0], 1);
            Vmath::Zero(nq, &JacVisc[2][0][2][1][0], 1);
            Vmath::Zero(nq, &JacVisc[2][0][2][2][0], 1);
            Vmath::Zero(nq, &JacVisc[2][0][2][3][0], 1);
            Vmath::Zero(nq, &JacVisc[2][0][2][4][0], 1);
            
            // row 4;
            Vmath::Zero(nq, &JacVisc[2][0][2][0][0], 1);
            Vmath::Vcopy(nq, &TwoOThreeMuVector[0], 1, &JacVisc[2][0][2][1][0], 1);
            Vmath::Neg(nq, &JacVisc[2][0][2][1][0], 1);
            Vmath::Zero(nq, &JacVisc[2][0][2][2][0], 1);
            Vmath::Zero(nq, &JacVisc[2][0][2][3][0], 1);
            
            Vmath::Vmul(nq,
                        &TwoOThreeMuVector[0], 1,
                        &vel[0][0], 1,
                        &JacVisc[2][0][0][4][0], 1);
            
            Vmath::Neg(nq, &JacVisc[2][0][2][4][0], 1);
            
            // row 5;
            Vmath::Zero(nq, &JacVisc[2][0][4][0][0], 1);
            Vmath::Zero(nq, &JacVisc[2][0][4][1][0], 1);
            Vmath::Zero(nq, &JacVisc[2][0][4][2][0], 1);
            Vmath::Zero(nq, &JacVisc[2][0][4][3][0], 1);
            Vmath::Zero(nq, &JacVisc[2][0][4][4][0], 1);
            
            
            /* =============================================================
             Avyz = JacVisc[2][1] =
             [ 0,        0,       0,          0,                        0]
             [ 0,        0,       0,          0,                        0]
             [ 0,        0,       0,         mu,                     mu*w]
             [ 0,        0, -2/3*mu,          0,                -2/3*u*mu]
             [ 0,        0,       0,          0,                        0]
             =============================================================*/
            
            // row 1;
            Vmath::Zero(nq, &JacVisc[2][1][0][0][0], 1);
            Vmath::Zero(nq, &JacVisc[2][1][0][1][0], 1);
            Vmath::Zero(nq, &JacVisc[2][1][0][2][0], 1);
            Vmath::Zero(nq, &JacVisc[2][1][0][3][0], 1);
            Vmath::Zero(nq, &JacVisc[2][1][0][4][0], 1);
            
            // row 2;
            Vmath::Zero(nq, &JacVisc[2][1][1][0][0], 1);
            Vmath::Zero(nq, &JacVisc[2][1][1][1][0], 1);
            Vmath::Zero(nq, &JacVisc[2][1][1][2][0], 1);
            Vmath::Zero(nq, &JacVisc[2][1][1][3][0], 1);
            Vmath::Zero(nq, &JacVisc[2][1][1][4][0], 1);
            
            // row 3;
            Vmath::Zero(nq, &JacVisc[2][1][2][0][0], 1);
            Vmath::Zero(nq, &JacVisc[2][1][2][1][0], 1);
            Vmath::Zero(nq, &JacVisc[2][1][2][2][0], 1);
            Vmath::Vcopy(nq, &muVec[0], 1, &JacVisc[2][1][2][3][0], 1);
            Vmath::Vcopy(nq, &MuVel[2][0], 1, &JacVisc[2][1][2][4][0], 1);
            
            // row 4;
            Vmath::Zero(nq, &JacVisc[2][1][3][0][0], 1);
            Vmath::Zero(nq, &JacVisc[2][1][3][1][0], 1);
            
            Vmath::Vcopy(nq,
                         &TwoOThreeMuVector[0], 1,
                         &JacVisc[2][1][2][1][0], 1);
            Vmath::Neg(nq, &JacVisc[2][1][2][1][0], 1);
            
            Vmath::Zero(nq, &JacVisc[2][1][2][2][0], 1);
            Vmath::Zero(nq, &JacVisc[2][1][2][3][0], 1);
            
            Vmath::Vmul(nq,
                        &TwoOThreeMuVector[0], 1,
                        &vel[1][0], 1,
                        &JacVisc[2][1][2][4][0], 1);
            
            Vmath::Neg(nq, &JacVisc[2][1][2][4][0], 1);
            
            // row 5;
            Vmath::Zero(nq, &JacVisc[2][1][4][0][0], 1);
            Vmath::Zero(nq, &JacVisc[2][1][4][1][0], 1);
            Vmath::Zero(nq, &JacVisc[2][1][4][2][0], 1);
            Vmath::Zero(nq, &JacVisc[2][1][4][3][0], 1);
            Vmath::Zero(nq, &JacVisc[2][1][4][4][0], 1);
            
            /* =============================================================
             Avyz = JacVisc[2][2] =
             [ 0,        0,       0,          0,                        0]
             [ 0,        0,       0,          0,                        0]
             [ 0,        0,       0,         mu,                     mu*w]
             [ 0,        0, -2/3*mu,          0,                -2/3*u*mu]
             [ 0,        0,       0,          0,                        0]
             =============================================================*/
            
            // row 1;
            Vmath::Zero(nq, &JacVisc[2][2][0][0][0], 1);
            Vmath::Zero(nq, &JacVisc[2][2][0][1][0], 1);
            Vmath::Zero(nq, &JacVisc[2][2][0][2][0], 1);
            Vmath::Zero(nq, &JacVisc[2][2][0][3][0], 1);
            Vmath::Vcopy(nq, &kPrho2Vec[0], 1, &JacVisc[2][2][0][4][0], 1);
            Vmath::Neg(nq, &JacVisc[2][2][0][4][0], 1);
            
            // row 2;
            Vmath::Zero(nq, &JacVisc[2][2][1][0][0], 1);
            Vmath::Vcopy(nq, &muVec[0], 1, &JacVisc[2][2][1][1][0], 1);
            Vmath::Zero(nq, &JacVisc[2][2][1][2][0], 1);
            Vmath::Zero(nq, &JacVisc[2][2][1][3][0], 1);
            Vmath::Vcopy(nq, &MuVel[0][0], 1, &JacVisc[2][2][1][4][0], 1);
            
            // row 3;
            Vmath::Zero(nq, &JacVisc[2][2][2][0][0], 1);
            Vmath::Zero(nq, &JacVisc[2][2][2][1][0], 1);
            Vmath::Vcopy(nq, &muVec[0], 1, &JacVisc[2][2][2][2][0], 1);
            Vmath::Zero(nq, &JacVisc[2][2][2][3][0], 1);
            Vmath::Vcopy(nq, &MuVel[1][0], 1, &JacVisc[2][2][3][4][0], 1);
            
            // row 4;
            Vmath::Zero(nq, &JacVisc[2][2][3][0][0], 1);
            Vmath::Zero(nq, &JacVisc[2][2][3][1][0], 1);
            Vmath::Zero(nq, &JacVisc[2][2][3][2][0], 1);
            
            Vmath::Vcopy(nq,
                         &FourOThreeMuVector[0], 1,
                         &JacVisc[2][2][3][3][0], 1);
            
            Vmath::Vmul(nq,
                        &FourOThreeMuVector[0], 1,
                        &vel[1][0], 1,
                        &JacVisc[2][2][3][4][0], 1);
            
            // row 5;
            Vmath::Zero(nq, &JacVisc[2][2][4][0][0], 1);
            Vmath::Zero(nq, &JacVisc[2][2][4][1][0], 1);
            Vmath::Zero(nq, &JacVisc[2][2][4][2][0], 1);
            Vmath::Zero(nq, &JacVisc[2][2][4][3][0], 1);
            
            Vmath::Vdiv(nq,
                        &ones[0], 1,
                        &m_baseflow[0][0], 1,
                        &JacVisc[2][2][4][4][0], 1);
            
            Vmath::Vmul(nq,
                        &GamMu_O_GamMinOnePrVector[0], 1,
                        &JacVisc[2][2][4][4][0], 1,
                        &JacVisc[2][2][4][4][0], 1);
            
            
        }
    }
    
    void AdjointCompressibleFlowSystem::GetDerivJacobian(
        Array<OneD, Array<OneD, Array<OneD, Array<OneD, NekDouble> > > > &Jac,
        Array<OneD, Array<OneD, Array<OneD, Array<OneD, NekDouble> > > > &JacDiv)
    {
        int i,j;
        
        int nq = JacDiv[0][0][0].num_elements();
        int nvar = m_fields.num_elements();
        
        
        if (m_spacedim == 2)
        {
            Array<OneD, NekDouble> tmp(nq,0.0);
            
            for (i = 0; i < nvar; ++i)
            {
                for (j = 0; j < nvar; ++j)
                {
                    m_fields[0]->PhysDeriv(Jac[0][i][j],
                                           JacDiv[0][i][j],
                                           tmp);
                    
                    m_fields[0]->PhysDeriv(Jac[1][i][j],
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
                    m_fields[0]->PhysDeriv(Jac[0][i][j],
                                           JacDiv[0][i][j],
                                           tmp1,
                                           tmp2);
                    
                    m_fields[0]->PhysDeriv(Jac[1][i][j],
                                           tmp1,
                                           JacDiv[1][i][j],
                                           tmp2);
                    
                    m_fields[0]->PhysDeriv(Jac[2][i][j],
                                           tmp1,
                                           tmp2,
                                           JacDiv[2][i][j]);
                }
            }
        }
    }
    
    void AdjointCompressibleFlowSystem::GetJacobianConvFluxPrim(
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
        
        Array<OneD, Array<OneD, NekDouble> > tmpx(nvar);
        Array<OneD, Array<OneD, NekDouble> > tmpy(nvar);
        
        for (i = 0; i < m_spacedim; ++i)
        {
            vel[i]          = Array<OneD, NekDouble>(nq);
        }
        
        GetVelocityVector(m_baseflow, vel);
        GetPressure      (m_baseflow, vel, pressure);
        
        for (i = 0; i < m_spacedim; ++i)
        {
            Vmath::Vvtvp(nq, vel[i], 1, vel[i], 1, velsq, 1, velsq, 1);
        }
        
        // determine H = E + p/rho;
        Vmath::Vdiv(nq, m_baseflow[nvar-1], 1, m_baseflow[0], 1, H, 1);
        
        Vmath::Vdiv(nq, pressure, 1, m_baseflow[0], 1, Htmp, 1);
        
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
            Vmath::Vcopy(nq, &m_baseflow[0][0], 1, &Jac[0][1][0][0], 1);
            // 2 * rho * u
            Vmath::Vmul(nq, &m_baseflow[0][0], 1, &vel[0][0], 1, &Jac[0][1][1][0], 1);
            Vmath::Smul(nq, 2.0, &Jac[0][1][1][0], 1, &Jac[0][1][1][0], 1);
            // rho * v
            Vmath::Vmul(nq, &m_baseflow[0][0], 1, &vel[1][0], 1, &Jac[0][1][2][0], 1);
            // rho * (H + u^2)
            Vmath::Zero(nq, &tmp[0], 1);
            Vmath::Vmul(nq, &vel[0][0], 1, &vel[0][0], 1, &tmp[0], 1);
            Vmath::Vadd(nq, &tmp[0], 1, &H[0], 1, &tmp[0], 1);
            Vmath::Vmul(nq, &m_baseflow[0][0], 1, &tmp[0], 1, &Jac[0][1][3][0], 1);
            
            // row three
            // 0
            Vmath::Vcopy(nq, &zeros[0], 1, &Jac[0][2][0][0], 1);
            // 0
            Vmath::Vcopy(nq, &zeros[0], 1, &Jac[0][2][1][0], 1);
            // rho * u
            Vmath::Vmul(nq, &m_baseflow[0][0], 1, &vel[0][0], 1, &Jac[0][2][2][0], 1);
            // rho * u * v
            Vmath::Vmul(nq, &m_baseflow[0][0], 1, &vel[0][0], 1, &Jac[0][2][3][0], 1);
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
            Vmath::Vmul(nq, &m_baseflow[0][0], 1, &vel[1][0], 1, &Jac[1][1][1][0], 1);
            // 0
            Vmath::Vcopy(nq, &zeros[0], 1, &Jac[1][1][2][0], 1);
            // rho * u * v
            Vmath::Vmul(nq, &m_baseflow[0][0], 1, &vel[0][0], 1, &Jac[1][1][3][0], 1);
            Vmath::Vmul(nq, &Jac[1][1][3][0], 1, &vel[1][0], 1, &Jac[1][1][3][0], 1);
            
            // row three
            // rho
            Vmath::Vcopy(nq, &m_baseflow[0][0], 1, &Jac[1][2][0][0], 1);
            // rho * u
            Vmath::Vmul(nq, &m_baseflow[0][0], 1, &vel[0][0], 1, &Jac[1][2][1][0], 1);
            // 2 * rho * v
            Vmath::Vmul(nq, &m_baseflow[0][0], 1, &vel[1][0], 1, &Jac[1][2][2][0], 1);
            Vmath::Smul(nq, 2.0, &Jac[1][2][2][0], 1, &Jac[1][2][2][0], 1);
            // rho * (H + v^2)
            Vmath::Zero(nq, &tmp[0], 1);
            Vmath::Vmul(nq, &vel[1][0], 1, &vel[1][0], 1, &tmp[0], 1);
            Vmath::Vadd(nq, &tmp[0], 1, &H[0], 1, &tmp[0], 1);
            Vmath::Vmul(nq, &m_baseflow[0][0], 1, &tmp[0], 1, &Jac[1][2][3][0], 1);
            
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
            // construction of the Ax Matrix
            /* =================================================================
             Jac[0][nvar][nvar][nq] = Ax^t =
             [   u,      u^2,          u*v,         u*w     (u*(u^2 + v^2))/2]
             [ rho,    rho*u,      2*rho*v,       rho*w     rho*u^2 + rho * H]
             [   0,        0,        rho*v,           0               rho*u*v]
             [   0,        0,            0,       rho*v               rho*u*v]
             [   0,        1,            0,           0     (gam*v)/(gam - 1)]
             */
            // =================================================================
            // row one
            Vmath::Vcopy(nq, &vel[0][0], 1, &Jac[0][0][0][0], 1);
            // u^2
            Vmath::Vmul(nq, &vel[0][0], 1, &vel[0][0], 1, &Jac[0][0][1][0], 1);
            // u * v
            Vmath::Vmul(nq, &vel[0][0], 1, &vel[1][0], 1, &Jac[0][0][2][0], 1);
            // u * w
            Vmath::Vmul(nq, &vel[0][0], 1, &vel[2][0], 1, &Jac[0][0][3][0], 1);
            // 1/2 * u * (u^2 + v^2)
            Vmath::Vmul(nq, &vel[0][0], 1, &velsq[0], 1, &Jac[0][0][4][0], 1);
            Vmath::Smul(nq, 0.5, &Jac[0][0][4][0], 1, &Jac[0][0][4][0], 1);
            
            // row two
            // rho
            Vmath::Vcopy(nq, &m_baseflow[0][0], 1, &Jac[0][1][0][0], 1);
            // 2 * rho * u
            Vmath::Vmul(nq, &m_baseflow[0][0], 1, &vel[0][0], 1, &Jac[0][1][1][0], 1);
            Vmath::Smul(nq, 2.0, &Jac[0][1][1][0], 1, &Jac[0][1][1][0], 1);
            // rho * v
            Vmath::Vmul(nq, &m_baseflow[0][0], 1, &vel[1][0], 1, &Jac[0][1][2][0], 1);
            // rho * w
            Vmath::Vmul(nq, &m_baseflow[0][0], 1, &vel[2][0], 1, &Jac[0][1][3][0], 1);
            // rho * (H + u^2)
            Vmath::Zero(nq, &tmp[0], 1);
            Vmath::Vmul(nq, &vel[0][0], 1, &vel[0][0], 1, &tmp[0], 1);
            Vmath::Vadd(nq, &tmp[0], 1, &H[0], 1, &tmp[0], 1);
            Vmath::Vmul(nq, &m_baseflow[0][0], 1, &tmp[0], 1, &Jac[0][1][4][0], 1);
            
            // row three
            // 0
            Vmath::Vcopy(nq, &zeros[0], 1, &Jac[0][2][0][0], 1);
            // 0
            Vmath::Vcopy(nq, &zeros[0], 1, &Jac[0][2][1][0], 1);
            // rho * u
            Vmath::Vmul(nq, &m_baseflow[0][0], 1, &vel[0][0], 1, &Jac[0][2][2][0], 1);
            // 0
            Vmath::Vcopy(nq, &zeros[0], 1, &Jac[0][2][3][0], 1);
            // rho * u * v
            Vmath::Vmul(nq, &m_baseflow[0][0], 1, &vel[0][0], 1, &Jac[0][2][4][0], 1);
            Vmath::Vmul(nq, &Jac[0][2][4][0], 1, &vel[1][0], 1, &Jac[0][2][4][0], 1);
            
            // row four
            // 0
            Vmath::Vcopy(nq, &zeros[0], 1, &Jac[0][3][0][0], 1);
            // 0
            Vmath::Vcopy(nq, &zeros[0], 1, &Jac[0][3][1][0], 1);
            // 0
            Vmath::Vcopy(nq, &zeros[0], 1, &Jac[0][3][2][0], 1);
            // rho * u
            Vmath::Vmul(nq, &m_baseflow[0][0], 1, &vel[0][0], 1, &Jac[0][3][3][0], 1);
            // rho * u * w
            Vmath::Vmul(nq, &m_baseflow[0][0], 1, &vel[0][0], 1, &Jac[0][2][4][0], 1);
            Vmath::Vmul(nq, &Jac[0][2][4][0], 1, &vel[2][0], 1, &Jac[0][2][4][0], 1);
            
            // row five
            // 0
            Vmath::Vcopy(nq, &zeros[0], 1, &Jac[0][4][0][0], 1);
            // 1
            Vmath::Vcopy(nq, &ones[0], 1, &Jac[0][4][1][0], 1);
            // 0
            Vmath::Vcopy(nq, &zeros[0], 1, &Jac[0][4][2][0], 1);
            // 0
            Vmath::Vcopy(nq, &zeros[0], 1, &Jac[0][4][3][0], 1);
            // u * gamma / (gamma - 1.0)
            Vmath::Vmul(nq, &vel[0][0], 1, &G_over_GminOne[0], 1, &Jac[0][4][4][0], 1);
            
            // construction of the Ay Matrix
            /* =================================================================
             Jac[1][nvar][nvar][nq] = Ay^t =
             [  v,      u*v,          v^2,       v*w     (v*(u^2 + v^2))/2]
             [  0,    rho*v,            0,         0               rho*u*v]
             [rho,    rho*u,      2*rho*v,     rho*w     rho*v^2 + rho * H]
             [  0,        0,            0,     rho*v               rho*u*v]
             [  0,        0,            1,         0     (gam*v)/(gam - 1)]
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
            // v * w
            Vmath::Vmul(nq, &vel[1][0], 1, &vel[2][0], 1, &Jac[1][0][3][0], 1);
            // 1/2 * v * (u^2 + v^2)
            Vmath::Vmul(nq, &vel[1][0], 1, &velsq[0], 1, &Jac[1][0][4][0], 1);
            Vmath::Smul(nq, 0.5, &Jac[1][0][4][0], 1, &Jac[1][0][4][0], 1);
            
            // row two
            // 0
            Vmath::Vcopy(nq, &zeros[0], 1, &Jac[1][1][0][0], 1);
            // rho * v
            Vmath::Vmul(nq, &m_baseflow[0][0], 1, &vel[1][0], 1, &Jac[1][1][1][0], 1);
            // 0
            Vmath::Vcopy(nq, &zeros[0], 1, &Jac[1][1][2][0], 1);
            // 0
            Vmath::Vcopy(nq, &zeros[0], 1, &Jac[1][1][3][0], 1);
            // rho * u * v
            Vmath::Vmul(nq, &m_baseflow[0][0], 1, &vel[0][0], 1, &Jac[1][1][4][0], 1);
            Vmath::Vmul(nq, &Jac[1][1][4][0], 1, &vel[1][0], 1, &Jac[1][1][4][0], 1);
            
            // row three
            // rho
            Vmath::Vcopy(nq, &m_baseflow[0][0], 1, &Jac[1][2][0][0], 1);
            // rho * u
            Vmath::Vmul(nq, &m_baseflow[0][0], 1, &vel[0][0], 1, &Jac[1][2][1][0], 1);
            // 2 * rho * v
            Vmath::Vmul(nq, &m_baseflow[0][0], 1, &vel[1][0], 1, &Jac[1][2][2][0], 1);
            Vmath::Smul(nq, 2.0, &Jac[1][2][2][0], 1, &Jac[1][2][2][0], 1);
            // rho * w
            Vmath::Vmul(nq, &m_baseflow[0][0], 1, &vel[2][0], 1, &Jac[1][2][3][0], 1);
            // rho * (H + v^2)
            Vmath::Zero(nq, &tmp[0], 1);
            Vmath::Vmul(nq, &vel[1][0], 1, &vel[1][0], 1, &tmp[0], 1);
            Vmath::Vadd(nq, &tmp[0], 1, &H[0], 1, &tmp[0], 1);
            Vmath::Vmul(nq, &m_baseflow[0][0], 1, &tmp[0], 1, &Jac[1][2][4][0], 1);
            
            // row four
            // 0
            Vmath::Vcopy(nq, &zeros[0], 1, &Jac[1][3][0][0], 1);
            // 0
            Vmath::Vcopy(nq, &zeros[0], 1, &Jac[1][3][1][0], 1);
            // 0
            Vmath::Vcopy(nq, &zeros[0], 1, &Jac[1][3][2][0], 1);
            // rho * v
            Vmath::Vmul(nq, &m_baseflow[0][0], 1, &vel[1][0], 1, &Jac[1][3][3][0], 1);
            // rho * v * w
            Vmath::Vmul(nq, &m_baseflow[0][0], 1, &vel[1][0], 1, &Jac[1][3][4][0], 1);
            Vmath::Vmul(nq, &Jac[1][3][4][0], 1, &vel[2][0], 1, &Jac[1][3][4][0], 1);
            
            // row five
            // 0
            Vmath::Vcopy(nq, &zeros[0], 1, &Jac[1][4][0][0], 1);
            // 0
            Vmath::Vcopy(nq, &zeros[0], 1, &Jac[1][4][1][0], 1);
            // 1
            Vmath::Vcopy(nq, &ones[0], 1, &Jac[1][4][2][0], 1);
            // 0
            Vmath::Vcopy(nq, &zeros[0], 1, &Jac[1][4][3][0], 1);
            // v * gamma / (gamma - 1.0)
            Vmath::Vmul(nq, &vel[1][0], 1, &G_over_GminOne[0], 1, &Jac[1][4][4][0], 1);
            
            // construction of the Az Matrix
            /* =================================================================
             Jac[2][nvar][nvar][nq] = Az^t =
             [  w,      u*w,          v*w,       w^2     (w*(u^2 + v^2+w^2))/2]
             [  0,    rho*w,            0,         0                   rho*u*v]
             [  0,        0,        rho*w,         0                   rho*u*v]
             [rho,    rho*u,        rho*v,   2*rho*w         rho*w^2 + rho * H]
             [  0,        0,            1,         0         (gam*w)/(gam - 1)]
             */
            // =================================================================
            
            Vmath::Zero(nq, &tmp[0], 1);
            // row one
            // v
            Vmath::Vcopy(nq, &vel[2][0], 1, &Jac[2][0][0][0], 1);
            // u * w
            Vmath::Vmul(nq, &vel[2][0], 1, &vel[0][0], 1, &Jac[2][0][1][0], 1);
            // v * w
            Vmath::Vmul(nq, &vel[1][0], 1, &vel[2][0], 1, &Jac[2][0][2][0], 1);
            // w^2
            Vmath::Vmul(nq, &vel[2][0], 1, &vel[2][0], 1, &Jac[2][0][3][0], 1);
            // 1/2 * w * (u^2 + v^2)
            Vmath::Vmul(nq, &vel[2][0], 1, &velsq[0], 1, &Jac[2][0][4][0], 1);
            Vmath::Smul(nq, 0.5, &Jac[2][0][4][0], 1, &Jac[2][0][4][0], 1);
            
            // row two
            // 0
            Vmath::Vcopy(nq, &zeros[0], 1, &Jac[2][1][0][0], 1);
            // rho * w
            Vmath::Vmul(nq, &m_baseflow[0][0], 1, &vel[2][0], 1, &Jac[2][1][1][0], 1);
            // 0
            Vmath::Vcopy(nq, &zeros[0], 1, &Jac[2][1][2][0], 1);
            // 0
            Vmath::Vcopy(nq, &zeros[0], 1, &Jac[2][1][3][0], 1);
            // rho * u * w
            Vmath::Vmul(nq, &m_baseflow[0][0], 1, &vel[0][0], 1, &Jac[2][1][4][0], 1);
            Vmath::Vmul(nq, &Jac[2][1][4][0], 1, &vel[2][0], 1, &Jac[2][1][4][0], 1);
            
            // row three
            // 0
            Vmath::Vcopy(nq, &zeros[0], 1, &Jac[2][2][0][0], 1);
            // 0
            Vmath::Vcopy(nq, &zeros[0], 1, &Jac[2][2][1][0], 1);
            // rho*w
            Vmath::Vmul(nq, &m_baseflow[0][0], 1, &vel[2][0], 1, &Jac[2][2][2][0], 1);
            // 0
            Vmath::Vcopy(nq, &zeros[0], 1, &Jac[2][2][3][0], 1);
            // rho*v*w
            Vmath::Vmul(nq, &m_baseflow[0][0], 1, &vel[1][0], 1, &Jac[2][2][4][0], 1);
            Vmath::Vmul(nq, &Jac[2][1][4][0], 1, &vel[2][0], 1, &Jac[2][2][4][0], 1);
            
            // row four
            // rho
            Vmath::Vcopy(nq, &m_baseflow[0][0], 1, &Jac[2][3][0][0], 1);
            // rho * u
            Vmath::Vmul(nq, &m_baseflow[0][0], 1, &vel[0][0], 1, &Jac[2][3][1][0], 1);
            // rho * v
            Vmath::Vmul(nq, &m_baseflow[0][0], 1, &vel[1][0], 1, &Jac[2][3][2][0], 1);
            // 2 * rho * w
            Vmath::Vmul(nq, &m_baseflow[0][0], 1, &vel[2][0], 1, &Jac[2][3][3][0], 1);
            Vmath::Smul(nq, 2.0, &Jac[2][3][3][0], 1, &Jac[2][3][3][0], 1);
            // rho * (H + w^2)
            Vmath::Zero(nq, &tmp[0], 1);
            Vmath::Vmul(nq, &vel[2][0], 1, &vel[2][0], 1, &tmp[0], 1);
            Vmath::Vadd(nq, &tmp[0], 1, &H[0], 1, &tmp[0], 1);
            Vmath::Vmul(nq, &m_baseflow[0][0], 1, &tmp[0], 1, &Jac[2][3][4][0], 1);
            
            // row five
            // 0
            Vmath::Vcopy(nq, &zeros[0], 1, &Jac[2][4][0][0], 1);
            // 0
            Vmath::Vcopy(nq, &zeros[0], 1, &Jac[2][4][1][0], 1);
            // 0
            Vmath::Vcopy(nq, &zeros[0], 1, &Jac[2][4][2][0], 1);
            // 1
            Vmath::Vcopy(nq, &ones[0], 1, &Jac[2][4][3][0], 1);
            // v * gamma / (gamma - 1.0)
            Vmath::Vmul(nq, &vel[2][0], 1, &G_over_GminOne[0], 1, &Jac[1][4][4][0], 1);

        }
    }
    
    void AdjointCompressibleFlowSystem::GetJacobianAddConvFlux(
        Array<OneD, Array<OneD, Array<OneD, Array<OneD, NekDouble> > > > &Jac)
    {
        int i;
        int nq = Jac[0][0][0].num_elements();
        
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
                
                GetVelocityVector(m_baseflow, vel);
                GetPressure      (m_baseflow, vel, Pressure);
                
                for (i = 0; i < m_spacedim; ++i)
                {
                    m_fields[0]->PhysDeriv(vel[i], dui_dx[i], dui_dy[i]);
                    
                }
                
                m_fields[0]->PhysDeriv(m_baseflow[0], drho_dxi[0], drho_dxi[1]);
                m_fields[0]->PhysDeriv(Pressure, dP_dxi[0], dP_dxi[1]);
                
                // Calculate d(P/rho^2)/dx = (ƒDer∂rho*dP/dx-2*P*drho/dx)/rho^3
                Vmath::Vdiv(nq, &ones[0], 1, &m_baseflow[0][0], 1, &OneORho[0], 1);
                m_fields[0]->PhysDeriv(OneORho, dOneORho_dxi[0], dOneORho_dxi[1]);
                
                Vmath::Vmul(nq, &m_baseflow[0][0], 1, &m_baseflow[0][0], 1, &PoRho2[0], 1);
                Vmath::Vdiv(nq, &Pressure[0], 1, &PoRho2[0], 1, &PoRho2[0], 1);
                m_fields[0]->PhysDeriv(PoRho2, dPoRho2_dxi[0], dPoRho2_dxi[1]);
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
                
                // gamma/((gamma-1)*Pr)*mu*d(1/rho)/dy)
                
                Vmath::Vmul(nq,
                            &GamMu_O_GamMinOnePrVector[0], 1,
                            &dOneORho_dxi[1][0], 1,
                            &Jac[1][3][3][0], 1);
                
                Vmath::Neg(nq, &Jac[1][3][3][0], 1);
            }
            if (m_spacedim == 3)
            {
                ASSERTL0("false", "3D adjoint solver not yet implemented");
            }
        }
    }
    
    void AdjointCompressibleFlowSystem::GetConservToPrimVariableInvMat(
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &dUdVInv)
    {
        int i;
        int nq = dUdVInv[0][0].num_elements();
        
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
        
        for (i = 0; i < m_spacedim; ++i)
        {
            vel[i]          = Array<OneD, NekDouble>(nq);
            velOverRho[i]   = Array<OneD, NekDouble>(nq);
            OneMinGamVel[i] = Array<OneD, NekDouble>(nq);
        }
        
        GetVelocityVector(m_baseflow, vel);
        GetPressure      (m_baseflow, vel, pressure);
        
        // determine u^2+v^2+w^2;
        for (i = 0; i < m_spacedim; ++i)
        {
            Vmath::Vvtvp(nq, vel[i], 1, vel[i], 1, velsq, 1, velsq, 1);
            Vmath::Vdiv(nq, vel[i], 1, m_baseflow[0], 1, velOverRho[i], 1);
            Vmath::Neg(nq, &velOverRho[i][0], 1);
            Vmath::Vmul(nq, &OneMinGamVec[0], 1,
                        &vel[i][0], 1,
                        &OneMinGamVel[i][0], 1);
        }
        
        Vmath::Smul(nq, 0.5, &velsq[0], 1, &HalfVelsq[0], 1);
        
        Vmath::Vmul(nq, &GamMinOneOverTwoVec[0], 1,
                    &velsq[0], 1,
                    &GamMinOneOverTwoVelsq[0], 1);
        
        Vmath::Vdiv(nq, &ones[0], 1, &m_baseflow[0][0], 1, &oneOverRho[0], 1);
        
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
    
    
    void AdjointCompressibleFlowSystem::GetConservToPrimVariableInvMatDiv(
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &dUdVInv,
        Array<OneD, Array<OneD, Array<OneD, Array<OneD, NekDouble> > > >
                                                        &dUdVInvdXi)
    {
        int i, j;
        int nq = dUdVInvdXi[0][0][0].num_elements();
        int nvar = m_fields.num_elements();
        
        Array<OneD, NekDouble> ones(nq,  1.0);
        Array<OneD, NekDouble> velsq(nq, 0.0);
        Array<OneD, NekDouble> HalfVelsq(nq, 0.0);
        Array<OneD, NekDouble> pressure(nq, 0.0);
        Array<OneD, NekDouble> oneOverRho(nq, 0.0);
        Array<OneD, Array<OneD, NekDouble> > vel(m_spacedim);
        Array<OneD, Array<OneD, NekDouble> > velOverRho(m_spacedim);
        
        for (i = 0; i < m_spacedim; ++i)
        {
            vel[i]            = Array<OneD, NekDouble>(nq, 0.0);
            velOverRho[i]     = Array<OneD, NekDouble>(nq, 0.0);
        }
        
        GetVelocityVector(m_baseflow, vel);
        GetPressure      (m_baseflow, vel, pressure);
        
        // determine u^2+v^2+w^2;
        for (i = 0; i < m_spacedim; ++i)
        {
            Vmath::Vvtvp(nq, vel[i], 1, vel[i], 1, velsq, 1, velsq, 1);
            Vmath::Vdiv(nq, vel[i], 1, m_baseflow[0], 1, velOverRho[i], 1);
            Vmath::Neg(nq, &velOverRho[i][0], 1);
        }
        
        Vmath::Smul(nq, 0.5, &velsq[0], 1, &HalfVelsq[0], 1);
        
        Vmath::Vdiv(nq, &ones[0], 1, &m_baseflow[0][0], 1, &oneOverRho[0], 1);
        
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
                    m_fields[0]->PhysDeriv(dUdVInv[i][j],
                                           dUdVInvdXi[0][i][j],
                                           dUdVInvdXi[1][i][j]);
                }
            }
        }
        // Needs to be checked still
        if (m_spacedim == 3)
        {
            for (i = 0; i < nvar; ++i)
            {
                for (j = 0; j < nvar; ++j)
                {
                    m_fields[0]->PhysDeriv(dUdVInv[i][j],
                                           dUdVInvdXi[0][i][j],
                                           dUdVInvdXi[1][i][j],
                                           dUdVInvdXi[2][i][j]);
                }
            }
        }
    }
    
    void AdjointCompressibleFlowSystem::GetFwdBwdBaseFlow(
                                Array<OneD, Array<OneD, NekDouble> > &FwdBase,
                                Array<OneD, Array<OneD, NekDouble> > &BwdBase)
    {
        int nConvectiveFields = FwdBase.num_elements();
        Array<OneD, Array<OneD, NekDouble> > directSol(nConvectiveFields);
        
        for (int i = 0; i < nConvectiveFields; ++i)
        {
            m_fields[i]->GetFwdBwdTracePhys(m_baseflow[i],
                                            FwdBase[i],
                                            BwdBase[i]);
        }
        
        // correct for boundary conditions // HACK so check for better implementation
        
        //int bcRegion = 0;
        int eMax     = 0;
        int cnt = 0;
        const Array<OneD, const int> &traceBndMap =
        m_fields[0]->GetTraceBndMap();
        
        int nBCregions = m_fields[0]->GetBndConditions().num_elements();
        
        for (int k = 0; k < nBCregions; ++k)
        {
            eMax = m_fields[0]->GetBndCondExpansions()[k]->GetExpSize();
            for (int e = 0; e < eMax; ++e)
            {
                int nBCEdgePts = m_fields[0]->GetBndCondExpansions()[k]->
                GetExp(e)->GetTotPoints();
                int id1 = m_fields[0]->GetBndCondExpansions()[k]->
                GetPhys_Offset(e);
                int id2 = m_fields[0]->GetTrace()->GetPhys_Offset(traceBndMap[cnt+e]);
                
                
                for (int i = 0; i < nConvectiveFields; ++i)
                {
                    Vmath::Vcopy(nBCEdgePts,
                                 &FwdBase[i][id2], 1,
                                 &BwdBase[i][id2], 1);
                }
            }
            if (m_Target == "Drag" || m_Target == "Lift")
            {
                for (int e = 0; e < eMax; ++e)
                {
                    int nBCEdgePts = m_fields[0]->GetBndCondExpansions()[k]->
                    GetExp(e)->GetTotPoints();
                    int id1 = m_fields[0]->GetBndCondExpansions()[k]->
                    GetPhys_Offset(e);
                    int id2 = m_fields[0]->GetTrace()->GetPhys_Offset(traceBndMap[cnt+e]);
                    
                    Array<OneD, NekDouble> tmp(nBCEdgePts, 0.0);
                    Array<OneD, NekDouble> tmp1(nBCEdgePts, 0.0);
                    Array<OneD, NekDouble> tmpfwdnew(nBCEdgePts, 0.0);
                    Array<OneD, NekDouble> tmpfwd(nBCEdgePts, 0.0);
                    
                    //==========================================================
                    // For Euler
                    
                    if (m_EqTypeStr=="AdjointEulerCFE" && m_EqTypeStr=="AdjointEulerADCFE")
                    {
                        for (int i = 0; i < m_spacedim; ++i)
                        {
                            Vmath::Vvtvp(nBCEdgePts,
                                         &FwdBase[1+i][id2], 1,
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
                                         &BwdBase[1+i][id2], 1,
                                         &BwdBase[1+i][id2], 1);
                        }
                    }
                    
                    if (m_EqTypeStr == "AdjointNavierStokesCFE")
                    {
                        Array<OneD, NekDouble> zeros(nBCEdgePts, 0.0);
                        
                        Vmath::Vcopy(nBCEdgePts,
                                     &FwdBase[0][id2], 1,
                                     &BwdBase[0][id2], 1);
                        
                        Vmath::Vcopy(nBCEdgePts,
                                     &FwdBase[1][id2], 1,
                                     &BwdBase[1][id2], 1);
                        
                        Vmath::Vcopy(nBCEdgePts,
                                     &FwdBase[2][id2], 1,
                                     &BwdBase[2][id2], 1);
                        
                        for (int i = 0; i < m_spacedim; i++)
                        {
                            Vmath::Neg(nBCEdgePts, &BwdBase[i+1][id2], 1);
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
                                         &BwdBase[i+1][id2], 1,
                                         &BwdBase[i+1][id2], 1,
                                         &qterm[0],         1,
                                         &qterm[0],         1);
                        }
                        
                        Vmath::Vdiv(nBCEdgePts,
                                    &qterm[0], 1,
                                    &BwdBase[0][id2], 1,
                                    &BwdBase[m_spacedim+1][id2], 1);
                        
                        Vmath::Vadd(nBCEdgePts,
                                    &TwallVec[0], 1,
                                    &BwdBase[0][id2], 1,
                                    &TwallVec[0], 1);
                        
                        Vmath::Vadd(nBCEdgePts,
                                    &TwallVec[0], 1,
                                    &BwdBase[m_spacedim+1][id2], 1,
                                    &BwdBase[m_spacedim+1][id2], 1);
                        
                    }
                }
            }
            
            cnt += m_fields[0]->GetBndCondExpansions()[k]->GetExpSize();
        }
    }
    
    void AdjointCompressibleFlowSystem::GetPressure(
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
    
    void AdjointCompressibleFlowSystem::GetPressure(
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
    
    bool AdjointCompressibleFlowSystem::v_PostIntegrate(int step)
    {
        NekDouble maxL2 = CalcSteadyState();
        return false;
    }
    
    NekDouble AdjointCompressibleFlowSystem::CalcSteadyState()
    {
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
                cout << "L2_max = "  << maxL2 << endl;
            }
        }
        else if (m_fields.num_elements() == 4)
        {
            if (m_comm->GetRank() == 0)
            {
                cout << "L2_max = "  << maxL2 << endl;
            }
        }
        else if (m_fields.num_elements() == 5)
        {
            if (m_comm->GetRank() == 0)
            {
                cout << "L2_max = "  << maxL2 << endl;
            }
        }
        
        std::ofstream myfile;
        myfile.open ("Convergence.txt", std::ios_base::app);
        
        myfile << maxL2 << endl;
        
        myfile.close();
        
        return maxL2;
    }
    
    void AdjointCompressibleFlowSystem::v_GenerateSummary(SolverUtils::SummaryList& s)
    {
        UnsteadySystem::v_GenerateSummary(s);
    }
    
    void AdjointCompressibleFlowSystem::GetVelocityVector(
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
}

