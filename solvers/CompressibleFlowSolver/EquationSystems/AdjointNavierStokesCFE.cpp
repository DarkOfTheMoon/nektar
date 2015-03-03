///////////////////////////////////////////////////////////////////////////////
//
// File AdjointNavierStokesCFE.cpp
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

#include <CompressibleFlowSolver/EquationSystems/AdjointNavierStokesCFE.h>

namespace Nektar
{
    string AdjointNavierStokesCFE::className = 
        SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
            "AdjointNavierStokesCFE", AdjointNavierStokesCFE::create, 
            "NavierStokes equations in conservative variables.");

    AdjointNavierStokesCFE::AdjointNavierStokesCFE(
            const LibUtilities::SessionReaderSharedPtr& pSession)
        : AdjointCompressibleFlowSystem(pSession)
    {
    }

    void AdjointNavierStokesCFE::v_InitObject()
    {

    
        AdjointCompressibleFlowSystem::v_InitObject();
    
        m_session->LoadParameter ("numCheck",       m_numCheck,            0.0);
        if(m_session->DefinesSolverInfo("PROBLEMTYPE"))
        {
            
            std::string ProblemTypeStr = m_session->GetSolverInfo("PROBLEMTYPE");
            int i;
            for(i = 0; i < (int) SIZE_ProblemType; ++i)
            {
                if(NoCaseStringCompare(ProblemTypeMap[i],ProblemTypeStr) == 0)
                {
                    m_problemType = (ProblemType)i;
                    break;
                }
            }
        }
        else
        {
            m_problemType = (ProblemType)0;
        }
        
        if (m_explicitAdvection)
        {
            m_ode.DefineOdeRhs     (&AdjointNavierStokesCFE::DoOdeRhs,        this);
            m_ode.DefineProjection (&AdjointNavierStokesCFE::DoOdeProjection, this);
        }
        else
        {
            ASSERTL0(false, "Implicit CFE not set up.");
        }
    }

    AdjointNavierStokesCFE::~AdjointNavierStokesCFE()
    {

    }

    void AdjointNavierStokesCFE::v_GenerateSummary(SolverUtils::SummaryList& s)
    {
        AdjointCompressibleFlowSystem::v_GenerateSummary(s);
        SolverUtils::AddSummaryItem(s, "Problem Type", ProblemTypeMap[m_problemType]);
    }

    void AdjointNavierStokesCFE::v_SetInitialConditions(
        NekDouble initialtime, 
        bool dumpInitialConditions,
        const int domain)
    {
        int nvariables = m_fields.num_elements();
        if (m_session->GetComm()->GetRank() == 0)
        {
            cout << endl;
            cout << "======================================================================="<< endl;
            cout << "======================= Initial Conditions ============================" << endl;
            cout << "======================================================================="<< endl;
            cout << endl;
        }
        // The unsteady adjoint solution starts with a zero initial condition for
        // force sensitivity analysis since the problem is driven by the boundary conditions. This might change once we start looking at different types of adjoint problems.
        
        if (m_numCheck == 0)
        {
            EquationSystem::v_SetInitialConditions(initialtime, false);
            
            //insert white noise in initial condition
            NekDouble Noise;
            int phystot = m_fields[0]->GetTotPoints();
            Array<OneD, NekDouble> noise(phystot);
            
            m_session->LoadParameter("Noise", Noise,0.0);
            int m_nConvectiveFields =  m_fields.num_elements();
            
            if(Noise > 0.0)
            {
                for(int i = 0; i < m_nConvectiveFields; i++)
                {
                    Vmath::FillWhiteNoise(phystot,Noise,noise,1,m_comm->GetColumnComm()->GetRank()+1);
                    Vmath::Vadd(phystot,m_fields[i]->GetPhys(),1,noise,1,m_fields[i]->UpdatePhys(),1);
                    m_fields[i]->FwdTrans_IterPerExp(m_fields[i]->GetPhys(),m_fields[i]->UpdateCoeffs());
                }
            }
            
            
            if (dumpInitialConditions)
            {
                // Dump initial conditions to file
                Checkpoint_Output(0);
            }
        }
        else
        {
            if (m_cnt == 0)
            {
                int nq = m_fields[0]->GetNpoints();                
                for (int i = 0; i < m_fields.num_elements(); i++)
                {
                    Vmath::Zero(nq, m_fields[i]->UpdatePhys(), 1);
                    m_fields[i]->SetPhysState(true);
                    Vmath::Zero(m_fields[i]->GetNcoeffs(),
                                m_fields[i]->UpdateCoeffs(), 1);
                    
                    if (m_session->GetComm()->GetRank() == 0)
                    {
                        cout << "  - Field "    << m_session->GetVariable(i)
                        << ": 0 (default)" << endl;
                    }
                }
            }
            else
            {
                
                std::string filename = m_session->GetFunctionFilename("InitialConditions", 0.0);
                std::string chkout = boost::lexical_cast<string>(m_cnt);
                for (int q = 0; q < nvariables; ++q)
                {
                    std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef;
                    std::vector<std::vector<NekDouble> > FieldData;
                    std::string fileVar = m_session->GetVariable(q);
                    
                    Array<OneD, NekDouble> vCoeffs(m_fields[0]->GetNcoeffs());
                    Vmath::Zero(vCoeffs.num_elements(),vCoeffs,1);
                    
                    int numexp = m_fields[0]->GetExpSize();
                    Array<OneD,int> ElementGIDs(numexp);
                    
                    // Define list of global element ids
                    for(int i = 0; i < numexp; ++i)
                    {
                        ElementGIDs[i] = m_fields[0]->GetExp(i)->GetGeom()->GetGlobalID();
                    }
                    
                    std::string filename_new = filename+"_"+chkout+".chk";
                    m_fld->Import(filename_new, FieldDef, FieldData,
                                  LibUtilities::NullFieldMetaDataMap,
                                  ElementGIDs);
                    
                    int idx = -1;
                    
                    // Loop over all the expansions
                    for (int i = 0; i < FieldDef.size(); ++i)
                    {
                        // Find the index of the required field in the
                        // expansion segment
                        for(int j = 0; j < FieldDef[i]->m_fields.size(); ++j)
                        {
                            if (FieldDef[i]->m_fields[j] == fileVar)
                            {
                                idx = j;
                            }
                        }
                        if (idx >= 0)
                        {
                            m_fields[0]->ExtractDataToCoeffs(
                                        FieldDef[i], FieldData[i],
                                        FieldDef[i]->m_fields[idx], vCoeffs);
                        }
                        else
                        {
                            cout << "Field " + fileVar + " not found." << endl;
                        }
                    }
                    
                    m_fields[0]->BwdTrans_IterPerExp(vCoeffs, m_fields[q]->UpdatePhys());
                }
                
                if (m_session->GetComm()->GetRank() == 0)
                {
                    for (int i = 0; i < m_fields.num_elements(); ++i)
                    {
                        std::string varName = m_session->GetVariable(i);
                        cout << "  - Field " << varName << ": "
                        << "from file " << filename+"_"+chkout+".chk"
                        << endl;
                    }
                }
            }
            
            if (dumpInitialConditions)
            {
                // Dump initial conditions to file
                Checkpoint_Output(0);
            }
        }
        
        if (m_session->GetComm()->GetRank() == 0)
        {
            cout << endl;
            cout << "======================================================================="<< endl;
            cout << endl;
        }
    }

    void AdjointNavierStokesCFE::DoOdeRhs(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
              Array<OneD,       Array<OneD, NekDouble> > &outarray,
        const NekDouble                                   time)
    {
        int i;
        int nvariables = inarray.num_elements();
        int npoints    = GetNpoints();
        
        Array<OneD, Array<OneD, NekDouble> > advVel(m_spacedim);
        Array<OneD, Array<OneD, NekDouble> > outarrayAdv(nvariables);
        Array<OneD, Array<OneD, NekDouble> > outarrayDiff(nvariables);
        
        // Advection term in physical rhs form
        for (i = 0; i < nvariables; ++i)
        {
            outarrayAdv[i] = Array<OneD, NekDouble>(npoints, 0.0);
            outarrayDiff[i] = Array<OneD, NekDouble>(npoints, 0.0);
        }
            
        // Advection term in physical rhs form
        
        m_advection->Advect(nvariables,
                            m_fields,
                            advVel,
                            inarray,
                            outarrayAdv,
                            time);
            
        // Diffusion term in physical rhs form
        //
        m_diffusion->Diffuse(nvariables, m_fields, inarray, outarrayDiff);
        
        for (i = 0; i < nvariables; ++i)
        {
            Vmath::Vadd(npoints,
                        outarrayAdv[i], 1,
                        outarrayDiff[i], 1,
                        outarray[i], 1);
        }
    }

    void AdjointNavierStokesCFE::DoOdeProjection(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
              Array<OneD,       Array<OneD, NekDouble> > &outarray,
        const NekDouble                                   time)
    {
        int i;
        int nvariables = inarray.num_elements();
        
        switch(m_projectionType)
        {
            case MultiRegions::eDiscontinuous:
            {
                // Just copy over array
                int npoints = GetNpoints();
                
                for(i = 0; i < nvariables; ++i)
                {
                    Vmath::Vcopy(npoints, inarray[i], 1, outarray[i], 1);
                }
                SetBoundaryConditions(outarray, time);
                break;
            }
            case MultiRegions::eGalerkin:
            case MultiRegions::eMixed_CG_Discontinuous:
            {
                ASSERTL0(false, "No Continuous Galerkin for full compressible "
                                "Navier-Stokes equations");
                break;
            }
            default:
                ASSERTL0(false, "Unknown projection scheme");
                break;
        }
    }
    
    void AdjointNavierStokesCFE::SetBoundaryConditions(
        Array<OneD, Array<OneD, NekDouble> > &inarray,
        NekDouble                             time)
    {
        std::string varName;
        int nvariables = m_fields.num_elements();
        int cnt        = 0;
        
        // loop over Boundary Regions
        for (int n = 0; n < m_fields[0]->GetBndConditions().num_elements(); ++n)
        {
            // Symmetric Boundary Condition
            if (m_fields[0]->GetBndConditions()[n]->GetUserDefined() == 
                SpatialDomains::eSymmetry)
            {
                SymmetryBC(n, cnt, inarray);
            }
            
            // Adjoint wall Condition
            if (m_fields[0]->GetBndConditions()[n]->GetUserDefined() ==
                SpatialDomains::eAdjointWall)
            {
                AdjointWallBC(n, cnt, inarray);
            }
            
            // Time Dependent Boundary Condition (specified in meshfile)
            if (m_fields[0]->GetBndConditions()[n]->GetUserDefined() 
                == SpatialDomains::eTimeDependent)
            {
                for (int i = 0; i < nvariables; ++i)
                {
                    varName = m_session->GetVariable(i);
                    m_fields[i]->EvaluateBoundaryConditions(time, varName);
                }
            }
    
            cnt += m_fields[0]->GetBndCondExpansions()[n]->GetExpSize();
        }
    }
}
