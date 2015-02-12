///////////////////////////////////////////////////////////////////////////////
//
// File PorousMedia.cpp
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
// Description: Porous Media class definition
//
///////////////////////////////////////////////////////////////////////////////

#include <PorousMediaSolver/EquationSystems/PorousMedia.h>
#include <LibUtilities/BasicUtils/Timer.h>
#include <cstdio>
#include <cstdlib>
#include <LibUtilities/Communication/Comm.h>
#include <SolverUtils/Filters/Filter.h>
#include <iomanip>

namespace Nektar
{

    /**
     * Constructor. Creates ...
     *
     * \param 
     * \param
     */
    PorousMedia::PorousMedia(const LibUtilities::SessionReaderSharedPtr& pSession):
        UnsteadySystem(pSession)
    {
    }

    void PorousMedia::v_InitObject()
    {
        int i,j;
        int numfields = m_fields.num_elements();
        std::string velids[] = {"u","v","w"};

        // Set up Velocity field to point to the first m_expdim of m_fields; 
        m_velocity = Array<OneD,int>(m_spacedim);
        
        for(i = 0; i < m_spacedim; ++i)
        {
            for(j = 0; j < numfields; ++j)
            {
                std::string var = m_boundaryConditions->GetVariable(j);
                if(NoCaseStringCompare(velids[i],var) == 0)
                {
                    m_velocity[i] = j;
                    break;
                }
                
                if(j == numfields)
                {
                    std::string error = "Failed to find field: " + var; 
                    ASSERTL0(false,error.c_str());
                }
            }
        }
        
        // Set up equation type enum using kEquationTypeStr
        for(i = 0; i < (int) eEquationTypeSize; ++i)
        {
            bool match;
            m_session->MatchSolverInfo("EQTYPE",kEquationTypeStr[i],match,false);
            if(match)
            {
                m_equationType = (EquationType)i; 
                break;
            }
        }
        ASSERTL0(i != eEquationTypeSize,"EQTYPE not found in SOLVERINFO section");


        // Determine how Darcy term is dealt with
        for (i = 0; i < (int)SIZE_DarcyTerm; ++i)
        {
            bool match;
            m_session->MatchSolverInfo("PERMEABILITYADVANCEMENT",DarcyTermMethodStr[i],match,false);
            if(match)
            {
                m_darcyType = (DarcyTermMethod)i; 
                break;
            }
        }
        ASSERTL0(i != SIZE_DarcyTerm,"PERMEABILITYADVANCEMENT not found in SOLVERINFO section");

        //m_session->MatchSolverInfo("PERMEABILITYADVANCEMENT", "Explicit",
        //                           m_explicitPermeability, true);
        
        // This probably should to into specific implementations 
        // Equation specific Setups 
        switch(m_equationType)
        {
        case eUnsteadyPorousMedia:
            m_session->LoadParameter("IO_InfoSteps", m_infosteps, 0);
            m_session->LoadParameter("SteadyStateSteps", m_steadyStateSteps, 0);
            m_session->LoadParameter("SteadyStateTol", m_steadyStateTol, 1e-6);
            
            // check to see if any user defined boundary condition is
            // indeed implemented
            for(int n = 0; n < m_fields[0]->GetBndConditions().num_elements(); ++n)
            {
                // Time Dependent Boundary Condition (if no user
                // defined then this is empty)
                ASSERTL0 (
                    m_fields[0]->GetBndConditions()[n]->GetUserDefined() ==
                    SpatialDomains::eNoUserDefined ||
                    m_fields[0]->GetBndConditions()[n]->GetUserDefined() ==
                    SpatialDomains::eWall_Forces ||
                    m_fields[0]->GetBndConditions()[n]->GetUserDefined() ==
                    SpatialDomains::eTimeDependent ||
                    m_fields[0]->GetBndConditions()[n]->GetUserDefined() ==
                    SpatialDomains::eRadiation ||
                    m_fields[0]->GetBndConditions()[n]->GetUserDefined() ==
                    SpatialDomains::eI ||
                    m_fields[0]->GetBndConditions()[n]->GetUserDefined() ==
                    SpatialDomains::eHighOutflow,
                    "Unknown USERDEFINEDTYPE boundary condition");
            }
            break;
        case eNoEquationType:
        default:
            ASSERTL0(false,"Unknown or undefined equation type");
        }

        // Initialise advection
        m_advObject = GetAdvectionTermFactory().CreateInstance("NoAdvection", m_session, m_graph);

        // Forcing terms
        m_forcing = SolverUtils::Forcing::Load(m_session, m_fields,
                                               v_GetForceDimension());

        //Kinematic viscosity
        m_session->LoadParameter("Kinvis", m_kinvis);

        m_darcyEvaluation=GetDarcyTermFactory().CreateInstance(DarcyTermMethodStr[m_darcyType],m_session,m_fields);

        //Setup permeability matrix
        m_darcyEvaluation->SetupPermeability();

        /*     //Set Darcy factor for implicit step
        m_darcy_fac = Array<OneD, NekDouble> (3*(m_spacedim-1));
        m_darcyEvaluation->GetImplicitDarcyFactor(m_darcy_fac);

        //temp - need to decide if we keep the implicit implementation
        m_perm_inv = Array<OneD, NekDouble> (3*(m_spacedim-1));
        m_darcyEvaluation->GetImplicitDarcyFactor(m_perm_inv);
        //m_darcyEvaluation->GetInversePermeability(m_perm_inv)
        //Kinematic viscosity
        m_session->LoadParameter("Kinvis", m_kinvis);

        // Load variable coefficients
        m_perm = Array<OneD, NekDouble> (3*(m_spacedim-1));
        m_spatialperm = Array<OneD, Array<OneD, NekDouble> > (m_spacedim);

        // Inverted Permeability Matrix
        //m_perm_inv = Array<OneD, NekDouble> (3*(m_spacedim-1));
        */
        m_extrapolation = MemoryManager<Extrapolate>::AllocateSharedPtr(
            m_session,
            m_fields,
            m_pressure,
            m_velocity);
        m_extrapolation->TimeIntegrationSteps(m_intScheme->GetIntegrationMethod(), m_intScheme);
    }

    PorousMedia::~PorousMedia(void)
    {

    }
    

    // Evaluation -N(V) for all fields except pressure using m_velocity

    void PorousMedia::EvaluateAdvectionTerms(const Array<OneD, const Array<OneD, NekDouble> > &inarray, 
                                                 Array<OneD, Array<OneD, NekDouble> > &outarray, 
                                                 Array<OneD, NekDouble> &wk)
    {
        int i;
        int nqtot      = m_fields[0]->GetTotPoints();
        int VelDim     = m_velocity.num_elements();
        Array<OneD, Array<OneD, NekDouble> > velocity(VelDim);
        Array<OneD, NekDouble > Deriv;
        for(i = 0; i < VelDim; ++i)
        {
            velocity[i] = inarray[m_velocity[i]]; 
        }
        // Set up Derivative work space; 
        if(wk.num_elements())
        {
            ASSERTL0(wk.num_elements() >= nqtot*VelDim,"Workspace is not sufficient");            
            Deriv = wk;
        }
        else
        {
            Deriv = Array<OneD, NekDouble> (nqtot*VelDim);
        }

        m_advObject->DoAdvection(m_fields,m_nConvectiveFields, 
                                 m_velocity,inarray,outarray,m_time,Deriv);
    }
    
    //time dependent boundary conditions updating    
    void PorousMedia::SetBoundaryConditions(NekDouble time)
    {
        int  nvariables = m_fields.num_elements();
	
        for (int i = 0; i < nvariables; ++i)
        {
            for(int n = 0; n < m_fields[i]->GetBndConditions().num_elements(); ++n)
            {	
                if(m_fields[i]->GetBndConditions()[n]->GetUserDefined() == SpatialDomains::eTimeDependent)
                {
                    m_fields[i]->EvaluateBoundaryConditions(time);
                }
            }
        }
    }
    

    /**
     * Add an additional forcing term programmatically.
     */
    void PorousMedia::AddForcing(const SolverUtils::ForcingSharedPtr& pForce)
    {
        m_forcing.push_back(pForce);
    }

} //end of namespace

