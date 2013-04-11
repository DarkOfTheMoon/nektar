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

        // Read the geometry and the expansion information
        //m_graph = SpatialDomains::MeshGraph::Read(pSession);
        int numdomains = m_graph->GetDomain().size();
        int numexpansions = numfields*numdomains;

            /*if(m_graph->GetDomain().size() > 1)
        {
            m_regions = Array<OneD, MultiRegions::ExpListSharedPtr> (numexpansions);
            const SpatialDomains::CompositeMap domain = (m_graph->GetDomain());
			
            for(i = 0 ; i < numexpansions; i++)
            {
                m_regions[i] = MemoryManager<MultiRegions::ContField3D>
                    ::AllocateSharedPtr(m_session,domain,m_graph,m_session->GetVariable(i%numfields),(i/numfields));
            }
			
            // Only needed for output: whole field
            m_outfields = Array<OneD, MultiRegions::ExpListSharedPtr> (nvariables);
            for(i = 0 ; i < m_outfields.num_elements(); i++)
            {
                m_outfields[i] = MemoryManager<MultiRegions::DisContField1D>::AllocateSharedPtr(m_session,m_graph,
                                                                                                m_session->GetVariable(i));
            }
            }*/


        
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
                if (m_fields[0]->GetBndConditions()[n]->GetUserDefined() != SpatialDomains::eNoUserDefined)
                {
                    if (m_fields[0]->GetBndConditions()[n]->GetUserDefined() != SpatialDomains::eTimeDependent)
                    {                     	     
                        if(m_fields[0]->GetBndConditions()[n]->GetUserDefined() != SpatialDomains::eI)
                        {  	 	 
                            ASSERTL0(false,"Unknown USERDEFINEDTYPE boundary condition");
                        }
                    }
                }
            }
            break;
        case eNoEquationType:
        default:
            ASSERTL0(false,"Unknown or undefined equation type");
        }

        //Kinematic viscosity
        m_session->LoadParameter("Kinvis", m_kinvis);

        // Load variable coefficients
        m_perm = Array<OneD, NekDouble> (3*(m_spacedim-1));
        m_spatialperm = Array<OneD, Array<OneD, NekDouble> > (m_spacedim);

        // Inverted Permeability Matrix
        m_perm_inv = Array<OneD, NekDouble> (3*(m_spacedim-1));

        //Setup of spatially varying anisotropic permeability
        if(m_session->DefinesFunction("SpatialAnisotropicPermeability"))
        {
            ASSERTL0(!m_explicitPermeability,"implicit spatially varying permeability not implemented");

            int nq = m_fields[0]->GetNpoints();
            
            if (m_spacedim == 2)
            {
                std::string varCoeffs[2] = {
                    "kxx",
                    "kyy",
                };

                Array<OneD, NekDouble> vTemp;
                for (int i = 0; i < m_spacedim; ++i)
                {
                    EvaluateFunction(varCoeffs[i], vTemp, "SpatialAnisotropicPermeability");
                    m_spatialperm[i] = Array<OneD, NekDouble>(nq);
                    Vmath::Sdiv(nq,1.0,vTemp,1,m_spatialperm[i],1);
                }
                Array<OneD,NekDouble> x0(nq);
                Array<OneD,NekDouble> x1(nq);
                Array<OneD,NekDouble> x2(nq);
                
                // Get the coordinates (assuming all fields have the same
                // discretisation)
                NekDouble scalefac = 10;
                m_fields[0]->GetCoords(x0,x1,x2);
                for(int j=0; j<nq; ++j)
                {
                    NekDouble x = x0[j];
                    NekDouble y = x1[j];
                    if(x > 0.375 && x <0.625)
                    {
                        //cout<<"x: "<<x0[j]<<" y: "<<x1[j]<<" z: "<<x2[j]<<endl;
                        m_spatialperm[0][j]=m_spatialperm[0][j]*scalefac;
                        m_spatialperm[1][j]=m_spatialperm[1][j]*scalefac;
                    }
                }

                // Transform variable coefficient and write out to file.
                m_fields[0]->FwdTrans_IterPerExp(m_spatialperm[i],
                                                 m_fields[0]->UpdateCoeffs());
                std::stringstream filename;
                filename << "AnisotropicPerm_" << varCoeffs[i];
                if (m_comm->GetSize() > 1)
                {
                    filename << "_P" << m_comm->GetRank();
                }
                filename << ".fld";
                WriteFld(filename.str());
                
            }
            if(m_spacedim == 3)
            {
                //std::string varName = "k";
                std::string varCoeffs[3] = {
                    "kxx",
                    "kyy",
                    "kzz",
                }; 

                //Explicit implementation
                Array<OneD, NekDouble> vTemp;
                for (int i = 0; i < m_spacedim; ++i)
                {
                    EvaluateFunction(varCoeffs[i], vTemp, "SpatialAnisotropicPermeability");
                    m_spatialperm[i] = Array<OneD, NekDouble>(nq);
                    Vmath::Sdiv(nq,1.0,vTemp,1,m_spatialperm[i],1);
                }
            }
        }
        else if (m_session->DefinesFunction("AnisotropicPermeability"))
        {
            if (m_spacedim == 2)
            {
                std::string varCoeffs[3] = {
                    "kxx",
                    "kyy",
                    "kxy"
                };
                for (int i = 0; i < (3*(m_spacedim-1)); ++i)
                {
                    ASSERTL0(m_session->DefinesFunction("AnisotropicPermeability", varCoeffs[i]),
                             "Function '" + varCoeffs[i] + "' not correctly defined.");
                    m_perm[i] = m_session->GetFunction("AnisotropicPermeability", varCoeffs[i])->Evaluate();
                }

                NekDouble detTemp = m_perm[0]*m_perm[1]-m_perm[2]*m_perm[2];
                
                // Check if permeability matrix is positive definite
                ASSERTL0(m_perm[0] > 0,"Permeability Matrix is not positive definite");
                ASSERTL0(detTemp > 0,"Permeability Matrix is not positive definite");
            
                m_perm_inv[0] = m_perm[1];
                m_perm_inv[1] = m_perm[0];
                m_perm_inv[2] = -m_perm[2];
                Vmath::Smul(3, 1/detTemp, m_perm_inv, 1, m_perm_inv, 1);
            }
            else
            {
                std::string varCoeffs[6] = {
                    "kxx",
                    "kyy",
                    "kzz",
                    "kxy",
                    "kxz",
                    "kyz"
                };                
                for (int i = 0; i < (3*(m_spacedim-1)); ++i)
                {
                    ASSERTL0(m_session->DefinesFunction("AnisotropicPermeability", varCoeffs[i]),
                             "Function '" + varCoeffs[i] + "' not correctly defined.");
                    m_perm[i] = m_session->GetFunction("AnisotropicPermeability", varCoeffs[i])->Evaluate();
                }

                NekDouble detTemp = m_perm[0]*(m_perm[1]*m_perm[2]-m_perm[5]*m_perm[5])
                    -m_perm[3]*(m_perm[2]*m_perm[3]-m_perm[4]*m_perm[5])
                    +m_perm[4]*(m_perm[3]*m_perm[5]-m_perm[1]*m_perm[4]);
                
                // Check if permeability matrix is positive definite
                ASSERTL0(m_perm[0] > 0,"Permeability Matrix is not positive definite");
                NekDouble pd_chk = m_perm[0]*m_perm[1]-m_perm[3]*m_perm[3];
                ASSERTL0(pd_chk > 0,"Permeability Matrix is not positive definite");
                ASSERTL0(detTemp > 0,"Permeability Matrix is not positive definite");

                m_perm_inv[0] = m_perm[1]*m_perm[2]-m_perm[5]*m_perm[5];
                m_perm_inv[1] = m_perm[0]*m_perm[2]-m_perm[4]*m_perm[4];
                m_perm_inv[2] = m_perm[0]*m_perm[1]-m_perm[3]*m_perm[3];
                m_perm_inv[3] = m_perm[4]*m_perm[5]-m_perm[2]*m_perm[3];
                m_perm_inv[4] = m_perm[3]*m_perm[5]-m_perm[1]*m_perm[4];
                m_perm_inv[5] = m_perm[3]*m_perm[4]-m_perm[0]*m_perm[5];
                
                Vmath::Smul(6, 1/detTemp, m_perm_inv, 1, m_perm_inv, 1);
            }
        }
        else if (m_session->DefinesParameter("Permeability"))
        {
            NekDouble kTemp;
            m_session->LoadParameter("Permeability", kTemp);
            
            for (int i = 0; i < m_spacedim; ++i)
            {
                m_perm[i] = kTemp;
            }
            for (int i = (m_spacedim+1); i < (3*(m_spacedim-1)); ++i)
            {
                m_perm[i] = 0;
            }

            NekDouble detTemp = m_perm[0]*(m_perm[1]*m_perm[2]-m_perm[5]*m_perm[5])
                -m_perm[3]*(m_perm[2]*m_perm[3]-m_perm[4]*m_perm[5])
                +m_perm[4]*(m_perm[3]*m_perm[5]-m_perm[1]*m_perm[4]);
            
            // Check if permeability matrix is positive definite
            ASSERTL0(m_perm[0] > 0,"Permeability Matrix is not positive definite");
            NekDouble pd_chk = m_perm[0]*m_perm[1]-m_perm[3]*m_perm[3];
            ASSERTL0(pd_chk > 0,"Permeability Matrix is not positive definite");
            ASSERTL0(detTemp > 0,"Permeability Matrix is not positive definite");
            
            m_perm_inv[0] = m_perm[1]*m_perm[2]-m_perm[5]*m_perm[5];
            m_perm_inv[1] = m_perm[0]*m_perm[2]-m_perm[4]*m_perm[4];
            m_perm_inv[2] = m_perm[0]*m_perm[1]-m_perm[3]*m_perm[3];
            m_perm_inv[3] = m_perm[4]*m_perm[5]-m_perm[2]*m_perm[3];
            m_perm_inv[4] = m_perm[3]*m_perm[5]-m_perm[1]*m_perm[4];
            m_perm_inv[5] = m_perm[3]*m_perm[4]-m_perm[0]*m_perm[5];
            
            Vmath::Smul(6, 1/detTemp, m_perm_inv, 1, m_perm_inv, 1);
        }
        else
        {
            ASSERTL0(0,"Permeability not defined");
        }

        std::string vConvectiveType = "NoAdvection";
        m_advObject = GetAdvectionTermFactory().CreateInstance(vConvectiveType, m_session, m_graph);
	
#if 0 // Not required if building on an UnsteadySystem rather than an EquationSystem
        // Set up filters
        LibUtilities::FilterMap::const_iterator x;
        LibUtilities::FilterMap f = m_session->GetFilters();
        for (x = f.begin(); x != f.end(); ++x)
        {
            m_filters.push_back(SolverUtils::GetFilterFactory().CreateInstance(x->first, m_session, x->second));
        }
#endif
    }

    PorousMedia::~PorousMedia(void)
    {

    }
    
    void PorousMedia::AdvanceInTime(int nsteps)
    {
        int i,n;
        static int nchk = 0;
		
        Timer timer;
	
        // Set up wrapper to fields data storage. 
        Array<OneD, Array<OneD, NekDouble> >  fields(m_nConvectiveFields);
        for(i = 0; i < m_nConvectiveFields; ++i)
        {
            fields[i]  = m_fields[i]->UpdatePhys();
        }
		
        // Initialise NS solver which is set up to use a GLM method
        // with calls to EvaluateAdvection_Permeability and
        // SolveUnsteadyStokesSystem
        m_integrationSoln = m_integrationScheme[m_intSteps-1]->InitializeScheme(m_timestep, fields, m_time, m_integrationOps);
		
        std::string   mdlname = m_session->GetSessionName() + ".mdl";
        std::ofstream mdlFile;

        
        std::vector<SolverUtils::FilterSharedPtr>::iterator x;
        for (x = m_filters.begin(); x != m_filters.end(); ++x)
        {
            (*x)->Initialise(m_fields, m_time);
        }

        //Time advance
        for(n = 0; n < nsteps; ++n)
        {

            timer.Start();

            // Advance velocity fields
            fields = m_integrationScheme[min(n,m_intSteps-1)]->TimeIntegrate(m_timestep, m_integrationSoln, m_integrationOps);
            
            m_time += m_timestep;
            
            timer.Stop();

            // Write out current time step
            if(m_infosteps && !((n+1)%m_infosteps) && m_comm->GetRank() == 0)
            {
                cout << "Step: " << n+1 << "  Time: " << m_time << " CPU-Time: " << timer.TimePerTest(1) << " s" << endl;
            }
            
            // dump data in m_fields->m_coeffs to file. 
            if(m_checksteps && n&&(!((n+1)%m_checksteps)))
            {
                for(i = 0; i < m_nConvectiveFields; ++i)
                {
                    m_fields[i]->SetPhys(fields[i]);
                    m_fields[i]->SetPhysState(true);
                }
                nchk++;
                Checkpoint_Output(nchk);
            }
            
            
            if(m_steadyStateSteps && n && (!((n+1)%m_steadyStateSteps)))
            {
                if(CalcSteadyState() == true)
                {
                    cout << "Reached Steady State to tolerance " << m_steadyStateTol << endl; 
                    break;
                }
            }

            // Transform data into coefficient space
            if (m_filters.size() > 0)
            {
                for (i = 0; i < m_nConvectiveFields; ++i)
                {
                    m_fields[i]->FwdTrans_IterPerExp(fields[i],
                                                     m_fields[i]->UpdateCoeffs());
                    m_fields[i]->SetPhysState(false);
                }
            }

            std::vector<SolverUtils::FilterSharedPtr>::iterator x;
            for (x = m_filters.begin(); x != m_filters.end(); ++x)
            {
                (*x)->Update(m_fields, m_time);
            }

        }
	
        for(i = 0; i < m_nConvectiveFields; ++i)
        {
            m_fields[i]->SetPhys(fields[i]);
            m_fields[i]->SetPhysState(true);
        }
        
        
        for (x = m_filters.begin(); x != m_filters.end(); ++x)
        {
            (*x)->Finalise(m_fields, m_time);
        }
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
    
    // Decide if at a steady state if the discrerte L2 sum of the
    // coefficients is the same as the previous step to within the
    // tolerance m_steadyStateTol;
    bool PorousMedia::CalcSteadyState(void)
    {
        static NekDouble previousL2 = 0.0;
        bool returnval = false;

        NekDouble L2 = 0.0;
        
        // calculate L2 discrete summation 
        int ncoeffs = m_fields[0]->GetNcoeffs(); 
        
        for(int i = 0; i < m_fields.num_elements(); ++i)
        {
            L2 += Vmath::Dot(ncoeffs,m_fields[i]->GetCoeffs(),1,m_fields[i]->GetCoeffs(),1);
        }
        
        if(fabs(L2-previousL2) < ncoeffs*m_steadyStateTol)
        {
            returnval = true;
        }

        previousL2 = L2;

        return returnval;
    }
} //end of namespace

