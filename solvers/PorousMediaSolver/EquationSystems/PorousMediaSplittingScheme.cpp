///////////////////////////////////////////////////////////////////////////////
//
// File PorousMediaSplittingScheme.cpp
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
// Description: Splitting Scheme for the Porous media
///////////////////////////////////////////////////////////////////////////////

#include <PorousMediaSolver/EquationSystems/PorousMediaSplittingScheme.h>

namespace Nektar
{
    string PorousMediaSplittingScheme::className = 
        SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction
        ("PorousMediaSplittingScheme", PorousMediaSplittingScheme::create);
    
     /**
     * Constructor. Creates ...
     *
     * \param 
     * \param
     */
    PorousMediaSplittingScheme::PorousMediaSplittingScheme(
            const LibUtilities::SessionReaderSharedPtr& pSession):
        PorousMedia(pSession)
    {
        
    }

    void PorousMediaSplittingScheme::v_InitObject()
    {
        UnsteadySystem::v_InitObject();
        PorousMedia::v_InitObject();
        // Set m_pressure to point to last field of m_fields; 
        if(NoCaseStringCompare(m_session->GetVariable(m_fields.num_elements()-1),"p") == 0)
        {
            m_nConvectiveFields = m_fields.num_elements()-1;
            m_pressure = m_fields[m_nConvectiveFields];
            m_pressureCalls = 1;
        }
        else
        {
            ASSERTL0(false,"Need to set up pressure field definition");
        }
        
        LibUtilities::TimeIntegrationMethod intMethod;
        std::string TimeIntStr = m_session->GetSolverInfo("TimeIntegrationMethod");
        int i;
        for(i = 0; i < (int) LibUtilities::SIZE_TimeIntegrationMethod; ++i)
        {
            if(NoCaseStringCompare(LibUtilities::TimeIntegrationMethodMap[i],TimeIntStr) == 0 )
            {
                intMethod = (LibUtilities::TimeIntegrationMethod)i; 
                break;
            }
        }
        
        ASSERTL0(i != (int) LibUtilities::SIZE_TimeIntegrationMethod, "Invalid time integration type.");

        if(!m_explicitPermeability && (m_session->GetSolverInfo("AdvectionForm") == "NoAdvection"))
        {
            switch(intMethod)
            {
                case LibUtilities::eIMEXOrder1:
                case LibUtilities::eIMEXOrder2:
                case LibUtilities::eIMEXOrder3:
                {
                    ASSERTL0(0,"Integration method not suitable: replace IMEX scheme with an purely implicit scheme");
                }
                break;
                default:
                break;
            }
        }

        switch(intMethod)
        {
            case LibUtilities::eBackwardEuler:
            case LibUtilities::eDIRKOrder2:
            case LibUtilities::eDIRKOrder3:
            {
                ASSERTL0( !m_explicitPermeability && (m_session->GetSolverInfo("AdvectionForm") == "NoAdvection"),
                          "Integration method not suitable: Options include IMEXOrder1, IMEXOrder2 or IMEXOrder3");
                m_intSteps = 1;
                m_integrationScheme = Array<OneD, LibUtilities::TimeIntegrationSchemeSharedPtr> (m_intSteps);
                LibUtilities::TimeIntegrationSchemeKey       IntKey0(intMethod);
                m_integrationScheme[0] = LibUtilities::TimeIntegrationSchemeManager()[IntKey0];
            }
            case LibUtilities::eIMEXOrder1: 
            {
                m_intSteps = 1;
                m_integrationScheme = Array<OneD, LibUtilities::TimeIntegrationSchemeSharedPtr> (m_intSteps);
                LibUtilities::TimeIntegrationSchemeKey       IntKey0(intMethod);
                m_integrationScheme[0] = LibUtilities::TimeIntegrationSchemeManager()[IntKey0];
            }
            break;
            case LibUtilities::eIMEXOrder2: 
            {
                m_intSteps = 2;
                m_integrationScheme = Array<OneD, LibUtilities::TimeIntegrationSchemeSharedPtr> (m_intSteps);
                LibUtilities::TimeIntegrationSchemeKey       IntKey0(LibUtilities::eIMEXOrder1);
                m_integrationScheme[0] = LibUtilities::TimeIntegrationSchemeManager()[IntKey0];
                LibUtilities::TimeIntegrationSchemeKey       IntKey1(intMethod);
                m_integrationScheme[1] = LibUtilities::TimeIntegrationSchemeManager()[IntKey1];
            }
            break;
            case LibUtilities::eIMEXOrder3: 
            {
                m_intSteps = 3;
                m_integrationScheme = Array<OneD, LibUtilities::TimeIntegrationSchemeSharedPtr> (m_intSteps);
                LibUtilities::TimeIntegrationSchemeKey       IntKey0(LibUtilities::eIMEXdirk_3_4_3);
                m_integrationScheme[0] = LibUtilities::TimeIntegrationSchemeManager()[IntKey0];
                LibUtilities::TimeIntegrationSchemeKey       IntKey1(LibUtilities::eIMEXdirk_3_4_3);
                m_integrationScheme[1] = LibUtilities::TimeIntegrationSchemeManager()[IntKey1];
                LibUtilities::TimeIntegrationSchemeKey       IntKey2(intMethod);
                m_integrationScheme[2] = LibUtilities::TimeIntegrationSchemeManager()[IntKey2];
            }
            break;
            default:
                ASSERTL0(0,"Integration method not suitable: Options include IMEXOrder1, IMEXOrder2 or IMEXOrder3");
            break;
        }        
        
        // set explicit time-intregration class operators
        m_integrationOps.DefineOdeRhs(&PorousMediaSplittingScheme::EvaluateAdvection_Permeability,this);

        // Count number of HBC conditions
        Array<OneD, const SpatialDomains::BoundaryConditionShPtr > PBndConds = m_pressure->GetBndConditions();
        Array<OneD, MultiRegions::ExpListSharedPtr>  PBndExp = m_pressure->GetBndCondExpansions();
        
        // Set up mapping from pressure boundary condition to pressure
        // element details.
        m_pressure->GetBoundaryToElmtMap(m_pressureBCtoElmtID,m_pressureBCtoTraceID);
        
        // Storage array for high order pressure BCs
        m_pressureHBCs = Array<OneD, Array<OneD, NekDouble> > (m_intSteps);
        
        int n,cnt;
        m_HBCnumber = 0;
        for(cnt = n = 0; n < PBndConds.num_elements(); ++n)
        {
            // High order boundary condition;
            if(PBndConds[n]->GetUserDefined() == SpatialDomains::eHigh)
            {
                cnt += PBndExp[n]->GetNcoeffs();
                m_HBCnumber += PBndExp[n]->GetExpSize();
            }
        }
        
        if (m_HBCnumber > 0) 
        {
            for(n = 0; n < m_intSteps; ++n)
            {
                m_pressureHBCs[n] = Array<OneD, NekDouble>(cnt);
            }
        }
        
        // set implicit time-intregration class operators
        m_integrationOps.DefineImplicitSolve(&PorousMediaSplittingScheme::SolveUnsteadyStokesSystem,this);
    }
    
    PorousMediaSplittingScheme::~PorousMediaSplittingScheme(void)
    {
        
    }
        
    
    void PorousMediaSplittingScheme::v_PrintSummary(std::ostream &out)
    {

        cout <<  "\tSolver Type     : Splitting Scheme" <<endl;
        
        
        if(m_session->DefinesSolverInfo("EvolutionOperator"))
        {
            cout << "\tEvolutionOp     : " << m_session->GetSolverInfo("EvolutionOperator")<< endl;
            
        }
        else
        {
            cout << "\tEvolutionOp     : " << endl;
        }
        if(m_session->DefinesSolverInfo("Driver"))
        {
            cout << "\tDriver          : " << m_session->GetSolverInfo("Driver")<< endl;
            
        }
        else{
            cout << "\tDriver          : "<< endl;
        }
        
        TimeParamSummary(out);
        cout << "\tTime integ.     : " << LibUtilities::TimeIntegrationMethodMap
            [m_integrationScheme[m_intSteps-1]->GetIntegrationMethod()] << endl;
    }
    
    void PorousMediaSplittingScheme::v_DoInitialise(void)
    {
        // Set initial condition using time t=0
        SetInitialConditions(0.0);
    }
    
    void PorousMediaSplittingScheme::v_DoSolve(void)
    {
        switch(m_equationType)
        {
            case eUnsteadyPorousMedia:
            {
                // Integrate from start time to end time
                AdvanceInTime(m_steps);
                break;
            }
            case eNoEquationType:
            default:
                ASSERTL0(false,"Unknown or undefined equation type for PorousMediaSplittingScheme solver");
        }
    }
    
    void PorousMediaSplittingScheme:: v_TransCoeffToPhys(void)
    {
        int nfields = m_fields.num_elements() - 1;
        for (int k=0 ; k < nfields; ++k)
        {
            //Backward Transformation in physical space for time evolution
            m_fields[k]->BwdTrans_IterPerExp(m_fields[k]->GetCoeffs(),
                                             m_fields[k]->UpdatePhys());
        }
    }
    
    void PorousMediaSplittingScheme:: v_TransPhysToCoeff(void)
    {
        
        int nfields = m_fields.num_elements() - 1;
        for (int k=0 ; k < nfields; ++k)
        {
            //Forward Transformation in physical space for time evolution
            m_fields[k]->FwdTrans_IterPerExp(m_fields[k]->GetPhys(),m_fields[k]->UpdateCoeffs());
        }
    }
    
    Array<OneD, bool> PorousMediaSplittingScheme::v_GetSystemSingularChecks()
    {
        int vVar = m_session->GetVariables().size();
        Array<OneD, bool> vChecks(vVar, false);
        vChecks[vVar-1] = true;
        return vChecks;
    }

    int PorousMediaSplittingScheme::v_GetForceDimension()
    {
        return m_session->GetVariables().size() - 1;
    }
    
    void PorousMediaSplittingScheme::EvaluateAdvection_Permeability(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray, 
        Array<OneD, Array<OneD, NekDouble> > &outarray, 
        const NekDouble time)
    {
        int nqtot        = m_fields[0]->GetTotPoints();
         
        // evaluate convection terms
        m_advObject->DoAdvection(m_fields, m_nConvectiveFields, m_velocity,inarray,outarray,m_time);

        if (m_explicitPermeability)
        {
            // add the force
            if(m_session->DefinesFunction("BodyForce"))
            {
                if(m_fields[0]->GetWaveSpace())
                {
                    for(int i = 0; i < m_nConvectiveFields; ++i)
                    {
                        m_forces[i]->SetWaveSpace(true);					
                        m_forces[i]->BwdTrans(m_forces[i]->GetCoeffs(),m_forces[i]->UpdatePhys());
                    }
                }
                for(int i = 0; i < m_nConvectiveFields; ++i)
                {
                    Vmath::Vadd(nqtot,outarray[i],1,(m_forces[i]->GetPhys()),1,outarray[i],1);
                }
            }
            // add permeability term for explicit permeability
            if (m_nConvectiveFields == 2)
            {
                if(m_session->DefinesFunction("SpatialAnisotropicPermeability"))
                {
                    Array <OneD, NekDouble> tmp(nqtot);
                    Vmath::Smul(nqtot,-m_kinvis,m_spatialperm[0],1,tmp,1);
                    Vmath::Vvtvp(nqtot,tmp,1,inarray[0],1,outarray[0],1,outarray[0],1);

                    //Get coordinate of quadrature points

                    Vmath::Smul(nqtot,-m_kinvis,m_spatialperm[1],1,tmp,1);
                    Vmath::Vvtvp(nqtot,tmp,1,inarray[1],1,outarray[1],1,outarray[1],1);
                }
                else
                {
                    Vmath::Svtvp(nqtot,-m_kinvis*m_perm_inv[0],inarray[0],1,outarray[0],1,outarray[0],1);
                    Vmath::Svtvp(nqtot,-m_kinvis*m_perm_inv[2],inarray[0],1,outarray[1],1,outarray[1],1);

                    Vmath::Svtvp(nqtot,-m_kinvis*m_perm_inv[2],inarray[1],1,outarray[0],1,outarray[0],1);
                    Vmath::Svtvp(nqtot,-m_kinvis*m_perm_inv[1],inarray[1],1,outarray[1],1,outarray[1],1);
                }
            }
            else
            {
                if(m_session->DefinesFunction("SpatialAnisotropicPermeability"))
                {
                    Array <OneD, NekDouble> tmp(nqtot);
                    Vmath::Smul(nqtot,-m_kinvis,m_spatialperm[0],1,tmp,1);
                    Vmath::Vvtvp(nqtot,tmp,1,inarray[0],1,outarray[0],1,outarray[0],1);

                    Vmath::Smul(nqtot,-m_kinvis,m_spatialperm[1],1,tmp,1);
                    Vmath::Vvtvp(nqtot,tmp,1,inarray[1],1,outarray[1],1,outarray[1],1);

                    Vmath::Smul(nqtot,-m_kinvis,m_spatialperm[2],1,tmp,1);
                    Vmath::Vvtvp(nqtot,tmp,1,inarray[2],1,outarray[2],1,outarray[2],1);
                }
                else
                {
                    Vmath::Svtvp(nqtot,-m_kinvis*m_perm_inv[0],inarray[0],1,outarray[0],1,outarray[0],1);
                    Vmath::Svtvp(nqtot,-m_kinvis*m_perm_inv[3],inarray[0],1,outarray[1],1,outarray[1],1);
                    Vmath::Svtvp(nqtot,-m_kinvis*m_perm_inv[4],inarray[0],1,outarray[2],1,outarray[2],1);
                    
                    Vmath::Svtvp(nqtot,-m_kinvis*m_perm_inv[3],inarray[1],1,outarray[0],1,outarray[0],1);
                    Vmath::Svtvp(nqtot,-m_kinvis*m_perm_inv[1],inarray[1],1,outarray[1],1,outarray[1],1);
                    Vmath::Svtvp(nqtot,-m_kinvis*m_perm_inv[5],inarray[1],1,outarray[2],1,outarray[2],1);
                    
                    Vmath::Svtvp(nqtot,-m_kinvis*m_perm_inv[4],inarray[2],1,outarray[0],1,outarray[0],1);
                    Vmath::Svtvp(nqtot,-m_kinvis*m_perm_inv[5],inarray[2],1,outarray[1],1,outarray[1],1);
                    Vmath::Svtvp(nqtot,-m_kinvis*m_perm_inv[2],inarray[2],1,outarray[2],1,outarray[2],1);
                }
            }
        }

        if(m_HBCnumber > 0)
        {
            // Set pressure BCs
            EvaluatePressureBCs(inarray, outarray);
        }
    }
    
    void PorousMediaSplittingScheme::SolveUnsteadyStokesSystem(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray, 
        Array<OneD, Array<OneD, NekDouble> > &outarray, 
        const NekDouble time, 
        const NekDouble aii_Dt)
    {
        int i,n;
        int phystot = m_fields[0]->GetTotPoints();
        int ncoeffs = m_fields[0]->GetNcoeffs();
        Array<OneD, Array< OneD, NekDouble> > F(m_nConvectiveFields);
        StdRegions::ConstFactorMap factors;
        factors[StdRegions::eFactorLambda] = 0.0;

        // Added here instead of in DefineOdeRhs
        // therefore not only IMEX time integration schemes can be used
        if( !m_explicitPermeability && (m_session->GetSolverInfo("AdvectionForm") == "NoAdvection") )
        {
            EvaluateAdvection_Permeability(inarray, outarray, time);
        }
                
        for(n = 0; n < m_nConvectiveFields; ++n)
        {
            F[n] = Array<OneD, NekDouble> (phystot);
        }
        
        SetBoundaryConditions(time);
        
        // Pressure Forcing = Divergence Velocity; 
        SetUpPressureForcing(inarray, F, aii_Dt);
                    
        // Solver Pressure Poisson Equation 
        m_pressure->HelmSolve(F[0], m_pressure->UpdateCoeffs(), NullFlagList, factors);
            
        // Viscous Term forcing
        SetUpViscousForcing(inarray, F, aii_Dt);

        if (m_explicitPermeability || m_session->DefinesFunction("SpatialAnisotropicPermeability"))
        {
            factors[StdRegions::eFactorLambda] = (1.0/aii_Dt/m_kinvis);
        }
        
        // Solve Helmholtz system and put in Physical space
        for(i = 0; i < m_nConvectiveFields; ++i)
        {
            // needs to be changed
            // is wrong for anisotropic permeability with off-diagonal entries! (isotropic works)
            if (!m_explicitPermeability)
            {
                factors[StdRegions::eFactorLambda] = (m_perm_inv[i])+(1.0/aii_Dt/m_kinvis);
            }    

            m_fields[i]->HelmSolve(F[i], m_fields[i]->UpdateCoeffs(), NullFlagList, factors, m_varperm);
            m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(),outarray[i]);
        }
    }
    
    void PorousMediaSplittingScheme::EvaluatePressureBCs(
        const Array<OneD, const Array<OneD, NekDouble> >  &fields, 
        const Array<OneD, const Array<OneD, NekDouble> >  &N, 
        const NekDouble Aii_Dt)
    {		
        Array<OneD, NekDouble> tmp;
        Array<OneD, const SpatialDomains::BoundaryConditionShPtr > PBndConds;
        Array<OneD, MultiRegions::ExpListSharedPtr>  PBndExp;
        int  n,cnt;
        int  nint    = min(m_pressureCalls++,m_intSteps);
        int  nlevels = m_pressureHBCs.num_elements();
        
        PBndConds   = m_pressure->GetBndConditions();
        PBndExp     = m_pressure->GetBndCondExpansions();
        
        // Reshuffle BC Storage vector
        tmp = m_pressureHBCs[nlevels-1];
        for(n = nlevels-1; n > 0; --n)
        {
            m_pressureHBCs[n] = m_pressureHBCs[n-1];
        }
        m_pressureHBCs[0] = tmp;
        
        // Calculate BCs at current level
        CalcPressureBCs(fields,N);
        
        // Copy High order values into storage array 
        for(cnt = n = 0; n < PBndConds.num_elements(); ++n)
        {
            // High order boundary condition;
            if(PBndConds[n]->GetUserDefined() == SpatialDomains::eHigh)
            {
                int nq = PBndExp[n]->GetNcoeffs();
                Vmath::Vcopy(nq,&(PBndExp[n]->GetCoeffs()[0]),1,&(m_pressureHBCs[0])[cnt],1);
                cnt += nq;
                }
        }
            
        // Extrapolate to n+1
        Vmath::Smul(cnt,kHighOrderBCsExtrapolation[nint-1][nint-1],
                    m_pressureHBCs[nint-1],1,m_pressureHBCs[nlevels-1],1);
        for(n = 0; n < nint-1; ++n)
        {
            Vmath::Svtvp(cnt,kHighOrderBCsExtrapolation[nint-1][n],
                         m_pressureHBCs[n],1,m_pressureHBCs[nlevels-1],1,
                         m_pressureHBCs[nlevels-1],1);
        }
        
        // Reset Values into Pressure BCs
        for(cnt = n = 0; n < PBndConds.num_elements(); ++n)
        {
            // High order boundary condition;
            if(PBndConds[n]->GetUserDefined() == SpatialDomains::eHigh)
            {
                int nq = PBndExp[n]->GetNcoeffs();
                Vmath::Vcopy(nq,&(m_pressureHBCs[nlevels-1])[cnt],1,&(PBndExp[n]->UpdateCoeffs()[0]),1);
                    cnt += nq;
            }
        }
    }
    
    void PorousMediaSplittingScheme::CalcPressureBCs(
        const Array<OneD, const Array<OneD, NekDouble> > &fields, 
        const Array<OneD, const Array<OneD, NekDouble> >  &N)
    {
        switch(m_velocity.num_elements())
        {
        case 1:
            ASSERTL0(false,"Porous media splitting scheme not designed to have just one velocity component");
            break;
        case 2:
            CalcPressureBCs2D(fields,N);
                break;
        case 3:
            CalcPressureBCs3D(fields,N);
            break;
        }
    }
    
    void PorousMediaSplittingScheme::CalcPressureBCs2D(
        const Array<OneD, const Array<OneD, NekDouble> > &fields, 
        const Array<OneD, const Array<OneD, NekDouble> >  &N)
    {
        int  i,n;
            
        Array<OneD, const SpatialDomains::BoundaryConditionShPtr > PBndConds;
        Array<OneD, MultiRegions::ExpListSharedPtr>  PBndExp;
        
        PBndConds = m_pressure->GetBndConditions();
        PBndExp   = m_pressure->GetBndCondExpansions();
        
        StdRegions::StdExpansionSharedPtr elmt;
        StdRegions::StdExpansion1DSharedPtr Pbc;
        
        Array<OneD, NekDouble> Pvals;
            
        int maxpts = 0,cnt;
        int elmtid,nq,offset, boundary;
        
        // find the maximum values of points 
        for(cnt = n = 0; n < PBndConds.num_elements(); ++n)
        {
            for(i = 0; i < PBndExp[n]->GetExpSize(); ++i)
            {
                maxpts = max(maxpts, m_fields[0]->GetExp(m_pressureBCtoElmtID[cnt++])->GetTotPoints());
            }
        }
            
        Array<OneD, const NekDouble> U,V,Nu,Nv;
        // Elemental tempspace: call once adn set others as offset 
        Array<OneD, NekDouble> Uy(5*maxpts); 
        Array<OneD, NekDouble> Vx = Uy + maxpts; 
        Array<OneD, NekDouble> Qx = Vx + maxpts;
        Array<OneD, NekDouble> Qy = Qx + maxpts; 
        Array<OneD, NekDouble> Q  = Qy + maxpts;
        
        for(cnt = n = 0; n < PBndConds.num_elements(); ++n)
        {            
            SpatialDomains::BndUserDefinedType type = PBndConds[n]->GetUserDefined(); 
                
            if(type == SpatialDomains::eHigh)
            {
                for(i = 0; i < PBndExp[n]->GetExpSize(); ++i,cnt++)
                {
                    // find element and edge of this expansion. 
                    // calculate curl x curl v;
                    elmtid = m_pressureBCtoElmtID[cnt];
                    elmt   = m_fields[0]->GetExp(elmtid);
                    nq     = elmt->GetTotPoints();
                    offset = m_fields[0]->GetPhys_Offset(elmtid);
                        
                    U = fields[m_velocity[0]] + offset;
                    V = fields[m_velocity[1]] + offset; 
                    
                    // Calculating vorticity Q = (dv/dx - du/dy)
                    elmt->PhysDeriv(MultiRegions::DirCartesianMap[0],V,Vx);
                    elmt->PhysDeriv(MultiRegions::DirCartesianMap[1],U,Uy);  
                    
                    Vmath::Vsub(nq,Vx,1,Uy,1,Q,1);
                        
                    // Calculate  Curl(Q) = Qy i - Qx j 
                    elmt->PhysDeriv(Q,Qx,Qy);
                    
                    Nu = N[0] + offset;
                    Nv = N[1] + offset; 

                    // Evaluate [N - kinvis Curlx Curl V].n
                    // x-component (stored in Qy)
                    //Vmath::Zero(nq,Qy,1);
                    Vmath::Svtvp(nq,-m_kinvis,Qy,1,Nu,1,Qy,1);
                        
                    // y-component (stored in Qx )
                    //Vmath::Zero(nq,Qx,1);
                    Vmath::Svtvp(nq,m_kinvis,Qx,1,Nv,1,Qx,1);
                    
                    if (!m_explicitPermeability)
                    {
                        Vmath::Svtvp(nq,-m_kinvis*m_perm_inv[0],U,1,Qy,1,Qy,1);
                        Vmath::Svtvp(nq,-m_kinvis*m_perm_inv[2],U,1,Qx,1,Qx,1);

                        Vmath::Svtvp(nq,-m_kinvis*m_perm_inv[2],V,1,Qy,1,Qy,1);
                        Vmath::Svtvp(nq,-m_kinvis*m_perm_inv[1],V,1,Qx,1,Qx,1);

                        if(m_session->DefinesFunction("BodyForce"))
                        {
                            if(m_fields[0]->GetWaveSpace())
                            {
                                for(int i = 0; i < m_nConvectiveFields; ++i)
                                {
                                    m_forces[i]->SetWaveSpace(true);					
                                    m_forces[i]->BwdTrans(m_forces[i]->GetCoeffs(),m_forces[i]->UpdatePhys());
                                }
                            }
                            Vmath::Vadd(nq,Qy,1,(m_forces[0]->GetPhys()),1,Qy,1);
                            Vmath::Vadd(nq,Qx,1,(m_forces[1]->GetPhys()),1,Qx,1);
                        }
                    }
                    
                    Pbc =  boost::dynamic_pointer_cast<StdRegions::StdExpansion1D> (PBndExp[n]->GetExp(i));
                        
                    boundary = m_pressureBCtoTraceID[cnt];
                    
                    // Get edge values and put into Uy, Vx
                    elmt->GetEdgePhysVals(boundary,Pbc,Qy,Uy);
                    elmt->GetEdgePhysVals(boundary,Pbc,Qx,Vx);

                    // calcuate (phi, dp/dn = [N-kinvis curl x curl v].n) 
                    Pvals = PBndExp[n]->UpdateCoeffs()+PBndExp[n]->GetCoeff_Offset(i);
                    Pbc->NormVectorIProductWRTBase(Uy,Vx,Pvals); 
                }
            }
                // setting if just standard BC no High order
            else if(type == SpatialDomains::eNoUserDefined || type == SpatialDomains::eTimeDependent) 
            {
                cnt += PBndExp[n]->GetExpSize();
            }
            else
            {
                ASSERTL0(false,"Unknown USERDEFINEDTYPE in pressure boundary condition");
            }
        }
    }
    
    void PorousMediaSplittingScheme::CalcPressureBCs3D(
            const Array<OneD, const Array<OneD, NekDouble> > &fields, 
            const Array<OneD, const Array<OneD, NekDouble> >  &N)
    {
        Array<OneD, const SpatialDomains::BoundaryConditionShPtr > PBndConds;
        Array<OneD, MultiRegions::ExpListSharedPtr>  PBndExp;
        
        PBndConds = m_pressure->GetBndConditions();
        PBndExp   = m_pressure->GetBndCondExpansions();
        
        int elmtid,nq,offset, boundary,cnt,n, i;
        int maxpts = 0;
        int phystot = m_fields[0]->GetTotPoints();
        
        Array<OneD, NekDouble> Pvals;
        
        // find the maximum values of points 
        for(cnt = n = 0; n < PBndConds.num_elements(); ++n)
        {
            for(int i = 0; i < PBndExp[n]->GetExpSize(); ++i)
            {
                maxpts = max(maxpts, m_fields[0]->GetExp(m_pressureBCtoElmtID[cnt++])->GetTotPoints());
            }
        }
        
        Array<OneD, const NekDouble> U,V,W,Nu,Nv,Nw;
        Array<OneD, NekDouble> Uy(maxpts);
        Array<OneD, NekDouble> Uz(maxpts);
        Array<OneD, NekDouble> Vx(maxpts);
        Array<OneD, NekDouble> Vz(maxpts);
        Array<OneD, NekDouble> Wx(maxpts);
        Array<OneD, NekDouble> Wy(maxpts);
        Array<OneD, NekDouble> Qx(maxpts);  
        Array<OneD, NekDouble> Qy(maxpts);
        Array<OneD, NekDouble> Qz(maxpts);
        
        StdRegions::StdExpansionSharedPtr elmt;
        StdRegions::StdExpansion2DSharedPtr Pbc;
        
        for(cnt = n = 0; n < PBndConds.num_elements(); ++n)
        {
            
            SpatialDomains::BndUserDefinedType type = PBndConds[n]->GetUserDefined();
                
            if(type == SpatialDomains::eHigh)
            {
                for(i = 0; i < PBndExp[n]->GetExpSize(); ++i,cnt++)
                {
                    // find element and face of this expansion. 
                    // calculate curl x curl v;
                    elmtid = m_pressureBCtoElmtID[cnt];
                    elmt   = m_fields[0]->GetExp(elmtid);
                    nq     = elmt->GetTotPoints();
                    offset = m_fields[0]->GetPhys_Offset(elmtid);
                    
                    U = fields[m_velocity[0]] + offset;
                    V = fields[m_velocity[1]] + offset; 
                    W = fields[m_velocity[2]] + offset;
                        
                    // Calculating vorticity Q = (dv/dx - du/dy)
                    elmt->PhysDeriv(MultiRegions::DirCartesianMap[1],U,Uy);
                    elmt->PhysDeriv(MultiRegions::DirCartesianMap[2],U,Uz);  
                    elmt->PhysDeriv(MultiRegions::DirCartesianMap[0],V,Vx);  
                    elmt->PhysDeriv(MultiRegions::DirCartesianMap[2],V,Vz);
                    elmt->PhysDeriv(MultiRegions::DirCartesianMap[0],W,Wx); 
                    elmt->PhysDeriv(MultiRegions::DirCartesianMap[1],W,Wy); 	
                    
                    Vmath::Vsub(nq,Wy,1,Vz,1,Qx,1);
                    Vmath::Vsub(nq,Uz,1,Wx,1,Qy,1);
                    Vmath::Vsub(nq,Vx,1,Uy,1,Qz,1);
                    
                    // Calculate  NxQ = Curl(Q) = (Qzy-Qyz) i + (Qxz-Qzx) j + (Qyx-Qxy) k
                    // NxQ = NxQ_x i + NxQ_y j + NxQ_z k
                    // Using the velocity derivatives memory space to
                    // store the vorticity derivatives.
                    // Qzy => Uy // Qyz => Uz // Qxz => Vx // Qzx => Vz // Qyx => Wx // Qxy => Wy 
                    elmt->PhysDeriv(MultiRegions::DirCartesianMap[1],Qz,Uy);
                    elmt->PhysDeriv(MultiRegions::DirCartesianMap[0],Qz,Vz);
                    elmt->PhysDeriv(MultiRegions::DirCartesianMap[2],Qy,Uz);
                    elmt->PhysDeriv(MultiRegions::DirCartesianMap[2],Qx,Vx);
                    elmt->PhysDeriv(MultiRegions::DirCartesianMap[0],Qy,Wx);
                    elmt->PhysDeriv(MultiRegions::DirCartesianMap[1],Qx,Wy);
                    
                    // Using the storage space associated with
                    // the 3 components of the vorticity to
                    // store the 3 components od the vorticity
                    // curl to save space Qx = Qzy-Qyz = Uy-Uz
                    // // Qy = Qxz-Qzx = Vx-Vz // Qz= Qyx-Qxy
                    // = Wx-Wy
                    Vmath::Vsub(nq,Uy,1,Uz,1,Qx,1);
                    Vmath::Vsub(nq,Vx,1,Vz,1,Qy,1);
                    Vmath::Vsub(nq,Wx,1,Wy,1,Qz,1);
                    
                    Nu = N[0] + offset;
                    Nv = N[1] + offset;
                    Nw = N[2] + offset;

                    // Evaluate [N - kinvis Curlx Curl V-k*V].n
                    // x-component (stored in Qx)
                    Vmath::Svtvp(nq,-m_kinvis,Qx,1,Nu,1,Qx,1);
                    // y-component (stored in Qy)
                    Vmath::Svtvp(nq,-m_kinvis,Qy,1,Nv,1,Qy,1);
                    // z-component (stored in Qz)
                    Vmath::Svtvp(nq,-m_kinvis,Qz,1,Nw,1,Qz,1);
                    
                    if (!m_explicitPermeability)
                    {
                        Vmath::Svtvp(nq,-m_kinvis*m_perm_inv[0],U,1,Qx,1,Qx,1);
                        Vmath::Svtvp(nq,-m_kinvis*m_perm_inv[3],U,1,Qy,1,Qy,1);
                        Vmath::Svtvp(nq,-m_kinvis*m_perm_inv[4],U,1,Qz,1,Qz,1);

                        Vmath::Svtvp(nq,-m_kinvis*m_perm_inv[3],V,1,Qx,1,Qx,1);
                        Vmath::Svtvp(nq,-m_kinvis*m_perm_inv[1],V,1,Qy,1,Qy,1);
                        Vmath::Svtvp(nq,-m_kinvis*m_perm_inv[5],V,1,Qz,1,Qz,1);

                        Vmath::Svtvp(nq,-m_kinvis*m_perm_inv[4],W,1,Qx,1,Qx,1);
                        Vmath::Svtvp(nq,-m_kinvis*m_perm_inv[5],W,1,Qy,1,Qy,1);
                        Vmath::Svtvp(nq,-m_kinvis*m_perm_inv[2],W,1,Qz,1,Qz,1);

                        if(m_session->DefinesFunction("BodyForce"))
                        {
                            if(m_fields[0]->GetWaveSpace())
                            {
                                for(int i = 0; i < m_nConvectiveFields; ++i)
                                {
                                    m_forces[i]->SetWaveSpace(true);					
                                    m_forces[i]->BwdTrans(m_forces[i]->GetCoeffs(),m_forces[i]->UpdatePhys());
                                }
                            }
                            Vmath::Vadd(nq,Qx,1,(m_forces[0]->GetPhys()),1,Qx,1);
                            Vmath::Vadd(nq,Qy,1,(m_forces[1]->GetPhys()),1,Qy,1);
                            Vmath::Vadd(nq,Qz,1,(m_forces[2]->GetPhys()),1,Qz,1);
                        }
                    }
                            
                    Pbc =  boost::dynamic_pointer_cast<StdRegions::StdExpansion2D> (PBndExp[n]->GetExp(i));
                    
                    boundary = m_pressureBCtoTraceID[cnt];
                    // Get face values and put into Uy, Vx and Wx
                    elmt->GetFacePhysVals(boundary,Pbc,Qx,Uy);
                    elmt->GetFacePhysVals(boundary,Pbc,Qy,Vx);
                    elmt->GetFacePhysVals(boundary,Pbc,Qz,Wx);
                    
                    // calcuate (phi, dp/dn = [N-kinvis curl x curl v].n) 
                    Pvals = PBndExp[n]->UpdateCoeffs()+PBndExp[n]->GetCoeff_Offset(i);
                    Pbc->NormVectorIProductWRTBase(Uy,Vx,Wx,Pvals); 
                }
            }
            // setting if just standard BC no High order
            else if(type == SpatialDomains::eNoUserDefined || type == SpatialDomains::eTimeDependent)
            {
                cnt += PBndExp[n]->GetExpSize();
            }
            else
            {
                    ASSERTL0(false,"Unknown USERDEFINEDTYPE in pressure boundary condition");
            }
        }
    }
    
        // Evaluate divergence of velocity field. 
    void PorousMediaSplittingScheme::SetUpPressureForcing(
        const Array<OneD, const Array<OneD, NekDouble> > &fields, 
        Array<OneD, Array<OneD, NekDouble> > &Forcing, 
        const NekDouble aii_Dt)
    {                
        int   i;
        int   physTot = m_fields[0]->GetTotPoints();
        Array<OneD, NekDouble> wk = Array<OneD, NekDouble>(physTot);
            
        Vmath::Zero(physTot,Forcing[0],1);

        for(i = 0; i < m_nConvectiveFields; ++i)
        {
            m_fields[i]->PhysDeriv(MultiRegions::DirCartesianMap[i],fields[i], wk);
            Vmath::Vadd(physTot,wk,1,Forcing[0],1,Forcing[0],1);
        }
        
        Vmath::Smul(physTot,1.0/aii_Dt,Forcing[0],1,Forcing[0],1);        
    }
    
    void PorousMediaSplittingScheme::SetUpViscousForcing(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray, 
        Array<OneD, Array<OneD, NekDouble> > &Forcing, 
        const NekDouble aii_Dt)
    {
        NekDouble aii_dtinv = 1.0/aii_Dt;
        int ncoeffs = m_fields[0]->GetNcoeffs();
        int phystot = m_fields[0]->GetTotPoints();
        
        // Grad p
#ifdef UseContCoeffs
        m_pressure->BwdTrans(m_pressure->GetContCoeffs(),m_pressure->UpdatePhys(),true);
            #else
        m_pressure->BwdTrans(m_pressure->GetCoeffs(),m_pressure->UpdatePhys());
#endif
        
        if(m_nConvectiveFields == 2)
        {
            m_pressure->PhysDeriv(m_pressure->GetPhys(), Forcing[0], Forcing[1]);
        }
        else
        {
            m_pressure->PhysDeriv(m_pressure->GetPhys(), Forcing[0], Forcing[1],Forcing[2]);
        }

        // Subtract inarray/(aii_dt) and divide by kinvis. Kinvis will
        // need to be updated for the convected fields.
        for(int i = 0; i < m_nConvectiveFields; ++i)
        {
            Vmath::Svtvp(phystot,-aii_dtinv,inarray[i],1,Forcing[i],1,Forcing[i],1);
            //Blas::Daxpy(phystot,-aii_dtinv,inarray[i],1,Forcing[i],1);
            //Blas::Dscal(phystot,1.0/m_kinvis,&(Forcing[i])[0],1);
        }

        if(!m_explicitPermeability)
        {
            // add the force
            if(m_session->DefinesFunction("BodyForce"))
            {
                if(m_fields[0]->GetWaveSpace())
                {
                    for(int i = 0; i < m_nConvectiveFields; ++i)
                    {
                        m_forces[i]->SetWaveSpace(true);					
                        m_forces[i]->BwdTrans(m_forces[i]->GetCoeffs(),m_forces[i]->UpdatePhys());
                    }
                }
                for(int i = 0; i < m_nConvectiveFields; ++i)
                {
                    Vmath::Vsub(phystot,Forcing[i],1,(m_forces[i]->GetPhys()),1,Forcing[i],1);
                }
            }
        }
        
        for(int i = 0; i < m_nConvectiveFields; ++i)
        {
            Vmath::Smul(phystot,1.0/m_kinvis,Forcing[i],1,Forcing[i],1);
        }

        // if (!m_explicitPermeability && m_session->DefinesFunction("AnisotropicPermeability"))
        // {
        //     Array<OneD, Array< OneD, NekDouble> > tempF(m_nConvectiveFields);
        //     for (int i = 0; i < m_nConvectiveFields; ++i)
        //     {
        //         tempF[i] = Array<OneD, NekDouble> (phystot);
        //     }

        //     if(m_nConvectiveFields == 2)
        //     {
        //         Vmath::Smul(phystot,m_perm[0],Forcing[0],1,tempF[0],1);
        //         Vmath::Smul(phystot,m_perm[2],Forcing[0],1,tempF[1],1);
                
        //         Vmath::Svtvp(phystot,m_perm[2],Forcing[1],1,tempF[0],1,tempF[0],1);
        //         Vmath::Svtvp(phystot,m_perm[1],Forcing[1],1,tempF[1],1,tempF[1],1);
        //     }
        //     else
        //     {
        //         Vmath::Smul(phystot,m_perm[0],Forcing[0],1,tempF[0],1);
        //         Vmath::Smul(phystot,m_perm[3],Forcing[0],1,tempF[1],1);
        //         Vmath::Smul(phystot,m_perm[4],Forcing[0],1,tempF[2],1);

        //         Vmath::Svtvp(phystot,m_perm[3],Forcing[1],1,tempF[0],1,tempF[0],1);
        //         Vmath::Svtvp(phystot,m_perm[1],Forcing[1],1,tempF[1],1,tempF[1],1);
        //         Vmath::Svtvp(phystot,m_perm[5],Forcing[1],1,tempF[2],1,tempF[2],1);

        //         Vmath::Svtvp(phystot,m_perm[4],Forcing[2],1,tempF[0],1,tempF[0],1);
        //         Vmath::Svtvp(phystot,m_perm[5],Forcing[2],1,tempF[1],1,tempF[1],1);
        //         Vmath::Svtvp(phystot,m_perm[2],Forcing[2],1,tempF[2],1,tempF[2],1);
        //     }

        //     for (int i = 0; i < m_nConvectiveFields; ++i)
        //     {
        //         Forcing[i] = tempF[i];
        //     }
        // }
    }
        
} //end of namespace
