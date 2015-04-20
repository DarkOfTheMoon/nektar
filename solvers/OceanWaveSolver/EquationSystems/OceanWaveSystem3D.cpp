///////////////////////////////////////////////////////////////////////////////
//
// File OceanWaveSystem3D.cpp
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
// Description: Unsteady diffusion solve routines
//
///////////////////////////////////////////////////////////////////////////////

#include <OceanWaveSolver/EquationSystems/OceanWaveSystem3D.h>
#include <iostream>
#include <iomanip>

namespace Nektar
{
    string OceanWaveSystem3D::className = GetEquationSystemFactory().RegisterCreatorFunction("OceanWaveSystem3D", OceanWaveSystem3D::create);

    OceanWaveSystem3D::OceanWaveSystem3D(
            const LibUtilities::SessionReaderSharedPtr& pSession)
        : UnsteadySystem(pSession)
    {
        m_rank = pSession->GetComm()->GetRank();
    }

    /**
     * @brief Initialisation object for the unsteady diffusion problem.
     */
    void OceanWaveSystem3D::v_InitObject()
    {
        UnsteadySystem::v_InitObject();

        //m_session->LoadParameter("wavefreq",   m_waveFreq, 0.0);
        //m_session->LoadParameter("epsilon",    m_epsilon,  0.0);

        //m_session->MatchSolverInfo(
        //    "SpectralVanishingViscosity", "True", m_useSpecVanVisc, false);

        if(m_useSpecVanVisc)
        {
            //m_session->LoadParameter("SVVCutoffRatio",m_sVVCutoffRatio,0.75);
            //m_session->LoadParameter("SVVDiffCoeff",m_sVVDiffCoeff,0.1);
        }


        std::cout << "hello" << std::endl;

        //load 3d graph:
        std::string graph3DName;
        m_session->LoadSolverInfo("Graph3D", graph3DName, "");

        std::cout << graph3DName << std::endl;

        ASSERTL0(graph3DName.size()>0, "No 3d graph file provided...");

        char graph3D_text[graph3DName.size()+1];
        graph3DName.copy(graph3D_text, graph3DName.size());
        graph3D_text[graph3DName.size()] = '\0';

        std::cout << "hello" << std::endl;

        int tmpNargs = 2;
        char * tmpPtr[2];
        tmpPtr[0] =  graph3D_text;
        tmpPtr[1] =  graph3D_text;

        LibUtilities::SessionReaderSharedPtr session3D = LibUtilities::SessionReader::CreateInstance(tmpNargs, tmpPtr);
        solver3D = EqSys::create(session3D);

        int npoints = m_fields[0]->GetPhys().num_elements();
        

        Array<OneD, NekDouble> dirs_x(npoints);
        Array<OneD, NekDouble> dirs_y(npoints);
        Array<OneD, NekDouble> dirs_z(npoints);
        m_fields[0]->GetCoords(dirs_x, dirs_y, dirs_z);

        std::cout << "hello" << std::endl;


        //init connection: 2D->3D
        boost::static_pointer_cast<EqSys>(solver3D)->InitConnection(dirs_x, dirs_y);

        std::cout << "hello" << std::endl;



        //SHOULD BE CHANGED:
        m_dt = m_timestep;

        m_Twave = 1;
        m_hd = 1;
        m_wwave = 2*M_PI/m_Twave;
        m_g = 9.81;
        m_Hwave = 0.02;
        m_kh = 4.0269; //todo, change to use disper.
        m_kwave = m_kh/m_hd;
        m_lwave = 2*M_PI/m_kwave;
        m_cwave = m_wwave/m_kwave;
        m_h0 = 1;
        
        


        //should be replaced with masks:
        

        bool has_iz2 = false;
        bool has_iz3 = false;

        //create iZ2, iZ3:
        NekDouble xmin = Vmath::Vmin(npoints, dirs_x, 1);
        NekDouble xmax = Vmath::Vmax(npoints, dirs_x, 1);
        int niZ2 = 0, niZ3 = 0;
        for (int i = 0; i < npoints; ++i)
        {
            if (dirs_x[i] <= -0.5) //xmin + 2*m_lwave)
            {
                //creating zone
                if (dirs_y[i] < 0)
                {
                    niZ2++;
                    has_iz2 = true;
                }

                //absorb zone
                if (dirs_y[i] > 0)
                {
                    niZ3++;
                    has_iz3 = true;
                }
                
            }
        }



        
        if (has_iz2)
            iZ2 = Array<OneD, int>(niZ2);

        if (has_iz3)
            iZ3 = Array<OneD, int>(niZ3);

        niZ2 = 0; 
        niZ3 = 0;
        for (int i = 0; i < npoints; ++i)
        {
            if (dirs_x[i] <= -0.5) //xmin + 2*m_lwave)
            {
                //creating zone
                if (dirs_y[i] < 0)
                {
                    iZ2[niZ2] = i;
                    niZ2++;
                }

                //absorb zone
                if (dirs_y[i] > 0)
                {
                    iZ3[niZ3] = i;
                    niZ3++;
                }
                
            }
        }
        if (has_iz2)
        {
            crm = Array<OneD, NekDouble>(niZ2);
            XiZ2 = Array<OneD, NekDouble>(niZ2);
            Vmath::Gathr(niZ2, dirs_x, iZ2, XiZ2);
            spongefunction4(XiZ2, 0, 0, 10, crm);

            m_XiZ2 = Array<OneD, NekDouble>(npoints);
            Array<OneD, NekDouble> Y(npoints);
            Array<OneD, NekDouble> Z(npoints);
            m_fields[0]->GetCoords(m_XiZ2,Y,Z);

            NekDouble xmaxZ2 = -0.5;

            Vmath::Sadd(npoints, xmaxZ2, m_XiZ2, 1, m_XiZ2, 1);
        }
        if (has_iz3)
        {
            cfr = Array<OneD, NekDouble>(niZ3);
            Array<OneD, NekDouble> XiZ3(niZ3);
            Vmath::Gathr(niZ3, dirs_x, iZ3, XiZ3);
            spongefunction4(XiZ3, 0, 5.0, 9, cfr);
        }
        

  

        m_homoInitialFwd = false;

        


        switch (m_projectionType)
        {
            case MultiRegions::eDiscontinuous:
            {
                std::string diffName;

                // Do not forwards transform initial condition
                m_homoInitialFwd = false;

                //std::cout << "END" << std::endl;

                m_session->LoadSolverInfo("DiffusionType", diffName, "LDG");
                //m_diffusion = SolverUtils::GetDiffusionFactory().
                //    CreateInstance(diffName, diffName);
                //m_diffusion->SetFluxVector(&OceanWaveSystem3D::
                //                           GetFluxVector, this);
                //std::cout << "END" << std::endl;
//
                //m_diffusion->InitObject(m_session, m_fields);
                break;
            }
        
            case MultiRegions::eGalerkin:
            case MultiRegions::eMixed_CG_Discontinuous:
            {
                // In case of Galerkin explicit diffusion gives an error
                if (m_explicitDiffusion)
                {
                    //ASSERTL0(false, "Explicit Galerkin diffusion not set up.");
                }
                // In case of Galerkin implicit diffusion: do nothing
            }
        }
        
        
        //if (m_explicitDiffusion)
        //{
            m_ode.DefineOdeRhs    (&OceanWaveSystem3D::DoOdeRhs,        this);
            m_ode.DefineProjection(&OceanWaveSystem3D::DoOdeProjection, this);

        //}
        //else
        //{
        //    m_ode.DefineImplicitSolve(
        //                        &OceanWaveSystem3D::DoImplicitSolve, this);
        //}
            std::cout << "final" << std::endl;
    }

    /**
     * @brief Unsteady diffusion problem destructor.
     */
    OceanWaveSystem3D::~OceanWaveSystem3D()
    {
    }

    void OceanWaveSystem3D::v_GenerateSummary(SummaryList& s)
    {
        UnsteadySystem::v_GenerateSummary(s);
    }
    
    
    /* @brief Compute the right-hand side for the unsteady diffusion problem.
     * 
     * @param inarray    Given fields.
     * @param outarray   Calculated solution.
     * @param time       Time.
     */

     void OceanWaveSystem3D::SetBoundaryConditions(
    Array<OneD, Array<OneD, NekDouble> > &inarray,
    NekDouble time)
  { 
  }


     void OceanWaveSystem3D::WallBoundary(
        int                                   bcRegion,
        int                                   cnt, 
        Array<OneD, Array<OneD, NekDouble> > &physarray)
    { 
        
    }



    void OceanWaveSystem3D::DoOdeRhs(
        const Array<OneD, const  Array<OneD, NekDouble> > &inarray,
              Array<OneD,        Array<OneD, NekDouble> > &outarray,
        const NekDouble time)
    {

        int npoints = inarray[0].num_elements();
        Vmath::Zero(outarray[0].num_elements(), outarray[0], 1);
        Vmath::Zero(outarray[1].num_elements(), outarray[1], 1);

        m_fields[0]->FwdTrans(inarray[0], m_fields[0]->UpdateCoeffs());
        boost::static_pointer_cast<EqSys>(solver3D)->Solve(m_fields[0]->GetCoeffs(), outarray[1]);

        Vmath::Vcopy(npoints, inarray[1], 1, outarray[0], 1);
        Vmath::Smul(npoints, -m_g, outarray[0], 1, outarray[0], 1);

        WaveForcing(time, inarray[1], inarray[0], outarray[1], outarray[0]);
    }


      void OceanWaveSystem3D::WallBoundary2D(int bcRegion, int cnt, Array<OneD, Array<OneD, NekDouble> > &physarray)
      { 
      }

    /**
     * @brief Compute the projection for the unsteady diffusion problem.
     * 
     * @param inarray    Given fields.
     * @param outarray   Calculated solution.
     * @param time       Time.
     */
    void OceanWaveSystem3D::DoOdeProjection(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
              Array<OneD,       Array<OneD, NekDouble> > &outarray,
        const NekDouble time)
    {

        int i;
        int nvariables = inarray.num_elements();

        //SetBoundaryConditions(outarray, time);


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
                //SetBoundaryConditions(outarray, time);
                break;
            }
            case MultiRegions::eGalerkin:
            case MultiRegions::eMixed_CG_Discontinuous:
            {
                Array<OneD, NekDouble> coeffs(m_fields[0]->GetNcoeffs());
                int npoints = GetNpoints();
                for(i = 0; i < nvariables; ++i)
                {  
                    m_fields[i]->FwdTrans(inarray[i], coeffs);
                    m_fields[i]->BwdTrans_IterPerExp(coeffs, outarray[i]);
                    //Vmath::Vcopy(npoints, inarray[i], 1, outarray[i], 1);
                }
                break;
            }
            default:
            {
                ASSERTL0(false, "Unknown projection scheme");
                break;
            }
        }

        //std::cout << "max: " << Vmath::Vmax(tmpPhi.num_elements(), tmpPhi, 1) << " min: " << Vmath::Vmin(tmpPhi.num_elements(), tmpPhi, 1) << std::endl;

    }
    
    /** 
     * @brief Implicit solution of the unsteady diffusion problem.
     */
    void OceanWaveSystem3D::DoImplicitSolve(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
              Array<OneD,       Array<OneD, NekDouble> > &outarray,
        const NekDouble time,
        const NekDouble lambda)
    {

    }
    
    /** 
     * @brief Return the flux vector for the unsteady diffusion problem.
     */
    void OceanWaveSystem3D::GetFluxVector(
        const int i, 
        const int j,
        const Array<OneD, Array<OneD, NekDouble> > &physfield,
              Array<OneD, Array<OneD, NekDouble> > &derivatives,
              Array<OneD, Array<OneD, NekDouble> > &flux)
    {

    }

    void OceanWaveSystem3D::spongefunction4(const Array<OneD, NekDouble> & x, NekDouble alpha, NekDouble p, int type, Array<OneD, NekDouble> & cr)
    {
        int n = x.num_elements();
        NekDouble xmin = -5; //Vmath::Vmin(n, x, 1);
        NekDouble xmax = -0.5; //Vmath::Vmax(n, x, 1);
        
        Array<OneD, NekDouble> x2(n);
        Vmath::Sadd(n, -xmin, x, 1, x2, 1);
        Vmath::Smul(n, 1.0/(xmax-xmin), x2, 1, x2, 1);


        switch(type)
        {
            case 9:
                Vmath::Vpow(n, x2, 1, p, cr, 1);
                Vmath::Neg(n, cr, 1);
                Vmath::Sadd(n, 1.0, cr, 1, cr, 1);
                break;
            case 10:
                Vmath::Vpow(n, x2, 1, 2.0, cr, 1);
                Vmath::Smul(n, 3.0, cr, 1, cr, 1);
                Vmath::Vpow(n, x2, 1, 3.0, x2, 1);
                Vmath::Smul(n, -2.0, x2, 1, x2, 1);
                Vmath::Vadd(n, cr, 1, x2, 1, cr, 1);
                break;
            default:
                ASSERTL0(false, "unknown sponge function...");
                break;

        }
        

    }

    //void OceanWaveSystem3D::disper(NekDouble g, NekDouble T, NekDouble h, NekDouble)

    /*

        function [kh] = disper(g,T,h,n)
        % determine the dimensionless dispersion parameter
        % for wave equations using the dispersion relation.
        w = 2*pi/T;
        kh = w^2*h/g;
        kh = khsolve(g,w,kh,h,n);

        function kh = khsolve(g,w,kh,h,n)
        % to be used within the function disper for iterative solution of dispersion relation
        for i = 1:n
            kh = sqrt(w^2/g*h.*kh.*coth(kh));
        end


    */

    
    void OceanWaveSystem3D::lineartravellingwave1D(NekDouble H,NekDouble c, NekDouble k, NekDouble z, NekDouble h, NekDouble w, NekDouble t, const Array<OneD, NekDouble> & x, Array<OneD, NekDouble> & eta, Array<OneD, NekDouble> & pp )
    {
        int n = x.num_elements();

        Array<OneD, NekDouble> xtime(n);
        Vmath::Smul(n, -k, x, 1, xtime, 1);
        Vmath::Sadd(n, w*t, xtime, 1, xtime, 1);


        Array<OneD, NekDouble> xcos(n);
        Array<OneD, NekDouble> xsin(n);
        Vmath::Vcos(n, xtime, 1, xcos, 1);
        Vmath::Vsin(n, xtime, 1, xsin, 1);

        //eta   =    H/2*cos(w*t-k*x);
        Vmath::Smul(n, H/2.0, xcos, 1, eta, 1);

        //pp    = -H*c/2*cosh(k*(z+h))./sinh(k*h).*sin(w*t-k*x);
        Vmath::Smul(n, -H*c/2.0*cosh(k*(z+h))/sinh(k*h), xsin, 1, pp, 1);
    }

    void OceanWaveSystem3D::WaveForcing(NekDouble time, const Array<OneD, NekDouble> & E, const Array<OneD, NekDouble> & P, Array<OneD, NekDouble> & rhsE, Array<OneD, NekDouble> & rhsP)
    {        

        int n = m_fields[0]->GetTotPoints();
        int niZ2 = iZ2.num_elements();
        int niZ3 = iZ3.num_elements();
        

        
        Array<OneD, NekDouble> rhsP2(n, 0.0);
        Array<OneD, NekDouble> rhsE2(n, 0.0);
        

        //secure with multithreading:
        if (iZ2.num_elements() > 0)
        {
            Array<OneD, NekDouble> E2(n,0.0);
            Array<OneD, NekDouble> P2(n,0.0);
            lineartravellingwave1D(m_Hwave, m_cwave, m_kwave, 0, m_hd, m_wwave, time, m_XiZ2, E2, P2);

            if (time < 5*m_Twave)
            {
                Vmath::Smul(n, time/(5*m_Twave), E2, 1, E2, 1);
                Vmath::Smul(n, time/(5*m_Twave), P2, 1, P2, 1);
            }

            //Create submatrices:
            Array<OneD, NekDouble> EiZ2(niZ2);
            Vmath::Gathr(niZ2, E, iZ2, EiZ2);
            Array<OneD, NekDouble> E2iZ2(niZ2);
            Vmath::Gathr(niZ2, E2, iZ2, E2iZ2);
            Array<OneD, NekDouble> PiZ2(niZ2);
            Vmath::Gathr(niZ2, P, iZ2, PiZ2);
            Array<OneD, NekDouble> P2iZ2(niZ2);
            Vmath::Gathr(niZ2, P2, iZ2, P2iZ2);

            Array<OneD, NekDouble> crm_neg(niZ2);
            Vmath::Sadd(niZ2, -1.0, crm, 1, crm_neg, 1);
            Vmath::Neg(niZ2, crm_neg, 1);

            //SPAWNING
            //modify rhsE:
            Vmath::Vsub(niZ2, E2iZ2, 1, EiZ2, 1, EiZ2, 1);
            Vmath::Vmul(niZ2, EiZ2, 1, crm_neg, 1, EiZ2, 1);
            Vmath::Smul(niZ2, 1.0/m_dt, EiZ2, 1, EiZ2, 1);
            Vmath::Scatr(niZ2, EiZ2, iZ2, rhsE2);
        
            //modify rhsP:
            Vmath::Vsub(niZ2, P2iZ2, 1, PiZ2, 1, PiZ2, 1);
            Vmath::Vmul(niZ2, PiZ2, 1, crm_neg, 1, PiZ2, 1);
            Vmath::Smul(niZ2, 1.0/m_dt, PiZ2, 1, PiZ2, 1);
            Vmath::Scatr(niZ2, PiZ2, iZ2, rhsP2);
        }

        if (iZ3.num_elements() > 0)
        {
            return;
            Array<OneD, NekDouble> EiZ3(niZ3,0.0);
            Vmath::Gathr(niZ3, E, iZ3, EiZ3);
            Array<OneD, NekDouble> PiZ3(niZ3,0.0);
            Vmath::Gathr(niZ3, P, iZ3, PiZ3);

            Array<OneD, NekDouble> cfr_neg(niZ3);
            Vmath::Sadd(niZ3, -1.0, cfr, 1, cfr_neg, 1);
            Vmath::Neg(niZ3, cfr_neg, 1);

            //DESPAWNING
            //modify rhsE:
            
            Vmath::Vmul(niZ3, EiZ3, 1, cfr_neg, 1, EiZ3, 1);
            Vmath::Smul(niZ3, -0.2/m_dt, EiZ3, 1, EiZ3, 1);
            Vmath::Scatr(niZ3, EiZ3, iZ3, rhsE2); 
            
            //modify rhsP:
            Vmath::Vmul(niZ3, PiZ3, 1, cfr_neg, 1, PiZ3, 1);
            Vmath::Smul(niZ3, -0.2/m_dt, PiZ3, 1, PiZ3, 1);
            Vmath::Scatr(niZ3, PiZ3, iZ3, rhsP2);
            
        }

        //final:
        Vmath::Vadd(n, rhsE, 1, rhsE2, 1, rhsE, 1);
        Vmath::Vadd(n, rhsP, 1, rhsP2, 1, rhsP, 1);

        
    }
    






        string EqSys::className = GetEquationSystemFactory().
        RegisterCreatorFunction("EqSys", EqSys::create);

    EqSys::EqSys(
            const LibUtilities::SessionReaderSharedPtr& pSession)
        : EquationSystem(pSession),
          m_factors()
    {
        m_factors[StdRegions::eFactorLambda] = 0.0;
        m_factors[StdRegions::eFactorTau] = 1.0;
    }

    void EqSys::v_InitObject()
    {
        EquationSystem::v_InitObject();
    }

    EqSys::~EqSys()
    {

    }

    void EqSys::v_GenerateSummary(SolverUtils::SummaryList& s)
    {
        EquationSystem::SessionSummary(s);
        SolverUtils::AddSummaryItem(s, "Lambda",
                                    m_factors[StdRegions::eFactorLambda]);
    }

    void EqSys::v_DoSolve()
    {
        for(int i = 0; i < m_fields.num_elements(); ++i)
        {
            // Zero field so initial conditions are zero
            Vmath::Zero(m_fields[i]->GetNcoeffs(),
                        m_fields[i]->UpdateCoeffs(), 1);
            m_fields[i]->HelmSolve(m_fields[i]->GetPhys(),
                                   m_fields[i]->UpdateCoeffs(),
                                   NullFlagList,
                                   m_factors);
            m_fields[i]->SetPhysState(false);
        }
    }

    Array<OneD, bool> EqSys::v_GetSystemSingularChecks()
    {
        return Array<OneD, bool>(m_session->GetVariables().size(), true);
    }

    void EqSys::Solve(const Array<OneD, NekDouble> & top, Array<OneD, NekDouble> & W)
    {
        int npoints = top.num_elements();
        for(int i = 0; i < m_fields.num_elements(); ++i)
        {
            //Update boundaries for 3d problem:

            //Vmath::Gathr(m_PhiToBnd.num_elements(), top, m_PhiToBnd, m_fields[i]->GetBndCondExpansions()[0]->UpdatePhys());
            m_fields[i]->GetBndCondExpansions()[0]->BwdTrans(top, m_fields[i]->GetBndCondExpansions()[0]->UpdatePhys());
            m_fields[i]->GetBndCondExpansions()[0]->FwdTrans(m_fields[i]->GetBndCondExpansions()[0]->GetPhys(), m_fields[i]->GetBndCondExpansions()[0]->UpdateCoeffs());

            // Zero field so initial conditions are zero
            Vmath::Zero(m_fields[i]->GetNcoeffs(),
                        m_fields[i]->UpdateCoeffs(), 1);
            m_fields[i]->HelmSolve(m_fields[i]->GetPhys(),
                                   m_fields[i]->UpdateCoeffs(),
                                   NullFlagList,
                                   m_factors);

            m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(),
                                              m_fields[i]->UpdatePhys());


            int n = m_fields[i]->GetPhys().num_elements();


            //Find derivative of function.
            Array<OneD, NekDouble> out_d2(n);
            m_fields[i]->PhysDeriv(2, m_fields[i]->GetPhys(), out_d2);

            //gather:
            Vmath::Gathr(m_PhiTo2D.num_elements(), out_d2, m_PhiTo2D, W);
        }
    }

    void EqSys::InitConnection(const Array<OneD, NekDouble> & X, const Array<OneD, NekDouble> & Y)
    {
        //int npoints  = m_fields[0]->GetPhys().num_elements();
        int npoints2 = X.num_elements();

        m_PhiTo2D = Array<OneD, int>(npoints2, 0);


        std::cout << m_fields[0]->GetBndCondExpansions()[0]->GetNcoeffs() << " " << m_fields[0]->GetBndCondExpansions()[0]->GetNpoints() << std::endl;
        std::cout << X.num_elements() << std::endl;


        int count3 = 0;
        int e, id1, id2;

        int npts = m_fields[0]->GetTotPoints();
        dirs_x = Array<OneD, NekDouble>(npts);
        dirs_y = Array<OneD, NekDouble>(npts);
        dirs_z = Array<OneD, NekDouble>(npts);
        m_fields[0]->GetCoords(dirs_x, dirs_y, dirs_z);

        int count = 0;
        for (int i = 0; i < npts; ++i)
        {
            if (dirs_z[i] > -0.01)
            {
                m_PhiTo2D[count] = i;
                count++;
            }
            
        }
        std::cout << count  << std::endl;


        int npts3 = m_fields[0]->GetBndCondExpansions()[0]->GetTotPoints();
        Array<OneD, NekDouble> dirs_xx(npts3);
        Array<OneD, NekDouble> dirs_yy(npts3);
        Array<OneD, NekDouble> dirs_zz(npts3);
        m_fields[0]->GetBndCondExpansions()[0]->GetCoords(dirs_xx, dirs_yy, dirs_zz);

        m_PhiToBnd = Array<OneD, int>(npts3, 0);

        count = 0;
        for (int i = 0; i < npts3; ++i)
        {
            for (int j = 0; j < X.num_elements(); ++j)
            {
                if ( fabs(dirs_xx[i]-X[j]) < 0.000001 && fabs(dirs_yy[i]-Y[j]) < 0.000001 )
                {
                    m_PhiToBnd[i] = j;
                    count++;
                    break;
                }
            }
        }
        std::cout << count  << std::endl;

        //ASSERTL0(false, "unknown sponge function...");
    }
    

}
