///////////////////////////////////////////////////////////////////////////////
//
// File OceanWaveSystem.cpp
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

#include <OceanWaveSolver/EquationSystems/OceanWaveSystem.h>
#include <iostream>
#include <iomanip>

namespace Nektar
{
    string OceanWaveSystem::className = GetEquationSystemFactory().RegisterCreatorFunction("OceanWaveSystem", OceanWaveSystem::create);

    OceanWaveSystem::OceanWaveSystem(
            const LibUtilities::SessionReaderSharedPtr& pSession)
        : UnsteadySystem(pSession)
    {
    }

    /**
     * @brief Initialisation object for the unsteady diffusion problem.
     */
    void OceanWaveSystem::v_InitObject()
    {
        UnsteadySystem::v_InitObject();

        m_session->LoadParameter("wavefreq",   m_waveFreq, 0.0);
        m_session->LoadParameter("epsilon",    m_epsilon,  0.0);

        m_session->MatchSolverInfo(
            "SpectralVanishingViscosity", "True", m_useSpecVanVisc, false);

        if(m_useSpecVanVisc)
        {
            m_session->LoadParameter("SVVCutoffRatio",m_sVVCutoffRatio,0.75);
            m_session->LoadParameter("SVVDiffCoeff",m_sVVDiffCoeff,0.1);
        }


        int npoints = m_fields[0]->GetTotPoints();
        int npoints2 = m_fields[0]->GetBndCondExpansions()[0]->GetTotPoints();

        m_PhiTo2D = Array<OneD, int>(npoints);
        m_PhiToTop = Array<OneD, int>(npoints);
        m_PhiTo3D = Array<OneD, int>(npoints2);


        //tmpPhi:
        tmpPhi = Array<OneD, NekDouble>(npoints);
        tmpRhsE = Array<OneD, NekDouble>(npoints);
        tmpRhsP = Array<OneD, NekDouble>(npoints);
        
        Array<OneD, NekDouble> bnd_x(npoints2);
        Array<OneD, NekDouble> bnd_y(npoints2);
        Array<OneD, NekDouble> bnd_z(npoints2);
        m_fields[0]->GetBndCondExpansions()[0]->GetCoords(bnd_x, bnd_y, bnd_z);

        Array<OneD, NekDouble> dirs_x(npoints);
        Array<OneD, NekDouble> dirs_y(npoints);
        Array<OneD, NekDouble> dirs_z(npoints);
        m_fields[0]->GetCoords(dirs_x, dirs_y, dirs_z);

        int nElms = m_fields[0]->GetNumElmts();
        int index = 0;
        int bndIndex = 0;
        int supercount = 0;
        for (int i = 0; i < nElms; ++i)
        {
            int nElmsSize = m_fields[0]->GetTotPoints(i);
            int nElmsSize2 = m_fields[0]->GetBndCondExpansions()[0]->GetTotPoints(i);

            int ratio = nElmsSize / nElmsSize2;
            for (int ki = 0; ki < nElmsSize2; ++ki)
            {
                int ki_index = bndIndex + ki;
                m_PhiTo3D[ki_index] = -1;

                float max_dist = 1000;
                for (int kj = 0; kj < nElmsSize; ++kj)
                {
                    
                    int kj_index = index + kj;

                    float dis_x = fabs(bnd_x[ki_index]-dirs_x[kj_index]); 
                    float dis_y = fabs(bnd_y[ki_index]-dirs_y[kj_index]);      

                    if (dirs_z[kj_index] < 0.8 && dirs_z[kj_index] > 0.2 && dis_x + dis_y < max_dist)
                    {
                        m_PhiTo3D[ki_index] = kj_index;
                        max_dist = dis_x + dis_y;
                    }
                }
                if (m_PhiTo3D[ki_index != -1])
                {
                    supercount++;
                }
            }
            for (int ki = 0; ki < nElmsSize; ++ki)
            {
                for (int kj = 0; kj < nElmsSize; ++kj)
                {
                    int ki_index = index + ki;
                    int kj_index = index + kj;

                    if (dirs_z[kj_index] < 0.8 && dirs_z[kj_index] > 0.2 && fabs(dirs_x[ki_index]-dirs_x[kj_index])<0.0001 && fabs(dirs_y[ki_index]-dirs_y[kj_index])<0.0001)
                    {
                        m_PhiTo2D[ki_index] = kj_index;
                    }
                    if (dirs_z[kj_index] > 0.99 && fabs(dirs_x[ki_index]-dirs_x[kj_index])<0.0001 && fabs(dirs_y[ki_index]-dirs_y[kj_index])<0.0001)
                    {
                        m_PhiToTop[ki_index] = kj_index;
                    }
                }
            }

            /*
            for (int k = 0; k < ratio; ++k)
            {
                for (int j = 0; j < nElmsSize2; ++j)
                {
                    m_PhiTo2D[index + k*nElmsSize2 + j] = index + nElmsSize2 + j;
                    if (i <= 2)
                    {
                        std::cout << index + k*nElmsSize2 + j << " " << index + nElmsSize2 + j << std::endl;
                    }
                }
                
            }
            */
            index += nElmsSize;
            bndIndex += nElmsSize2;
        }

        std::cout << "finishing setup..." << std::endl;



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
        
        
        
        
        
        
        
        
        

        std::cout << "timestep: " << m_timestep << std::endl;

        

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
                }

                //absorb zone
                if (dirs_y[i] > 0)
                {
                    niZ3++;
                }
                
            }
        }

        
        iZ2 = Array<OneD, int>(niZ2);
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

        std::cout << "finishing mappings..." << std::endl;
        

        XiZ2 = Array<OneD, NekDouble>(niZ2);
        Vmath::Gathr(niZ2, dirs_x, iZ2, XiZ2);

        Array<OneD, NekDouble> XiZ3(niZ3);
        Vmath::Gathr(niZ3, dirs_x, iZ3, XiZ3);

        std::cout << "hi" << std::endl;

        crm = Array<OneD, NekDouble>(niZ2);
        cfr = Array<OneD, NekDouble>(niZ3);
        spongefunction4(XiZ2, 0, 0, 10, crm);
        spongefunction4(XiZ3, 0, 5.0, 9, cfr);
        
        std::cout << "hi" << std::endl;   


        switch (m_projectionType)
        {
            case MultiRegions::eDiscontinuous:
            {
                std::string diffName;

                // Do not forwards transform initial condition
                //m_homoInitialFwd = false;

                //std::cout << "END" << std::endl;

                m_session->LoadSolverInfo("DiffusionType", diffName, "LDG");
                //m_diffusion = SolverUtils::GetDiffusionFactory().
                //    CreateInstance(diffName, diffName);
                //m_diffusion->SetFluxVector(&OceanWaveSystem::
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
            m_ode.DefineOdeRhs    (&OceanWaveSystem::DoOdeRhs,        this);
            m_ode.DefineProjection(&OceanWaveSystem::DoOdeProjection, this);

        //}
        //else
        //{
        //    m_ode.DefineImplicitSolve(
        //                        &OceanWaveSystem::DoImplicitSolve, this);
        //}
    }

    /**
     * @brief Unsteady diffusion problem destructor.
     */
    OceanWaveSystem::~OceanWaveSystem()
    {
    }

    void OceanWaveSystem::v_GenerateSummary(SummaryList& s)
    {
        UnsteadySystem::v_GenerateSummary(s);
        if(m_useSpecVanVisc)
        {
            stringstream ss;
            ss << "SVV (cut off = " << m_sVVCutoffRatio
               << ", coeff = "      << m_sVVDiffCoeff << ")";
            AddSummaryItem(s, "Smoothing", ss.str());
        }
    }
    
    
    /* @brief Compute the right-hand side for the unsteady diffusion problem.
     * 
     * @param inarray    Given fields.
     * @param outarray   Calculated solution.
     * @param time       Time.
     */

     void OceanWaveSystem::SetBoundaryConditions(
    Array<OneD, Array<OneD, NekDouble> > &inarray,
    NekDouble time)
  { 
  }


     void OceanWaveSystem::WallBoundary(
        int                                   bcRegion,
        int                                   cnt, 
        Array<OneD, Array<OneD, NekDouble> > &physarray)
    { 
        
    }



    void OceanWaveSystem::DoOdeRhs(
        const Array<OneD, const  Array<OneD, NekDouble> > &inarray,
              Array<OneD,        Array<OneD, NekDouble> > &outarray,
        const NekDouble time)
    {
        //std::cout << "here" << std::endl;
        std::cout << time << std::endl;

        int npoints = inarray[0].num_elements();

        Array<OneD, NekDouble> d(npoints, m_h0);
        Array<OneD, NekDouble> Dzsigma(npoints);

        Vmath::Sdiv(npoints, 1.0, d, 1, Dzsigma, 1);


        int npoints2 = m_PhiTo3D.num_elements();

        Array<OneD, NekDouble> phi2(npoints);
        Array<OneD, NekDouble> phi3(npoints, 1.0);

        Vmath::Vmul(npoints, Dzsigma, 1, Dzsigma, 1, phi2, 1);

        //Update boundaries for 3d problem:
        for (int i = 0; i < npoints2; ++i)
        {
            m_fields[0]->GetBndCondExpansions()[0]->UpdatePhys()[i] = inarray[1][m_PhiTo3D[i]];
        }
        m_fields[0]->GetBndCondExpansions()[0]->FwdTrans(m_fields[0]->GetBndCondExpansions()[0]->GetPhys(), 
                                                    m_fields[0]->GetBndCondExpansions()[0]->UpdateCoeffs());

        
        //setup new conditions for laplace solver:
        m_factors[StdRegions::eFactorLambda] = 0.0;
        m_factors[StdRegions::eFactorTau] = 1.0;
        StdRegions::VarCoeffMap m_varcoeff;

        //should be changed in the non-linear case:
        m_varcoeff[StdRegions::eVarCoeffD22] = phi3;
        m_varcoeff[StdRegions::eVarCoeffD11] = phi3;
        m_varcoeff[StdRegions::eVarCoeffD00] = phi3;

        
    
        // RHS computation using the new advection base class
        Vmath::Zero(m_fields[0]->GetNcoeffs(),m_fields[0]->UpdateCoeffs(),1);

        //m_fields[0]->ClearGlobalLinSysManager();
        m_fields[0]->HelmSolve(m_fields[0]->GetPhys(),
                                   m_fields[0]->UpdateCoeffs(),
                                   NullFlagList,
                                   m_factors,
                                   m_varcoeff);

        m_fields[0]->BwdTrans(m_fields[0]->GetCoeffs(),
                                              tmpPhi);


        //Find derivative of function.
        Array<OneD, NekDouble> out_d0(npoints);
        Array<OneD, NekDouble> out_d1(npoints);
        Array<OneD, NekDouble> out_d2(npoints);
        m_fields[0]->PhysDeriv(tmpPhi, out_d0, out_d1, out_d2);
        

        //rhsE update:        
        for (int i = 0; i < npoints; ++i)
        {
            NekDouble eta = 0; //inarray[2][m_PhiToTop[i]];
            //std::cout << outarray[1][i] << " ";
            outarray[2][i] = out_d2[m_PhiToTop[i]] / ( m_h0 );
        }
        

        //rhsP update:
        for (int i = 0; i < npoints; ++i)
        {
            outarray[1][i] = -m_g * inarray[2][m_PhiToTop[i]];
        }
        
        
        //zero changes to 3d problem:
        Vmath::Zero(npoints, outarray[0], 1);


        if (inarray.num_elements() > 5)
        {
            Vmath::Vcopy(npoints, outarray[2], 1, tmpRhsE, 1);
            Vmath::Vcopy(npoints, outarray[1], 1, tmpRhsP, 1);
        }


        //wave generation:
        WaveForcing(time, inarray[2], inarray[1], outarray[2], outarray[1]);

        //display out_d2:
        //Vmath::Vcopy(npoints, out_d2, 1, tmpPhi, 1);

        for (int i = 3; i < inarray.num_elements(); ++i)
        {
            Vmath::Zero(npoints, outarray[i], 1);
        }

        
    }


      void OceanWaveSystem::WallBoundary2D(int bcRegion, int cnt, Array<OneD, Array<OneD, NekDouble> > &physarray)
      { 
      }

    /**
     * @brief Compute the projection for the unsteady diffusion problem.
     * 
     * @param inarray    Given fields.
     * @param outarray   Calculated solution.
     * @param time       Time.
     */
    void OceanWaveSystem::DoOdeProjection(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
              Array<OneD,       Array<OneD, NekDouble> > &outarray,
        const NekDouble time)
    {

        int i;
        int nvariables = inarray.num_elements();

        //SetBoundaryConditions(outarray, time);


        //Find derivative of function.
        int n = tmpPhi.num_elements();
        Array<OneD, NekDouble> out_d0(n);
        Array<OneD, NekDouble> out_d1(n);
        Array<OneD, NekDouble> out_d2(n);
        m_fields[0]->PhysDeriv(tmpPhi, out_d0, out_d1, out_d2);

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
                Vmath::Vcopy(npoints, tmpPhi, 1, outarray[0], 1);
                //SetBoundaryConditions(outarray, time);
                break;
            }
            case MultiRegions::eGalerkin:
            case MultiRegions::eMixed_CG_Discontinuous:
            {
                Array<OneD, NekDouble> coeffs(m_fields[0]->GetNcoeffs());
                for(i = 0; i < nvariables; ++i)
                {  
                    m_fields[i]->FwdTrans(inarray[i], coeffs);
                    m_fields[i]->BwdTrans(coeffs, outarray[i]);
                }
                m_fields[0]->FwdTrans(tmpPhi, coeffs);
                m_fields[0]->BwdTrans(coeffs, outarray[0]);
                if (inarray.num_elements() > 3)
                {
                    m_fields[3]->FwdTrans(out_d2, coeffs);
                    m_fields[3]->BwdTrans(coeffs, outarray[3]);
                }
                
                if (inarray.num_elements() > 5)
                {
                    /*
                    m_fields[4]->FwdTrans(tmpRhsE, coeffs);
                    m_fields[4]->BwdTrans(coeffs, outarray[4]);

                    m_fields[5]->FwdTrans(tmpRhsP, coeffs);
                    m_fields[5]->BwdTrans(coeffs, outarray[5]);
                    */

                    m_fields[4]->FwdTrans(out_d0, coeffs);
                    m_fields[4]->BwdTrans(coeffs, outarray[4]);

                    m_fields[5]->FwdTrans(out_d1, coeffs);
                    m_fields[5]->BwdTrans(coeffs, outarray[5]);
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
    void OceanWaveSystem::DoImplicitSolve(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
              Array<OneD,       Array<OneD, NekDouble> > &outarray,
        const NekDouble time,
        const NekDouble lambda)
    {

    }
    
    /** 
     * @brief Return the flux vector for the unsteady diffusion problem.
     */
    void OceanWaveSystem::GetFluxVector(
        const int i, 
        const int j,
        const Array<OneD, Array<OneD, NekDouble> > &physfield,
              Array<OneD, Array<OneD, NekDouble> > &derivatives,
              Array<OneD, Array<OneD, NekDouble> > &flux)
    {
    }

    void OceanWaveSystem::spongefunction4(const Array<OneD, NekDouble> & x, NekDouble alpha, NekDouble p, int type, Array<OneD, NekDouble> & cr)
    {
        int n = x.num_elements();
        NekDouble xmin = Vmath::Vmin(n,x,1);
        NekDouble xmax = Vmath::Vmax(n,x,1);
        
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

    //void OceanWaveSystem::disper(NekDouble g, NekDouble T, NekDouble h, NekDouble)

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

    
    void OceanWaveSystem::lineartravellingwave1D(NekDouble H,NekDouble c, NekDouble k, NekDouble z, NekDouble h, NekDouble w, NekDouble t, const Array<OneD, NekDouble> & x, Array<OneD, NekDouble> & eta, Array<OneD, NekDouble> & pp )
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

    void OceanWaveSystem::WaveForcing(NekDouble time, const Array<OneD, NekDouble> & E, const Array<OneD, NekDouble> & P, Array<OneD, NekDouble> & rhsE, Array<OneD, NekDouble> & rhsP)
    {        
        int n = m_fields[0]->GetTotPoints();
        int niZ2 = iZ2.num_elements();
        int niZ3 = iZ3.num_elements();
        Array<OneD, NekDouble> X(n);
        Array<OneD, NekDouble> Y(n);
        Array<OneD, NekDouble> Z(n);
        m_fields[0]->GetCoords(X,Y,Z);

        NekDouble xmaxZ2 = Vmath::Vmax(niZ2, XiZ2, 1);

        Vmath::Sadd(n, -xmaxZ2, X, 1, X, 1);

        std::cout << Vmath::Vmax(n,X,1) << std::endl;
    
        Array<OneD, NekDouble> E2(n);
        Array<OneD, NekDouble> P2(n);
        lineartravellingwave1D(m_Hwave, m_cwave, m_kwave, 0, m_hd, m_wwave, time, X, E2, P2);


        if (time < 5*m_Twave)
        {
            Vmath::Smul(n, time/(5*m_Twave), E2, 1, E2, 1);
            Vmath::Smul(n, time/(5*m_Twave), P2, 1, P2, 1);
        }


        Array<OneD, NekDouble> crm_neg(niZ2);
        Vmath::Sadd(niZ2, -1.0, crm, 1, crm_neg, 1);
        Vmath::Neg(niZ2, crm_neg, 1);

        Array<OneD, NekDouble> cfr_neg(niZ3);
        Vmath::Sadd(niZ3, -1.0, cfr, 1, cfr_neg, 1);
        Vmath::Neg(niZ3, cfr_neg, 1);

        
        //Create submatrices:
        Array<OneD, NekDouble> EiZ2(niZ2);
        Vmath::Gathr(niZ2, E, iZ2, EiZ2);
        Array<OneD, NekDouble> E2iZ2(niZ2);
        Vmath::Gathr(niZ2, E2, iZ2, E2iZ2);
        Array<OneD, NekDouble> EiZ3(niZ3);
        Vmath::Gathr(niZ3, E, iZ3, EiZ3);

        Array<OneD, NekDouble> PiZ2(niZ2);
        Vmath::Gathr(niZ2, P, iZ2, PiZ2);
        Array<OneD, NekDouble> P2iZ2(niZ2);
        Vmath::Gathr(niZ2, P2, iZ2, P2iZ2);
        Array<OneD, NekDouble> PiZ3(niZ3);
        Vmath::Gathr(niZ3, P, iZ3, PiZ3);

        
        Array<OneD, NekDouble> rhsP2(n, 0.0);
        Array<OneD, NekDouble> rhsE2(n, 0.0);
        
        //modify rhsE:
        Vmath::Vsub(niZ2, E2iZ2, 1, EiZ2, 1, EiZ2, 1);
        Vmath::Vmul(niZ2, EiZ2, 1, crm_neg, 1, EiZ2, 1);
        Vmath::Smul(niZ2, 1.0/m_dt, EiZ2, 1, EiZ2, 1);
        Vmath::Scatr(niZ2, EiZ2, iZ2, rhsE2);

        Vmath::Vmul(niZ3, EiZ3, 1, cfr_neg, 1, EiZ3, 1);
        Vmath::Smul(niZ3, -1.0/m_dt, EiZ3, 1, EiZ3, 1);
        Vmath::Scatr(niZ3, EiZ3, iZ3, rhsE2);

        //modify rhsP:
        Vmath::Vsub(niZ2, P2iZ2, 1, PiZ2, 1, PiZ2, 1);
        Vmath::Vmul(niZ2, PiZ2, 1, crm_neg, 1, PiZ2, 1);
        Vmath::Smul(niZ2, 1.0/m_dt, PiZ2, 1, PiZ2, 1);
        Vmath::Scatr(niZ2, PiZ2, iZ2, rhsP2);

        Vmath::Vmul(niZ3, PiZ3, 1, cfr_neg, 1, PiZ3, 1);
        Vmath::Smul(niZ3, -1.0/m_dt, PiZ3, 1, PiZ3, 1);
        Vmath::Scatr(niZ3, PiZ3, iZ3, rhsP2);

        //final:
        Vmath::Vadd(n, rhsE, 1, rhsE2, 1, rhsE, 1);
        Vmath::Vadd(n, rhsP, 1, rhsP2, 1, rhsP, 1);

        
    }
    


}
