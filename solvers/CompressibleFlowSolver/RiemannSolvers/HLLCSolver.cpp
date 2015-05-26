///////////////////////////////////////////////////////////////////////////////
//
// File: HLLCSolver.cpp
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
// Description: HLLC Riemann solver.
//
///////////////////////////////////////////////////////////////////////////////

#include <CompressibleFlowSolver/RiemannSolvers/HLLCSolver.h>

namespace Nektar
{
    std::string HLLCSolver::solverName =
        SolverUtils::GetRiemannSolverFactory().RegisterCreatorFunction(
            "HLLC", HLLCSolver::create, "HLLC Riemann solver");
    
    HLLCSolver::HLLCSolver() : CompressibleSolver()
    {
        
    }
    
    /**
     * @brief HLLC Riemann solver
     *
     * @param rhoL      Density left state.
     * @param rhoR      Density right state.
     * @param rhouL     x-momentum component left state.
     * @param rhouR     x-momentum component right state.
     * @param rhovL     y-momentum component left state.
     * @param rhovR     y-momentum component right state.
     * @param rhowL     z-momentum component left state.
     * @param rhowR     z-momentum component right state.
     * @param EL        Energy left state.
     * @param ER        Energy right state.
     * @param rhof      Computed Riemann flux for density.
     * @param rhouf     Computed Riemann flux for x-momentum component
     * @param rhovf     Computed Riemann flux for y-momentum component
     * @param rhowf     Computed Riemann flux for z-momentum component
     * @param Ef        Computed Riemann flux for energy.
     */
    void HLLCSolver::v_PointSolve(
        NekDouble  rhoL, NekDouble  rhouL, NekDouble  rhovL, NekDouble  rhowL, NekDouble  EL,
        NekDouble  rhoR, NekDouble  rhouR, NekDouble  rhovR, NekDouble  rhowR, NekDouble  ER,
        NekDouble &rhof, NekDouble &rhouf, NekDouble &rhovf, NekDouble &rhowf, NekDouble &Ef)
    {
        static NekDouble gamma = m_params["gamma"]();
        
        // Left and Right velocities
        NekDouble uL = rhouL / rhoL;
        NekDouble vL = rhovL / rhoL;
        NekDouble wL = rhowL / rhoL;
        NekDouble uR = rhouR / rhoR;
        NekDouble vR = rhovR / rhoR;
        NekDouble wR = rhowR / rhoR;
        
        // Left and right pressure, sound speed and enthalpy.
        NekDouble pL = (gamma - 1.0) *
            (EL - 0.5 * (rhouL * uL + rhovL * vL + rhowL * wL));
        NekDouble pR = (gamma - 1.0) *
            (ER - 0.5 * (rhouR * uR + rhovR * vR + rhowR * wR));
        NekDouble cL = sqrt(gamma * pL / rhoL);
        NekDouble cR = sqrt(gamma * pR / rhoR);
        NekDouble hL = (EL + pL) / rhoL;
        NekDouble hR = (ER + pR) / rhoR;
        
        // Square root of rhoL and rhoR.
        NekDouble srL  = sqrt(rhoL);
        NekDouble srR  = sqrt(rhoR);
        NekDouble srLR = srL + srR;
        
        // Velocity Roe averages
        NekDouble uRoe   = (srL * uL + srR * uR) / srLR;
        NekDouble vRoe   = (srL * vL + srR * vR) / srLR;
        NekDouble wRoe   = (srL * wL + srR * wR) / srLR;
        NekDouble hRoe   = (srL * hL + srR * hR) / srLR;
        NekDouble cRoe   = sqrt((gamma - 1.0)*(hRoe - 0.5 *
                                               (uRoe * uRoe + vRoe * vRoe + wRoe * wRoe)));
        
        // Maximum wave speeds
        NekDouble SL = std::min(uL-cL, uRoe-cRoe);
        NekDouble SR = std::max(uR+cR, uRoe+cRoe);
        
        // HLLC Riemann fluxes (positive case)
        if (SL >= 0)
        {
            rhof  = rhouL;
            rhouf = rhouL * uL + pL;
            rhovf = rhouL * vL;
            rhowf = rhouL * wL;
            Ef    = uL * (EL + pL);
        }
        // HLLC Riemann fluxes (negative case)
        else if (SR <= 0)
        {
            rhof  = rhouR;
            rhouf = rhouR * uR + pR;
            rhovf = rhouR * vR;
            rhowf = rhouR * wR;
            Ef    = uR * (ER + pR);
        }
        // HLLC Riemann fluxes (general case (SL < 0 | SR > 0)
        else
        {
            NekDouble SM = (pR - pL + rhouL * (SL - uL) - rhouR * (SR - uR)) /
            (rhoL * (SL - uL) - rhoR * (SR - uR));
            NekDouble rhoML  = rhoL * (SL - uL) / (SL - SM);
            NekDouble rhouML = rhoML * SM;
            NekDouble rhovML = rhoML * vL;
            NekDouble rhowML = rhoML * wL;
            NekDouble EML    = rhoML * (EL / rhoL +
                                        (SM - uL) * (SM + pL / (rhoL * (SL - uL))));
            
            NekDouble rhoMR  = rhoR * (SR - uR) / (SR - SM);
            NekDouble rhouMR = rhoMR * SM;
            NekDouble rhovMR = rhoMR * vR;
            NekDouble rhowMR = rhoMR * wR;
            NekDouble EMR    = rhoMR * (ER / rhoR +
                                        (SM - uR) * (SM + pR / (rhoR * (SR - uR))));
            
            if (SL < 0.0 && SM >= 0.0)
            {
                rhof  = rhouL + SL * (rhoML - rhoL);
                rhouf = rhouL * uL + pL + SL * (rhouML - rhouL);
                rhovf = rhouL * vL + SL * (rhovML - rhovL);
                rhowf = rhouL * wL + SL * (rhowML - rhowL);
                Ef    = uL * (EL + pL) + SL * (EML - EL);
            }
            else if(SM < 0.0 && SR > 0.0)
            {
                rhof  = rhouR + SR * (rhoMR - rhoR);
                rhouf = rhouR * uR + pR + SR * (rhouMR - rhouR);
                rhovf = rhouR * vR + SR * (rhovMR - rhovR);
                rhowf = rhouR * wR + SR * (rhowMR - rhowR);
                Ef    = uR * (ER + pR) + SR * (EMR - ER);
            }
        }
    }

    void HLLCSolver::v_PointSolveVisc(
        NekDouble  rhoL, NekDouble  rhouL, NekDouble  rhovL, NekDouble  rhowL, NekDouble  EL, NekDouble  EpsL,
        NekDouble  rhoR, NekDouble  rhouR, NekDouble  rhovR, NekDouble  rhowR, NekDouble  ER, NekDouble  EpsR,
        NekDouble &rhof, NekDouble &rhouf, NekDouble &rhovf, NekDouble &rhowf, NekDouble &Ef, NekDouble &Epsf)
    {
        static NekDouble gamma = m_params["gamma"]();
        
        // Left and Right velocities
        NekDouble uL = rhouL / rhoL;
        NekDouble vL = rhovL / rhoL;
        NekDouble wL = rhowL / rhoL;
        NekDouble uR = rhouR / rhoR;
        NekDouble vR = rhovR / rhoR;
        NekDouble wR = rhowR / rhoR;
        
        // Left and right pressure, sound speed and enthalpy.
        NekDouble pL = (gamma - 1.0) *
            (EL - 0.5 * (rhouL * uL + rhovL * vL + rhowL * wL));
        NekDouble pR = (gamma - 1.0) *
            (ER - 0.5 * (rhouR * uR + rhovR * vR + rhowR * wR));
        NekDouble cL = sqrt(gamma * pL / rhoL);
        NekDouble cR = sqrt(gamma * pR / rhoR);
        NekDouble hL = (EL + pL) / rhoL;
        NekDouble hR = (ER + pR) / rhoR;
        
        // Square root of rhoL and rhoR.
        NekDouble srL  = sqrt(rhoL);
        NekDouble srR  = sqrt(rhoR);
        NekDouble srLR = srL + srR;
        
        // Velocity Roe averages
        NekDouble uRoe   = (srL * uL + srR * uR) / srLR;
        NekDouble vRoe   = (srL * vL + srR * vR) / srLR;
        NekDouble wRoe   = (srL * wL + srR * wR) / srLR;
        NekDouble hRoe   = (srL * hL + srR * hR) / srLR;
        NekDouble cRoe   = sqrt((gamma - 1.0)*(hRoe - 0.5 *
                            (uRoe * uRoe + vRoe * vRoe + wRoe * wRoe)));
        
        // Maximum wave speeds
        NekDouble SL = std::min(uL-cL, uRoe-cRoe);
        NekDouble SR = std::max(uR+cR, uRoe+cRoe);
        
        // HLLC Riemann fluxes (positive case)
        if (SL >= 0)
        {
            rhof  = rhouL;
            rhouf = rhouL * uL + pL;
            rhovf = rhouL * vL;
            rhowf = rhouL * wL;
            Ef    = uL * (EL + pL);
            Epsf  = 0.0;
        }
        // HLLC Riemann fluxes (negative case)
        else if (SR <= 0)
        {
            rhof  = rhouR;
            rhouf = rhouR * uR + pR;
            rhovf = rhouR * vR;
            rhowf = rhouR * wR;
            Ef    = uR * (ER + pR);
            Epsf  = 0.0;
        }
        // HLLC Riemann fluxes (general case (SL < 0 | SR > 0)
        else
        {
            NekDouble SM = (pR - pL + rhouL * (SL - uL) - rhouR * (SR - uR)) /
            (rhoL * (SL - uL) - rhoR * (SR - uR));
            NekDouble rhoML  = rhoL * (SL - uL) / (SL - SM);
            NekDouble rhouML = rhoML * SM;
            NekDouble rhovML = rhoML * vL;
            NekDouble rhowML = rhoML * wL;
            NekDouble EML    = rhoML * (EL / rhoL +
                                        (SM - uL) * (SM + pL / (rhoL * (SL - uL))));
            NekDouble EpsML  = 0;
            
            NekDouble rhoMR  = rhoR * (SR - uR) / (SR - SM);
            NekDouble rhouMR = rhoMR * SM;
            NekDouble rhovMR = rhoMR * vR;
            NekDouble rhowMR = rhoMR * wR;
            NekDouble EMR    = rhoMR * (ER / rhoR +
                                        (SM - uR) * (SM + pR / (rhoR * (SR - uR))));
            NekDouble EpsMR    = 0;
            
            if (SL < 0.0 && SM >= 0.0)
            {
                rhof  = rhouL + SL * (rhoML - rhoL);
                rhouf = rhouL * uL + pL + SL * (rhouML - rhouL);
                rhovf = rhouL * vL + SL * (rhovML - rhovL);
                rhowf = rhouL * wL + SL * (rhowML - rhowL);
                Ef    = uL * (EL + pL) + SL * (EML - EL);
                Epsf  = 0.0;
            }
            else if(SM < 0.0 && SR > 0.0)
            {
                rhof  = rhouR + SR * (rhoMR - rhoR);
                rhouf = rhouR * uR + pR + SR * (rhouMR - rhouR);
                rhovf = rhouR * vR + SR * (rhovMR - rhovR);
                rhowf = rhouR * wR + SR * (rhowMR - rhowR);
                Ef    = uR * (ER + pR) + SR * (EMR - ER);
                Epsf  = 0.0;
            }
        }
    }
    
    void HLLCSolver::v_PointAdjointSolve(
        NekDouble  rhoL, NekDouble  rhouL, NekDouble  rhovL, NekDouble  rhowL, NekDouble  EL,
        NekDouble  rhoR, NekDouble  rhouR, NekDouble  rhovR, NekDouble  rhowR, NekDouble  ER,
        NekDouble  rhoLdir, NekDouble  rhouLdir, NekDouble  rhovLdir, NekDouble  rhowLdir, NekDouble  ELdir,
        NekDouble  rhoRdir, NekDouble  rhouRdir, NekDouble  rhovRdir, NekDouble  rhowRdir, NekDouble  ERdir,
        NekDouble &rhof, NekDouble &rhouf, NekDouble &rhovf, NekDouble &rhowf, NekDouble &Ef)
    {
        static NekDouble gamma = m_params["gamma"]();
        
        // Left and Right velocities
        NekDouble uLdir = rhouLdir / rhoLdir;
        NekDouble vLdir = rhovLdir / rhoLdir;
        NekDouble wLdir = rhowLdir / rhoLdir;
        NekDouble uRdir = rhouRdir / rhoRdir;
        NekDouble vRdir = rhovRdir / rhoRdir;
        NekDouble wRdir = rhowRdir / rhoRdir;
        
        // Left and right pressures
        NekDouble pLdir = (gamma - 1.0) *
        (ELdir - 0.5 * (rhouLdir*uLdir + rhovLdir*vLdir + rhowLdir*wLdir));
        
        NekDouble pRdir = (gamma - 1.0) *
        (ERdir - 0.5 * (rhouRdir*uRdir + rhovRdir*vRdir + rhowRdir*wRdir));
        
        // Left and right speeds of sound
        NekDouble cLdir = sqrt(gamma * pLdir / rhoLdir);
        NekDouble cRdir = sqrt(gamma * pRdir / rhoRdir);
        
        // Left and right entalpies
        NekDouble hLdir = (ELdir + pLdir) / rhoLdir;
        NekDouble hRdir = (ERdir + pRdir) / rhoRdir;
        
        // Square root of rhoL and rhoR.
        NekDouble srLdir  = sqrt(rhoLdir);
        NekDouble srRdir  = sqrt(rhoRdir);
        NekDouble srLRdir = srLdir + srRdir;
        
        // Velocity Roe averages
        NekDouble uRoe   = (srLdir * uLdir + srRdir * uRdir) / srLRdir;
        NekDouble vRoe   = (srLdir * vLdir + srRdir * vRdir) / srLRdir;
        NekDouble wRoe   = (srLdir * wLdir + srRdir * wRdir) / srLRdir;
        NekDouble hRoe   = (srLdir * hLdir + srRdir * hRdir) / srLRdir;
        NekDouble cRoe   = sqrt((gamma - 1.0)*(hRoe - 0.5 *
                            (uRoe * uRoe + vRoe * vRoe + wRoe * wRoe)));
        
        // Maximum wave speeds
        NekDouble SL = std::min(uLdir-cLdir, uRoe-cRoe);
        NekDouble SR = std::max(uRdir+cRdir, uRoe+cRoe);
        
        NekDouble zrhoL_flux  = 0.0, zrhoR_flux  = 0.0,
                  zrhouL_flux = 0.0, zrhouR_flux = 0.0,
                  zrhovL_flux = 0.0, zrhovR_flux = 0.0,
                  zrhowL_flux = 0.0, zrhowR_flux = 0.0,
                  zEL_flux    = 0.0, zER_flux    = 0.0;
        
        NekDouble vsqLdir = pow(uLdir,2)+pow(vLdir,2);
        NekDouble vsqRdir = pow(uRdir,2)+pow(vRdir,2);
        
        // HLLC Riemann fluxes (positive case)
        if (SL >= 0)
        {
            zrhoL_flux =     (0.5*(gamma-1)*vsqLdir-pow(uLdir,2))*rhouL
                            -uLdir*vLdir*rhovL
                            -uLdir*(hLdir - vsqLdir*0.5*(gamma-1))*EL;
            
            zrhouL_flux =    rhoL
                            -(uLdir*(gamma-3))*rhouL
                            +vLdir*rhovL
                            +(hLdir+pow(uLdir,2)*(1-gamma))*EL;
            
            zrhovL_flux =    -(vLdir*(gamma-1))*rhouL
                            +uLdir*rhovL
                            -uLdir*vLdir*(gamma-1)*EL;
            
            zrhowL_flux =   0.0;
            
            zEL_flux    =    (gamma-1)*rhouL + (gamma*uLdir)*EL;
            
            
            rhof  = zrhoL_flux;
            rhouf = zrhouL_flux;
            rhovf = zrhovL_flux;
            rhowf = zrhowL_flux;
            Ef    = zEL_flux;
        }
        // HLLC Riemann fluxes (negative case)
        else if (SR <= 0)
        {
            zrhoR_flux =     (0.5*(gamma-1)*vsqRdir-pow(uRdir,2))*rhouR
                            -uRdir*vRdir*rhovR
                            -uRdir*(hRdir - vsqRdir*0.5*(gamma-1))*ER;
            
            zrhouR_flux =    rhoR
                            -(uRdir*(gamma-3))*rhouR
                            +vRdir*rhovR
                            +(hRdir+pow(uRdir,2)*(1-gamma))*ER;
            
            zrhovR_flux =    -(vRdir*(gamma-1))*rhouR
                            +uRdir*rhovR
                            -uRdir*vRdir*(gamma-1)*ER;
            
            zrhowL_flux =   0.0;
            
            zER_flux    =    (gamma-1)*rhouR + (gamma*uRdir)*ER;
            
            
            rhof  = zrhoR_flux;
            rhouf = zrhouR_flux;
            rhovf = zrhovR_flux;
            rhowf = zrhowL_flux;
            Ef    = zER_flux;
        }
        // HLLC Riemann fluxes (general case (SL < 0 | SR > 0)
        else
        {
            
            NekDouble SM = (pRdir - pLdir + rhouLdir * (SL - uLdir) - rhouRdir * (SR - uRdir)) /
            (rhoL * (SL - uLdir) - rhoRdir * (SR - uRdir));
            
            
            NekDouble rhoML  = rhoLdir * (SL - uLdir) / (SL - SM);
            NekDouble rhouML = rhoML * SM;
            NekDouble rhovML = rhoML * vLdir;
            NekDouble rhowML = rhoML * wLdir;
            NekDouble EML    = rhoML * (ELdir / rhoLdir +
                    (SM - uLdir) * (SM + pLdir / (rhoLdir * (SL - uLdir))));
        
            
            NekDouble rhoMR  = rhoR * (SR - uRdir) / (SR - SM);
            NekDouble rhouMR = rhoMR * SM;
            NekDouble rhovMR = rhoMR * vRdir;
            NekDouble rhowMR = rhoMR * wRdir;
            NekDouble EMR    = rhoMR * (ERdir / rhoRdir +
                    (SM - uRdir) * (SM + pRdir / (rhoRdir * (SR - uRdir))));
            
            if (rhoRdir == 0.0)
            {
                SM  = 0.0;
                EMR = 0.0;
            }
            
            if (SL < 0.0 && SM >= 0.0)
            {
                zrhoL_flux =     (0.5*(gamma-1)*vsqLdir-pow(uLdir,2))*rhouL
                                -uLdir*vLdir*rhovL
                                -uLdir*(hLdir - vsqLdir*0.5*(gamma-1))*EL;
                
                zrhouL_flux =    rhoL
                                -(uLdir*(gamma-3))*rhouL
                                +vLdir*rhovL
                                +(hLdir+pow(uLdir,2)*(1-gamma))*EL;
                
                zrhovL_flux =    -(vLdir*(gamma-1))*rhouL
                                +uLdir*rhovL
                                -uLdir*vLdir*(gamma-1)*EL;
                
                zrhowL_flux =   0.0;
                
                zEL_flux    =    (gamma-1)*rhouL + (gamma*uLdir)*EL;

                rhof  = zrhoL_flux + SL * (rhoML - rhoL);
                rhouf = zrhouL_flux + SL * (rhouML - rhouL);
                rhovf = zrhovL_flux + SL * (rhovML - rhovL);
                rhowf = zrhowL_flux + SL * (rhowML - rhowL);
                Ef    = zEL_flux + SL * (EML - EL);
            }
            else if(SM < 0.0 && SR > 0.0)
            {
                zrhoR_flux =     (0.5*(gamma-1)*vsqRdir-pow(uRdir,2))*rhouR
                                -uRdir*vRdir*rhovR
                                -uRdir*(hRdir - vsqRdir*0.5*(gamma-1))*ER;
                
                zrhouR_flux =    rhoR
                                -(uRdir*(gamma-3))*rhouR
                                +vRdir*rhovR
                                +(hRdir+pow(uRdir,2)*(1-gamma))*ER;
                
                zrhovR_flux =    -(vRdir*(gamma-1))*rhouR
                                +uRdir*rhovR
                                -uRdir*vRdir*(gamma-1)*ER;
                
                zrhowL_flux =   0.0;
                
                zER_flux    =    (gamma-1)*rhouR + (gamma*uRdir)*ER;
                
                
                rhof  = zrhoR_flux + SR * (rhoMR - rhoR);
                rhouf = zrhouR_flux + SR * (rhouMR - rhouR);
                rhovf = zrhovR_flux + SR * (rhovMR - rhovR);
                rhowf = zrhowL_flux + SR * (rhowMR - rhowR);
                Ef    = zER_flux + SR * (EMR - ER);
            }
        }
    }
    
    void HLLCSolver::v_PointAdjointNSSolve(
        double  rhoL, double  rhouL, double  rhovL, double  rhowL, double  EL,
        double  rhoR, double  rhouR, double  rhovR, double  rhowR, double  ER,
        double  rhoLdir, double  rhouLdir, double  rhovLdir, double  rhowLdir, double  ELdir,
        double  rhoRdir, double  rhouRdir, double  rhovRdir, double  rhowRdir, double  ERdir,
        double  DrhoLdirDX, double  DrhouLdirDX, double  DrhovLdirDX, double  DrhowLdirDX, double  DPLdirDX,
        double  DrhoRdirDX, double  DrhouRdirDX, double  DrhovRdirDX, double  DrhowRdirDX, double  DPRdirDX,
        double  DrhoLdirDY, double  DrhouLdirDY, double  DrhovLdirDY, double  DrhowLdirDY, double  DPLdirDY,
        double  DrhoRdirDY, double  DrhouRdirDY, double  DrhovRdirDY, double  DrhowRdirDY, double  DPRdirDY,
        double  DrhoLdirDZ, double  DrhouLdirDZ, double  DrhovLdirDZ, double  DrhowLdirDZ, double  DPLdirDZ,
        double  DrhoRdirDZ, double  DrhouRdirDZ, double  DrhovRdirDZ, double  DrhowRdirDZ, double  DPRdirDZ,
        double &rhof, double &rhouf, double &rhovf, double &rhowf, double &Ef)
    {
        
    }
}
