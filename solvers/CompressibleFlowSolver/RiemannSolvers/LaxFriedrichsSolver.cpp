///////////////////////////////////////////////////////////////////////////////
//
// File: LaxFriedrichsSolver.cpp
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
// Description: Lax-Friedrichs Riemann solver.
//
///////////////////////////////////////////////////////////////////////////////

#include <CompressibleFlowSolver/RiemannSolvers/LaxFriedrichsSolver.h>

namespace Nektar
{
    std::string LaxFriedrichsSolver::solverName =
        SolverUtils::GetRiemannSolverFactory().RegisterCreatorFunction(
            "LaxFriedrichs", LaxFriedrichsSolver::create,
            "Lax-Friedrichs Riemann solver");

    LaxFriedrichsSolver::LaxFriedrichsSolver() : CompressibleSolver()
    {
        
    }
    
    /**
     * @brief Lax-Friedrichs Riemann solver
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
    void LaxFriedrichsSolver::v_PointSolve(
        NekDouble  rhoL, NekDouble  rhouL, NekDouble  rhovL, NekDouble  rhowL, NekDouble  EL,
        NekDouble  rhoR, NekDouble  rhouR, NekDouble  rhovR, NekDouble  rhowR, NekDouble  ER,
        NekDouble &rhof, NekDouble &rhouf, NekDouble &rhovf, NekDouble &rhowf, NekDouble &Ef)
    {
        static NekDouble gamma = m_params["gamma"]();
        
        // Left and right velocities
        NekDouble uL = rhouL / rhoL;
        NekDouble vL = rhovL / rhoL;
        NekDouble wL = rhowL / rhoL;
        NekDouble uR = rhouR / rhoR;
        NekDouble vR = rhovR / rhoR;
        NekDouble wR = rhowR / rhoR;
        
        // Left and right pressures
        NekDouble pL = (gamma - 1.0) *
        (EL - 0.5 * (rhouL*uL + rhovL*vL + rhowL*wL));
        
        NekDouble pR = (gamma - 1.0) *
        (ER - 0.5 * (rhouR*uR + rhovR*vR + rhowR*wR));
        
        // Left and right speeds of sound
        NekDouble cL = sqrt(gamma * pL / rhoL);
        NekDouble cR = sqrt(gamma * pR / rhoR);
        
        // Left and right entalpies
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
        
        // Minimum and maximum wave speeds
        NekDouble S    = std::max(uRoe+cRoe, std::max(uR+cR, uL+cL));
        NekDouble sign = 1.0;
        
        /*if(S == -uL+cL)
        {
            sign = -1.0;
        }*/
        
        // Lax-Friedrichs Riemann rho flux
        rhof  = 0.5 * ((rhouL + rhouR) - sign * S * (rhoR -rhoL));
        
        // Lax-Friedrichs Riemann rhou flux
        rhouf = 0.5 * ((rhoL * uL * uL + pL + rhoR * uR * uR + pR) -
                       sign * S * (rhouR - rhouL));
        
        // Lax-Friedrichs Riemann rhov flux
        rhovf = 0.5 * ((rhoL * uL * vL + rhoR * uR * vR) -
                       sign * S * (rhovR - rhovL));
        
        // Lax-Friedrichs Riemann rhow flux
        rhowf = 0.5 * ((rhoL * uL * wL + rhoR * uR * wR) -
                       sign * S * (rhowR - rhowL));
        
        // Lax-Friedrichs Riemann E flux
        Ef    = 0.5 * ((uL * (EL + pL) + uR * (ER + pR)) -
                       sign * S * (ER - EL));
    }
    
    void LaxFriedrichsSolver::v_PointSolveVisc(
        NekDouble  rhoL, NekDouble  rhouL, NekDouble  rhovL, NekDouble  rhowL, NekDouble  EL, NekDouble  EpsL,
        NekDouble  rhoR, NekDouble  rhouR, NekDouble  rhovR, NekDouble  rhowR, NekDouble  ER, NekDouble  EpsR,
        NekDouble &rhof, NekDouble &rhouf, NekDouble &rhovf, NekDouble &rhowf, NekDouble &Ef, NekDouble &Epsf)
    {
        static NekDouble gamma = m_params["gamma"]();
        
        // Left and right velocities
        NekDouble uL = rhouL / rhoL;
        NekDouble vL = rhovL / rhoL;
        NekDouble wL = rhowL / rhoL;
        NekDouble uR = rhouR / rhoR;
        NekDouble vR = rhovR / rhoR;
        NekDouble wR = rhowR / rhoR;
        
        // Left and right pressures
        NekDouble pL = (gamma - 1.0) *
        (EL - 0.5 * (rhouL*uL + rhovL*vL + rhowL*wL));
        
        NekDouble pR = (gamma - 1.0) *
        (ER - 0.5 * (rhouR*uR + rhovR*vR + rhowR*wR));
        
        // Left and right speeds of sound
        NekDouble cL = sqrt(gamma * pL / rhoL);
        NekDouble cR = sqrt(gamma * pR / rhoR);
        
        // Left and right entalpies
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
        
        // Minimum and maximum wave speeds
        NekDouble S    = std::max(uRoe+cRoe, std::max(uR+cR, uL+cL));
        NekDouble sign = 1.0;
        /*
        if(S == -uL+cL)
        {
            sign = -1.0;
        }*/
        
        // Lax-Friedrichs Riemann rho flux
        rhof  = 0.5 * ((rhouL + rhouR) - sign * S * (rhoR -rhoL));
        
        // Lax-Friedrichs Riemann rhou flux
        rhouf = 0.5 * ((rhoL * uL * uL + pL + rhoR * uR * uR + pR) -
                       sign * S * (rhouR - rhouL));
        
        // Lax-Friedrichs Riemann rhov flux
        rhovf = 0.5 * ((rhoL * uL * vL + rhoR * uR * vR) -
                       sign * S * (rhovR - rhovL));
        
        // Lax-Friedrichs Riemann rhow flux
        rhowf = 0.5 * ((rhoL * uL * wL + rhoR * uR * wR) -
                       sign * S * (rhowR - rhowL));
        
        // Lax-Friedrichs Riemann E flux
        Ef    = 0.5 * ((uL * (EL + pL) + uR * (ER + pR)) -
                       sign * S * (ER - EL));
        
        Epsf  = 0.5 * ((EpsL + EpsR) -
                       sign * S * (EpsR - EpsL));
    }
    
    void LaxFriedrichsSolver::v_PointAdjointSolve(
        NekDouble  rhoL, NekDouble  rhouL, NekDouble  rhovL, NekDouble  rhowL, NekDouble  EL,
        NekDouble  rhoR, NekDouble  rhouR, NekDouble  rhovR, NekDouble  rhowR, NekDouble  ER,
        NekDouble  rhoLdir, NekDouble  rhouLdir, NekDouble  rhovLdir, NekDouble  rhowLdir, NekDouble  ELdir,
        NekDouble  rhoRdir, NekDouble  rhouRdir, NekDouble  rhovRdir, NekDouble  rhowRdir, NekDouble  ERdir,
        NekDouble &rhof, NekDouble &rhouf, NekDouble &rhovf, NekDouble &rhowf, NekDouble &Ef)
    {
        static NekDouble gamma = m_params["gamma"]();
        
        /*rhof  = 0.5*(rhoL+rhoR);
        rhouf = 0.5*(rhouL+rhouR);
        rhovf = 0.5*(rhovL+rhovR);
        rhowf = 0.5*(rhowL+rhowR);
        Ef    = 0.5*(EL+ER);*/
        
        // Left and right velocities
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
        
        // Minimum and maximum wave speeds
        NekDouble S    = std::min(-uRoe-cRoe, std::min(-uRdir-cRdir, -uLdir-cLdir));
        S = abs(S);
        
        NekDouble sign = 1.0;
        
        /*if(S == -uLdir+cLdir)
        {
            sign = -1.0;
        }*/
        
        NekDouble zrhoL_flux  = 0.0, zrhoR_flux  = 0.0,
                  zrhouL_flux = 0.0, zrhouR_flux = 0.0,
                  zrhovL_flux = 0.0, zrhovR_flux = 0.0,
                  zrhowL_flux = 0.0, zrhowR_flux = 0.0,
                  zEL_flux    = 0.0, zER_flux    = 0.0;
        
        NekDouble vsqLdir = pow(uLdir,2)+pow(vLdir,2)+pow(wLdir,2);
        NekDouble vsqRdir = pow(uRdir,2)+pow(vRdir,2)+pow(wRdir,2);
        // =====================================================================
        
        zrhoL_flux =     (0.5*(gamma-1)*vsqLdir-pow(uLdir,2))*rhouL
                         -uLdir*vLdir*rhovL
                         -uLdir*wLdir*rhowL
                         -uLdir*(hLdir - vsqLdir*0.5*(gamma-1))*EL;
        
        zrhoR_flux =     (0.5*(gamma-1)*vsqRdir-pow(uRdir,2))*rhouR
                         -uRdir*vRdir*rhovR
                         -uRdir*wRdir*rhowR
                         -uRdir*(hRdir - vsqRdir*0.5*(gamma-1))*ER;

        // =====================================================================
        
        zrhouL_flux =    rhoL
                         -(uLdir*(gamma-3))*rhouL
                         +vLdir*rhovL
                         +wLdir*rhowL
                         +(hLdir+pow(uLdir,2)*(1-gamma))*EL;
        
        zrhouR_flux=     rhoR
                         -(uRdir*(gamma-3))*rhouR
                         +vRdir*rhovR
                         +wRdir*rhowR
                         +(hRdir+pow(uRdir,2)*(1-gamma))*ER;

        // =====================================================================
        
        zrhovL_flux =    -(vLdir*(gamma-1))*rhouL
                         +uLdir*rhovL
                         -uLdir*vLdir*(gamma-1)*EL;
        
        
        zrhovR_flux =    -(vRdir*(gamma-1))*rhouR
                         +uRdir*rhovR
                         -uRdir*vRdir*(gamma-1)*ER;

        // =====================================================================
        
        zrhowL_flux =    -(wLdir*(gamma-1))*rhouL
                         +uLdir*rhowL
                         -uLdir*wLdir*(gamma-1)*EL;
        
        
        zrhowR_flux =    -(wRdir*(gamma-1))*rhouR
                         +uRdir*rhowR
                         -uRdir*wRdir*(gamma-1)*ER;

        // =====================================================================
        
        zEL_flux    =    (gamma-1)*rhouL + (gamma*uLdir)*EL;
        zER_flux    =    (gamma-1)*rhouR + (gamma*uRdir)*ER;
        

        // =====================================================================
        // Lax-Friedrichs Riemann rho flux
        rhof  = 0.5 * (-(zrhoL_flux + zrhoR_flux) - sign * S * (rhoR - rhoL));// - sign * S * (rhoR - rhoL));
        
        // Lax-Friedrichs Riemann rhou flux
        rhouf = 0.5 * (-(zrhouL_flux + zrhouR_flux)  - sign *  S * (rhouR - rhouL));// - sign * S * (rhouR - rhouL));
        
        // Lax-Friedrichs Riemann rhov flux
        rhovf = 0.5 * (-(zrhovL_flux + zrhovR_flux)  - sign *  S * (rhovR - rhovL));// - sign * S * (rhovR - rhovL));
        
        // Lax-Friedrichs Riemann rhow flux
        rhowf = 0.5 * (-(zrhowL_flux + zrhowR_flux)  - sign *  S * (rhowR - rhowL));// - sign * S * (rhowR - rhowL));
        
        // Lax-Friedrichs Riemann E flux
        Ef    = 0.5 * (-(zEL_flux + zER_flux)  - sign *   S * (ER - EL));// - sign * S * (ER - EL));
    }
}
