///////////////////////////////////////////////////////////////////////////////
//
// File: RoeSolver.cpp
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
// Description: Roe Riemann solver.
//
///////////////////////////////////////////////////////////////////////////////

#include <CompressibleFlowSolver/RiemannSolvers/RoeSolver.h>

namespace Nektar
{
    std::string RoeSolver::solverName =
        SolverUtils::GetRiemannSolverFactory().RegisterCreatorFunction(
            "Roe",
            RoeSolver::create,
            "Roe Riemann solver");

    RoeSolver::RoeSolver() : CompressibleSolver()
    {

    }

    /**
     * @brief Roe Riemann solver.
     *
     * Stated equations numbers are from:
     *
     *   "Riemann Solvers and Numerical Methods for Fluid Dynamics: A Practical
     *   Introduction", E. F. Toro (3rd edition, 2009).
     *
     * We follow the algorithm prescribed following equation 11.70.
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
    void RoeSolver::v_PointSolve(
        double  rhoL, double  rhouL, double  rhovL, double  rhowL, double  EL,
        double  rhoR, double  rhouR, double  rhovR, double  rhowR, double  ER,
        double &rhof, double &rhouf, double &rhovf, double &rhowf, double &Ef)
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
            (EL - 0.5 * (rhouL * uL + rhovL * vL + rhowL * wL));
        NekDouble pR = (gamma - 1.0) *
            (ER - 0.5 * (rhouR * uR + rhovR * vR + rhowR * wR));
        
        // Left and right enthalpy
        NekDouble hL = (EL + pL) / rhoL;
        NekDouble hR = (ER + pR) / rhoR;

        // Square root of rhoL and rhoR.
        NekDouble srL  = sqrt(rhoL);
        NekDouble srR  = sqrt(rhoR);
        NekDouble srLR = srL + srR;
        
        // Velocity, enthalpy and sound speed Roe averages (equation 11.60).
        NekDouble uRoe   = (srL * uL + srR * uR) / srLR;
        NekDouble vRoe   = (srL * vL + srR * vR) / srLR;
        NekDouble wRoe   = (srL * wL + srR * wR) / srLR;
        NekDouble hRoe   = (srL * hL + srR * hR) / srLR;
        NekDouble URoe   = (uRoe * uRoe + vRoe * vRoe + wRoe * wRoe);
        NekDouble cRoe   = sqrt((gamma - 1.0)*(hRoe - 0.5 * URoe));
        
        // Compute eigenvectors (equation 11.59).
        NekDouble k[5][5] = {
            {1, uRoe - cRoe, vRoe, wRoe, hRoe - uRoe * cRoe},
            {1, uRoe,        vRoe, wRoe, 0.5 * URoe},
            {0, 0,           1,    0,    vRoe},
            {0, 0,           0,    1,    wRoe},
            {1, uRoe+cRoe,  vRoe,  wRoe, hRoe + uRoe*cRoe}
        };

        // Calculate jumps \Delta u_i (defined preceding equation 11.67).
        NekDouble jump[5] = {
            rhoR  - rhoL,
            rhouR - rhouL,
            rhovR - rhovL,
            rhowR - rhowL,
            ER    - EL
        };

        // Define \Delta u_5 (equation 11.70).
        NekDouble jumpbar = jump[4] - (jump[2]-vRoe*jump[0])*vRoe -
            (jump[3]-wRoe*jump[0])*wRoe;
        
        // Compute wave amplitudes (equations 11.68, 11.69).
        NekDouble alpha[5];
        
        alpha[1] = ((gamma-1.0)/(cRoe*cRoe)) * (jump[0]*(hRoe - uRoe*uRoe)
                                            + uRoe*jump[1] - jumpbar);
        
        alpha[0] = (jump[0]*(uRoe + cRoe) - jump[1] - cRoe*alpha[1])/(2.0*cRoe);
        alpha[4] = jump[0] - (alpha[0] + alpha[1]);
        alpha[2] = jump[2] - vRoe * jump[0];
        alpha[3] = jump[3] - wRoe * jump[0];
        
        // Compute average of left and right fluxes needed for equation 11.29.
        rhof  = 0.5*(rhoL*uL + rhoR*uR);
        rhouf = 0.5*(pL + rhoL*uL*uL + pR + rhoR*uR*uR);
        rhovf = 0.5*(rhoL*uL*vL + rhoR*uR*vR);
        rhowf = 0.5*(rhoL*uL*wL + rhoR*uR*wR);
        Ef    = 0.5*(uL*(EL + pL) + uR*(ER + pR));

        // Compute eigenvalues \lambda_i (equation 11.58).
        NekDouble uRoeAbs = fabs(uRoe);
        NekDouble lambda[5] = {
            fabs(uRoe - cRoe),
            uRoeAbs,
            uRoeAbs,
            uRoeAbs,
            fabs(uRoe + cRoe)
        };

        // Finally perform summation (11.29).
        for (int i = 0; i < 5; ++i)
        {
            NekDouble ahat = 0.5*alpha[i]*lambda[i];
            rhof  -= ahat*k[i][0];
            rhouf -= ahat*k[i][1];
            rhovf -= ahat*k[i][2];
            rhowf -= ahat*k[i][3];
            Ef    -= ahat*k[i][4];
        }
    }
    void RoeSolver::v_PointAdjointSolve(
                                        double  rhoL, double  rhouL, double  rhovL, double  rhowL, double  EL,
                                        double  rhoR, double  rhouR, double  rhovR, double  rhowR, double  ER,
                                        double  rhoLdir, double  rhouLdir, double  rhovLdir, double  rhowLdir, double  ELdir,
                                        double  rhoRdir, double  rhouRdir, double  rhovRdir, double  rhowRdir, double  ERdir,
                                        double &rhof, double &rhouf, double &rhovf, double &rhowf, double &Ef)
    {
        static NekDouble gamma = m_params["gamma"]();
        
        // Left and right velocities
        NekDouble uLdir = rhouLdir / rhoLdir;
        NekDouble vLdir = rhovLdir / rhoLdir;
        NekDouble wLdir = rhowLdir / rhoLdir;
        NekDouble uRdir = rhouRdir / rhoRdir;
        NekDouble vRdir = rhovRdir / rhoRdir;
        NekDouble wRdir = rhowRdir / rhoRdir;
        
        // Left and right pressures
        NekDouble pLdir = (gamma - 1.0) *
        (ELdir - 0.5 * (rhouLdir * uLdir + rhovLdir * vLdir + rhowLdir * wLdir));
        NekDouble pRdir = (gamma - 1.0) *
        (ERdir - 0.5 * (rhouRdir * uRdir + rhovRdir * vRdir + rhowRdir * wRdir));
        
        // Left and right enthalpy
        NekDouble hLdir = (ELdir + pLdir) / rhoLdir;
        NekDouble hRdir = (ERdir + pRdir) / rhoRdir;
        
        // Square root of rhoL and rhoR.
        NekDouble srL  = sqrt(rhoLdir);
        NekDouble srR  = sqrt(rhoRdir);
        NekDouble srLR = srL + srR;
        
        // Velocity, enthalpy and sound speed Roe averages (equation 11.60).
        NekDouble uRoe   = (srL * uLdir + srR * uRdir) / srLR;
        NekDouble vRoe   = (srL * vLdir + srR * vRdir) / srLR;
        NekDouble wRoe   = (srL * wLdir + srR * wRdir) / srLR;
        NekDouble hRoe   = (srL * hLdir + srR * hRdir) / srLR;
        NekDouble URoe   = (uRoe * uRoe + vRoe * vRoe + wRoe * wRoe);
        NekDouble cRoe   = sqrt((gamma - 1.0)*(hRoe - 0.5 * URoe));
        
        NekDouble k[5][5] = {
            {0.5*URoe + (cRoe*uRoe)/(gamma-1.0), -uRoe-cRoe/(gamma-1.0), -vRoe, -wRoe, 1},
            {uRoe*uRoe-hRoe, -uRoe, 0, 0, 1},
            {-vRoe, 0, 1, 0, 0},
            {-wRoe, 0, 0, 1, 0},
            {0.5*URoe - (cRoe*uRoe)/(gamma-1.0), -uRoe+cRoe/(gamma-1.0), -vRoe, -wRoe, 1},
        };
        
        // Calculate jumps \Delta u_i (defined preceding equation 11.67).
        NekDouble jump[5] = {
            rhoR  - rhoL,
            rhouR - rhouL,
            rhovR - rhovL,
            rhowR - rhowL,
            ER    - EL
        };
        
        // wave strength
        NekDouble alpha[5];
        
        // =====================================================================
        alpha[0] = ((gamma - 1)/(2*(cRoe*cRoe)))*jump[0]
        + (-((cRoe - uRoe)*(gamma - 1))/(2*(cRoe*cRoe)))*jump[1]
        + ((vRoe*(gamma - 1))/(2*(cRoe*cRoe)))*jump[2]
        + ((wRoe*(gamma - 1))/(2*(cRoe*cRoe)))*jump[3]
        + (((gamma - 1)*(hRoe - cRoe*uRoe))/(2*(cRoe*cRoe)))*jump[4];
        
        alpha[1] = (-(gamma - 1)/(cRoe*cRoe))*jump[0]
        + (-(uRoe*(gamma - 1))/(cRoe*cRoe))*jump[1]
        + (-(vRoe*(gamma - 1))/(cRoe*cRoe))*jump[2]
        + (-(wRoe*(gamma - 1))/(cRoe*cRoe))*jump[3]
        + (-(URoe*(gamma - 1))/(2*(cRoe*cRoe)))*jump[4];
        
        alpha[2] = ((vRoe*(gamma - 1))/(cRoe*cRoe))*jump[0]
        + ((uRoe*vRoe*(gamma - 1))/(cRoe*cRoe))*jump[1]
        + (((vRoe*vRoe)*(gamma - 1))/((cRoe*cRoe)) + 1)*jump[2]
        + ((vRoe*wRoe*(gamma - 1))/(cRoe*cRoe))*jump[3]
        + ((hRoe*vRoe*(gamma - 1))/(cRoe*cRoe))*jump[4];
        
        alpha[3] = ((wRoe*(gamma - 1))/(cRoe*cRoe))*jump[0]
        + ((uRoe*wRoe*(gamma - 1))/(cRoe*cRoe))*jump[1]
        + ((vRoe*wRoe*(gamma - 1))/(cRoe*cRoe))*jump[2]
        + (((wRoe*wRoe)*(gamma - 1))/((cRoe*cRoe)) + 1)*jump[3]
        + ((hRoe*wRoe*(gamma - 1))/(cRoe*cRoe))*jump[4];
        
        alpha[4] = ((gamma - 1)/(2*(cRoe*cRoe)))*jump[0]
        + (((cRoe + uRoe)*(gamma - 1))/(2*(cRoe*cRoe)))*jump[1]
        + ((vRoe*(gamma - 1))/(2*(cRoe*cRoe)))*jump[2]
        + ((wRoe*(gamma - 1))/(2*(cRoe*cRoe)))*jump[3]
        + (((gamma - 1)*(hRoe + cRoe*uRoe))/(2*(cRoe*cRoe)))*jump[4];
        //======================================================================
        
        NekDouble   zrhoL_flux  = 0.0, zrhoR_flux  = 0.0,
        zrhouL_flux = 0.0, zrhouR_flux = 0.0,
        zrhovL_flux = 0.0, zrhovR_flux = 0.0,
        zrhowL_flux = 0.0, zrhowR_flux = 0.0,
        zEL_flux    = 0.0, zER_flux    = 0.0;
        
        NekDouble vsqLdir = pow(uLdir,2)+pow(vLdir,2)+pow(wLdir,2);
        NekDouble vsqRdir = pow(uRdir,2)+pow(vRdir,2)+pow(wRdir,2);
        
        // HLLC Riemann fluxes (positive case)
        
        zrhoL_flux =     (0.5*(gamma-1)*vsqLdir-uLdir*uLdir)*rhouL
        -uLdir*vLdir*rhovL
        -uLdir*wLdir*rhowL
        -uLdir*(hLdir - vsqLdir*0.5*(gamma-1))*EL;
        
        zrhouL_flux =   rhoL
        -(uLdir*(gamma-3))*rhouL
        +vLdir*rhovL
        +wLdir*rhowL
        +(hLdir+uLdir*uLdir*(1-gamma))*EL;
        
        zrhovL_flux =   -(vLdir*(gamma-1))*rhouL
        +uLdir*rhovL
        -uLdir*vLdir*(gamma-1)*EL;
        
        zrhowL_flux =    -(wLdir*(gamma-1))*rhouL
        +uLdir*rhowL
        -uLdir*wLdir*(gamma-1)*EL;
        
        zEL_flux    =    (gamma-1)*rhouL + (gamma*uLdir)*EL;
        
        //==================================
        
        zrhoR_flux =     (0.5*(gamma-1)*vsqRdir-pow(uRdir,2))*rhouR
        -uRdir*vRdir*rhovR
        -uRdir*wRdir*rhowR
        -uRdir*(hRdir - vsqRdir*0.5*(gamma-1))*ER;
        
        zrhouR_flux =    rhoR
        -(uRdir*(gamma-3))*rhouR
        +vRdir*rhovR
        +wRdir*rhowR
        +(hRdir+pow(uRdir,2)*(1-gamma))*ER;
        
        zrhovR_flux =   -(vRdir*(gamma-1))*rhouR
        +uRdir*rhovR
        -uRdir*vRdir*(gamma-1)*ER;
        
        zrhowR_flux =   -(wRdir*(gamma-1))*rhouR
        +uRdir*rhowR
        -uRdir*wRdir*(gamma-1)*ER;
        
        zER_flux    =    (gamma-1)*rhouR + (gamma*uRdir)*ER;
        
        //======================================================================
        
        rhof  = -0.5*(zrhoR_flux+zrhoL_flux);
        rhouf = -0.5*(zrhouR_flux+zrhouL_flux);
        rhovf = -0.5*(zrhovR_flux+zrhovL_flux);
        rhowf = -0.5*(zrhowR_flux+zrhowL_flux);
        Ef    = -0.5*(zER_flux+zEL_flux);
        
        // Compute eigenvalues \lambda_i (equation 11.58).
        NekDouble uRoeAbs = fabs(uRoe);
        NekDouble lambda[5] = {
            fabs(-uRoe - cRoe),
            uRoeAbs,
            uRoeAbs,
            uRoeAbs,
            fabs(uRoe + cRoe)
        };
        
        // Finally perform summation (11.29).
        for (int i = 0; i < 5; ++i)
        {
            NekDouble ahat = 0.5*alpha[i]*lambda[i];
            rhof  -= ahat*k[i][0];
            rhouf -= ahat*k[i][1];
            rhovf -= ahat*k[i][2];
            rhowf -= ahat*k[i][3];
            Ef    -= ahat*k[i][4];
        }
    }
    
    void RoeSolver::v_PointAdjointNSSolve(
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
        static NekDouble gamma = m_params["gamma"]();
        static NekDouble mu = m_params["mu"]();

        std::cout << mu << std::endl;
        static NekDouble R  = 287.058;
        static NekDouble k  = 0.0257;
        
        NekDouble FadjL = 0.0;
        NekDouble FadjR = 0.0;
        
        NekDouble uLdir = rhouLdir/rhoLdir;
        NekDouble uRdir = rhouRdir/rhoRdir;
        
        NekDouble vLdir = rhovLdir/rhoLdir;
        NekDouble vRdir = rhovRdir/rhoRdir;
        
        NekDouble vsqL = uLdir * uLdir + vLdir * vLdir;
        NekDouble vsqR = uRdir * uRdir + vRdir * vRdir;
        
        NekDouble PLDir = (gamma - 1.0) * (ELdir - 0.5*vsqL);
        NekDouble PRDir = (gamma - 1.0) * (ERdir - 0.5*vsqR);
        
        FadjL = -(0.5 * (gamma-1)/(R*rhoLdir*rhoLdir)*k*vsqL*DPLdirDX)*EL
        +(mu * vLdir/rhoLdir*(DrhouLdirDY+DrhovLdirDX))*EL
        +(2.0/3.0 * mu * uLdir*(2.0 * DrhouLdirDX - DrhovLdirDY))*EL
        +(k/(R*rhoLdir*rhoLdir*rhoLdir)*(2.0*DrhoLdirDX*PLDir+DPLdirDX*rhoLdir))*EL;
        
        FadjR = -(0.5 * (gamma-1)/(R*rhoRdir*rhoRdir)*k*vsqR*DPRdirDX)*ER
        +(mu * vRdir/rhoRdir*(DrhouRdirDY+DrhovRdirDX))*ER
        +(2.0/3.0 * mu * uRdir*(2.0 * DrhouRdirDX - DrhovRdirDY))*ER
        +(k/(R*rhoRdir*rhoRdir*rhoRdir)*(2.0*DrhoRdirDX*PRDir+DPRdirDX*rhoRdir))*ER;
        
        rhof = 0.5*(FadjL+FadjR);
        
        FadjL = 0.0;
        FadjR = 0.0;
        
        FadjL = -(2.0/3.0*mu/rhoLdir*(2.0*DrhouLdirDX-DrhovLdirDY))*EL
        +(DrhoLdirDX*k*uLdir*(gamma-1)/(R*rhoLdir*rhoLdir))*EL;
        
        FadjR = -(2.0/3.0*mu/rhoRdir*(2.0*DrhouRdirDX-DrhovRdirDY))*ER
        +(DrhoRdirDX*k*uRdir*(gamma-1)/(R*rhoRdir*rhoRdir))*ER;
        
        rhouf = 0.5*(FadjL+FadjR);
        
        FadjL = 0.0;
        FadjR = 0.0;
        
        FadjL = -(mu/rhoLdir*(DrhouLdirDY+DrhovLdirDX))*EL
        +(DrhoLdirDX*k/R*(gamma-1)*vLdir/(rhoLdir*rhoLdir))*EL;

        FadjR = -(mu/rhoRdir*(DrhouRdirDY+DrhovRdirDX))*ER
        +(DrhoRdirDX*k/R*(gamma-1)*vRdir/(rhoRdir*rhoRdir))*EL;
        
        rhovf = 0.5*(FadjL+FadjR);
        
        FadjL = 0.0;
        FadjR = 0.0;
        
        FadjL = -(DrhoLdirDX*k/R*(gamma-1)/(rhoLdir*rhoLdir))*EL;
        
        FadjR = -(DrhoRdirDX*k/R*(gamma-1)/(rhoRdir*rhoRdir))*ER;
        
        
        Ef = 0.5*(FadjL+FadjR);
        
        FadjL = 0.0;
        FadjR = 0.0;
    }
}
