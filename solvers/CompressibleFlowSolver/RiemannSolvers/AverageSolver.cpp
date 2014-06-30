///////////////////////////////////////////////////////////////////////////////
//
// File: AverageSolver.cpp
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
// Description: Average Riemann solver.
//
///////////////////////////////////////////////////////////////////////////////

#include <CompressibleFlowSolver/RiemannSolvers/AverageSolver.h>

using namespace std;

namespace Nektar
{
    std::string AverageSolver::solverName =
        SolverUtils::GetRiemannSolverFactory().RegisterCreatorFunction(
            "Average",
            AverageSolver::create,
            "Average Riemann solver");

    AverageSolver::AverageSolver() : CompressibleSolver()
    {
        m_pointSolve = false;
    }

    /**
     * @brief Average Riemann solver.
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
    void AverageSolver::v_ArraySolve(
        const Array<OneD, const Array<OneD, NekDouble> > &Fwd,
        const Array<OneD, const Array<OneD, NekDouble> > &Bwd,
              Array<OneD,       Array<OneD, NekDouble> > &flux)
    {
        static NekDouble gamma = m_params["gamma"]();
        
        int expDim = Fwd.num_elements()-2;
        int i, j;
        
        for (j = 0; j < Fwd[0].num_elements(); ++j)
        {
            NekDouble tmp1 = 0.0, tmp2 = 0.0;
            Array<OneD, NekDouble> Ufwd(expDim);
            Array<OneD, NekDouble> Ubwd(expDim);
            
            for (i = 0; i < expDim; ++i)
            {
                Ufwd[i] = Fwd[i+1][j]/Fwd[0][j];
                Ubwd[i] = Bwd[i+1][j]/Bwd[0][j];
                tmp1   += Ufwd[i]*Fwd[i+1][j];
                tmp2   += Ubwd[i]*Bwd[i+1][j];
            }
            
            NekDouble Pfwd = (gamma - 1.0) * (Fwd[expDim+1][j] - 0.5 * tmp1);
            NekDouble Pbwd = (gamma - 1.0) * (Bwd[expDim+1][j] - 0.5 * tmp2);
            
            // Compute the average flux
            flux[0][j] = 0.5 * (Fwd[1][j] + Bwd[1][j]);
            flux[expDim+1][j] = 0.5 * (Ufwd[0] * (Fwd[expDim+1][j] + Pfwd) + 
                                       Ubwd[0] * (Bwd[expDim+1][j] + Pbwd));
            
            for (i = 0; i < expDim; ++i)
            {
                flux[i+1][j] = 0.5 * (Fwd[0][j] * Ufwd[0] * Ufwd[i] + 
                                      Bwd[0][j] * Ubwd[0] * Ubwd[i]);
            }

            // Add in pressure contribution to u field
            flux[1][j] += 0.5 * (Pfwd + Pbwd);
        }
    }
    
    void AverageSolver::v_ArrayAdjointSolve(
                const Array<OneD, const Array<OneD, NekDouble> > &Fwd,
                const Array<OneD, const Array<OneD, NekDouble> > &Bwd,
                const Array<OneD, const Array<OneD, NekDouble> > &FwdDir,
                const Array<OneD, const Array<OneD, NekDouble> > &BwdDir,
                Array<OneD,       Array<OneD, NekDouble> > &flux)
    {
        static NekDouble gamma = m_params["gamma"]();
    
        int expDim = Fwd.num_elements()-2;
        int nvariables = Fwd.num_elements();
        int i, j;
        
        for (j = 0; j < Fwd[0].num_elements(); ++j)
        {
            NekDouble tmp1 = 0.0, tmp2 = 0.0, vsqFwd = 0.0, vsqBwd = 0.0,
            FadjL = 0.0, FadjR = 0.0;
            
            Array<OneD, NekDouble> Ufwd(expDim);
            Array<OneD, NekDouble> Ubwd(expDim);
            
            Array<OneD, NekDouble> UfwdDir(expDim);
            Array<OneD, NekDouble> UbwdDir(expDim);
            
            for (i = 0; i < expDim; ++i)
            {
                UfwdDir[i] = FwdDir[i+1][j]/FwdDir[0][j];
                UbwdDir[i] = BwdDir[i+1][j]/BwdDir[0][j];
                
                tmp1    +=  UfwdDir[i]*FwdDir[i+1][j];
                tmp2    +=  UbwdDir[i]*BwdDir[i+1][j];
                vsqFwd  +=  UfwdDir[i]*UfwdDir[i];
                vsqBwd  +=  UbwdDir[i]*UbwdDir[i];
                
            }
            
            NekDouble PfwdDir = (gamma - 1.0) * (FwdDir[expDim+1][j] - 0.5*tmp1);
            NekDouble PbwdDir = (gamma - 1.0) * (BwdDir[expDim+1][j] - 0.5*tmp2);
            
            NekDouble HfwdDir = (FwdDir[expDim+1][j]+PfwdDir)/FwdDir[0][j];
            NekDouble HbwdDir = (BwdDir[expDim+1][j]+PbwdDir)/BwdDir[0][j];
        
            if (expDim == 1)
            {
                ASSERTL0(false, "1D adjoint solver not yet implemented");
            }
            if (expDim == 2)
            {
               FadjL = (0.5*(gamma-1)*vsqFwd-pow(UfwdDir[0],2))*Fwd[1][j]
                       -UfwdDir[0]*UfwdDir[1]*Fwd[2][j]
                       -UfwdDir[0]*(HfwdDir - vsqFwd*0.5*(gamma-1))*Fwd[3][j];
                
               FadjR = (0.5*(gamma-1)*vsqBwd-pow(UbwdDir[0],2))*Bwd[1][j]
                       -UbwdDir[0]*UbwdDir[1]*Bwd[2][j]
                       -UbwdDir[0]*(HbwdDir - vsqBwd*0.5*(gamma-1))*Bwd[3][j];
                
               flux[0][j] = 0.5*(FadjL+FadjR);
                
               FadjL = 0.0;
               FadjR = 0.0;
                
               FadjL =  Fwd[0][j]
                        -(UfwdDir[0]*(gamma-3))*Fwd[1][j]
                        +UfwdDir[1]*Fwd[2][j]
                        +(HfwdDir+pow(UfwdDir[0],2)*(1-gamma))*Fwd[3][j];
                
               FadjR =  Bwd[0][j]
                        -(UbwdDir[0]*(gamma-3))*Bwd[1][j]
                        +UbwdDir[1]*Bwd[2][j]
                        +(HbwdDir+pow(UbwdDir[0],2)*(1-gamma))*Bwd[3][j];
               
               flux[1][j] = 0.5*(FadjL+FadjR);
               
               FadjL = 0.0;
               FadjR = 0.0;
                
               FadjL =  -(UfwdDir[1]*(gamma-1))*Fwd[1][j]
                        +UfwdDir[0]*Fwd[2][j]
                        -UfwdDir[0]*UfwdDir[1]*(gamma-1)*Fwd[3][j];
               
               FadjR =  -(UbwdDir[1]*(gamma-1))*Bwd[1][j]
                        +UbwdDir[0]*Bwd[2][j]
                        -UbwdDir[0]*UbwdDir[1]*(gamma-1)*Bwd[3][j];
                
               flux[2][j] = 0.5*(FadjL+FadjR);
               
               FadjL = 0.0;
               FadjR = 0.0;
                
               FadjL = (gamma-1)*Fwd[1][j] + (gamma*UfwdDir[0])*Fwd[3][j];
               FadjR = (gamma-1)*Bwd[1][j] + (gamma*UbwdDir[0])*Bwd[3][j];
               
               flux[3][j] = 0.5*(FadjL+FadjR);

               FadjL = 0.0;
               FadjR = 0.0;
                
            }
            if (expDim == 3)
            {
               FadjL = (-pow(UfwdDir[0],2)+(gamma-1)/2*vsqFwd)*Fwd[1][j]
                            -UfwdDir[0]*UfwdDir[1]*Fwd[2][j]
                            -UfwdDir[0]*UfwdDir[2]*Fwd[3][j]
                            -(HfwdDir- vsqFwd*(gamma-1)/2)*Fwd[4][j];
                
               FadjR = (-pow(UbwdDir[0],2)+(gamma-1)/2*vsqBwd)*Bwd[1][j]
                            -UbwdDir[0]*UbwdDir[1]*Bwd[2][j]
                            -UbwdDir[0]*UbwdDir[2]*Bwd[3][j]
                            -(HbwdDir- vsqBwd*(gamma-1)/2)*Bwd[3][j];
               
               flux[0][j] = 0.5*(FadjL+FadjR);
                
               FadjL =  Fwd[0][j]
                              -(UfwdDir[0]*(gamma-3))*Fwd[1][j]
                              +UfwdDir[1]*Fwd[2][j]
                              +UfwdDir[2]*Fwd[3][j]
                              -(HfwdDir- UfwdDir[0]*(gamma-1))*Fwd[4][j];
                                    
               FadjR =  Bwd[0][j]
                              -(UbwdDir[0]*(gamma-3))*Bwd[1][j]
                              +UbwdDir[1]*Bwd[2][j]
                              +UbwdDir[2]*Bwd[2][j]
                              -(HbwdDir- UbwdDir[0]*(gamma-1))*Bwd[3][j];
               
               flux[1][j] = 0.5*(FadjL+FadjR);
                
               FadjL =  -(UfwdDir[1]*(gamma-1))*Fwd[1][j]
                              +UfwdDir[0]*Fwd[2][j]
                              +UfwdDir[0]*UfwdDir[1]*(gamma-1)*Fwd[3][j];
                         
               FadjR =  -(UbwdDir[1]*(gamma-1))*Bwd[1][j]
                              +UbwdDir[0]*Bwd[2][j]
                              +UbwdDir[0]*UbwdDir[1]*(gamma-1)*Bwd[3][j];
               
               flux[2][j] = 0.5*(FadjL+FadjR);
                          
               FadjL =  -(UfwdDir[2]*(gamma-1))*Fwd[1][j]
                              +UfwdDir[0]*Fwd[2][j]
                              +UfwdDir[0]*UfwdDir[2]*(gamma-1)*Fwd[3][j];
                                                 
               FadjR =  -(UbwdDir[2]*(gamma-1))*Bwd[1][j]
                              +UbwdDir[0]*Bwd[2][j]
                              +UbwdDir[0]*UbwdDir[2]*(gamma-1)*Bwd[3][j];
               
               flux[3][j] = 0.5*(FadjL+FadjR);
                          
               FadjL = (gamma-1)*Fwd[1][j] + (gamma*UfwdDir[0])*Fwd[4][j];
               FadjR = (gamma-1)*Bwd[1][j] + (gamma*UbwdDir[0])*Bwd[4][j];
                                
               flux[4][j] = 0.5*(FadjL+FadjR);
            }
        }
    }
    
    void AverageSolver::v_ArrayAdjointNSSolve(
            const Array<OneD, const Array<OneD, NekDouble> > &Fwd,
            const Array<OneD, const Array<OneD, NekDouble> > &Bwd,
            const Array<OneD, const Array<OneD, NekDouble> > &FwdDir,
            const Array<OneD, const Array<OneD, NekDouble> > &BwdDir,
            Array<OneD, Array<OneD, Array<OneD, NekDouble > > > &FwdDirDIFF,
            Array<OneD, Array<OneD, Array<OneD, NekDouble > > > &BwdDirDIFF,
                 Array<OneD,       Array<OneD, NekDouble> > &flux)
    {

        static NekDouble gamma = m_params["gamma"]();
        
        int expDim = Fwd.num_elements()-2;
        int nvariables = Fwd.num_elements();
        int i, j;
        
        static NekDouble R  = 287.05;
        static NekDouble k  = 0.0257;
        static NekDouble mu = 0.4169;
        
        //NekDouble R = 287.05; // gas constant
        //NekDouble k = 0.0257; // thermal conductivity
        
        for (j = 0; j < Fwd[0].num_elements(); ++j)
        {
            NekDouble tmp1 = 0.0, tmp2 = 0.0, vsqFwd = 0.0, vsqBwd = 0.0,
            FadjL = 0.0, FadjR = 0.0;
            
            Array<OneD, NekDouble> Ufwd(expDim);
            Array<OneD, NekDouble> Ubwd(expDim);
            
            Array<OneD, NekDouble> UfwdDir(expDim);
            Array<OneD, NekDouble> UbwdDir(expDim);
            
            for (i = 0; i < expDim; ++i)
            {
                UfwdDir[i] = FwdDir[i+1][j]/FwdDir[0][j];
                UbwdDir[i] = BwdDir[i+1][j]/BwdDir[0][j];
                
                tmp1    +=  UfwdDir[i]*FwdDir[i+1][j];
                tmp2    +=  UbwdDir[i]*BwdDir[i+1][j];
                vsqFwd  +=  UfwdDir[i]*UfwdDir[i];
                vsqBwd  +=  UbwdDir[i]*UbwdDir[i];
                
            }
            
            NekDouble PfwdDir = (gamma - 1.0) * (FwdDir[expDim+1][j] - 0.5*tmp1);
            NekDouble PbwdDir = (gamma - 1.0) * (BwdDir[expDim+1][j] - 0.5*tmp2);
            
            NekDouble HfwdDir = (FwdDir[expDim+1][j]+PfwdDir)/FwdDir[0][j];
            NekDouble HbwdDir = (BwdDir[expDim+1][j]+PbwdDir)/BwdDir[0][j];
            
            if (expDim == 1)
            {
                ASSERTL0(false, "1D adjoint solver not yet implemented");
            }
            if (expDim == 2)
            {
                /*FadjL = (0.5*(gamma-1)*vsqFwd-pow(UfwdDir[0],2))*Fwd[1][j]
                -UfwdDir[0]*UfwdDir[1]*Fwd[2][j]
                -UfwdDir[0]*(HfwdDir - vsqFwd*0.5*(gamma-1))*Fwd[3][j]*/
/*extra*/       FadjL = -(0.5 * (gamma-1)/(R*FwdDir[0][j]*FwdDir[0][j])*k*vsqFwd*FwdDirDIFF[0][3][j])*Fwd[3][j]
                +(mu * UfwdDir[1]/FwdDir[0][j]*(FwdDirDIFF[1][1][j]+FwdDirDIFF[0][2][j]))*Fwd[3][j]
                +(2.0/3.0 * mu * UfwdDir[0]*(2.0 * FwdDirDIFF[0][1][j] - FwdDirDIFF[1][2][j]))*Fwd[3][j]
                +(k/(R*FwdDir[0][j]*FwdDir[0][j]*FwdDir[0][j])*(2.0*FwdDirDIFF[0][0][j]*PfwdDir-FwdDirDIFF[0][3][j]*FwdDir[0][j]))*Fwd[3][j];
                
                /*FadjR = (0.5*(gamma-1)*vsqBwd-pow(UbwdDir[0],2))*Bwd[1][j]
                -UbwdDir[0]*UbwdDir[1]*Bwd[2][j]
                -UbwdDir[0]*(HbwdDir - vsqBwd*0.5*(gamma-1))*Bwd[3][j]*/
/*extra*/       FadjR = -(0.5 * (gamma-1)/(R*BwdDir[0][j]*BwdDir[0][j])*k*vsqBwd*BwdDirDIFF[0][3][j])*Bwd[3][j]
                +(mu * UbwdDir[1]/BwdDir[0][j]*(BwdDirDIFF[1][1][j]+BwdDirDIFF[0][2][j]))*Bwd[3][j]
                +(2.0/3.0 * mu * UbwdDir[0]*(2.0 * BwdDirDIFF[0][1][j] - BwdDirDIFF[1][2][j]))*Bwd[3][j]
                +(k/(R*BwdDir[0][j]*BwdDir[0][j]*BwdDir[0][j])*(2.0*BwdDirDIFF[0][0][j]*PbwdDir-BwdDirDIFF[0][3][j]*BwdDir[0][j]))*Bwd[3][j];

                
                flux[0][j] = 0.5*(FadjL+FadjR);
                
                FadjL = 0.0;
                FadjR = 0.0;
                
                /*FadjL =  Fwd[0][j]
                -(UfwdDir[0]*(gamma-3))*Fwd[1][j]
                +UfwdDir[1]*Fwd[2][j]
                +(HfwdDir+pow(UfwdDir[0],2)*(1-gamma))*Fwd[3][j]*/
/*extra*/       FadjL = -(2.0/3.0*mu/FwdDir[0][j]*(2.0*FwdDirDIFF[0][1][j]-FwdDirDIFF[1][2][j]))*Fwd[3][j]
                +(FwdDirDIFF[0][0][j]*k*UfwdDir[0]*(gamma-1)/(R*FwdDir[0][j]*FwdDir[0][j]))*Fwd[3][j];
                
                /*FadjR =  Bwd[0][j]
                -(UbwdDir[0]*(gamma-3))*Bwd[1][j]
                +UbwdDir[1]*Bwd[2][j]
                +(HbwdDir+pow(UbwdDir[0],2)*(1-gamma))*Bwd[3][j]*/
/*extra*/       FadjR =  -(2.0/3.0*mu/BwdDir[0][j]*(2.0*BwdDirDIFF[0][1][j]-BwdDirDIFF[1][2][j]))*Bwd[3][j]
                +(BwdDirDIFF[0][0][j]*k*UbwdDir[0]*(gamma-1)/(R*BwdDir[0][j]*BwdDir[0][j]))*Bwd[3][j];
                
                flux[1][j] = 0.5*(FadjL+FadjR);
                
                FadjL = 0.0;
                FadjR = 0.0;
                
                /*FadjL =  -(UfwdDir[1]*(gamma-1))*Fwd[1][j]
                +UfwdDir[0]*Fwd[2][j]
                -UfwdDir[0]*UfwdDir[1]*(gamma-1)*Fwd[3][j]*/
/*extra*/       FadjL =  -(mu/FwdDir[0][j]*(FwdDirDIFF[1][1][j]+FwdDirDIFF[0][2][j]))*Fwd[3][j]
                +(FwdDirDIFF[0][0][j]*k/R*(gamma-1)*UfwdDir[1]/(FwdDir[0][j]*FwdDir[0][j]))*Fwd[3][j];
                
                /*FadjR =  -(UbwdDir[1]*(gamma-1))*Bwd[1][j]
                +UbwdDir[0]*Bwd[2][j]
                -UbwdDir[0]*UbwdDir[1]*(gamma-1)*Bwd[3][j]*/
/*extra*/       FadjR = - (mu/BwdDir[0][j]*(BwdDirDIFF[1][1][j]+BwdDirDIFF[0][2][j]))*Bwd[3][j]
                +(BwdDirDIFF[0][0][j]*k/R*(gamma-1)*UbwdDir[1]/(BwdDir[0][j]*BwdDir[0][j]))*Fwd[3][j];
                
                flux[2][j] = 0.5*(FadjL+FadjR);
                
                FadjL = 0.0;
                FadjR = 0.0;
                
                /*FadjL = (gamma-1)*Fwd[1][j] + (gamma*UfwdDir[0])*Fwd[3][j]*/
                FadjL = -(FwdDirDIFF[0][0][j]*k/R*(gamma-1)/(FwdDir[0][j]*FwdDir[0][j]))*Fwd[3][j];
                
                /*FadjR = (gamma-1)*Bwd[1][j] + (gamma*UbwdDir[0])*Bwd[3][j]*/
                FadjR = -(BwdDirDIFF[0][0][j]*k/R*(gamma-1)/(BwdDir[0][j]*BwdDir[0][j]))*Bwd[3][j];
                
                
                flux[3][j] = 0.5*(FadjL+FadjR);
                
                FadjL = 0.0;
                FadjR = 0.0;
                
            }
        }
        
        /*cout << "FLUXES ADDED " << endl;
        cout << Vmath::Vmax(Fwd[0].num_elements(), flux[0], 1) << endl;
        cout << Vmath::Vmax(Fwd[0].num_elements(), flux[1], 1) << endl;
        cout << Vmath::Vmax(Fwd[0].num_elements(), flux[2], 1) << endl;
        cout << Vmath::Vmax(Fwd[0].num_elements(), flux[3], 1) << endl;
        cout << endl;*/
    }
}
