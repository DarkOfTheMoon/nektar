///////////////////////////////////////////////////////////////////////////////
//
// File AdjointCompressibleFlowSystem.h
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
// Description: Auxiliary functions for the compressible flow system
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_EQUATIONSYSTEMS_ADJOINTCOMPRESSIBLEFLOWSYSTEM_H
#define NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_EQUATIONSYSTEMS_ADJOINTCOMPRESSIBLEFLOWSYSTEM_H

#include <SolverUtils/UnsteadySystem.h>
#include <SolverUtils/RiemannSolvers/RiemannSolver.h>
#include <SolverUtils/AdvectionSystem.h>
#include <SolverUtils/Diffusion/Diffusion.h>
#include <StdRegions/StdQuadExp.h>
#include <StdRegions/StdHexExp.h>


#define EPSILON 0.000001

#define CROSS(dest, v1, v2){                 \
          dest[0] = v1[1] * v2[2] - v1[2] * v2[1]; \
          dest[1] = v1[2] * v2[0] - v1[0] * v2[2]; \
          dest[2] = v1[0] * v2[1] - v1[1] * v2[0];}

#define DOT(v1, v2) (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2])

#define SUB(dest, v1, v2){       \
          dest[0] = v1[0] - v2[0]; \
          dest[1] = v1[1] - v2[1]; \
          dest[2] = v1[2] - v2[2];}

namespace Nektar
{  
    /**
     * 
     */
    class AdjointCompressibleFlowSystem: public SolverUtils::UnsteadySystem
    {
    public:

        friend class MemoryManager<AdjointCompressibleFlowSystem>;

        /// Creates an instance of this class
        static SolverUtils::EquationSystemSharedPtr create(
            const LibUtilities::SessionReaderSharedPtr& pSession)
        {
            return MemoryManager<AdjointCompressibleFlowSystem>::AllocateSharedPtr(pSession);
        }
        /// Name of class
        static std::string className;
      
        virtual ~AdjointCompressibleFlowSystem();
        
        /// Function to calculate the stability limit for DG/CG.
        NekDouble GetStabilityLimit(int n);
        
        /// Function to calculate the stability limit for DG/CG 
        /// (a vector of them).
        Array<OneD, NekDouble> GetStabilityLimitVector(
            const Array<OneD,int> &ExpOrder);
      
    protected:
        SolverUtils::RiemannSolverSharedPtr m_riemannSolver;
        SolverUtils::RiemannSolverSharedPtr m_riemannSolverLDG;
        SolverUtils::AdvectionSharedPtr     m_advection;
        SolverUtils::DiffusionSharedPtr     m_diffusion;
        Array<OneD, Array<OneD, NekDouble> >m_vecLocs;
        NekDouble                           m_gamma;
        NekDouble                           m_pInf;
        NekDouble                           m_rhoInf;
        NekDouble                           m_uInf;
        NekDouble                           m_vInf;
        NekDouble                           m_wInf;
        NekDouble                           m_gasConstant;
        NekDouble                           m_Twall;
        std::string                         m_ViscosityType;
        std::string                         m_shockCaptureType;
	    std::string                         m_EqTypeStr;
        NekDouble                           m_mu;
        NekDouble                           m_thermalConductivity;
        NekDouble                           m_Cp;
        NekDouble                           m_Prandtl;
        NekDouble                           m_rhoInfBase;
        NekDouble                           m_uInfBase;
        NekDouble                           m_vInfBase;
        NekDouble                           m_pInfBase;
        NekDouble                           m_UInf;
        NekDouble                           m_alphaInfBase;
        NekDouble                           m_Skappa;
        NekDouble                           m_Kappa;
        NekDouble                           m_mu0;
        NekDouble                           m_muMAX;
        NekDouble                           m_muMAXLIM;
        NekDouble                           m_muMIN;
        NekDouble                           m_muMINLIM;
        NekDouble                           m_Lref;
        NekDouble                           m_alpha;
        int                                 m_numCheck;
        int                                 m_finalCheck;
        int                                 m_cnt;
        StdRegions::StdQuadExpSharedPtr     m_OrthoQuadExp;
        StdRegions::StdHexExpSharedPtr      m_OrthoHexExp;
        NekDouble                           m_Fx;
        NekDouble                           m_Fy;
        NekDouble                           m_Fz;
        
        // L2 error file
        std::ofstream m_errFile;
        
        // Check for steady state at step interval
        int m_steadyStateSteps;
        
        // Tolerance to which steady state should be evaluated at
        NekDouble m_steadyStateTol;
        
        Array<OneD, Array<OneD, NekDouble> > m_un;
        
        // Storage for the forward solution
        Array<OneD, Array<OneD, NekDouble> > m_baseflow;
        
        Array<OneD, MultiRegions::ExpListSharedPtr> m_basefields;
        
        // Storage for the jacobians
        Array<OneD, Array<OneD, Array<OneD, Array<OneD, NekDouble> > > >
                                                                      m_dVdUdXi;
        
        Array<OneD, Array<OneD, Array<OneD, Array<OneD, NekDouble> > > >
                                                                      m_JacPrim;
        Array<OneD, Array<OneD, Array<OneD, Array<OneD, NekDouble> > > >
                                                                      m_JacCons;
        
        Array<OneD, Array<OneD, Array<OneD, Array<OneD, NekDouble> > > >
                                                                   m_JacDivPrim;
        Array<OneD, Array<OneD, Array<OneD, Array<OneD, NekDouble> > > >
                                                                   m_JacAddPrim;
        Array<OneD, Array<OneD, Array<OneD, Array<OneD, NekDouble> > > >
                                                                m_JacAddDivPrim;
        Array<OneD, Array<OneD, Array<OneD, Array<OneD, NekDouble> > > >
                                                                   m_JacAddCons;
        Array<OneD, Array<OneD, Array<OneD, Array<OneD, NekDouble> > > >
                                                                   m_JacDivCons;
        
        Array<OneD, Array<OneD, Array<OneD, Array<OneD, NekDouble> > > >
                                                                m_JacAddDivCons;
        
        Array<OneD, Array<OneD, Array<OneD, Array<OneD,
                                            Array<OneD, NekDouble> > > > >
                                                                      m_JacVisc;
        
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > > m_dVdU;
        
        AdjointCompressibleFlowSystem(
            const LibUtilities::SessionReaderSharedPtr& pSession);

        virtual void v_InitObject();
      
        /// Print a summary of time stepping parameters.
        virtual void v_GenerateSummary(SolverUtils::SummaryList& s);
    
        
        //======================================================================
        void GetConservToPrimVariableMat(
           const Array<OneD, Array<OneD, NekDouble> > &inarray,
                 Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &dUdV);
        
        void GetConservToPrimVariableInvMat(
           Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &dUdVInv);
        
        void GetConservToPrimVariableInvMatDiv(
           Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &dUdVInv,
           Array<OneD, Array<OneD, Array<OneD, Array<OneD, NekDouble> > > >
                                                                   &dUdVInvdXi);
       
        //======================================================================
        void GetJacobianConvFluxPrim(
          Array<OneD, Array<OneD, Array<OneD, Array<OneD, NekDouble> > > > &Jac);
        void GetJacobianConvFlux(
        Array<OneD, Array<OneD, Array<OneD, Array<OneD, NekDouble> > > > &Jac);
        void GetAdjointDerivJacVector(
          const Array<OneD, Array<OneD, NekDouble> > &inarray,
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &outarray);
        
        void GetAdjointViscousFluxVector(
         const Array<OneD, Array<OneD, NekDouble> > &inarray,
               Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &z_derivatives,
               Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &outarray);
        
        void GetJacobianAddConvFlux(
           Array<OneD, Array<OneD, Array<OneD, Array<OneD, NekDouble> > > > &Jac);
        
        void GetJacobianViscousFluxPrim(
        Array<OneD, Array<OneD, Array<OneD, Array<OneD,
                                   Array<OneD, NekDouble> > > > >        JacVisc);
        
        void GetAdjointAddConvFluxVector(
              const Array<OneD, Array<OneD, NekDouble> > &inarray,
                    Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &outarray);
        
        void GetAdjointDerivAddJacVector(
             const Array<OneD, Array<OneD, NekDouble> > &inarray,
                   Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &outarray);
        
        void GetDerivJacobian(
        Array<OneD, Array<OneD, Array<OneD, Array<OneD, NekDouble> > > > &Jac,
        Array<OneD, Array<OneD, Array<OneD, Array<OneD, NekDouble> > > > &JacDiv);
        //======================================================================
        
        void GetAdjointFluxVector(
            const Array<OneD, Array<OneD, NekDouble> > &inarray,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &outarray);
        void GetAdjointFluxVectorAlternative(
            const Array<OneD, Array<OneD, NekDouble> > &inarray,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &outarray);

        void GetFwdBwdBaseFlow(
                  Array<OneD, Array<OneD, NekDouble> > &FwdDir,
                  Array<OneD, Array<OneD, NekDouble> > &BwdDir);
        void GetFwdBwdDerivBaseFlow(
                 Array<OneD, Array<OneD, NekDouble> > &FwdDir,
                 Array<OneD, Array<OneD, NekDouble> > &BwdDir,
                 Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &FwdDirDIFF,
                 Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &BwdDirDIFF);
        void GetViscousFluxVectorDeAlias(
            const Array<OneD, Array<OneD, NekDouble> >         &physfield,
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &derivatives,
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &viscousTensor);
        void GetPressure(
            const Array<OneD, const Array<OneD,       NekDouble> >&physfield,
            Array<OneD,                         NekDouble>  &pressure);
        void GetPressure(
            const Array<OneD, const Array<OneD,       NekDouble> >&physfield,
            const Array<OneD, const Array<OneD,       NekDouble> >&velocity,
                  Array<OneD,                         NekDouble>  &pressure);
        void GetVelocityVector(
            const Array<OneD,       Array<OneD,       NekDouble> >&physfield,
                  Array<OneD,       Array<OneD,       NekDouble> >&velocity);
        void SymmetryBC(
            int                                                 bcRegion,
            int                                                 cnt,
            Array<OneD, Array<OneD, NekDouble> >               &physarray);
        void AdjointWallBC(
            int                                                 bcRegion,
            int                                                 cnt,
            Array<OneD, Array<OneD, NekDouble> >               &physarray);
        void GetSensor(
        const Array<OneD, const Array<OneD,       NekDouble> > &physarray,
              Array<OneD,                         NekDouble>   &Sensor,
              Array<OneD,                         NekDouble>   &SensorKappa);
        void GetArtificialDynamicViscosity(
        const Array<OneD,  Array<OneD, NekDouble> > &physfield,
              Array<OneD,                    NekDouble  > &mu_var);
        void GetSoundSpeed(
        const Array<OneD,       Array<OneD,       NekDouble> >&physfield,
              Array<OneD,                         NekDouble>  &pressure,
              Array<OneD,                         NekDouble>  &soundspeed);
        void GetAbsoluteVelocity(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
              Array<OneD,                   NekDouble>   &Vtot);
        
        void GetJacobianConvFluxPrimitiveVar(
                                             Array<OneD, Array<OneD, Array<OneD, Array<OneD, NekDouble> > > > &Jac);
        
        void GetDerivJacVectorFromPrimitiveVar(
                                               const Array<OneD, Array<OneD, NekDouble> > &inarray,
                                               Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &outarray);
        
        void GetConservToPrimVariableInvMatDiv(
                                               Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &dUdVInv,
                                               Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &dUdVInvdX,
                                               Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &dUdVInvdY);
        
        void GetElementDimensions(
                                  Array<OneD,                   NekDouble >  &hmin);
        
        void v_ExtraFldOutput(
                std::vector<Array<OneD, NekDouble> > &fieldcoeffs,
                              std::vector<std::string>             &variables);
  
        
        void SetCommonBC(const std::string &userDefStr,
                         const int n,
                         const NekDouble time,
                         int &cnt,
                         Array<OneD, Array<OneD, NekDouble> > &inarray);
        
        virtual bool v_PostIntegrate(int step);
        bool CalcSteadyState(bool output);

        
        void CheckForRestart(NekDouble &time);
        
        virtual void v_SetInitialConditions(
            NekDouble initialtime = 0.0,
            bool dumpInitialConditions = true,
            const int domain = 0)
        {
        }
        NekDouble GetGasConstant()
        {
            return m_gasConstant;
        }
        
        NekDouble GetGamma()
        {
            return m_gamma;
        }
      
        const Array<OneD, const Array<OneD, NekDouble> > &GetVecLocs()
        {
            return m_vecLocs;
        }
        
        const Array<OneD, const Array<OneD, NekDouble> > &GetNormals()
        {
            
            /*int nDim = m_traceNormals.num_elements();
            int nTracePts = m_traceNormals[0].num_elements();
            
            for (int i = 0; i < nDim; ++i)
            {
                Vmath::Neg(nTracePts, m_traceNormals[i], 1);
            }*/
            
            return m_traceNormals;
        }
        
        virtual void v_DoInitialise(void);
    };
}
#endif
