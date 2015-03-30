///////////////////////////////////////////////////////////////////////////////
//
// File OceanWaveSolver.h
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

#ifndef NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_OceanWaveSystem_H
#define NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_OceanWaveSystem_H

#include <SolverUtils/UnsteadySystem.h>
#include <SolverUtils/Diffusion/Diffusion.h>

using namespace Nektar::SolverUtils;

namespace Nektar
{
    class OceanWaveSystem : public UnsteadySystem
    {
    public:
        friend class MemoryManager<OceanWaveSystem>;

        /// Creates an instance of this class
        static EquationSystemSharedPtr create(
                const LibUtilities::SessionReaderSharedPtr& pSession) {
            EquationSystemSharedPtr p = MemoryManager<OceanWaveSystem>::AllocateSharedPtr(pSession);
            p->InitObject();
            return p;
        } 
        /// Name of class
        static std::string className;

        /// Destructor
        virtual ~OceanWaveSystem();

    protected:
        bool m_useSpecVanVisc;
        NekDouble m_sVVCutoffRatio;   // cut off ratio from which to start decayhing modes
        NekDouble m_sVVDiffCoeff;     // Diffusion coefficient of SVV modes
        SolverUtils::DiffusionSharedPtr         m_diffusion;        
        SolverUtils::RiemannSolverSharedPtr     m_riemannSolver;
        
        virtual void v_GenerateSummary(SummaryList &s);
        
        OceanWaveSystem(
                const LibUtilities::SessionReaderSharedPtr& pSession);
        
        virtual void v_InitObject();
        
        void GetFluxVector(
            const int i, 
            const int j,
            const Array<OneD, Array<OneD, NekDouble> > &physfield,
                  Array<OneD, Array<OneD, NekDouble> > &derivatives,
                  Array<OneD, Array<OneD, NekDouble> > &flux);
        void DoOdeRhs(
            const Array<OneD, const Array<OneD, NekDouble> >&inarray,
                  Array<OneD,       Array<OneD, NekDouble> >&outarray,
            const NekDouble time);
        void DoOdeProjection(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                  Array<OneD,       Array<OneD, NekDouble> > &outarray,
            const NekDouble time);
        virtual void DoImplicitSolve(
            const Array<OneD, const Array<OneD, NekDouble> >&inarray,
                  Array<OneD,       Array<OneD, NekDouble> >&outarray,
            NekDouble time,
            NekDouble lambda);

        void spongefunction4(const Array<OneD, NekDouble> & x, NekDouble alpha, NekDouble p, int type, Array<OneD, NekDouble> & cr);
        void lineartravellingwave1D(NekDouble H,NekDouble c, NekDouble k, NekDouble z, NekDouble h, NekDouble w, NekDouble t, const Array<OneD, NekDouble> & x, Array<OneD, NekDouble> & eta, Array<OneD, NekDouble> & pp );
        void WallBoundary2D(int bcRegion, int cnt, Array<OneD, Array<OneD, NekDouble> > &physarray);
        void WaveForcing(NekDouble time, 
            const Array<OneD, NekDouble> & P, 
            const Array<OneD, NekDouble> & E, 
            Array<OneD, NekDouble> & rhsP, 
            Array<OneD, NekDouble> & rhsE);
    

        void WallBoundary(
        int                                   bcRegion,
        int                                   cnt, 
        Array<OneD, Array<OneD, NekDouble> > &physarray);

        void SetBoundaryConditions(
    Array<OneD, Array<OneD, NekDouble> > &inarray,
    NekDouble time);

    private:
        NekDouble m_waveFreq;
        NekDouble m_epsilon;
        Array<OneD, int> m_PhiTo2D;
        Array<OneD, int> m_PhiToTop;
        Array<OneD, int> m_PhiTo3D;
        StdRegions::ConstFactorMap m_factors;

        Array<OneD, int> iZ2;
        Array<OneD, int> iZ3;

        //
        Array<OneD, NekDouble> crm;
        Array<OneD, NekDouble> cfr;

        Array<OneD, NekDouble> XiZ2;
    };
}

#endif
