///////////////////////////////////////////////////////////////////////////////
//
// File: HLLCSolver.h
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

#ifndef NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_RIEMANNSOLVER_HLLCSOLVER
#define NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_RIEMANNSOLVER_HLLCSOLVER

#include <CompressibleFlowSolver/RiemannSolvers/CompressibleSolver.h>

namespace Nektar
{
    class HLLCSolver : public CompressibleSolver
    {
    public:
        static RiemannSolverSharedPtr create()
        {
            return RiemannSolverSharedPtr(new HLLCSolver());
        }
        
        static std::string solverName;
        
    protected:
        HLLCSolver();
        
        virtual void v_PointSolve(
            NekDouble  rhoL, NekDouble  rhouL, NekDouble  rhovL, NekDouble  rhowL, NekDouble  EL,
            NekDouble  rhoR, NekDouble  rhouR, NekDouble  rhovR, NekDouble  rhowR, NekDouble  ER,
            NekDouble &rhof, NekDouble &rhouf, NekDouble &rhovf, NekDouble &rhowf, NekDouble &Ef);
        
        virtual void v_PointAdjointSolve(
            NekDouble  rhoL, NekDouble  rhouL, NekDouble  rhovL, NekDouble  rhowL, NekDouble  EL,
            NekDouble  rhoR, NekDouble  rhouR, NekDouble  rhovR, NekDouble  rhowR, NekDouble  ER,
            NekDouble  rhoLdir, NekDouble  rhouLdir, NekDouble  rhovLdir, NekDouble  rhowLdir, NekDouble  ELdir,
            NekDouble  rhoRdir, NekDouble  rhouRdir, NekDouble  rhovRdir, NekDouble  rhowRdir, NekDouble  ERdir,
            NekDouble &rhof, NekDouble &rhouf, NekDouble &rhovf, NekDouble &rhowf, NekDouble &Ef);
        
        virtual void v_PointSolveVisc(
            NekDouble  rhoL, NekDouble  rhouL, NekDouble  rhovL, NekDouble  rhowL, NekDouble  EL, NekDouble  EpsL,
            NekDouble  rhoR, NekDouble  rhouR, NekDouble  rhovR, NekDouble  rhowR, NekDouble  ER, NekDouble  EpsR,
            NekDouble &rhof, NekDouble &rhouf, NekDouble &rhovf, NekDouble &rhowf, NekDouble &Ef, NekDouble &Epsf);
        
        virtual void v_PointAdjointNSSolve(
                                           double  rhoL, double  rhouL, double  rhovL, double  rhowL, double  EL,
                                           double  rhoR, double  rhouR, double  rhovR, double  rhowR, double  ER,
                                           double  rhoLdir, double  rhouLdir, double  rhovLdir, double  rhowLdir, double  ELdir,
                                           double  rhoRdir, double  rhouRdir, double  rhovRdir, double  rhowRdir, double  ERdir,
                                           double  DrhoLdirDX, double  DrhouLdirDX, double  DrhovLdirDX, double  DrhowLdirDX, double  DELdirDX,
                                           double  DrhoRdirDX, double  DrhouRdirDX, double  DrhovRdirDX, double  DrhowRdirDX, double  DERdirDX,
                                           double  DrhoLdirDY, double  DrhouLdirDY, double  DrhovLdirDY, double  DrhowLdirDY, double  DELdirDY,
                                           double  DrhoRdirDY, double  DrhouRdirDY, double  DrhovRdirDY, double  DrhowRdirDY, double  DERdirDY,
                                           double  DrhoLdirDZ, double  DrhouLdirDZ, double  DrhovLdirDZ, double  DrhowLdirDZ, double  DELdirDZ,
                                           double  DrhoRdirDZ, double  DrhouRdirDZ, double  DrhovRdirDZ, double  DrhowRdirDZ, double  DERdirDZ,
                                           double &rhof, double &rhouf, double &rhovf, double &rhowf, double &Ef);
    };
}

#endif
