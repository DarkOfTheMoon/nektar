///////////////////////////////////////////////////////////////////////////////
//
// File: RiemannSolver.h
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
// Description: Abstract base class for Riemann solvers with factory.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_RIEMANNSOLVER
#define NEKTAR_SOLVERUTILS_RIEMANNSOLVER

#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/LinearAlgebra/NekTypeDefs.hpp>
#include <SolverUtils/SolverUtilsDeclspec.h>

#include <string>

#include <boost/function.hpp>
#include <boost/call_traits.hpp>

namespace Nektar
{
    template <typename Dim, typename DataType>
    class Array;

    namespace SolverUtils
    {
        typedef boost::function<
            const Array<OneD, const NekDouble>& ()> RSScalarFuncType;
        typedef boost::function<
            const Array<OneD, const Array<OneD, NekDouble> >& ()> RSVecFuncType;
        typedef boost::function<
            NekDouble ()> RSParamFuncType;
        
        typedef boost::function<void (
        Array<OneD, Array<OneD, NekDouble> >&,
        Array<OneD, Array<OneD, NekDouble> >&)>
        FwdBwdDirectSolutionVecCB;
        
        typedef boost::function<void (
        Array<OneD, Array<OneD, NekDouble> >&,
        Array<OneD, Array<OneD, NekDouble> >&,
        Array<OneD, Array<OneD, Array<OneD, NekDouble > > >&,
        Array<OneD, Array<OneD, Array<OneD, NekDouble > > >&)>
        FwdBwdDIFFDirectSolutionVecCB;

        class RiemannSolver
        {
        public:
            SOLVER_UTILS_EXPORT void Solve(
                const Array<OneD, const Array<OneD, NekDouble> > &Fwd,
                const Array<OneD, const Array<OneD, NekDouble> > &Bwd,
                      Array<OneD,       Array<OneD, NekDouble> > &flux);

            template<typename FuncPointerT, typename ObjectPointerT>
            void SetScalar(std::string    name,
                           FuncPointerT   func,
                           ObjectPointerT obj)
            {
                m_scalars[name] = boost::bind(func, obj);
            }

            void SetScalar(std::string name, RSScalarFuncType fp)
            {
                m_scalars[name] = fp;
            }

            template<typename FuncPointerT, typename ObjectPointerT>
            void SetVector(std::string    name,
                           FuncPointerT   func,
                           ObjectPointerT obj)
            {
                m_vectors[name] = boost::bind(func, obj);
            }

            void SetVector(std::string name, RSVecFuncType fp)
            {
                m_vectors[name] = fp;
            }

            template<typename FuncPointerT, typename ObjectPointerT>
            void SetParam(std::string    name,
                          FuncPointerT   func,
                          ObjectPointerT obj)
            {
                m_params[name] = boost::bind(func, obj);
            }

            void SetParam(std::string name, RSParamFuncType fp)
            {
                m_params[name] = fp;
            }
            
            template<typename FuncPointerT, typename ObjectPointerT>
            void SetFwdBwdDirectSolution(FuncPointerT func, ObjectPointerT obj)
            {
                m_FwdBwdDirectSolution = boost::bind(func, obj, _1, _2);
            }
            
            template<typename FuncPointerT, typename ObjectPointerT>
            void SetFwdBwdDIFFDirectSolution(FuncPointerT func, ObjectPointerT obj)
            {
                m_FwdBwdDIFFDirectSolution = boost::bind(func, obj, _1, _2, _3, _4);
            }
            
            template<typename FuncPointerT, typename ObjectPointerT>
            void SetAuxiliary(std::string    name,
                              FuncPointerT   func,
                              ObjectPointerT obj)
            {
                m_auxiliary[name] = boost::bind(func, obj);
            }

            std::map<std::string, RSScalarFuncType> &GetScalars()
            {
                return m_scalars;
            }

            std::map<std::string, RSVecFuncType>    &GetVectors()
            {
                return m_vectors;
            }

            std::map<std::string, RSParamFuncType>  &GetParams()
            {
                return m_params;
            }
            
            inline void SetFwdBwdDirectSolution(FwdBwdDirectSolutionVecCB FwdBwdDirectSol)
            {
                m_FwdBwdDirectSolution = FwdBwdDirectSol;
            }
            
            inline void SetFwdBwdDIFFDirectSolution(FwdBwdDIFFDirectSolutionVecCB FwdBwdDirectSol)
            {
                m_FwdBwdDIFFDirectSolution = FwdBwdDirectSol;
            }

        protected:
            /// Indicates whether the Riemann solver requires a rotation to be
            /// applied to the velocity fields.
            bool                                    m_requiresRotation;
            /// Map of scalar function types.
            std::map<std::string, RSScalarFuncType> m_scalars;
            /// Map of vector function types.
            std::map<std::string, RSVecFuncType>    m_vectors;
            /// Map of parameter function types.
            std::map<std::string, RSParamFuncType > m_params;
            /// Map of scalar function types.
            std::map<std::string, RSScalarFuncType> m_auxiliary;
            /// Rotation matrices for each trace quadrature point.
            Array<OneD, Array<OneD, NekDouble> >    m_rotMat;
            /// Rotation storage
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > m_rotStorage;
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > m_rotStoragetmp;
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > m_rotStorageDirSol;
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > m_rotStorageDirDIFFXSol;
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > m_rotStorageDirDIFFYSol;
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > m_rotStorageDirDIFFZSol;
            
            SOLVER_UTILS_EXPORT RiemannSolver();

            virtual void v_Solve(
                const Array<OneD, const Array<OneD, NekDouble> > &Fwd,
                const Array<OneD, const Array<OneD, NekDouble> > &Bwd,
                      Array<OneD,       Array<OneD, NekDouble> > &flux) = 0;
            
            virtual void v_AdjointSolve(
                const Array<OneD, const Array<OneD, NekDouble> > &FwdDir,
                const Array<OneD, const Array<OneD, NekDouble> > &BwdDir,
                const Array<OneD, const Array<OneD, NekDouble> > &Fwd,
                const Array<OneD, const Array<OneD, NekDouble> > &Bwd,
                      Array<OneD,       Array<OneD, NekDouble> > &flux) = 0;
            
            virtual void v_AdjointNSSolve(
                const Array<OneD, const Array<OneD, NekDouble> > &FwdDir,
                const Array<OneD, const Array<OneD, NekDouble> > &BwdDir,
                const Array<OneD, const Array<OneD, NekDouble> > &Fwd,
                const Array<OneD, const Array<OneD, NekDouble> > &Bwd,
                Array<OneD, Array<OneD, Array<OneD, NekDouble > > > &FwdDIFF,
                Array<OneD, Array<OneD, Array<OneD, NekDouble > > > &BwdDIFF,
                      Array<OneD,       Array<OneD, NekDouble> > &flux) = 0;

            void GenerateRotationMatrices(
                const Array<OneD, const Array<OneD, NekDouble> > &normals);
            void FromToRotation(
                Array<OneD, const NekDouble> &from,
                Array<OneD, const NekDouble> &to,
                NekDouble                    *mat);
            void rotateToNormal  (
                const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                const Array<OneD, const Array<OneD, NekDouble> > &normals,
                      Array<OneD,       Array<OneD, NekDouble> > &outarray);
            void rotateFromNormal(
                const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                const Array<OneD, const Array<OneD, NekDouble> > &normals,
                      Array<OneD,       Array<OneD, NekDouble> > &outarray);
            bool CheckScalars(std::string name);
            bool CheckVectors(std::string name);
            bool CheckParams (std::string name);
            
            FwdBwdDirectSolutionVecCB         m_FwdBwdDirectSolution;
            FwdBwdDIFFDirectSolutionVecCB     m_FwdBwdDIFFDirectSolution;
        };

        /// A shared pointer to an EquationSystem object
        typedef boost::shared_ptr<RiemannSolver> RiemannSolverSharedPtr;
        /// Datatype of the NekFactory used to instantiate classes derived
        /// from the RiemannSolver class.
        typedef LibUtilities::NekFactory<std::string, RiemannSolver>
            RiemannSolverFactory;
        SOLVER_UTILS_EXPORT RiemannSolverFactory& GetRiemannSolverFactory();
    }
}

#endif
