///////////////////////////////////////////////////////////////////////////////
//
// File GlobalLinSysIterativeKokkos.h
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
// Description: GlobalLinSysIterativeKokkos header
//
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_LIB_MULTIREGIONS_GlobalLinSysIterativeKOKKOS_H
#define NEKTAR_LIB_MULTIREGIONS_GlobalLinSysIterativeKOKKOS_H

#include <MultiRegions/MultiRegionsDeclspec.h>
#include <MultiRegions/GlobalLinSys.h>
#include <MultiRegions/Preconditioner.h>

#include <boost/circular_buffer.hpp>

#ifdef NEKTAR_USE_KOKKOS
    #include <Kokkos_Core.hpp>
#endif

namespace Nektar
{
namespace MultiRegions
{
// Forward declarations
class ExpList;

/// A global linear system.
class GlobalLinSysIterativeKokkos : virtual public GlobalLinSys
{
public:
    /// Constructor for full direct matrix solve.
    MULTI_REGIONS_EXPORT GlobalLinSysIterativeKokkos(
            const GlobalLinSysKey                &pKey,
            const boost::weak_ptr<ExpList>       &pExpList,
            const boost::shared_ptr<AssemblyMap> &pLocToGloMap);

    MULTI_REGIONS_EXPORT virtual ~GlobalLinSysIterativeKokkos();

protected:
    /// Global to universal unique map
    Array<OneD, int>                            m_map;

    /// maximum iterations
    int                                         m_maxiter;

    /// Tolerance of iterative solver.
    NekDouble                                   m_tolerance;

    /// dot product of rhs to normalise stopping criterion
    NekDouble                                   m_rhs_magnitude;

    /// cnt to how many times rhs_magnitude is called
    NekDouble                                   m_rhs_mag_sm;

    PreconditionerSharedPtr                     m_precon;

    MultiRegions::PreconditionerType            m_precontype;

    int                                         m_totalIterations;

    /// Provide verbose output and root if parallel.
    bool                                        m_root;
    bool                                        m_verbose;

    /// Storage for solutions to previous linear problems
    boost::circular_buffer<Array<OneD, NekDouble> > m_prevLinSol;

    /// Total counter of previous solutions
    int m_numPrevSols;

    /// Actual iterative solve
    void DoConjugateGradient(
            const int pNumRows,
            const Array<OneD,const NekDouble> &pInput,
                  Array<OneD,      NekDouble> &pOutput,
            const AssemblyMapSharedPtr &locToGloMap,
            const int pNumDir);


    void Set_Rhs_Magnitude(const NekVector<NekDouble> &pIn);

    virtual void v_UniqueMap() = 0;

private:

    /// Solve the matrix system
    virtual void v_SolveLinearSystem(
            const int pNumRows,
            const Array<OneD,const NekDouble> &pInput,
                  Array<OneD,      NekDouble> &pOutput,
            const AssemblyMapSharedPtr &locToGloMap,
            const int pNumDir);

    virtual void v_DoMatrixMultiply(
            const Array<OneD, NekDouble>& pInput,
                  Array<OneD, NekDouble>& pOutput) = 0;
};
}
}

#endif
