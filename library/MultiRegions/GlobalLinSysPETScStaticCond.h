///////////////////////////////////////////////////////////////////////////////
//
// File GlobalLinSys.h
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
// Description: GlobalLinSysStaticCond header
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_MULTIREGIONS_GLOBALLINSYSPETSCSTATICCOND_H
#define NEKTAR_LIB_MULTIREGIONS_GLOBALLINSYSPETSCSTATICCOND_H

#include <MultiRegions/GlobalLinSysPETSc.h>
#include <MultiRegions/GlobalLinSysStaticCond.h>
#include <MultiRegions/MultiRegionsDeclspec.h>

namespace Nektar
{
    namespace MultiRegions
    {
        // Forward declarations
        class ExpList;
        class GlobalLinSysPETScStaticCond;

        typedef boost::shared_ptr<GlobalLinSysPETScStaticCond>
                                        GlobalLinSysPETScStaticCondSharedPtr;

        /// A global linear system.
        class GlobalLinSysPETScStaticCond : virtual public GlobalLinSysPETSc,
                                            virtual public GlobalLinSysStaticCond
        {
        public:
            /// Creates an instance of this class
            static GlobalLinSysSharedPtr create(
                const GlobalLinSysKey &pLinSysKey,
                const boost::weak_ptr<ExpList> &pExpList,
                const boost::shared_ptr<AssemblyMap> &pLocToGloMap)
            {
                GlobalLinSysSharedPtr p = MemoryManager<
                    GlobalLinSysPETScStaticCond>::AllocateSharedPtr(
                        pLinSysKey, pExpList, pLocToGloMap);
                p->InitObject();
                return p;
            }

            /// Name of class
            MULTI_REGIONS_EXPORT static std::string className;
            static std::string className2;

            /// Constructor for full direct matrix solve.
            MULTI_REGIONS_EXPORT GlobalLinSysPETScStaticCond(
                const GlobalLinSysKey                &mkey,
                const boost::weak_ptr<ExpList>       &pExpList,
                const boost::shared_ptr<AssemblyMap> &locToGloMap);

            /// Constructor for full direct matrix solve.
            MULTI_REGIONS_EXPORT GlobalLinSysPETScStaticCond(
                const GlobalLinSysKey                &mkey,
                const boost::weak_ptr<ExpList>       &pExpList,
                const DNekScalBlkMatSharedPtr         pSchurCompl,
                const DNekScalBlkMatSharedPtr         pBinvD,
                const DNekScalBlkMatSharedPtr         pC,
                const DNekScalBlkMatSharedPtr         pInvD,
                const boost::shared_ptr<AssemblyMap> &locToGloMap,
                const PreconditionerSharedPtr         pPrecon =
                                                     PreconditionerSharedPtr());

            MULTI_REGIONS_EXPORT virtual ~GlobalLinSysPETScStaticCond();

        protected:
            DNekScalBlkMatSharedPtr                  m_S1Blk;

            virtual void v_InitObject();

            /// Assemble the Schur complement matrix.
            virtual void v_AssembleSchurComplement(
                boost::shared_ptr<AssemblyMap> locToGloMap);
            virtual void v_DoMatrixMultiply(
                const Array<OneD, const NekDouble> &input,
                      Array<OneD,       NekDouble> &output);
            virtual DNekScalBlkMatSharedPtr v_GetStaticCondBlock(unsigned int n);
            virtual DNekScalBlkMatSharedPtr v_PreSolve(
                int                     scLevel,
                NekVector<NekDouble>   &F_GlobBnd);
            virtual void v_BasisTransform(
                Array<OneD, NekDouble>& pInOut,
                int                     offset);
            virtual void v_BasisInvTransform(
                Array<OneD, NekDouble>& pInOut);

            virtual GlobalLinSysStaticCondSharedPtr v_Recurse(
                const GlobalLinSysKey                &mkey,
                const boost::weak_ptr<ExpList>       &pExpList,
                const DNekScalBlkMatSharedPtr         pSchurCompl,
                const DNekScalBlkMatSharedPtr         pBinvD,
                const DNekScalBlkMatSharedPtr         pC,
                const DNekScalBlkMatSharedPtr         pInvD,
                const boost::shared_ptr<AssemblyMap> &locToGloMap);

        };
    }
}

#endif
