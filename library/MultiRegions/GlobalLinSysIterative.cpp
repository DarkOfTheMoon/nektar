///////////////////////////////////////////////////////////////////////////////
//
// File GlobalLinSysIterative.cpp
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
// Description: GlobalLinSysIterative definition
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/GlobalLinSysIterative.h>
#include <MultiRegions/LocalToGlobalC0ContMap.h>

namespace Nektar
{
    namespace MultiRegions
    {
        /**
         * @class GlobalLinSysIterative
         *
         * Solves a linear system using direct methods.
         */

        /// Constructor for full direct matrix solve.
        GlobalLinSysIterative::GlobalLinSysIterative(const GlobalLinSysKey &mkey,
                     const boost::shared_ptr<LocalMatrixSystem> &pLocMatSys,
                     const boost::shared_ptr<LocalToGlobalBaseMap>
                                                            &locToGloMap)
                : GlobalLinSys(mkey, pLocMatSys, locToGloMap)
        {
        }

        GlobalLinSysIterative::~GlobalLinSysIterative()
        {
        }

        /// Solve the linear system for given input and output vectors.
        void GlobalLinSysIterative::Solve( const Array<OneD,const NekDouble> &in,
                          Array<OneD,      NekDouble> &out)
        {
        }

        /// Solve the linear system for given input and output vectors
        /// using a specified local to global map.
        void GlobalLinSysIterative::Solve( const Array<OneD, const NekDouble> &in,
                          Array<OneD,       NekDouble> &out,
                    const LocalToGlobalBaseMapSharedPtr &locToGloMap,
                    const Array<OneD, const NekDouble> &dirForcing)
        {
        }
    }
}