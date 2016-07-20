///////////////////////////////////////////////////////////////////////////////
//
// File: ExtrapOrder0BC.cpp
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
// Description: Extrapolation of order 0 boundary condition
//
///////////////////////////////////////////////////////////////////////////////

#include "ExtrapOrder0BC.h"

using namespace std;

namespace Nektar
{

std::string ExtrapOrder0BC::className = GetCFSBndCondFactory().
    RegisterCreatorFunction("ExtrapOrder0",
                            ExtrapOrder0BC::create,
                            "Extrapolation of order 0 boundary condition.");

ExtrapOrder0BC::ExtrapOrder0BC(const LibUtilities::SessionReaderSharedPtr& pSession,
           const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
           const Array<OneD, Array<OneD, NekDouble> >& pTraceNormals,
           const int pSpaceDim,
           const int bcRegion)
    : CFSBndCond(pSession, pFields, pTraceNormals, pSpaceDim, bcRegion)
{
}

void ExtrapOrder0BC::v_Apply(
        int                                                 bcRegion,
        int                                                 cnt,
        Array<OneD, Array<OneD, NekDouble> >               &Fwd,
        Array<OneD, Array<OneD, NekDouble> >               &physarray,
        const NekDouble                                    &time)
{
    int i, j;
    int e, pnt;
    int id1, id2, nBCEdgePts;
    int nVariables = physarray.num_elements();
    int nDimensions = m_spacedim;

    const Array<OneD, const int> &traceBndMap
        = m_fields[0]->GetTraceBndMap();

    int eMax;

    eMax = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExpSize();

    // Loop on bcRegions
    for (e = 0; e < eMax; ++e)
    {
        nBCEdgePts = m_fields[0]->GetBndCondExpansions()[bcRegion]->
            GetExp(e)->GetTotPoints();
        id1 = m_fields[0]->GetBndCondExpansions()[bcRegion]->
            GetPhys_Offset(e) ;
        id2 = m_fields[0]->GetTrace()->GetPhys_Offset(traceBndMap[cnt+e]);

        // Loop on points of bcRegion 'e'
        for (i = 0; i < nBCEdgePts; i++)
        {
            pnt = id2+i;

            // Setting up bcs for density
            (m_fields[0]->GetBndCondExpansions()[bcRegion]->
                UpdatePhys())[id1+i] = Fwd[0][pnt];

            // Setting up bcs for velocities
            for (j = 1; j <=nDimensions; ++j)
            {
                (m_fields[j]->GetBndCondExpansions()[bcRegion]->
                 UpdatePhys())[id1+i] = Fwd[j][pnt];
            }

            // Setting up bcs for energy
            (m_fields[nVariables-1]->GetBndCondExpansions()[bcRegion]->
                UpdatePhys())[id1+i] = Fwd[nVariables-1][pnt];
        }
    }
}

}