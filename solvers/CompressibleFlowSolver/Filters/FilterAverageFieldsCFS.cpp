///////////////////////////////////////////////////////////////////////////////
//
// File FilterAverageField.cpp
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
// Description: Average solution fields during time-stepping.
//
///////////////////////////////////////////////////////////////////////////////
#include <CompressibleFlowSolver/Filters/FilterAverageFieldsCFS.h>

namespace Nektar
{

std::string FilterAverageFieldsCFS::className =
    SolverUtils::GetFilterFactory().RegisterCreatorFunction(
                "AverageFieldsCFS", FilterAverageFieldsCFS::create);

FilterAverageFieldsCFS::FilterAverageFieldsCFS(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const ParamMap &pParams) :
    FilterAverageFields(pSession, pParams)
{
    m_session->LoadParameter("Gamma", m_gamma, 1.4);
}

FilterAverageFieldsCFS::~FilterAverageFieldsCFS()
{
}

void FilterAverageFieldsCFS::v_ProcessSample(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time)
{
    int nfield = pFields.num_elements();
    Array<OneD, Array<OneD, NekDouble> > velocity(nfield - 2);
    int ncoeff = pFields[0]->GetNcoeffs();
    
    for(int n = 0; n < nfield - 2; ++n)
    {
        velocity[n] = Array<OneD, NekDouble>(ncoeff, 0.0);
    }
    
    /**************************************************************************/
    //Change for Compressible Flow Solver
    //Rho calculation
    Vmath::Vadd(ncoeff,
                pFields[0]->GetCoeffs(),1,
                m_outFields[0],1,
                m_outFields[0],1);
    
    //Velocity components calculation
    for(int n = 0; n < nfield - 2; ++n)
    {
        Vmath::Vdiv(ncoeff,
                    pFields[n+1]->GetCoeffs(),1,
                    pFields[0]->GetCoeffs(),1,
                    velocity[n],1);
        
        Vmath::Vadd(ncoeff,
                    velocity[n], 1,
                    m_outFields[n+1], 1,
                    m_outFields[n+1], 1);
    }
    
    Array<OneD, NekDouble> tmp(ncoeff, 0.0);
    Array<OneD, NekDouble> pressure(ncoeff, 0.0);
    
    //Calculate pressure
    for (int n = 0; n < nfield - 2; n++)
    {
        Vmath::Vmul(ncoeff,
                    velocity[n], 1,
                    velocity[n], 1,
                    tmp,1);
        
        
        Vmath::Smul(ncoeff, 0.5,
                    tmp, 1,
                    tmp, 1);
        
        Vmath::Vadd(ncoeff,
                    pressure, 1,
                    tmp, 1,
                    pressure, 1);
    }
    
    Vmath::Vmul(ncoeff,
                pressure, 1,
                pFields[0]->GetCoeffs(), 1,
                pressure,1);
    
    Vmath::Vsub(ncoeff,
                pFields[nfield - 1]->GetCoeffs(), 1,
                pressure, 1,
                pressure,1);
    
    NekDouble gammaMinusOne    = m_gamma - 1.0;
    
    Vmath::Smul(ncoeff, gammaMinusOne,
                pressure, 1,
                pressure, 1);
    
    Vmath::Vadd(ncoeff,
                pressure, 1,
                m_outFields[nfield - 1], 1,
                m_outFields[nfield - 1], 1);
    
    
    /**************************************************************************/
}
    
}
