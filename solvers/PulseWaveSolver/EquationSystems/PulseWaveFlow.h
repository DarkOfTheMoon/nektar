///////////////////////////////////////////////////////////////////////////////
//
// File PulseWaveFlow.h
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
// Description: PulseWaveFlow header
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_PULSEWAVEFLOW_H
#define NEKTAR_PULSEWAVEFLOW_H

#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <MultiRegions/ExpList.h>

namespace Nektar
{
    class PulseWaveFlow;
    typedef boost::shared_ptr<PulseWaveFlow>  PulseWaveFlowSharedPtr;
    
    static PulseWaveFlowSharedPtr NullPulseWaveFlowSharedPtr;

    typedef LibUtilities::NekFactory< std::string, 
        PulseWaveFlow, 
      Array<OneD, MultiRegions::ExpListSharedPtr>&, 
      const LibUtilities::SessionReaderSharedPtr& > FlowFactory;
    FlowFactory& GetFlowFactory();
    
    class PulseWaveFlow
    {
    public:
PulseWaveFlow(Array<OneD, MultiRegions::ExpListSharedPtr> pVessel,
		const LibUtilities::SessionReaderSharedPtr pSession);

        virtual ~PulseWaveFlow();

        inline void DoBoundary(int omega, int n);

    protected:
        virtual void v_DoBoundary(int omega, int n) = 0;

        Array<OneD, MultiRegions::ExpListSharedPtr> m_vessels;
	LibUtilities::SessionReaderSharedPtr m_session;

    private:
    };

    /**
     *
     */
    inline void PulseWaveFlow::DoBoundary(int omega, int n)
    {
        v_DoBoundary(omega, n);
    }



}
#endif