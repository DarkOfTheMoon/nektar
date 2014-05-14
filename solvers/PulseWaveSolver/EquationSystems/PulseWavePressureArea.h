///////////////////////////////////////////////////////////////////////////////
//
// File PulseWavePressureArea.h
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
// Description: PulseWavePressureArea header
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_PULSEWAVEPRESSUREAREA_H
#define NEKTAR_PULSEWAVEPRESSUREAREA_H

#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <MultiRegions/ExpList.h>

namespace Nektar
{
    class PulseWavePressureArea;
    typedef boost::shared_ptr<PulseWavePressureArea>  PulseWavePressureAreaSharedPtr;
    
    static PulseWavePressureAreaSharedPtr NullPulseWavePressureAreaSharedPtr;

    typedef LibUtilities::NekFactory< std::string, 
        PulseWavePressureArea, 
        Array<OneD, MultiRegions::ExpListSharedPtr>&, 
        const LibUtilities::SessionReaderSharedPtr&,
        const int& > PressureAreaFactory;
    PressureAreaFactory& GetPressureAreaFactory();
   
    class PulseWavePressureArea
    {
    public:
        PulseWavePressureArea(Array<OneD, MultiRegions::ExpListSharedPtr> &pVessel,
                              const LibUtilities::SessionReaderSharedPtr &pSession,
                              const int &nDomains);

        virtual ~PulseWavePressureArea();

        inline void ReadParameters(int omega, int nqTrace, MultiRegions::ExpListSharedPtr &field);
        inline void GetPacons(int omega, Array<OneD, Array<OneD, NekDouble> > &pacons_trace, int nqTrace);
        inline void getPressure(NekDouble A, Array<OneD, NekDouble> pacons, NekDouble &pressure);
        inline void getdpda(NekDouble A, NekDouble u, Array<OneD, NekDouble> pacons, NekDouble &dpda);
        inline void getAfromPressure(NekDouble pressure, Array<OneD, NekDouble> pacons, NekDouble &A);
        inline void getW(NekDouble A, NekDouble u, Array<OneD, NekDouble> pacons, NekDouble &W1, NekDouble &W2);
        inline void getAfromChars(NekDouble W_1, NekDouble W_2, Array<OneD, NekDouble> pacons, NekDouble &A);
        inline void getufromChars(NekDouble W_1, NekDouble W_2, Array<OneD, NekDouble> pacons, NekDouble &u);
        inline void getlambda1(NekDouble A, NekDouble u, Array<OneD, NekDouble> pacons, NekDouble &lambda_1);
        inline void getlambda2(NekDouble A, NekDouble u, Array<OneD, NekDouble> pacons, NekDouble &lambda_2);
        inline void getdpdu(NekDouble A, NekDouble u, Array<OneD, NekDouble> pacons, NekDouble &dpdu);
        inline void getdW1da(NekDouble A, NekDouble u, Array<OneD, NekDouble> pacons, NekDouble &dW1da);
        inline void getdW1du(NekDouble A, NekDouble u, Array<OneD, NekDouble> pacons, NekDouble &dW1du);
        inline void getdW2da(NekDouble A, NekDouble u, Array<OneD, NekDouble> pacons, NekDouble &dW2da);
        inline void getdW2du(NekDouble A, NekDouble u, Array<OneD, NekDouble> pacons, NekDouble &dW2du);
        inline void getCharIntegral(NekDouble A, NekDouble u, Array<OneD, NekDouble> pacons, NekDouble &intergralTerm);
        inline void solveAxb(int rows, const Array<OneD, NekDouble> &matrix_buf, const Array<OneD, NekDouble> &b_buf, Array<OneD, NekDouble> &x);

    protected:
        virtual void v_ReadParameters(int omega, int nqTrace, MultiRegions::ExpListSharedPtr &field)=0;
        virtual void v_GetPacons(int omega, Array<OneD, Array<OneD, NekDouble> > &pacons_trace, int nqTrace)=0;
        virtual void v_getPressure(NekDouble A, Array<OneD, NekDouble> pacons, NekDouble &pressure)=0;
        virtual void v_getdpda(NekDouble A, NekDouble u, Array<OneD, NekDouble> pacons, NekDouble &dpda)=0;
        virtual void v_getAfromPressure(NekDouble pressure, Array<OneD, NekDouble> pacons, NekDouble &A)=0;
        virtual void v_getW(NekDouble A, NekDouble u, Array<OneD, NekDouble> pacons, NekDouble &W1, NekDouble &W2)=0;
        virtual void v_getAfromChars(NekDouble W_1, NekDouble W_2, Array<OneD, NekDouble> pacons, NekDouble &A)=0;
        virtual void v_getufromChars(NekDouble W_1, NekDouble W_2, Array<OneD, NekDouble> pacons, NekDouble &u)=0;
        virtual void v_getlambda1(NekDouble A, NekDouble u, Array<OneD, NekDouble> pacons, NekDouble &lambda_1)=0;
        virtual void v_getlambda2(NekDouble A, NekDouble u, Array<OneD, NekDouble> pacons, NekDouble &lambda_2)=0;
        virtual void v_getdpdu(NekDouble A, NekDouble u, Array<OneD, NekDouble> pacons, NekDouble &dpdu)=0;
        virtual void v_getdW1da(NekDouble A, NekDouble u, Array<OneD, NekDouble> pacons, NekDouble &dW1da)=0;
        virtual void v_getdW1du(NekDouble A, NekDouble u, Array<OneD, NekDouble> pacons, NekDouble &dW1du)=0;
        virtual void v_getdW2da(NekDouble A, NekDouble u, Array<OneD, NekDouble> pacons, NekDouble &dW2da)=0;
        virtual void v_getdW2du(NekDouble A, NekDouble u, Array<OneD, NekDouble> pacons, NekDouble &dW2du)=0;
        virtual void v_getCharIntegral(NekDouble A, NekDouble u, Array<OneD, NekDouble> pacons, NekDouble &intergralTerm)=0;
        virtual void v_solveAxb(int rows, const Array<OneD, NekDouble> &matrix_buf, const Array<OneD, NekDouble> &b_buf, Array<OneD, NekDouble> &x)=0;

        Array<OneD, MultiRegions::ExpListSharedPtr> m_vessels;
	LibUtilities::SessionReaderSharedPtr m_session;
        const int m_domains;

        NekDouble m_rho;
        NekDouble m_time;

        void EvaluateFunction(
            Array<OneD, MultiRegions::ExpListSharedPtr> pFields,
            LibUtilities::SessionReaderSharedPtr pSession,
            std::string pFieldName, 
            Array<OneD, NekDouble>& pArray,
            const std::string& pFunctionName,
            NekDouble pTime = NekDouble(0),
            const int domain = 0,
            const int nq=0);

    private:
    };

    inline void PulseWavePressureArea::ReadParameters(int omega, int nqTrace, MultiRegions::ExpListSharedPtr &field)
    {
        v_ReadParameters(omega, nqTrace, field);
    }

    inline void PulseWavePressureArea::GetPacons(int omega, Array<OneD, Array<OneD, NekDouble> > &pacons_trace, int nqTrace)
    {
        v_GetPacons(omega, pacons_trace, nqTrace);
    }

    inline void PulseWavePressureArea::getPressure(NekDouble A, Array<OneD, NekDouble> pacons, NekDouble &pressure)
    {
        v_getPressure(A, pacons, pressure);
    };
    inline void PulseWavePressureArea::getdpda(NekDouble A, NekDouble u, Array<OneD, NekDouble> pacons, NekDouble &dpda)
    {
        v_getdpda(A, u, pacons, dpda);
    };
    inline void PulseWavePressureArea::getAfromPressure(NekDouble pressure, Array<OneD, NekDouble> pacons, NekDouble &A)
    {
        v_getAfromPressure(pressure, pacons, A);
    };
    inline void PulseWavePressureArea::getW(NekDouble A, NekDouble u, Array<OneD, NekDouble> pacons, NekDouble &W1, NekDouble &W2)
    {
        v_getW(A, u, pacons, W1, W2);
    };
    inline void PulseWavePressureArea::getAfromChars(NekDouble W_1, NekDouble W_2, Array<OneD, NekDouble> pacons, NekDouble &A)
    {
        v_getAfromChars(W_1, W_2, pacons, A);
    };
    inline void PulseWavePressureArea::getufromChars(NekDouble W_1, NekDouble W_2, Array<OneD, NekDouble> pacons, NekDouble &u)
    {
        v_getufromChars(W_1, W_2, pacons, u);
    };
    inline void PulseWavePressureArea::getlambda1(NekDouble A, NekDouble u, Array<OneD, NekDouble> pacons, NekDouble &lambda_1)
    {
        v_getlambda1(A, u, pacons, lambda_1);
    };
    inline void PulseWavePressureArea::getlambda2(NekDouble A, NekDouble u, Array<OneD, NekDouble> pacons, NekDouble &lambda_2)
    {
        v_getlambda2(A, u, pacons, lambda_2);
    };
    inline void PulseWavePressureArea::getdpdu(NekDouble A, NekDouble u, Array<OneD, NekDouble> pacons, NekDouble &dpdu)
    {
        v_getdpdu(A, u, pacons, dpdu);
    };
    inline void PulseWavePressureArea::getdW1da(NekDouble A, NekDouble u, Array<OneD, NekDouble> pacons, NekDouble &dW1da)
    {
        v_getdW1da(A, u, pacons, dW1da);
    };
    inline void PulseWavePressureArea::getdW1du(NekDouble A, NekDouble u, Array<OneD, NekDouble> pacons, NekDouble &dW1du)
    {
        v_getdW1du(A, u, pacons, dW1du);
    };
    inline void PulseWavePressureArea::getdW2da(NekDouble A, NekDouble u, Array<OneD, NekDouble> pacons, NekDouble &dW2da)
    {
        v_getdW2da(A, u, pacons, dW2da);
    };
    inline void PulseWavePressureArea::getdW2du(NekDouble A, NekDouble u, Array<OneD, NekDouble> pacons, NekDouble &dW2du)
    {
        v_getdW2du(A, u, pacons, dW2du);
    };
    inline void PulseWavePressureArea::getCharIntegral(NekDouble A, NekDouble u, Array<OneD, NekDouble> pacons, NekDouble &intergralTerm)
    {
        v_getCharIntegral(A, u, pacons, intergralTerm);
    };
    inline void PulseWavePressureArea::solveAxb(int rows, const Array<OneD, NekDouble> &matrix_buf, const Array<OneD, NekDouble> &b_buf, Array<OneD, NekDouble> &x)
    {
        v_solveAxb(rows, matrix_buf, b_buf, x);
    };
}
#endif
