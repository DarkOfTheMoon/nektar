///////////////////////////////////////////////////////////////////////////////
//
// File LymphaticPressureArea.h
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
// Description: ArterialPressureArea header
//
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_ARTERIALPRESSUREAREA_H
#define NEKTAR_ARTERIALPRESSUREAREA_H

#include <string>
#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <PulseWaveSolver/EquationSystems/PulseWavePressureArea.h>

namespace Nektar
{
    // Forward declarations
    class ArterialPressureArea;

    /// Pointer to a PulseWaveOutflow object.
    typedef boost::shared_ptr<ArterialPressureArea> ArterialPressureAreaSharedPtr;
    
    /// A global linear system.
    class ArterialPressureArea : public PulseWavePressureArea
    {
    public:
        /// Creates an instance of this class
      static PulseWavePressureAreaSharedPtr create(Array<OneD, MultiRegions::ExpListSharedPtr>& pVessel, 
                                               const LibUtilities::SessionReaderSharedPtr& pSession)
        {
            return MemoryManager<ArterialPressureArea>::AllocateSharedPtr(pVessel,pSession);
        }

        /// Name of class
        static std::string className;
        
        ArterialPressureArea(Array<OneD, MultiRegions::ExpListSharedPtr> pVessel, 
                 const LibUtilities::SessionReaderSharedPtr pSession); 

        virtual ~ArterialPressureArea();
    protected:

        virtual void v_ReadParameters(int nDomains, int nqTrace);

        virtual void v_GetPacons(int omega, Array<OneD, Array<OneD, NekDouble> > &pacons_trace, int nqTrace);

        /// Pressure-area relationship definition
        virtual void v_getPressure(
            NekDouble A, 
            Array<OneD, NekDouble> pacons, 
            NekDouble &pressure);
        
        /// Definition of dp/da
        virtual void v_getdpda(
            NekDouble A, 
            NekDouble u, 
            Array<OneD, NekDouble> pacons, 
            NekDouble &dpda);
        
        /// Calculates A from pressure
        virtual void v_getAfromPressure(
            NekDouble pressure, 
            Array<OneD, NekDouble> pacons, 
            NekDouble &A);
        
        /// Calculates the two characteristic variables
        virtual void v_getW(
            NekDouble A, 
            NekDouble u, 
            Array<OneD, NekDouble> pacons, 
            NekDouble &W1, 
            NekDouble &W2);
        
        /// Calculates A from W_1 and W_2
        virtual void v_getAfromChars(
            NekDouble W_1, 
            NekDouble W_2, 
            Array<OneD, NekDouble> pacons, 
            NekDouble &A);
        
        /// Definition of u in terms of W_1 and W_2
        virtual void v_getufromChars(
            NekDouble W_1, 
            NekDouble W_2, 
            Array<OneD, NekDouble> pacons, 
            NekDouble &u);
        
        /// lambda_1 relationship definition
        virtual void v_getlambda1(
            NekDouble A, 
            NekDouble u, 
            Array<OneD, NekDouble> pacons, 
            NekDouble &lambda_1);
        
        /// lambda_2 relationship definition
        virtual void v_getlambda2(
            NekDouble A, 
            NekDouble u, 
            Array<OneD, NekDouble> pacons, 
            NekDouble &lambda_2);
        
        /// Definition of dp/du
        virtual void v_getdpdu(
            NekDouble A, 
            NekDouble u, 
            Array<OneD, NekDouble> pacons, 
            NekDouble &dpdu);
        
        /// Definition of dW1/da
        virtual void v_getdW1da(
            NekDouble A, 
            NekDouble u, 
            Array<OneD, NekDouble> pacons, 
            NekDouble &dW1da);
        
        /// Definition of dW1/du
        virtual void v_getdW1du(
            NekDouble A, 
            NekDouble u, 
            Array<OneD, NekDouble> pacons, 
            NekDouble &dW1du);
        
        /// Definition of dW2/dA
        virtual void v_getdW2da(
            NekDouble A, 
            NekDouble u, 
            Array<OneD, NekDouble> pacons, 
            NekDouble &dW2da);
        
        /// Definition of dW2/du
        virtual void v_getdW2du(
            NekDouble A, 
            NekDouble u, 
            Array<OneD, NekDouble> pacons, 
            NekDouble &dW2du);
        
        /// Evaluates the integral \int^A_{A_{0}}{ \sqrt{\frac{1}{\rho A} \frac{\partial p\left(A,t\right)}{\partial A}} dA}
        virtual void v_getCharIntegral(
            NekDouble A, 
            NekDouble u, 
            Array<OneD, NekDouble> pacons, 
            NekDouble &intergralTerm);
        
        //void solveAxb(int rows, Array<OneD, double> A_buf, Array<OneD, double> b_buf);
        virtual void v_solveAxb(
            int rows, 
            const Array<OneD, NekDouble> &matrix_buf, 
            const Array<OneD, NekDouble> &b_buf, 
            Array<OneD, NekDouble> &x);

        Array<OneD, Array<OneD, NekDouble> > m_beta;
        Array<OneD, Array<OneD, NekDouble> > m_beta_trace;

    private:

    };
}

#endif
