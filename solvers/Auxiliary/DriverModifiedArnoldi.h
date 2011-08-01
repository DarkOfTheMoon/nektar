///////////////////////////////////////////////////////////////////////////////
//
// File DriverModifiedArnoldi.h
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
// Description: Driver class for eigenvalue analysis using the modified Arnoldi
//              method.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_DRIVERMODIFIEDARNOLDI_H
#define NEKTAR_SOLVERS_DRIVERMODIFIEDARNOLDI_H

#include <Auxiliary/DriverArnoldi.hpp>

namespace Nektar
{
    class DriverModifiedArnoldi: public DriverArnoldi
    {
    public:
        friend class MemoryManager<DriverModifiedArnoldi>;
        
        /// Creates an instance of this class
        static DriverSharedPtr create(LibUtilities::SessionReaderSharedPtr& pSession) {
            DriverSharedPtr p = MemoryManager<DriverModifiedArnoldi>::AllocateSharedPtr(pSession);
            p->InitObject();
            return p;
        }
	
        ///Name of the class
        static std::string className;
        
	protected:
        
        Array<OneD, MultiRegions::ExpListSharedPtr> m_fields;
        Array<OneD, MultiRegions::ExpListSharedPtr> m_forces;
        
        /// Constructor
        DriverModifiedArnoldi(LibUtilities::SessionReaderSharedPtr pSession);

        /// Destructor
        virtual ~DriverModifiedArnoldi();

        /// Virtual function for initialisation implementation.
        virtual void v_InitObject();

        /// Virtual function for solve implementation.
        virtual void v_Execute();

    private:
        void OutputEv(FILE* pFile, const int nev, Array<OneD, NekDouble> &workl, int* ipntr, NekDouble period, bool TimeSteppingAlgorithm);

        /// Generates a new vector in the sequence by applying the linear operator.
        void EV_update(Array<OneD, NekDouble> &src,
                       Array<OneD, NekDouble> &tgt);

        /// Generates the upper Hessenberg matrix H and computes its eigenvalues.
        void EV_small(Array<OneD, Array<OneD, NekDouble> > &Kseq,
                      const int ntot,
                      const Array<OneD, NekDouble> &alpha,
                      const int kdim,
                      Array<OneD, NekDouble> &zvec,
                      Array<OneD, NekDouble> &wr,
                      Array<OneD, NekDouble> &wi,
                      NekDouble &resnorm);

        /// Tests for convergence of eigenvalues of H.
        int EV_test(const int itrn,
                    const int kdim,
                    Array<OneD, NekDouble> &zvec,
                    Array<OneD, NekDouble> &wr,
                    Array<OneD, NekDouble> &wi,
                    const NekDouble resnorm,
                    const int nvec,
                    ofstream &evlout,
                    NekDouble &resid0);


        /// Sorts a sequence of eigenvectors/eigenvalues by magnitude.
        void EV_sort(Array<OneD, NekDouble> &evec,
                     Array<OneD, NekDouble> &wr,
                     Array<OneD, NekDouble> &wi,
                     Array<OneD, NekDouble> &test,
                     const int dim);

    };
	
} //end of namespace

#endif //NEKTAR_SOLVERS_AUXILIARY_ADRBASE_H
