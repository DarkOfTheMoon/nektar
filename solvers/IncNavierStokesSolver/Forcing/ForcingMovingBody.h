///////////////////////////////////////////////////////////////////////////////
//
// File: ForcingMovingBody.h
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
// Description: Moving Body (Wavyness and acceleration)
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_FORCINGMOVINGBODY
#define NEKTAR_SOLVERUTILS_FORCINGMOVINGBODY

#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/FFT/NektarFFT.h>
#include <SolverUtils/SolverUtilsDeclspec.h>
#include <SolverUtils/Forcing/Forcing.h>
#include <IncNavierStokesSolver/Filters/FilterMovingBody.h>
#include <GlobalMapping/Mapping.h>

namespace Nektar
{

class ForcingMovingBody : public SolverUtils::Forcing
{
    public:

        friend class MemoryManager<ForcingMovingBody>;

        /// Creates an instance of this class
        static SolverUtils::ForcingSharedPtr create(
                const LibUtilities::SessionReaderSharedPtr& pSession,
                const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
                const unsigned int& pNumForcingFields,
                const TiXmlElement* pForce)
        {
            SolverUtils::ForcingSharedPtr p =
                                    MemoryManager<ForcingMovingBody>::
                                            AllocateSharedPtr(pSession);
            p->InitObject(pFields, pNumForcingFields, pForce);
            return p;
        }

        ///Name of the class
        static std::string className;

    protected:
        // Mapping object
        GlobalMapping::MappingSharedPtr               m_mapping;

        virtual void v_InitObject(
            const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
            const unsigned int&                         pNumForcingFields,
            const TiXmlElement*                         pForce);

        virtual void v_Apply(
            const Array<OneD, MultiRegions::ExpListSharedPtr>& fields,
            const Array<OneD, Array<OneD, NekDouble> >& inarray,
                  Array<OneD, Array<OneD, NekDouble> >& outarray,
            const NekDouble&                            time);

    private:

        ForcingMovingBody(
            const LibUtilities::SessionReaderSharedPtr& pSession);

        void CheckIsFromFile(const TiXmlElement* pForce);

        void InitialiseCableModel(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields);

        void InitialiseFilter(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
            const TiXmlElement* pForce);

        void ModalDecompositionMethod(
            const Array<OneD, NekDouble> &Forces,
                  Array<OneD, Array<OneD, NekDouble> > &vib);

        void FiniteElementMethod(
                  Array<OneD, NekDouble> &Forces,
                  Array<OneD, Array<OneD, NekDouble> > &vibs,
                  Array<OneD, Array<OneD, NekDouble> > &rots);

        void StructDynModel(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
                  Array<OneD, Array<OneD, NekDouble> > &HydroForces,
                  NekDouble time );

        void SetDynEqCoeffMatrix();

        void SetStiffnessMatrix();

        void RollOver(Array<OneD, Array<OneD, NekDouble> > &input);

        int m_movingBodyCalls;     ///< number of times the movbody have been called
        int m_np;                  ///< number of planes per processors
        int m_vdim;                ///< vibration dimension

        NekDouble m_structrho;     ///< mass of the cable per unit length
        NekDouble m_structdamp;    ///< damping ratio of the cable
        NekDouble m_length;        ///< length ratio of the cable
        NekDouble m_timestep;      ///< time step
        ///
        LibUtilities::NektarFFTSharedPtr m_FFT;
        ///
        FilterMovingBodySharedPtr m_MovBodyfilter;
        /// storage the loading on the cable
        Array<OneD, Array<OneD, NekDouble> > m_q;
        /// storage for the cable's motion(x,y) variables:vibrations
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > > m_vib;
        /// storage for the cable's motion(x,y) variables:rotations
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > > m_rot;

        /// fictitious velocity storage
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > > m_fV;
        /// fictitious acceleration storage
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > > m_fA;
        /// matrices in Modal Decompostion method
        Array<OneD, DNekMatSharedPtr> m_CoeffMat_A;
        Array<OneD, DNekMatSharedPtr> m_CoeffMat_B;
        /// matrices in Finite Element method
        DNekMatSharedPtr m_CoeffMat_K;
		DNekMatSharedPtr m_CoeffMat_M;
		DNekMatSharedPtr m_CoeffMat_D;
        /// [0] is displacements, [1] is velocities, [2] is accelerations
        Array<OneD, std::string> m_funcName;
        /// motion direction: [0] is 'x' and [1] is 'y'
        Array<OneD, std::string> m_motion;
        /// do determine if the the body motion come from an extern file
        Array<OneD, bool>        m_IsFromFile;
        /// Store the derivatives of motion variables in x-direction
        Array<OneD, Array< OneD, NekDouble> > m_zta;
        /// Store the derivatives of motion variables in y-direction
        Array<OneD, Array< OneD, NekDouble> > m_eta;
};

}

#endif
