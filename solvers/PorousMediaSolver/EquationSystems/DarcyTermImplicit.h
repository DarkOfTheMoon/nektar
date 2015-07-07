///////////////////////////////////////////////////////////////////////////////
//
// File: DarcyTermImplicit.h
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
// Description: Abstract base class for DarcyTermImplicit.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_DARCYTERMIMPLICIT_H
#define NEKTAR_SOLVERS_DARCYTERMIMPLICIT_H

//#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/Memory/NekMemoryManager.hpp>
//#include <LibUtilities/BasicUtils/SessionReader.h>
#include <MultiRegions/ExpList.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <PorousMediaSolver/EquationSystems/DarcyTerm.h>

namespace Nektar
{
    //--------------------------
    // Implicit Darcy Term class
    // -------------------------
    
    //class DarcyTermImplicitIsotropic;
    
    //typedef boost::shared_ptr<DarcyTermImplicitIsotropic> DarcyTermImplicitSharedPtr;
    
    class DarcyTermImplicitIsotropic : public DarcyTerm
    {
    public:
        friend class MemoryManager<DarcyTermImplicitIsotropic>;

        /// Creates an instance of this class
        static DarcyTermSharedPtr create(
            const LibUtilities::SessionReaderSharedPtr &pSession,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields)
        {
            DarcyTermSharedPtr p = MemoryManager<DarcyTermImplicitIsotropic>::AllocateSharedPtr(pSession,pFields);
            return p;
        }

        /// Name of class
        static std::string className;

        DarcyTermImplicitIsotropic(
            const LibUtilities::SessionReaderSharedPtr &pSession,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields);

        virtual ~DarcyTermImplicitIsotropic();
        
    protected:
        virtual void v_EvaluateDarcyTerm(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray, 
            Array<OneD, Array<OneD, NekDouble> > &outarray,
            NekDouble kinvis);

        virtual void v_AddDarcyPressureTerm(
            int nq,
            NekDouble kinvis,
            Array<OneD, NekDouble> &Q, 
            Array<OneD, const NekDouble> &Vel,
            int i);

        virtual void v_SetupPermeability();

        virtual void v_GetImplicitDarcyFactor(
            Array<OneD, NekDouble> &permCoeff);
    };

    //--------------------------
    // Anisotropic Implicit Darcy Term class
    // -------------------------
    
    //class DarcyTermImplicitAnisotropic;
    
    //typedef boost::shared_ptr<DarcyTermImplicitAnisotropic> DarcyTermImplicitAnisotropicSharedPtr;
    
    class DarcyTermImplicitAnisotropic : public DarcyTerm
    {
    public:
        friend class MemoryManager<DarcyTermImplicitAnisotropic>;

        /// Creates an instance of this class
        static DarcyTermSharedPtr create(
            const LibUtilities::SessionReaderSharedPtr &pSession,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields)
        {
            DarcyTermSharedPtr p = MemoryManager<DarcyTermImplicitAnisotropic>::AllocateSharedPtr(pSession,pFields);
            return p;
        }

        /// Name of class
        static std::string className;

        DarcyTermImplicitAnisotropic(
            const LibUtilities::SessionReaderSharedPtr &pSession,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields);

        virtual ~DarcyTermImplicitAnisotropic();
        
    protected:
        virtual void v_EvaluateDarcyTerm(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray, 
            Array<OneD, Array<OneD, NekDouble> > &outarray,
            NekDouble kinvis);

        virtual void v_AddDarcyPressureTerm(
            int nq,
            NekDouble kinvis,
            Array<OneD, NekDouble> &Q, 
            Array<OneD, const NekDouble> &Vel,
            int i);

        virtual void v_SetupPermeability();

        virtual void v_GetImplicitDarcyFactor(
            Array<OneD, NekDouble> &permCoeff);
    };
    
}

#endif

