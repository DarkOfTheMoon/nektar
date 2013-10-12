///////////////////////////////////////////////////////////////////////////////
//
// File: DarcyTerm.h
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
// Description: Abstract base class for Darcy Term.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_DARCYTERM_H
#define NEKTAR_SOLVERS_DARCYTERM_H

#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <MultiRegions/ExpList.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>

namespace Nektar
{
    // Forward declaration
    class DarcyTerm;

    typedef boost::shared_ptr<DarcyTerm> DarcyTermSharedPtr;
    
    typedef LibUtilities::NekFactory< std::string, DarcyTerm,
        const LibUtilities::SessionReaderSharedPtr& ,
        Array<OneD, MultiRegions::ExpListSharedPtr>& > DarcyTermFactory; 

    DarcyTermFactory& GetDarcyTermFactory();


    class DarcyTerm
    {
    public:
        DarcyTerm(        
            const LibUtilities::SessionReaderSharedPtr pSession,
            Array<OneD, MultiRegions::ExpListSharedPtr> pFields);
        
        virtual ~DarcyTerm();

        inline void EvaluateDarcyTerm(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray, 
            Array<OneD, Array<OneD, NekDouble> > &outarray,
            NekDouble kinvis);

        inline void SetupPermeability();

        inline void GetImplicitDarcyFactor(
            Array<OneD, NekDouble> &permCoeff);
        
    protected:
        
        virtual void v_EvaluateDarcyTerm(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray, 
            Array<OneD, Array<OneD, NekDouble> > &outarray,
            NekDouble kinvis)=0;

        virtual void v_SetupPermeability()=0;

        virtual void v_GetImplicitDarcyFactor(
            Array<OneD, NekDouble> &permCoeff)=0;

        Array<OneD, NekDouble>  m_perm;  ///< Permeability matrix
        Array<OneD, Array<OneD, NekDouble> >  m_spatialperm;  ///< Permeability matrix
        Array<OneD, NekDouble>  m_perm_inv;  ///< inverse Permeability matrix

        
        LibUtilities::SessionReaderSharedPtr m_session;

        Array<OneD, MultiRegions::ExpListSharedPtr> m_fields;

        LibUtilities::CommSharedPtr m_comm;

    private:
    };

    /**
     *
     */
    inline void DarcyTerm::EvaluateDarcyTerm(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray, 
        Array<OneD, Array<OneD, NekDouble> > &outarray,
        NekDouble kinvis)
    {
        v_EvaluateDarcyTerm(inarray,outarray,kinvis);
    }

    /**
     *
     */
    inline void DarcyTerm::SetupPermeability()
    {
        v_SetupPermeability();
    }

    /**
     *
     */
    inline void DarcyTerm::GetImplicitDarcyFactor(
        Array<OneD, NekDouble> &permCoeff)
    {
        v_GetImplicitDarcyFactor(permCoeff);
    }


}

#endif

