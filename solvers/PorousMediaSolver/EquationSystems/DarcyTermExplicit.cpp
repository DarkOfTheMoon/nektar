///////////////////////////////////////////////////////////////////////////////
//
// File: DarcyTermExplcit.cpp
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
// Description: Abstract base class for DarcyTermExplicit.
//
///////////////////////////////////////////////////////////////////////////////

#include <PorousMediaSolver/EquationSystems/DarcyTermExplicit.h>
#include <LibUtilities/Communication/Comm.h>

namespace Nektar
{
    /**
     * Registers the class with the Factory.
     */
    std::string DarcyTermExplicit::className = GetDarcyTermFactory().RegisterCreatorFunction(
        "Explicit Term",
        DarcyTermExplicit::create,
        "Explicit Term");

    DarcyTermExplicit::DarcyTermExplicit(
        const LibUtilities::SessionReaderSharedPtr pSession,
        Array<OneD, MultiRegions::ExpListSharedPtr> pFields)
        : DarcyTerm(pSession,pFields)
    {
    }

    DarcyTermExplicit::~DarcyTermExplicit()
    {
    }
    
    /** 
     * 
     */
    void DarcyTermExplicit::v_EvaluateDarcyTerm(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray, 
        Array<OneD, Array<OneD, NekDouble> > &outarray,
        NekDouble kinvis)
    {
        int nqtot = m_fields[0]->GetTotPoints();
        int nDim = m_fields.num_elements()-1;

        if (nDim == 2)
        {
            Vmath::Svtvp(nqtot,-kinvis*m_perm_inv[0],inarray[0],1,outarray[0],1,outarray[0],1);
            Vmath::Svtvp(nqtot,-kinvis*m_perm_inv[2],inarray[0],1,outarray[1],1,outarray[1],1);
            
            Vmath::Svtvp(nqtot,-kinvis*m_perm_inv[2],inarray[1],1,outarray[0],1,outarray[0],1);
            Vmath::Svtvp(nqtot,-kinvis*m_perm_inv[1],inarray[1],1,outarray[1],1,outarray[1],1);
        }
        else
        {
            Vmath::Svtvp(nqtot,-kinvis*m_perm_inv[0],inarray[0],1,outarray[0],1,outarray[0],1);
            Vmath::Svtvp(nqtot,-kinvis*m_perm_inv[3],inarray[0],1,outarray[1],1,outarray[1],1);
            Vmath::Svtvp(nqtot,-kinvis*m_perm_inv[4],inarray[0],1,outarray[2],1,outarray[2],1);
            
            Vmath::Svtvp(nqtot,-kinvis*m_perm_inv[3],inarray[1],1,outarray[0],1,outarray[0],1);
            Vmath::Svtvp(nqtot,-kinvis*m_perm_inv[1],inarray[1],1,outarray[1],1,outarray[1],1);
            Vmath::Svtvp(nqtot,-kinvis*m_perm_inv[5],inarray[1],1,outarray[2],1,outarray[2],1);
        
            Vmath::Svtvp(nqtot,-kinvis*m_perm_inv[4],inarray[2],1,outarray[0],1,outarray[0],1);
            Vmath::Svtvp(nqtot,-kinvis*m_perm_inv[5],inarray[2],1,outarray[1],1,outarray[1],1);
            Vmath::Svtvp(nqtot,-kinvis*m_perm_inv[2],inarray[2],1,outarray[2],1,outarray[2],1);
        }
    }


    /**
     * Registers the class with the Factory.
     */
    std::string DarcyTermExplicitSpatial::className = GetDarcyTermFactory().RegisterCreatorFunction(
        "Explicit spatially varying Darcy term",
        DarcyTermExplicit::create,
        "Explicit spatially varying Darcy term");

    DarcyTermExplicitSpatial::DarcyTermExplicitSpatial(
        const LibUtilities::SessionReaderSharedPtr pSession,
        Array<OneD, MultiRegions::ExpListSharedPtr> pFields)
        : DarcyTerm(pSession,pFields)
    {
    }

    DarcyTermExplicitSpatial::~DarcyTermExplicitSpatial()
    {
    }
    
    /** 
     * 
     */
    void DarcyTermExplicitSpatial::v_EvaluateDarcyTerm(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray, 
        Array<OneD, Array<OneD, NekDouble> > &outarray,
        NekDouble kinvis)
    {
        Array <OneD, NekDouble> tmp(nqtot);
        int nqtot = m_fields[0]->GetTotPoints();
        int nDim = m_fields.num_elements()-1;
        
        for(int i=0; i<nDim; ++i)
        {
            Vmath::Smul(nqtot,-kinvis,m_spatialperm[i],1,tmp,1);
            Vmath::Vvtvp(nqtot,tmp,1,inarray[i],1,outarray[i],1,outarray[i],1);
        }
    }

}

