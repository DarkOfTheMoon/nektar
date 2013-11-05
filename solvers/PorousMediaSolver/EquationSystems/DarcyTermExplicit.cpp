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
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields)
        : DarcyTerm(pSession,pFields)
    {
    }

    DarcyTermExplicit::~DarcyTermExplicit()
    {
    }

    /** 
     * 
     */
    void DarcyTermExplicit::v_SetupPermeability()
    {
        int nDim = m_fields.num_elements()-1; 
        switch(nDim)
        {
            case 2:
            {
                std::string varCoeffs[3] = {
                    "kxx",
                    "kyy",
                    "kxy"
                };
                for (int i = 0; i < (3*(nDim-1)); ++i)
                {
                    ASSERTL0(m_session->DefinesFunction("AnisotropicPermeability", varCoeffs[i]),
                             "Function '" + varCoeffs[i] + "' not correctly defined.");
                    m_perm[i] = m_session->GetFunction("AnisotropicPermeability", varCoeffs[i])->Evaluate();
                }
                
                NekDouble detTemp = m_perm[0]*m_perm[1]-m_perm[2]*m_perm[2];
                
                // Check if permeability matrix is positive definite
                ASSERTL0(m_perm[0] > 0,"Permeability Matrix is not positive definite");
                ASSERTL0(detTemp > 0,"Permeability Matrix is not positive definite");
                
                m_perm_inv[0] = m_perm[1];
                m_perm_inv[1] = m_perm[0];
                m_perm_inv[2] = -m_perm[2];
                Vmath::Smul(3, 1/detTemp, m_perm_inv, 1, m_perm_inv, 1);
            }
            break;
            case 3:
            {
                std::string varCoeffs[6] = {
                    "kxx",
                    "kyy",
                    "kzz",
                    "kxy",
                    "kxz",
                    "kyz"
                };                
                for (int i = 0; i < (3*(nDim-1)); ++i)
                {
                    ASSERTL0(m_session->DefinesFunction("AnisotropicPermeability", varCoeffs[i]),
                             "Function '" + varCoeffs[i] + "' not correctly defined.");
                    m_perm[i] = m_session->GetFunction("AnisotropicPermeability", varCoeffs[i])->Evaluate();
                }

                NekDouble detTemp = m_perm[0]*(m_perm[1]*m_perm[2]-m_perm[5]*m_perm[5])
                    -m_perm[3]*(m_perm[2]*m_perm[3]-m_perm[4]*m_perm[5])
                    +m_perm[4]*(m_perm[3]*m_perm[5]-m_perm[1]*m_perm[4]);
                
                // Check if permeability matrix is positive definite
                ASSERTL0(m_perm[0] > 0,"Permeability Matrix is not positive definite");
                NekDouble pd_chk = m_perm[0]*m_perm[1]-m_perm[3]*m_perm[3];
                ASSERTL0(pd_chk > 0,"Permeability Matrix is not positive definite");
                ASSERTL0(detTemp > 0,"Permeability Matrix is not positive definite");

                m_perm_inv[0] = m_perm[1]*m_perm[2]-m_perm[5]*m_perm[5];
                m_perm_inv[1] = m_perm[0]*m_perm[2]-m_perm[4]*m_perm[4];
                m_perm_inv[2] = m_perm[0]*m_perm[1]-m_perm[3]*m_perm[3];
                m_perm_inv[3] = m_perm[4]*m_perm[5]-m_perm[2]*m_perm[3];
                m_perm_inv[4] = m_perm[3]*m_perm[5]-m_perm[1]*m_perm[4];
                m_perm_inv[5] = m_perm[3]*m_perm[4]-m_perm[0]*m_perm[5];
                
                Vmath::Smul(6, 1/detTemp, m_perm_inv, 1, m_perm_inv, 1);
            }
            break;
            default:
                ASSERTL0(0,"Dimension not supported");
                break;
        }
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

    void DarcyTermExplicit::v_GetImplicitDarcyFactor(
        Array<OneD, NekDouble> &permCoeff)
    {
    }

    /**
     * Registers the class with the Factory.
     */
    std::string DarcyTermExplicitSpatial::className = GetDarcyTermFactory().RegisterCreatorFunction(
        "Explicit spatially varying Darcy term",
        DarcyTermExplicit::create,
        "Explicit spatially varying Darcy term");

    DarcyTermExplicitSpatial::DarcyTermExplicitSpatial(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields)
        : DarcyTerm(pSession,pFields)
    {
    }

    DarcyTermExplicitSpatial::~DarcyTermExplicitSpatial()
    {
    }
    

    /** 
     * 
     */
    void DarcyTermExplicitSpatial::v_SetupPermeability()
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
        int nqtot = m_fields[0]->GetTotPoints();
        Array <OneD, NekDouble> tmp(nqtot);
        int nDim = m_fields.num_elements()-1;
        
        for(int i=0; i<nDim; ++i)
        {
            Vmath::Smul(nqtot,-kinvis,m_spatialperm[i],1,tmp,1);
            Vmath::Vvtvp(nqtot,tmp,1,inarray[i],1,outarray[i],1,outarray[i],1);
        }
    }

    void DarcyTermExplicitSpatial::v_GetImplicitDarcyFactor(
        Array<OneD, NekDouble> &permCoeff)
    {
    }


}

