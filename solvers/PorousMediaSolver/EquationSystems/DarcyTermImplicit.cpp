///////////////////////////////////////////////////////////////////////////////
//
// File: DarcyTermImplicit.cpp
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

#include <PorousMediaSolver/EquationSystems/DarcyTermImplicit.h>
#include <LibUtilities/Communication/Comm.h>

namespace Nektar
{
    /**
     * Registers the class with the Factory.
     */
    std::string DarcyTermImplicitIsotropic::className = GetDarcyTermFactory().RegisterCreatorFunction(
        "ImplicitIsotropic",
        DarcyTermImplicitIsotropic::create,
        "Implicit Term Isotropic");

    DarcyTermImplicitIsotropic::DarcyTermImplicitIsotropic(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields)
        : DarcyTerm(pSession,pFields)
    {
    }

    /** 
     * 
     */
    DarcyTermImplicitIsotropic::~DarcyTermImplicitIsotropic()
    {
    }

    /** 
     * 
     */
    void DarcyTermImplicitIsotropic::v_SetupPermeability()
    {
        int nDim = m_fields.num_elements()-1;
        NekDouble kTemp;
        m_session->LoadParameter("Permeability", kTemp);

        m_perm = Array<OneD, NekDouble> (6);
        m_perm_inv = Array<OneD, NekDouble> (6);
            
        for (int i = 0; i < nDim; ++i)
        {
            m_perm[i] = kTemp;
        }

        for (int i = (nDim+1); i < (3*(nDim-1)); ++i)
        {
            m_perm[i] = 0;
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
    
    /** 
     * 
     */
    void DarcyTermImplicitIsotropic::v_EvaluateDarcyTerm(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray, 
        Array<OneD, Array<OneD, NekDouble> > &outarray,
        NekDouble kinvis)
    {
    }

    void DarcyTermImplicitIsotropic::v_GetImplicitDarcyFactor(
        Array<OneD, NekDouble> &permCoeff)
    {
        int nDim = m_fields.num_elements()-1;
        for (int i=0; i<nDim; ++i)
        {
            permCoeff[i]=m_perm_inv[i];
        }
    }

    /**
     * Additional contribution to high order pressure boundary condition
     * due to darcy term
     */
    void DarcyTermImplicitIsotropic::v_AddDarcyPressureTerm(
        int nq,
        NekDouble kinvis,
        Array<OneD, Array<OneD, const NekDouble> > &Vel,
        Array<OneD, Array<OneD, NekDouble> > &Q)

    {
        int nVelfields = m_fields.num_elements()-1;

        switch(nVelfields)
        {
            case 2:
            {
                Vmath::Svtvp(nq,-kinvis*m_perm_inv[0],Vel[0],1,Q[0],1,Q[0],1);
                Vmath::Svtvp(nq,-kinvis*m_perm_inv[2],Vel[0],1,Q[1],1,Q[1],1);

                Vmath::Svtvp(nq,-kinvis*m_perm_inv[2],Vel[1],1,Q[0],1,Q[0],1);
                Vmath::Svtvp(nq,-kinvis*m_perm_inv[1],Vel[1],1,Q[1],1,Q[1],1);

                //m_forces[0]->BwdTrans(m_forces[0]->GetCoeffs(),m_forces[0]->UpdatePhys());
                //m_forces[1]->BwdTrans(m_forces[1]->GetCoeffs(),m_forces[1]->UpdatePhys());

                //Vmath::Vadd(nq,Q[0],1,(m_forces[0]->GetPhys()),1,Q[0],1);
                //Vmath::Vadd(nq,Q[1],1,(m_forces[1]->GetPhys()),1,Q[1],1);
            }
            break;
            case 3:
            {
                Vmath::Svtvp(nq,-kinvis*m_perm_inv[0],Vel[0],1,Q[0],1,Q[0],1);
                Vmath::Svtvp(nq,-kinvis*m_perm_inv[3],Vel[0],1,Q[1],1,Q[1],1);
                Vmath::Svtvp(nq,-kinvis*m_perm_inv[4],Vel[0],1,Q[2],1,Q[2],1);
                
                Vmath::Svtvp(nq,-kinvis*m_perm_inv[3],Vel[1],1,Q[0],1,Q[0],1);
                Vmath::Svtvp(nq,-kinvis*m_perm_inv[1],Vel[1],1,Q[1],1,Q[1],1);
                Vmath::Svtvp(nq,-kinvis*m_perm_inv[5],Vel[1],1,Q[2],1,Q[2],1);

                Vmath::Svtvp(nq,-kinvis*m_perm_inv[4],Vel[2],1,Q[0],1,Q[0],1);
                Vmath::Svtvp(nq,-kinvis*m_perm_inv[5],Vel[2],1,Q[1],1,Q[1],1);
                Vmath::Svtvp(nq,-kinvis*m_perm_inv[2],Vel[2],1,Q[2],1,Q[2],1);

                //m_forces[0]->BwdTrans(m_forces[0]->GetCoeffs(),m_forces[0]->UpdatePhys());
                //m_forces[1]->BwdTrans(m_forces[1]->GetCoeffs(),m_forces[1]->UpdatePhys());
                //m_forces[2]->BwdTrans(m_forces[2]->GetCoeffs(),m_forces[2]->UpdatePhys());

                //Vmath::Vadd(nq,Q[0],1,(m_forces[0]->GetPhys()),1,Q[0],1);
                //Vmath::Vadd(nq,Q[1],1,(m_forces[1]->GetPhys()),1,Q[1],1);
                //Vmath::Vadd(nq,Q[2],1,(m_forces[2]->GetPhys()),1,Q[2],1);
            }
            break;
            default:
                ASSERTL0(0,"Dimension not supported");
                break;
        }
    }


    /**
     * Registers the class with the Factory.
     */
    std::string DarcyTermImplicitAnisotropic::className = GetDarcyTermFactory().RegisterCreatorFunction(
        "ImplicitAnisotropic",
        DarcyTermImplicitAnisotropic::create,
        "Implicit Term Anisotropic");

    DarcyTermImplicitAnisotropic::DarcyTermImplicitAnisotropic(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields)
        : DarcyTerm(pSession,pFields)
    {
    }

    /** 
     * 
     */
    DarcyTermImplicitAnisotropic::~DarcyTermImplicitAnisotropic()
    {
    }

    /** 
     * 
     */
    void DarcyTermImplicitAnisotropic::v_SetupPermeability()
    {
        int nDim = m_fields.num_elements()-1;
        switch(nDim)
        {
            case 2:
            {
                m_perm = Array<OneD, NekDouble> (3);
                m_perm_inv = Array<OneD, NekDouble> (3);

                std::string varCoeffs[3] = {
                    "kxx",
                    "kyy",
                    "kxy"
                };
                for (int i = 0; i < 3; ++i)
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
                m_perm = Array<OneD, NekDouble> (6);
                m_perm_inv = Array<OneD, NekDouble> (6);

                std::string varCoeffs[6] = {
                    "kxx",
                    "kyy",
                    "kzz",
                    "kxy",
                    "kxz",
                    "kyz"
                };
                
                for (int i = 0; i < 6; ++i)
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
    void DarcyTermImplicitAnisotropic::v_EvaluateDarcyTerm(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray, 
        Array<OneD, Array<OneD, NekDouble> > &outarray,
        NekDouble kinvis)
    {
    }

    /**
     * Additional contribution to high order pressure boundary condition
     * due to darcy term
     */
    void DarcyTermImplicitAnisotropic::v_AddDarcyPressureTerm(
        int nq,
        NekDouble kinvis,
        Array<OneD, Array<OneD, const NekDouble> > &Vel,
        Array<OneD, Array<OneD, NekDouble> > &Q)

    {
        int nVelfields = m_fields.num_elements()-1;

        switch(nVelfields)
        {
            case 2:
            {
                Vmath::Svtvp(nq,-kinvis*m_perm_inv[0],Vel[0],1,Q[0],1,Q[0],1);
                Vmath::Svtvp(nq,-kinvis*m_perm_inv[2],Vel[0],1,Q[1],1,Q[1],1);

                Vmath::Svtvp(nq,-kinvis*m_perm_inv[2],Vel[1],1,Q[0],1,Q[0],1);
                Vmath::Svtvp(nq,-kinvis*m_perm_inv[1],Vel[1],1,Q[1],1,Q[1],1);

                //m_forces[0]->BwdTrans(m_forces[0]->GetCoeffs(),m_forces[0]->UpdatePhys());
                //m_forces[1]->BwdTrans(m_forces[1]->GetCoeffs(),m_forces[1]->UpdatePhys());

                //Vmath::Vadd(nq,Q[0],1,(m_forces[0]->GetPhys()),1,Q[0],1);
                //Vmath::Vadd(nq,Q[1],1,(m_forces[1]->GetPhys()),1,Q[1],1);
            }
            break;
            case 3:
            {
                Vmath::Svtvp(nq,-kinvis*m_perm_inv[0],Vel[0],1,Q[0],1,Q[0],1);
                Vmath::Svtvp(nq,-kinvis*m_perm_inv[3],Vel[0],1,Q[1],1,Q[1],1);
                Vmath::Svtvp(nq,-kinvis*m_perm_inv[4],Vel[0],1,Q[2],1,Q[2],1);
                
                Vmath::Svtvp(nq,-kinvis*m_perm_inv[3],Vel[1],1,Q[0],1,Q[0],1);
                Vmath::Svtvp(nq,-kinvis*m_perm_inv[1],Vel[1],1,Q[1],1,Q[1],1);
                Vmath::Svtvp(nq,-kinvis*m_perm_inv[5],Vel[1],1,Q[2],1,Q[2],1);

                Vmath::Svtvp(nq,-kinvis*m_perm_inv[4],Vel[2],1,Q[0],1,Q[0],1);
                Vmath::Svtvp(nq,-kinvis*m_perm_inv[5],Vel[2],1,Q[1],1,Q[1],1);
                Vmath::Svtvp(nq,-kinvis*m_perm_inv[2],Vel[2],1,Q[2],1,Q[2],1);

                //m_forces[0]->BwdTrans(m_forces[0]->GetCoeffs(),m_forces[0]->UpdatePhys());
                //m_forces[1]->BwdTrans(m_forces[1]->GetCoeffs(),m_forces[1]->UpdatePhys());
                //m_forces[2]->BwdTrans(m_forces[2]->GetCoeffs(),m_forces[2]->UpdatePhys());

                //Vmath::Vadd(nq,Q[0],1,(m_forces[0]->GetPhys()),1,Q[0],1);
                //Vmath::Vadd(nq,Q[1],1,(m_forces[1]->GetPhys()),1,Q[1],1);
                //Vmath::Vadd(nq,Q[2],1,(m_forces[2]->GetPhys()),1,Q[2],1);
            }
            break;
            default:
                ASSERTL0(0,"Dimension not supported");
                break;
        }
    }

    void DarcyTermImplicitAnisotropic::v_GetImplicitDarcyFactor(
        Array<OneD, NekDouble> &permCoeff)
    {
        int nDim = m_fields.num_elements()-1;
        for (int i=0; i<nDim; ++i)
        {
            permCoeff[i]=m_perm_inv[i];
        }
    }

}

