///////////////////////////////////////////////////////////////////////////////
//
// File AdjointAdvectionWeakDG.cpp
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
// Description: Evaluation of the adjoint advective term
//
///////////////////////////////////////////////////////////////////////////////

#include <CompressibleFlowSolver/AdvectionTerms/AdjointAdvectionWeakDG.h>
#include <StdRegions/StdSegExp.h>
#include <MultiRegions/ContField1D.h>
#include <MultiRegions/ContField2D.h>
#include <MultiRegions/ContField3D.h>
#include <MultiRegions/DisContField1D.h>
#include <MultiRegions/DisContField2D.h>
#include <MultiRegions/ContField3DHomogeneous1D.h>
#include <MultiRegions/DisContField3DHomogeneous1D.h>
#include <MultiRegions/ContField3DHomogeneous2D.h>
#include <MultiRegions/DisContField3DHomogeneous2D.h>

namespace Nektar
{

string AdjointAdvectionWeakDG::className = SolverUtils
        ::GetAdvectionFactory().RegisterCreatorFunction("AdjointWeakDG",
                                        AdjointAdvectionWeakDG::create);
    

/**
 *
 */
AdjointAdvectionWeakDG::AdjointAdvectionWeakDG()
    : Advection()
{
}


void AdjointAdvectionWeakDG::v_InitObject(
       LibUtilities::SessionReaderSharedPtr        pSession,
       Array<OneD, MultiRegions::ExpListSharedPtr> pFields)
{

    Advection::v_InitObject(pSession, pFields);

    m_session            = pSession;
    m_boundaryConditions = MemoryManager<SpatialDomains::BoundaryConditions>
                    ::AllocateSharedPtr(m_session, pFields[0]->GetGraph());
    m_spacedim           = pFields[0]->GetGraph()->GetSpaceDimension();
    m_expdim             = pFields[0]->GetGraph()->GetMeshDimension();
    m_CoeffState         = MultiRegions::eLocal;

    //Setting parameters for homogeneous problems
    m_HomoDirec          = 0;
    m_useFFT             = false;
    m_HomogeneousType    = eNotHomogeneous;
    m_SingleMode         = false;
    m_HalfMode           = false;
    m_MultipleModes      = false;
    m_homogen_dealiasing = false;
    if(m_session->DefinesSolverInfo("HOMOGENEOUS"))
    {
        std::string HomoStr = m_session->GetSolverInfo("HOMOGENEOUS");
        m_spacedim          = 3;

        if((HomoStr == "HOMOGENEOUS1D")||(HomoStr == "Homogeneous1D")||
           (HomoStr == "1D")||(HomoStr == "Homo1D"))
        {
            m_HomogeneousType = eHomogeneous1D;
            m_LhomZ           = m_session->GetParameter("LZ");
            m_HomoDirec       = 1;

            m_homogen_dealiasing = pSession->DefinesSolverInfo("DEALIASING");

            if(m_session->DefinesSolverInfo("ModeType"))
            {
                m_session->MatchSolverInfo("ModeType","SingleMode",m_SingleMode,false);
                m_session->MatchSolverInfo("ModeType","HalfMode",m_HalfMode,false);
                m_session->MatchSolverInfo("ModeType","MultipleModes",m_MultipleModes,false);
            }
            if(m_session->DefinesSolverInfo("ModeType"))
            {
                if(m_SingleMode)
                {
                    m_npointsZ=2;
                }
                else if(m_HalfMode)
                {
                    m_npointsZ=1;
                }
                else if(m_MultipleModes)
                {
                    m_npointsZ        = m_session->GetParameter("HomModesZ");
                }
                else
                {
                    ASSERTL0(false, "SolverInfo ModeType not valid");
                }
            }
            else
            {
                m_session->LoadParameter("HomModesZ",m_npointsZ);

            }
        }

        if((HomoStr == "HOMOGENEOUS2D")||(HomoStr == "Homogeneous2D")||
           (HomoStr == "2D")||(HomoStr == "Homo2D"))
        {
            m_HomogeneousType = eHomogeneous2D;
            m_session->LoadParameter("HomModesY", m_npointsY);
            m_session->LoadParameter("LY",        m_LhomY);
            m_session->LoadParameter("HomModesZ", m_npointsZ);
            m_session->LoadParameter("LZ",        m_LhomZ);
            m_HomoDirec       = 2;
        }

        if((HomoStr == "HOMOGENEOUS3D")||(HomoStr == "Homogeneous3D")||
           (HomoStr == "3D")||(HomoStr == "Homo3D"))
        {
            m_HomogeneousType = eHomogeneous3D;
            m_session->LoadParameter("HomModesX",m_npointsX);
            m_session->LoadParameter("LX",       m_LhomX   );
            m_session->LoadParameter("HomModesY",m_npointsY);
            m_session->LoadParameter("LY",       m_LhomY   );
            m_session->LoadParameter("HomModesZ",m_npointsZ);
            m_session->LoadParameter("LZ",       m_LhomZ   );
            
            // Give the final chk number.
            m_session->LoadParameter("chkfinal",   m_chkf, 0.0);
            m_HomoDirec       = 3;
        }

        if(m_session->DefinesSolverInfo("USEFFT"))
        {
            m_useFFT = true;
        }
    }
    else
    {
        m_npointsZ = 1; // set to default value so can use to identify 2d or 3D (homogeneous) expansions
    }

    if(m_session->DefinesSolverInfo("PROJECTION"))
    {
        std::string ProjectStr
        = m_session->GetSolverInfo("PROJECTION");

        if((ProjectStr == "Continuous")||(ProjectStr == "Galerkin")||
           (ProjectStr == "CONTINUOUS")||(ProjectStr == "GALERKIN"))
        {
            m_projectionType = MultiRegions::eGalerkin;
        }
        else if(ProjectStr == "DisContinuous")
        {
            m_projectionType = MultiRegions::eDiscontinuous;
        }
        else
        {
            ASSERTL0(false,"PROJECTION value not recognised");
        }
    }
    else
    {
        cerr << "Projection type not specified in SOLVERINFO,"
        "defaulting to continuous Galerkin" << endl;
        m_projectionType = MultiRegions::eGalerkin;
    }

    
    
    SpatialDomains::MeshGraphSharedPtr graphShrPtr =
    SpatialDomains::MeshGraph::Read(m_session);
    
    ASSERTL0(m_session->DefinesFunction("BaseFlow"),
             "Base flow must be defined for linearised forms.");
    string basefile = m_session->GetFunctionFilename("BaseFlow", 0);
    int nvar = m_session->GetVariables().size();
    m_baseflow = Array<OneD, Array<OneD, NekDouble> >(nvar);
    for (int i = 0; i < nvar; ++i)
    {
        m_baseflow[i] = Array<OneD, NekDouble>(pFields[i]->GetTotPoints(), 0.0);
    }
    

    m_slices=1;
        
        //BaseFlow from file
    if (m_session->GetFunctionType("BaseFlow", m_session->GetVariable(0))
        == LibUtilities::eFunctionTypeFile)
    {
        //ImportFldBase(basefile,pFields);
        
    }
    //analytic base flow
    else
    {
        int nq = pFields[0]->GetNpoints();
        Array<OneD,NekDouble> x0(nq);
        Array<OneD,NekDouble> x1(nq);
        Array<OneD,NekDouble> x2(nq);
        
        // get the coordinates (assuming all fields have the same
        // discretisation)
        pFields[0]->GetCoords(x0,x1,x2);
        for(unsigned int i = 0 ; i < m_baseflow.num_elements(); i++)
        {
            LibUtilities::EquationSharedPtr ifunc
            = m_session->GetFunction("BaseFlow", i);
            
            ifunc->Evaluate(x0,x1,x2,m_baseflow[i]);
        }
        
    }

}


AdjointAdvectionWeakDG::~AdjointAdvectionWeakDG()
{
}


void AdjointAdvectionWeakDG::v_Advect(
    const int nConvectiveFields,
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    const Array<OneD, Array<OneD, NekDouble> >        &advVel,
    const Array<OneD, Array<OneD, NekDouble> >        &inarray,
    Array<OneD, Array<OneD, NekDouble> >              &outarray,
    const NekDouble                                   &time)
{
    int nDim            = fields[0]->GetCoordim(0);
    int nPointsTot      = fields[0]->GetTotPoints();
    int nCoeffs         = fields[0]->GetNcoeffs();
    int nTracePointsTot = fields[0]->GetTrace()->GetTotPoints();
    int i, j;
   
    Array<OneD, Array<OneD, NekDouble> >        outarray_tmp(nConvectiveFields);
    Array<OneD, Array<OneD, NekDouble> >        outarray_tmp2(nConvectiveFields);
    
    Array<OneD, Array<OneD, NekDouble> > tmp(nConvectiveFields);
    Array<OneD, Array<OneD, Array<OneD, NekDouble> > > FluxVector(
                                                    nConvectiveFields);
    Array<OneD, Array<OneD, Array<OneD, NekDouble> > > JacDivVector(
                                                    nConvectiveFields);
        // Allocate storage for flux vector F(u).
    
    ASSERTL1(m_riemann,
             "Riemann solver must be provided for AdvectionWeakDG.");
    
    for (i = 0; i < nConvectiveFields; ++i)
    {
        FluxVector[i] =
        Array<OneD, Array<OneD, NekDouble> >(m_spaceDim);
        
        //JacDivVector[i] =
        //Array<OneD, Array<OneD, NekDouble> >(m_spaceDim);
        
        for (j = 0; j < m_spaceDim; ++j)
        {
            FluxVector[i][j] =
            Array<OneD, NekDouble>(nPointsTot, 0.0);
            //JacDivVector[i][j] =
            //Array<OneD, NekDouble>(nPointsTot, 0.0);
        }
    }

    m_AdjointFluxVector(inarray, FluxVector);
    
    //m_JacTransposeDivVector(inarray, JacDivVector);

    // Get the advection part (without numerical flux)
    for(i = 0; i < nConvectiveFields; ++i)
    {
        tmp[i] = Array<OneD, NekDouble>(nCoeffs, 0.0);
        
        outarray_tmp[i] = Array<OneD, NekDouble>(nCoeffs, 0.0);
        outarray_tmp2[i] = Array<OneD, NekDouble>(nCoeffs, 0.0);
        
        for (j = 0; j < nDim; ++j)
        {
            fields[i]->IProductWRTDerivBase(j,
                                            FluxVector[i][j],
                                            outarray[i]);
            
            // since the adjoint for the convective part is NOT given in
            // conservative form, an additional term where the derivatives
            // of the transposed jacobians are multiplied with the basis hence,
            
            //fields[i]->IProductWRTBase(JacDivVector[i][j],
            //                           outarray_tmp[i]);
            
            /*Vmath::Vadd(nCoeffs,
                        outarray[i], 1,
                        outarray_tmp[i], 1,
                        outarray_tmp2[i], 1);*/
            
            Vmath::Vadd(nCoeffs,
                        outarray[i], 1,
                        tmp[i], 1,
                        tmp[i], 1);
        }
    }
    
    // Store forwards/backwards space along trace space
    Array<OneD, Array<OneD, NekDouble> > Fwd    (nConvectiveFields);
    Array<OneD, Array<OneD, NekDouble> > Bwd    (nConvectiveFields);
    Array<OneD, Array<OneD, NekDouble> > numflux(nConvectiveFields);
    Array<OneD, Array<OneD, NekDouble> > numflux2(nConvectiveFields);
    
    for(i = 0; i < nConvectiveFields; ++i)
    {
        Fwd[i]      = Array<OneD, NekDouble>(nTracePointsTot, 0.0);
        Bwd[i]      = Array<OneD, NekDouble>(nTracePointsTot, 0.0);
        numflux[i]  = Array<OneD, NekDouble>(nTracePointsTot, 0.0);
        numflux2[i] = Array<OneD, NekDouble>(nTracePointsTot, 0.0);
        fields[i]->GetFwdBwdTracePhys(inarray[i], Fwd[i], Bwd[i]);
    }
    
    m_riemann->Solve(m_spaceDim, Fwd, Bwd, numflux);
    
    // Evaulate <\phi, \hat{F}\cdot n> - OutField[i]
    for(i = 0; i < nConvectiveFields; ++i)
    {
        Vmath::Neg                      (nCoeffs, tmp[i], 1);
        fields[i]->AddTraceIntegral     (numflux[i], tmp[i]);
        fields[i]->MultiplyByElmtInvMass(tmp[i], tmp[i]);
        fields[i]->BwdTrans             (tmp[i], outarray[i]);
    }
}

    
void AdjointAdvectionWeakDG::v_SetBaseFlow(
                const Array<OneD, Array<OneD, NekDouble> >    &inarray)
{
    int nvar = inarray.num_elements();
    int npts = inarray[0].num_elements();
    m_baseflow = Array<OneD, Array<OneD, NekDouble> >(nvar);
    for (int i = 0; i < nvar; ++i)
    {
        m_baseflow[i] = Array<OneD, NekDouble>(npts, 0.0);
        Vmath::Vcopy(npts, inarray[i], 1, m_baseflow[i], 1);
    }
}
    
void AdjointAdvectionWeakDG::v_GetBaseFlow(
       Array<OneD, Array<OneD, NekDouble> >    &baseflow)
{
    for (int i = 0; i < baseflow.num_elements(); ++i)
    {
        baseflow[i] = m_baseflow[i];
    }
}

void AdjointAdvectionWeakDG::v_ImportFldBase(std::string pInfile,
                Array<OneD, MultiRegions::ExpListSharedPtr>& pFields)
{
    std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef;
    std::vector<std::vector<NekDouble> > FieldData;
        
    int nqtot = pFields[0]->GetTotPoints();
    Array<OneD, NekDouble> tmp_coeff(pFields[0]->GetNcoeffs(), 0.0);
        
        //Get Homogeneous
    LibUtilities::FieldIOSharedPtr fld =
    MemoryManager<LibUtilities::FieldIO>::AllocateSharedPtr(
                                        m_session->GetComm());
    fld->Import(pInfile, FieldDef, FieldData);
        
    int nvar = m_session->GetVariables().size();
    int s;
        
    if(m_session->DefinesSolverInfo("HOMOGENEOUS"))
    {
        std::string HomoStr = m_session->GetSolverInfo("HOMOGENEOUS");
    }
        
        // copy FieldData into m_fields
    for(int j = 0; j < nvar; ++j)
    {
        for(int i = 0; i < FieldDef.size(); ++i)
        {
            if ((m_session->DefinesSolverInfo("HOMOGENEOUS") &&
                    (m_session->GetSolverInfo("HOMOGENEOUS")=="HOMOGENEOUS1D" ||
                    m_session->GetSolverInfo("HOMOGENEOUS")=="1D" ||
                    m_session->GetSolverInfo("HOMOGENEOUS")=="Homo1D")) &&
                nvar==4)
            {
                    // w-component must be ignored and set to zero.
                if (j != nvar - 2)
                {
                        // p component (it is 4th variable of the 3D and corresponds 3nd variable of 2D)
                    s = (j == nvar - 1) ? 2 : j;
                        
                        //extraction of the 2D
                    pFields[j]->ExtractDataToCoeffs(
                                            FieldDef[i],
                                            FieldData[i],
                                            FieldDef[i]->m_fields[s],
                                            tmp_coeff);
                }
                    
                // Put zero on higher modes
                int ncplane = (pFields[0]->GetNcoeffs()) / m_npointsZ;
                    
                if (m_npointsZ > 2)
                {
                        Vmath::Zero(ncplane*(m_npointsZ-2),
                                    &tmp_coeff[2*ncplane],1);
                }
            }
            //2D cases and Homogeneous1D Base Flows
            else
            {
                bool flag = FieldDef[i]->m_fields[j] ==
                m_session->GetVariable(j);
                    
                ASSERTL0(flag, (std::string("Order of ") + pInfile
                        + std::string(" data and that defined in "
                                    "m_boundaryconditions differs")).c_str());
                    
                pFields[j]->ExtractDataToCoeffs(FieldDef[i], FieldData[i],
                                                    FieldDef[i]->m_fields[j],
                                                    tmp_coeff);
            }
        }
        
        if(m_SingleMode || m_HalfMode)
        {
                pFields[j]->GetPlane(0)->BwdTrans(tmp_coeff, m_baseflow[j]);
                
            if(m_SingleMode)
            {
                    //copy the bwd into the second plane for single Mode Analysis
                    int ncplane=(pFields[0]->GetNpoints())/m_npointsZ;
                    Vmath::Vcopy(ncplane,&m_baseflow[j][0],1,&m_baseflow[j][ncplane],1);
            }
        }
        else
        {
            pFields[j]->BwdTrans(tmp_coeff, m_baseflow[j]);
        
        }
    }
}
/**
 * Import field from infile and load into \a m_fields. This routine will
 * also perform a \a BwdTrans to ensure data is in both the physical and
 * coefficient storage.
 * @param   infile          Filename to read.
 */
} //end of namespace
        
