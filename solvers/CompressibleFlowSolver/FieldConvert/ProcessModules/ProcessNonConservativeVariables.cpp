////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessNonConservativeVariables.cpp
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  License for the specific language governing rights and limitations under
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description: Computes non-conservative variables
//
////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <iostream>
using namespace std;

#include "ProcessNonConservativeVariables.h"

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/ParseUtils.hpp>

namespace Nektar
{
namespace FieldUtils
{

ModuleKey ProcessNonConservativeVariables::className =
    GetModuleFactory().RegisterCreatorFunction(
                                               ModuleKey(eProcessModule, "nonconservativevariables"),
                                               ProcessNonConservativeVariables::create, "Redefine conservative variables into non-conservative variables");

    ProcessNonConservativeVariables::ProcessNonConservativeVariables(FieldSharedPtr f) : ProcessModule(f)
    {
    }

    ProcessNonConservativeVariables::~ProcessNonConservativeVariables()
    {
    }

    void ProcessNonConservativeVariables::Process(po::variables_map &vm)
    {
        if (m_f->m_verbose)
        {
            cout << "ProcessNonConservativeVariables: Process non-conservative variables..." << endl;
        }

        int i, j, s;
        int expdim    = m_f->m_graph->GetMeshDimension();
        int spacedim  = expdim;
        if ((m_f->m_fielddef[0]->m_numHomogeneousDir) == 1 ||
            (m_f->m_fielddef[0]->m_numHomogeneousDir) == 2)
        {
            spacedim = 3;
        }
        int nfields = m_f->m_fielddef[0]->m_fields.size();

        int npoints = m_f->m_exp[0]->GetNpoints();
        Array<OneD, Array<OneD, NekDouble> > outfield(nfields);

        int nstrips;
        m_f->m_session->LoadParameter("Strip_Z",nstrips,1);
        ASSERTL0(nstrips == 1,"Routine not set up for strips");

        NekDouble gamma;    
        m_f->m_session->LoadParameter("Gamma", gamma, 1.4);
    
        NekDouble gammaMinusOne    = gamma - 1.0;
    
        m_f->m_exp.resize(nfields*nstrips);
    
        for (i = 0; i < nfields; ++i)
        {
            outfield[i] = Array<OneD, NekDouble>(npoints);
        }
    
        Array<OneD, NekDouble> tmp(npoints, 0.0);
    
        // Keep rho
        Vmath::Vcopy(npoints,m_f->m_exp[0]->GetPhys(),1,outfield[0],1);
    
        // Calculate velocity
        for (i = 1; i < nfields-1; ++i)
        {
            Vmath::Vdiv(npoints,
                        m_f->m_exp[i]->GetPhys(), 1,
                        m_f->m_exp[0]->GetPhys(), 1,
                        outfield[i], 1);
        }
        
        //Calculate pressure
        for (i = 0; i < spacedim; i++)
        {
            Vmath::Vmul(npoints,
                        m_f->m_exp[i + 1]->GetPhys(), 1,
                        m_f->m_exp[i + 1]->GetPhys(), 1,
                        tmp,1);
        
        
            Vmath::Smul(npoints, 0.5,
                        tmp, 1,tmp, 1);
            
            Vmath::Vadd(npoints,
                        outfield[nfields-1], 1,
                        tmp, 1,
                        outfield[nfields-1], 1);
        }
        
        Vmath::Vdiv(npoints,
                    outfield[nfields-1], 1,
                    m_f->m_exp[0]->GetPhys(), 1,
                    outfield[nfields-1],1);
    
        Vmath::Vsub(npoints,
                    m_f->m_exp[spacedim + 1]->GetPhys(), 1,
                    outfield[nfields-1], 1,
                    outfield[nfields-1], 1);
        
        Vmath::Smul(npoints, gammaMinusOne,
                    outfield[nfields-1], 1,
                    outfield[nfields-1], 1);
        
        for (i = 0; i < nfields; ++i)
        {
            m_f->m_exp[i]->SetPhys(outfield[i]);
            m_f->m_exp[i]->FwdTrans_IterPerExp(outfield[i],
                                               m_f->m_exp[i]->UpdateCoeffs());
        }
    
        vector<string > outname;
        outname.push_back("rho");
        if (spacedim == 1)
        {
            outname.push_back("u");
            outname.push_back("p");
        }

        if (spacedim == 2)
        {
            outname.push_back("u");
            outname.push_back("v");
            outname.push_back("p");
        }
    
        if (spacedim == 3)
        {
            outname.push_back("u");
            outname.push_back("v");
            outname.push_back("w");
            outname.push_back("p");

        }
    
        std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef
            = m_f->m_exp[0]->GetFieldDefinitions();
        std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());
    
    
        for (j = 0; j <  nfields; ++j)
        {
            for (i = 0; i < FieldDef.size(); ++i)
            {
                FieldDef[i]->m_fields.push_back(outname[j]);
                m_f->m_exp[j]->AppendFieldData(FieldDef[i], FieldData[i]);
            }
        }
        m_f->m_fielddef = FieldDef;
        m_f->m_data     = FieldData;
    }
}
    
}
