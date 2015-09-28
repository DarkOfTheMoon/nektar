////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessVorticity.cpp
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
//  Description: Computes vorticity field.
//
////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <iostream>
using namespace std;

#include "ProcessChangeVariables.h"

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/ParseUtils.hpp>

namespace Nektar
{
namespace Utilities
{

ModuleKey ProcessChangeVariables::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "changevariables"),
        ProcessChangeVariables::create, "Redefine CFS variables");

ProcessChangeVariables::ProcessChangeVariables(FieldSharedPtr f) : ProcessModule(f)
{
}

ProcessChangeVariables::~ProcessChangeVariables()
{
}

void ProcessChangeVariables::Process(po::variables_map &vm)
{
    if (m_f->m_verbose)
    {
        cout << "ProcessChangeVariables: Reddefine CFS variables..." << endl;
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
    int newfields = spacedim + 1;

    int npoints = m_f->m_exp[0]->GetNpoints();
    Array<OneD, Array<OneD, NekDouble> > grad(nfields*nfields);
    Array<OneD, Array<OneD, NekDouble> > outfield(newfields);

    int nstrips;
    NekDouble gamma;
    
    m_f->m_session->LoadParameter("Strip_Z",nstrips,1);
    m_f->m_session->LoadParameter("Gamma", gamma, 1.4);
    
    NekDouble gammaMinusOne    = gamma - 1.0;
    
    m_f->m_exp.resize(nfields*nstrips);
    
    for (i = 0; i < newfields; ++i)
    {
        outfield[i] = Array<OneD, NekDouble>(npoints);
    }
    
    vector<MultiRegions::ExpListSharedPtr> Exp(nstrips*newfields);
    Array<OneD, NekDouble> tmp(npoints, 0.0);
    
    for(s = 0; s < nstrips; ++s) //homogeneous strip varient
    {
        // Calculate velocity
        for (i = 0; i < newfields-1; ++i)
        {
            Vmath::Vdiv(npoints,
                        m_f->m_exp[s*nfields+i+1]->GetPhys(), 1,
                        m_f->m_exp[s*nfields]->GetPhys(), 1,
                        outfield[i], 1);
        }
        
        //Calculate pressure
        for (i = 0; i < spacedim; i++)
        {
            Vmath::Vmul(npoints,
                        m_f->m_exp[s*nfields + i + 1]->GetPhys(), 1,
                        m_f->m_exp[s*nfields + i + 1]->GetPhys(), 1,
                        tmp,1);
        
        
            Vmath::Smul(npoints, 0.5,
                        tmp, 1,
                        tmp, 1);
            
            Vmath::Vadd(npoints,
                        outfield[newfields-1], 1,
                        tmp, 1,
                        outfield[newfields-1], 1);
        }
        
        Vmath::Vdiv(npoints,
                    outfield[newfields-1], 1,
                    m_f->m_exp[s*nfields]->GetPhys(), 1,
                    outfield[newfields-1],1);
        
        Vmath::Vsub(npoints,
                    m_f->m_exp[s*nfields + spacedim + 1]->GetPhys(), 1,
                    outfield[newfields-1], 1,
                    outfield[newfields-1],1);
        
        Vmath::Smul(npoints, gammaMinusOne,
                    outfield[newfields-1], 1,
                    outfield[newfields-1], 1);
        
        for (i = 0; i < newfields; ++i)
        {
            int n = s*newfields + i;
            Exp[n] = m_f->AppendExpList(m_f->m_fielddef[0]->m_numHomogeneousDir);
            Exp[n]->UpdatePhys() = outfield[i];
            Exp[n]->FwdTrans_IterPerExp(outfield[i],
                                        Exp[n]->UpdateCoeffs());
        }
    }
    
    vector<MultiRegions::ExpListSharedPtr>::iterator it;
    for(s = 0; s < nstrips; ++s)
    {
        for(i = 0; i < newfields; ++i)
        {
            it = m_f->m_exp.begin()+s*(nfields+newfields)+nfields+i;
            m_f->m_exp.insert(it, Exp[s*newfields+i]);
        }
    }
    
    vector<string > outname;
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
    
    
    for(s = 0; s < nstrips; ++s) //homogeneous strip varient
    {
        for (j = 0; j <  newfields; ++j)
        {
            for (i = 0; i < FieldDef.size()/nstrips; ++i)
            {
                int n = s * FieldDef.size()/nstrips + i;
                FieldDef[n]->m_fields.push_back(outname[j]);
                m_f->m_exp[s*(nfields + newfields)+ nfields + j]->AppendFieldData(FieldDef[n], FieldData[n]);
            }
        }
    }

    m_f->m_fielddef = FieldDef;
    m_f->m_data     = FieldData;
}
    
}
}
