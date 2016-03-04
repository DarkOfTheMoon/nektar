////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessJacobianEnergy.cpp
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
//  Description: Compute energy of Jacobian.
//
////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <iostream>
using namespace std;

#include "ProcessGrabData.h"

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/ParseUtils.hpp>

namespace Nektar
{
namespace Utilities
{

ModuleKey ProcessGrabData::className =
        GetModuleFactory().RegisterCreatorFunction(
                ModuleKey(eProcessModule, "ProcessGrabData"),
                ProcessGrabData::create,
                "data.");

ProcessGrabData::ProcessGrabData(FieldSharedPtr f) :
    ProcessModule(f)
{
 m_config["totxt"]            = ConfigOption(false, "NotSet",
                                    "Txt file to write to");
}

ProcessGrabData::~ProcessGrabData()
{
}

void ProcessGrabData::Process(po::variables_map &vm)
{
    if (m_f->m_verbose)
    {
        cout << "starting file" << endl;
    }

    vector<Array<OneD, NekDouble> > fields;

    for(int i = 0; i < m_f->m_exp.size(); i++)
    {
        fields.push_back(m_f->m_exp[i]->UpdatePhys());
    }

    Array<OneD, NekDouble> x(m_f->m_exp[0]->GetTotPoints()),
                           y(m_f->m_exp[0]->GetTotPoints()),
                           z(m_f->m_exp[0]->GetTotPoints());
    m_f->m_exp[0]->GetCoords(x,y,z);

    ofstream file;
    string totxt = m_config["totxt"].as<string>();
    file.open(totxt.c_str());

    for(int i = 0; i < fields[0].num_elements(); i++)
    {
        file << x[i] << " " << y[i] << " " << z[i] << " ";
        for(int j = 0; j < m_f->m_exp.size(); j++)
        {
            file << fields[j][i] << " ";
        }
        file << endl;
    }

    file.close();

    exit(-1);

}

}
}
