////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessCyl.cpp
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
//  Description: create mesh from cad using mesh utils
//
////////////////////////////////////////////////////////////////////////////////

#include <NekMeshUtils/MeshElements/MeshElements.h>

#include <LocalRegions/SegExp.h>
#include <LocalRegions/QuadExp.h>
#include <LocalRegions/TriExp.h>
#include <LocalRegions/NodalTriExp.h>

#include <LibUtilities/BasicUtils/SharedArray.hpp>

#include "ProcessCyl.h"

using namespace std;
using namespace Nektar::MeshUtils;

namespace Nektar
{
namespace Utilities
{

ModuleKey ProcessCyl::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "cyl"),
        ProcessCyl::create);

/**
 * @brief Default constructor.
 */
ProcessCyl::ProcessCyl(MeshSharedPtr m) : ProcessModule(m)
{
    m_config["surf"] = ConfigOption(false, "-1",
        "Tag identifying surface to process.");
    m_config["r"] = ConfigOption(false, "0.0",
        "Radius of cylinder.");
    m_config["N"] = ConfigOption(false, "7",
        "Number of points along edge.");
}


/**
 * @brief Destructor.
 */
ProcessCyl::~ProcessCyl()
{
}


void ProcessCyl::Process()
{
    set<pair<int,int> > tmp;
    int                 surfTag         = m_config["surf"].as<int>();
    int                 prismedge[2][3] = {{0,5,4}, {2,6,7}};
    Node                zax(0,0,0,1);
    int                 dim             = m_mesh->m_expDim;

    for (int i = 0; i < m_mesh->m_element[dim].size(); ++i)
    {
        ElementSharedPtr el = m_mesh->m_element[dim][i];
        int nSurf = dim == 3 ? el->GetFaceCount() : el->GetEdgeCount();

        for (int j = 0; j < nSurf; ++j)
        {
            int bl = el->GetBoundaryLink(j);
            if (bl == -1)
            {
                continue;
            }

            ElementSharedPtr bEl  = m_mesh->m_element[dim - 1][bl];
            vector<int>      tags = bEl->GetTagList();

            if (find(tags.begin(), tags.end(), surfTag) ==
                tags.end())
            {
                continue;
            }

            ASSERTL0(j == 1 || j == 3, "rofl");

            // Check all edge interior points.
            for (int k = 0; k < 3; ++k)
            {
                EdgeSharedPtr edge = el->GetEdge(prismedge[(j-1)/2][k]);
                GenerateEdgeNodes(edge);
            }
        }
    }
}


void ProcessCyl::GenerateEdgeNodes(EdgeSharedPtr edge)
{
    NodeSharedPtr n1 = edge->m_n1;
    NodeSharedPtr n2 = edge->m_n2;

    int nq    = m_config["N"].as<int>();
    double r  = m_config["r"].as<double>();
    double t1 = atan2(n1->m_y, n1->m_x);
    double t2 = atan2(n2->m_y, n2->m_x);
    double dt;
    double dz;

    if (t1 < -M_PI/2.0 && t2 > 0.0)
    {
        t1 += 2*M_PI;
    }
    if (t2 < -M_PI/2.0 && t1 > 0.0)
    {
        t2 += 2*M_PI;
    }

    dt = (t2-t1) / (nq-1);
    dz = (n2->m_z - n1->m_z) / (nq-1);

    edge->m_edgeNodes.resize(nq-2);
    Node dx = (*n2-*n1) * (1.0/(nq-1.0));
    for (int i = 1; i < nq-1; ++i)
    {
        edge->m_edgeNodes[i-1] = NodeSharedPtr(
            new Node(0, r*cos(t1 + i*dt),
                        r*sin(t1 + i*dt),
                        n1->m_z + i*dz));
    }
    edge->m_curveType = LibUtilities::ePolyEvenlySpaced;
}

}
}
