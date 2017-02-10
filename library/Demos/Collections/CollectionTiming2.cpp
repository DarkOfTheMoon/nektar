///////////////////////////////////////////////////////////////////////////////
//
// File: CollectionTiming.cpp
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
// Description: Small demo to run timings of various operators.
//
///////////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <cstdlib>
#include <iomanip>

#include <boost/timer/timer.hpp>

#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/Communication/Comm.h>
#include <LibUtilities/BasicUtils/Timer.h>
#include <MultiRegions/ExpList3D.h>
#include <Collections/Collection.h>
#include <SpatialDomains/MeshGraph.h>

using namespace std;
using namespace Nektar;

using boost::timer::cpu_timer;
using boost::timer::cpu_times;
using boost::timer::nanosecond_type;
using boost::timer::format;

MultiRegions::ExpListSharedPtr SetupExpList(
    int                                  N,
    LibUtilities::SessionReaderSharedPtr session,
    SpatialDomains::MeshGraphSharedPtr   graph,
    Collections::ImplementationType      impType)
{
    graph->SetExpansionsToPolyOrder(N);

    MultiRegions::ExpListSharedPtr expList =
        MemoryManager<MultiRegions::ExpList3D>::AllocateSharedPtr(
            session, graph);

    expList->CreateCollections(impType);

    return expList;
}


int main(int argc, char *argv[])
{
    LibUtilities::SessionReader::RegisterCmdLineFlag(
        "data", "d", "Print in data format");
    LibUtilities::SessionReaderSharedPtr session
        = LibUtilities::SessionReader::CreateInstance(argc, argv);
    LibUtilities::CommSharedPtr vComm = vSession->GetComm();

    bool fmt = session->DefinesCmdLineArgument("data");

    MultiRegions::ExpListSharedPtr expList;

    cpu_timer timer;

    int Ntest, order;
    session->LoadParameter("Ntest", Ntest, 1000);
    session->LoadParameter("order", order, 7);

    string sl = fmt ? "# " : "";

    // Read in mesh
    SpatialDomains::MeshGraphSharedPtr graph =
        SpatialDomains::MeshGraph::Read(session);

    if (vComm->GetRank() == 0)
    {
        cout << "Testing BwdTrans " << Ntest << " " << order << endl;
    }

    // BwdTrans operator
    Collections::ImplementationType impType = Collections::eSumFac;

    expList = SetupExpList(order, session, graph, impType);
    Array<OneD, NekDouble> input (expList->GetNcoeffs());
    Array<OneD, NekDouble> output(expList->GetNpoints());

    Timer t;
    vComm->Block();
    t.Start();
    for (int i = 0; i < Ntest; ++i)
    {
        expList->BwdTrans(input, output);
    }
    vComm->Block();
    t.Stop();

    int nElmt = expList->GetExpSize();
    int nM = order + 1;
    int nQ = order + 2;

    // flops: 3 matrix-matrix multiplications
    long long flop = (nElmt*(nM*nM*nM*nQ + nQ*nQ*nM*nM + nQ*nQ*nQ*nM));
    NekDouble gflop = flop / 1024.0 / 1024.0 / 1024.0;
    NekDouble elapsed = t.TimePerTest(1);

    vComm->AllReduce(gflop, LibUtilities::ReduceSum);

    if (vComm->GetRank() == 0)
    {
        cout << "Time: " << elapsed << " GFLOP/s: " << gflop << endl;
    }
}
