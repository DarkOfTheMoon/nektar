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

#include <libxsmm.h>
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

void HexFlops(Collections::ImplementationType  impType,
              int                              order,
              int                              nElmt,
              NekDouble                       &gflop,
              NekDouble                       &matSize)
{
    const int nM = order + 1;
    const int nQ = order + 2;

    if (impType == Collections::eIterPerExp ||
        impType == Collections::eSumFac)
    {
        gflop = 2.0 * nElmt * (nM * nM * nM * nQ + nQ * nQ * nM * nM +
                               nQ * nQ * nQ * nM);
        matSize =
            ((nQ * nM) + (nM * nM * nM) + (nQ * nM * nM)) + // m_funcs[0]
            ((nQ * nM) + (nM * nQ) + (nQ * nQ)) * nM      + // m_funcs[1]
            ((nQ * nQ * nM) + (nM * nQ) + (nQ * nQ * nQ));  // m_funcs[2]
        matSize *= nElmt;
    }
    else if (impType == Collections::eStdMat)
    {
        gflop = 2.0 * nElmt * nM * nM * nM * nQ * nQ * nQ;
        matSize = nQ*nQ*nQ * nM*nM*nM * nElmt;
    }

    gflop *= 1e-9;
    matSize *= sizeof(NekDouble);
    matSize /= 1024.0 * 1024.0 * 1024.0;
}

void TetFlops(Collections::ImplementationType  impType,
              int                              order,
              int                              nElmt,
              NekDouble                       &gflop,
              NekDouble                       &matSize)
{
    const int nM = order + 1;
    const int nQ = order + 2;

    if (impType == Collections::eIterPerExp ||
        impType == Collections::eSumFac)
    {
        gflop = 2.0 * nM * nQ * nQ * nQ;

        for (int i = 0; i < nM; ++i)
        {
            for (int j = 0; j < nM - i; ++j)
            {
                gflop += 2.0 * nQ * (nM - i - j);
            }

            gflop += 2.0 * (nM - i) * nQ * nQ;
        }

        gflop *= nElmt;
    }
    else if (impType == Collections::eStdMat)
    {
        gflop = 2.0 * nElmt * nQ * nQ * nQ *
            LibUtilities::StdTetData::getNumberOfCoefficients(nM, nM, nM);
    }

    gflop *= 1e-9;
    matSize = 0;
}

int main(int argc, char *argv[])
{
    LibUtilities::SessionReader::RegisterCmdLineFlag(
        "data", "d", "Print in data format");
    LibUtilities::SessionReaderSharedPtr session
        = LibUtilities::SessionReader::CreateInstance(argc, argv);
    LibUtilities::CommSharedPtr vComm = session->GetComm();

    libxsmm_init();

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

    // BwdTrans operator
    Collections::ImplementationType impType = Collections::eIterPerExp;

    expList = SetupExpList(order + 1, session, graph, impType);
    Array<OneD, NekDouble> input (expList->GetNcoeffs(), 0.0);
    Array<OneD, NekDouble> output(expList->GetNpoints());

    // Do one BwdTrans operator
    expList->BwdTrans(input, output);

    // Loop over Ntest iterations
    Timer t;
    t.Start();
    for (int i = 0; i < Ntest; ++i)
    {
        expList->BwdTrans(input, output);
    }
    t.Stop();

    // Record average time for this rank
    NekDouble elapsed = t.TimePerTest(Ntest);

    // Reduce across all ranks
    vComm->AllReduce(elapsed, LibUtilities::ReduceSum);
    elapsed /= vComm->GetSize();

    // Calculate # flops for single operator
    int nElmt = expList->GetExpSize();
    NekDouble gflop, matSize;

    if (expList->GetExp(0)->DetShapeType() == LibUtilities::eHexahedron)
    {
        HexFlops(impType, order, nElmt, gflop, matSize);
    }
    else if (expList->GetExp(0)->DetShapeType() == LibUtilities::eTetrahedron)
    {
        TetFlops(impType, order, nElmt, gflop, matSize);
    }

    // Calculate total gflop
    vComm->AllReduce(gflop,   LibUtilities::ReduceSum);
    vComm->AllReduce(matSize, LibUtilities::ReduceSum);

    // Calculate gflop/sec
    NekDouble gflops = gflop / elapsed, bandwidth = matSize / elapsed;

    if (vComm->GetRank() == 0)
    {
        cout << order << " " << elapsed << " " << gflops //<< " " << bandwidth
             << endl;
    }

    vComm->Finalise();

    return 0;
}
