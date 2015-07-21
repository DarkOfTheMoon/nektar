///////////////////////////////////////////////////////////////////////////////
//
// File IncNavierStokesSolver.cpp
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
// Description: Incompressible Navier Stokes solver
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Driver.h>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/Thread.h>
#include <LibUtilities/Communication/Comm.h>
#include <ctime>

using namespace Nektar;
using namespace Nektar::SolverUtils;

class MPIJob;

class MPIJob : public Nektar::Thread::ThreadJob
{
    private:
    int m_cycles;
    int m_rank;
    LibUtilities::SessionReaderSharedPtr m_session;
    public:
    MPIJob( int cycles, int rank, LibUtilities::SessionReaderSharedPtr session) :
        m_cycles(cycles), m_rank(rank), m_session(session)
    {
        //empty
    }

    void Run()
    {
        cout << m_rank << ": Starting WORK part" << endl;
        m_session->DoMPIWait(m_cycles);
        cout << m_rank << ": Ended WORK part" << endl;
    }
};

int main(int argc, char *argv[])
{
    LibUtilities::SessionReaderSharedPtr session;
    string vDriverModule;
    DriverSharedPtr drv;

  
    try
    {
        // Create session reader.
        session = LibUtilities::SessionReader::CreateInstance(argc, argv);



        Nektar::Thread::ThreadManagerSharedPtr tm;
        tm = Nektar::Thread::GetThreadMaster().GetInstance(Nektar::Thread::ThreadMaster::SessionJob);
        tm->Wait();

        int size = 5e6;
        int sleeps = 10;
        int cycles = 1e9;

        Nektar::LibUtilities::CommSharedPtr comm = session->GetComm();
        ASSERTL0(comm->GetSize() > 1, "Need more than 1 MPI process");
        int rank = comm->GetRank();

        tm->QueueJob( new MPIJob(cycles, rank, session) );

        if (rank == 0)
        {
            std::vector<unsigned int> *vec = new std::vector<unsigned int>(size);
            for (int i=1; i<comm->GetSize(); i++)
            {
                cout << "0: Receiving data from " << i << endl;
                comm->Recv(i,*vec);
                cout << "0: Received data from " << i << endl;
            }
        } else
        {
            std::vector<unsigned int> *vec = new std::vector<unsigned int>(size,99);
            timespec ts;
            sleeps *= rank;
            ts.tv_sec  = sleeps;
            ts.tv_nsec = 0;
            cout << rank << ": Sleeping (" << sleeps << " seconds)" << endl;
            nanosleep(&ts, NULL);
            cout << rank << ": Finished sleeping" << endl;
            cout << rank << ": Sending data (" << size << " ints)" << endl;
            comm->Send(0,*vec);
            cout << rank << ": Sent" << endl;
        }
        
        // Finalise communications
        tm->Wait();
        session->Finalise();
    }
    catch (const std::runtime_error&)
    {
        return 1;
    }
    catch (const std::string& eStr)
    {
        cout << "Error: " << eStr << endl;
    }
    
    return 0;
}
