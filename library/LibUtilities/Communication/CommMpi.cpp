///////////////////////////////////////////////////////////////////////////////
//
// File CommMpi.cpp
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
// Description: MPI communication implementation
//
///////////////////////////////////////////////////////////////////////////////

#ifdef NEKTAR_USING_PETSC
#include "petscsys.h"
#endif

#include <iostream>
using namespace std;

#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/nonblocking.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/queue.hpp>
#include <boost/serialization/deque.hpp>
#include <boost/serialization/vector.hpp>
namespace mpi = boost::mpi;

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Communication/CommMpi.h>


namespace Nektar
{
namespace LibUtilities
{
std::string CommMpi::className = GetCommFactory().RegisterCreatorFunction(
    "ParallelMPI", CommMpi::create, "Parallel communication using MPI.");
int CommMpi::nSpares = 1;

/**
 *
 */
CommMpi::CommMpi(int narg, char *arg[]) : Comm(narg, arg)
{
    m_isLogging = false;
    m_isRecovering = false;

    int init = 0;
    MPI_Initialized(&init);
    ASSERTL0(!init, "MPI has already been initialised.");

    int retval = MPI_Init(&narg, &arg);
    if (retval != MPI_SUCCESS)
    {
        ASSERTL0(false, "Failed to initialise MPI");
    }

    int worldSize;
    int worldRank;
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);

    // Set MPI to call our error handler on failrue
    MPI_Errhandler errh;
    MPI_Comm_create_errhandler(HandleMpiError, &errh);
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, errh);

    MPI_Comm_dup(MPI_COMM_WORLD, &m_agreecomm);
    MPI_Comm_set_errhandler(m_agreecomm, MPI_ERRORS_RETURN);

    // Decide if we are a spare
    int spare = (worldRank > worldSize - nSpares - 1)? MPI_UNDEFINED : 1;

    // Create a communicator without the spares
    MPI_Comm_split( MPI_COMM_WORLD, spare, worldRank, &m_comm );

    // If we are a spare, sit and wait
    if ( MPI_COMM_NULL == m_comm )
    {
        std::cout << "Im a spare...rank " << worldRank << std::endl;
        do
        {
            // Always ready to complete
            int completed = 1;
            int x = MPIX_Comm_agree( m_agreecomm, &completed );
            std::cout << "Return value from comm agree is " << x << std::endl;
            if( completed )
            {
                std::cout << "Spare process invoking Finalize" << std::endl;
                MPI_Finalize();
                exit(0);
            }
            std::cout << "Spare process about to enroll" << std::endl;
            EnrolSpare();
            std::cout << "Completed enroll" << std::endl;
        } while ( MPI_COMM_NULL == m_comm );
        MPI_Comm_size(m_comm, &m_size);
        MPI_Comm_rank(m_comm, &m_rank);
    }
    else
    {
        std::cout << "Active process: rank " << worldRank << std::endl;
        MPI_Comm_size(m_comm, &m_size);
        MPI_Comm_rank(m_comm, &m_rank);
    }

#ifdef NEKTAR_USING_PETSC
    PetscInitializeNoArguments();
#endif

    m_type = "Parallel MPI";
}

/**
 *
 */
CommMpi::CommMpi(MPI_Comm pComm) : Comm()
{
    m_comm = pComm;
    MPI_Comm_size(m_comm, &m_size);
    MPI_Comm_rank(m_comm, &m_rank);

    m_type = "Parallel MPI";
    m_isRecovering = false;
    m_isLogging = false;
}

/**
 *
 */
CommMpi::~CommMpi()
{
    MPI_Comm_free(&m_comm);
}


/**
 *
 */
void CommMpi::v_Finalise()
{
#ifdef NEKTAR_USING_PETSC
    PetscFinalize();
#endif
    int flag;
    MPI_Finalized(&flag);
    if (!flag)
    {
        if ( MPI_COMM_NULL != m_comm )
        {
            std::cout << "Non-spare process invoked Finalize" << std::endl;
            int completed = 1;
            MPIX_Comm_agree(MPI_COMM_WORLD, &completed);
        }
        MPI_Comm_free(&m_comm);
        MPI_Finalize();
    }
}

/**
 *
 */
void* CommMpi::v_GetComm()
{
    return (void*)(m_comm);
}

/**
 *
 */
int CommMpi::v_GetRank()
{
    return m_rank;
}

/**
 *
 */
bool CommMpi::v_TreatAsRankZero(void)
{
    if (m_rank == 0)
    {
        return true;
    }
    else
    {
        return false;
    }
    return true;
}

/**
 *
 */
void CommMpi::v_Block()
{
    if (m_isRecovering)
    {
        return;
    }

    int rc = MPI_Barrier(m_comm);
    if (rc != MPI_SUCCESS)
    {
        throw 1;
    }
}

/**
 *
 */
double CommMpi::v_Wtime()
{
    return MPI_Wtime();
}

/**
 *
 */
void CommMpi::v_Send(void *buf, int count, CommDataType dt, int dest)
{
    if (m_isRecovering)
    {
        return;
    }

    if (MPISYNC)
    {
        MPI_Ssend(buf, count, dt, dest, 0, m_comm);
    }
    else
    {
        MPI_Send(buf, count, dt, dest, 0, m_comm);
    }
}

/**
 *
 */
void CommMpi::v_Recv(void *buf, int count, CommDataType dt, int source)
{
    if (m_isRecovering)
    {
        cout << "IMPLEMENTATION NEEDED" << endl;
    }
cout << "Recv" << endl;
    MPI_Recv(buf, count, dt, source, 0, m_comm, MPI_STATUS_IGNORE);
    // ASSERTL0(status.MPI_ERROR == MPI_SUCCESS,
    //         "MPI error receiving data.");
    int size;
    MPI_Type_size(dt, &size);
    std::vector<char> x(count * size);
    x.assign((char*)buf, (char*)buf + count*size);
    m_data.push(x);
    cout << "Received " << count * size << " bytes." << endl;
}

/**
 *
 */
void CommMpi::v_SendRecv(void *sendbuf, int sendcount, CommDataType sendtype,
                         int dest, void *recvbuf, int recvcount,
                         CommDataType recvtype, int source)
{
    MPI_Status status;
    int retval = MPI_Sendrecv(sendbuf, sendcount, sendtype, dest, 0, recvbuf,
                              recvcount, recvtype, source, 0, m_comm, &status);
cout << "SendRecv" << endl;
    ASSERTL0(retval == MPI_SUCCESS,
             "MPI error performing send-receive of data.");
}

/**
*
*/
void CommMpi::v_SendRecvReplace(void *buf, int count, CommDataType dt,
                                int pSendProc, int pRecvProc)
{
    MPI_Status status;
    int retval = MPI_Sendrecv_replace(buf, count, dt, pRecvProc, 0, pSendProc,
                                      0, m_comm, &status);

    ASSERTL0(retval == MPI_SUCCESS,
             "MPI error performing Send-Receive-Replace of data.");
}

/**
 *
 */
void CommMpi::v_AllReduce(void *buf, int count, CommDataType dt,
                          enum ReduceOperator pOp)
{
    if (GetSize() == 1)
    {
        return;
    }

    int dtsize;
    MPI_Type_size(dt, &dtsize);

    if (m_isRecovering)
    {
        ASSERTL0(!m_data.empty(), "QUEUE IS EMPTY!!");

        std::vector<char> x = m_data.front();
        m_data.pop();

        memcpy(buf, &x[0], count*dtsize);
        return;
    }

    MPI_Op vOp;
    switch (pOp)
    {
        case ReduceMax:
            vOp = MPI_MAX;
            break;
        case ReduceMin:
            vOp = MPI_MIN;
            break;
        case ReduceSum:
        default:
            vOp = MPI_SUM;
            break;
    }
    int retval = MPI_Allreduce(MPI_IN_PLACE, buf, count, dt, vOp, m_comm);

    ASSERTL0(retval == MPI_SUCCESS, "MPI error performing All-reduce.");

    if (m_isLogging)
    {
        std::vector<char> x;
        x.assign((char*)buf, (char*)buf+count*dtsize);
        m_data.push(x);
        cout << "AllReduce: Appended " << dtsize << " bytes of data." << endl;
    }
}

/**
 *
 */
void CommMpi::v_AlltoAll(void *sendbuf, int sendcount, CommDataType sendtype,
                         void *recvbuf, int recvcount, CommDataType recvtype)
{
    int retval = MPI_Alltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount,
                              recvtype, m_comm);
cout << "AllToAll" << endl;
    ASSERTL0(retval == MPI_SUCCESS, "MPI error performing All-to-All.");
}

/**
 *
 */
void CommMpi::v_AlltoAllv(void *sendbuf, int sendcounts[], int sdispls[],
                          CommDataType sendtype, void *recvbuf,
                          int recvcounts[], int rdispls[],
                          CommDataType recvtype)
{
    int retval = MPI_Alltoallv(sendbuf, sendcounts, sdispls, sendtype, recvbuf,
                               recvcounts, rdispls, recvtype, m_comm);

    ASSERTL0(retval == MPI_SUCCESS, "MPI error performing All-to-All-v.");
}

void CommMpi::v_Bcast(void *buffer, int count, CommDataType dt, int root)
{
    int dtsize;
    MPI_Type_size(dt, &dtsize);

    if (m_isRecovering)
    {
        ASSERTL0(!m_data.empty(), "QUEUE IS EMPTY!!");

        std::vector<char> x = m_data.front();
        m_data.pop();

        memcpy(buffer, &x[0], count*dtsize);
        return;
    }

    int retval = MPI_Bcast(buffer, count, dt, root, m_comm);
    ASSERTL0(retval == MPI_SUCCESS, "MPI error performing Bcast-v.");

    if (m_isLogging)
    {
        std::vector<char> x;
        x.assign((char*)buffer, (char*)buffer+count*dtsize);
        m_data.push(x);
        cout << "BCast: Appended " << dtsize << " bytes of data." << endl;
    }
}

void CommMpi::v_Exscan(Array<OneD, unsigned long long> &pData,
                       const enum ReduceOperator pOp,
                       Array<OneD, unsigned long long> &ans)
{
    int n = pData.num_elements();
    ASSERTL0(n == ans.num_elements(), "Array sizes differ in Exscan");

    MPI_Op vOp;
    switch (pOp)
    {
        case ReduceMax:
            vOp = MPI_MAX;
            break;
        case ReduceMin:
            vOp = MPI_MIN;
            break;
        case ReduceSum:
        default:
            vOp = MPI_SUM;
            break;
    }

    int retval = MPI_Exscan(pData.get(), ans.get(), n, MPI_UNSIGNED_LONG_LONG,
                            vOp, m_comm);
    ASSERTL0(retval == MPI_SUCCESS, "MPI error performing Exscan-v.");
}

void CommMpi::v_Gather(void *sendbuf, int sendcount, CommDataType sendtype,
                       void *recvbuf, int recvcount, CommDataType recvtype,
                       int root)
{
    int retval = MPI_Gather(sendbuf, sendcount, sendtype, recvbuf, recvcount,
                            recvtype, root, m_comm);

    ASSERTL0(retval == MPI_SUCCESS, "MPI error performing Gather.");
}

void CommMpi::v_Scatter(void *sendbuf, int sendcount, CommDataType sendtype,
                        void *recvbuf, int recvcount, CommDataType recvtype,
                        int root)
{
    int retval = MPI_Scatter(sendbuf, sendcount, sendtype, recvbuf, recvcount,
                             recvtype, root, m_comm);
    ASSERTL0(retval == MPI_SUCCESS, "MPI error performing Scatter.");
}

/**
 * Processes are considered as a grid of size pRows*pColumns. Comm
 * objects are created corresponding to the rows and columns of this
 * grid. The row and column to which this process belongs is stored in
 * #m_commRow and #m_commColumn.
 */
void CommMpi::v_SplitComm(int pRows, int pColumns)
{
    cout << "START SPLITCOMM" << endl;

    if (m_isRecovering)
    {
        cout << "Recovering row and column comm" << endl;
        m_commRow = m_derivedComm.front();
        m_derivedComm.pop_front();
        m_commColumn = m_derivedComm.front();
        m_derivedComm.pop_front();
        ASSERTL1(m_commRow->GetSize() == pColumns, "Row size does not match.");
        ASSERTL1(m_commColumn->GetSize() == pRows, "Column size does not match.");
        cout << "END SPLITCOMM (Recovering)" << endl;
        return;
    }

    ASSERTL0(pRows * pColumns == m_size,
             "Rows/Columns do not match comm size.");

    MPI_Comm newComm;

    // Compute row and column in grid.
    int myCol = m_rank % pColumns;
    int myRow = (m_rank - myCol) / pColumns;

    // Split Comm into rows - all processes with same myRow are put in
    // the same communicator. The rank within this communicator is the
    // column index.
    MPI_Comm_split(m_comm, myRow, myCol, &newComm);
    CommMpiSharedPtr commRowMpi = boost::shared_ptr<CommMpi>(new CommMpi(newComm));
    m_commRow = commRowMpi;
cout << "Original colour is: " << myRow << endl;
    // Split Comm into columns - all processes with same myCol are put
    // in the same communicator. The rank within this communicator is
    // the row index.
    MPI_Comm_split(m_comm, myCol, myRow, &newComm);
    CommMpiSharedPtr commColumnMpi = boost::shared_ptr<CommMpi>(new CommMpi(newComm));
    m_commColumn = commColumnMpi;
cout << "Original colour is: " << myCol << endl;
    if (m_isLogging)
    {
        cout << "Logging split comm" << endl;
        m_derivedCommFlag.push(myRow);
        m_derivedComm.push_back(commRowMpi);
        m_commRow->BeginTransactionLog();
        m_derivedCommFlag.push(myCol);
        m_derivedComm.push_back(commColumnMpi);
        m_commColumn->BeginTransactionLog();
    }
    cout << "END SPLITCOMM" << endl;
}

/**
 * Create a new communicator if the flag is non-zero.
 */
CommSharedPtr CommMpi::v_CommCreateIf(int flag)
{
    cout << "START CREATECOMMIF" << endl;
    CommMpiSharedPtr c;
    if (m_isRecovering)
    {
        c = m_derivedComm.front();
        m_derivedComm.pop_front();
    }
    else
    {
        MPI_Comm newComm;
        // color == MPI_UNDEF => not in the new communicator
        // key == 0 on all => use rank to order them. OpenMPI, at least,
        // implies this is faster than ordering them ourselves.
        MPI_Comm_split(m_comm, flag ? 0 : MPI_UNDEFINED, 0, &newComm);

        if (newComm == MPI_COMM_NULL)
        {
            // flag == 0 => get back MPI_COMM_NULL, return a null ptr instead.
            c = boost::shared_ptr<CommMpi>();
        }
        else
        {
            // Return a real communicator
            c = boost::shared_ptr<CommMpi>(new CommMpi(newComm));
        }

        if (m_isLogging)
        {
            if (c.get())
            {
                c->BeginTransactionLog();
            }
            m_derivedCommFlag.push(flag ? 0 : MPI_UNDEFINED);
            m_derivedComm.push_back(c);
        }
    }
    cout << "END CREATECOMMIF" << endl;
    return c;
}

int CommMpi::v_EnrolSpare()
{
    // Unlock spares so they join us
    if (MPI_COMM_NULL != m_comm)
    {
        int completed = 0;
        int x = MPIX_Comm_agree(m_agreecomm, &completed);
        std::cout << "Return value from comm agree is " << x << std::endl;
    }

    // =============================
    // First remove the dead process from our world (active + spare)
    MPI_Comm scomm, newcomm;
    int rc, flag, ssize, srank, oldsize, oldrank, dsize, drank;
    // Create a shunk comm from those live processes (scomm)
    MPIX_Comm_shrink(MPI_COMM_WORLD, &scomm);
    // Make sure we handle any new errors on scomm
    MPI_Comm_set_errhandler( scomm, MPI_ERRORS_RETURN );
    // Get size and rank
    MPI_Comm_size(scomm, &ssize);
    MPI_Comm_rank(scomm, &srank);

    // Keep trying until we succeed without further failures
    do
    {
        // I am a surviving rank, work out which ranks failed
        if (MPI_COMM_NULL != m_comm)
        {
            // Get our old rank and size
            MPI_Comm_size(m_comm, &oldsize);
            MPI_Comm_rank(m_comm, &oldrank);

            // First check we have enough spares left to replace all those
            // which have failed
            if ( oldsize > ssize )
            {
                MPI_Abort(MPI_COMM_WORLD, MPI_ERR_PROC_FAILED);
            }

            MPI_Group cgrp, sgrp, dgrp;

            // Get the group of dead processes
            MPI_Comm_group(m_comm, &cgrp);
            MPI_Comm_group(scomm,  &sgrp);
            MPI_Group_difference(cgrp, sgrp, &dgrp);
            MPI_Group_size(dgrp, &dsize);

            // Let rank 0 in the shrunk comm determine new assignments
            if ( 0 == srank )
            {
                // Give every spare rank a new assignment
                // [ C C C C C C C C C C ]
                // [ C C C D C D C S S S ]
                for(int i = 0; i < ssize - (oldsize - dsize); i++) {
                    // Assign spares to cover the dead ranks
                    // For each rank in the dead-group, we find the corresponding
                    // rank in the full communicator
                    if( i < dsize )
                    {
                        MPI_Group_translate_ranks(dgrp, 1, &i, cgrp, &drank);
                    }
                    // Any additional spares we have will not be in the newly
                    // created communicator created shortly.
                    else
                    {
                        drank=-1; /* still a spare */
                    }

                    // send their new assignment to all spares
                    MPI_Send(&drank, 1, MPI_INT, i + oldsize - dsize, 1, scomm);
                }
            }

            MPI_Group_free(&cgrp);
            MPI_Group_free(&sgrp);
            MPI_Group_free(&dgrp);
        }
        // I am a spare waiting for my assignment
        else
        {
            MPI_Recv(&oldrank, 1, MPI_INT, 0, 1, scomm, MPI_STATUS_IGNORE);
            std::cout << "Spare received assignment: " << oldrank << std::endl;

            m_isRecovering = true;
        }

        // Remove dead process and reassign spare processes to these ranks
        // oldrank contains the original rank for those processes which are
        // still alive, and contains the new assignment for spares which are
        // needed to replace dead processes.
        //
        // Spares which stay spare are assigned MPI_UNDEFINED and therefore
        // are not in the new communicator.
        rc = MPI_Comm_split(scomm, oldrank < 0 ? MPI_UNDEFINED : 1, oldrank, &newcomm);

        flag = MPIX_Comm_agree(scomm, &flag);
        MPI_Comm_free(&scomm);

        if( MPI_SUCCESS != flag ) {
            if( MPI_SUCCESS == rc )
            {
                MPI_Comm_free(&newcomm);
            }
        }
    } while ( MPI_SUCCESS != flag );

    // Replace the original comm
    if (MPI_COMM_NULL != m_comm)
    {
        MPI_Comm_free(&m_comm);
    }
    m_comm = newcomm;

    // Update rank and size
    MPI_Comm_rank(m_comm, &m_rank);
    MPI_Comm_size(m_comm, &m_size);

    RestoreState();

    return MPI_SUCCESS;
}

void CommMpi::BackupState()
{
    mpi::communicator c(m_comm, mpi::comm_attach);
    mpi::request reqs[4];
    int rank = c.rank();
    int size = c.size();

    if (size > 1)
    {
        int recv_rank = (rank + size - 1) % size;
        int send_rank = (rank + 1) % size;
        cout << "Backup: Sending " << m_data.size() << " items in queue." << endl;
        reqs[0] = c.isend(send_rank, 0, m_data);
        reqs[1] = c.isend(send_rank, 1, m_derivedCommFlag);
        reqs[2] = c.irecv(recv_rank, 0, m_dataBackup);
        reqs[3] = c.irecv(recv_rank, 1, m_derivedCommFlagBackup);
        cout << "Backup: Sent to " << send_rank << endl;
        cout << "Backup: Waiting for data from " << recv_rank << endl;
        mpi::wait_all(reqs, reqs + 4);
        cout << "Backup: Received " << m_dataBackup.size() << " items." << endl;
        cout << "Backup: Received " << m_derivedCommFlagBackup.size() << " derived comm flags." << endl;
    }
    else
    {
        cout << "Backup: Not backing up as comm of size 1" << endl;
    }
}

void CommMpi::RestoreState()
{
    cout << "Start Restore" << endl;
    if (m_size > 1)
    {
        cout << " -> More than 1 process" << endl;
        mpi::communicator c(m_comm, mpi::comm_attach);
        int send_rank = (m_rank + m_size - 1) % m_size;
        int recv_rank = (m_rank + 1) % m_size;

        // First work out who needs to send recovery data
        cout << " -> Determining recovery participants" << endl;
        int restore = m_isRecovering ? 1 : 0;
        int sendRecoveryData = 0;
        mpi::request reqs[2];
        // Retrieve recovery state of my buddy process
        reqs[0] = c.irecv(send_rank, 0, sendRecoveryData);
        // Send my recovery state
        reqs[1] = c.isend(recv_rank, 0, restore);
        mpi::wait_all(reqs, reqs + 2);

        cout << " -> I am rank: " << m_rank << endl;
        cout << " -> Send recovery data: " << sendRecoveryData << endl;
        cout << " -> Recovering? " << m_isRecovering << endl;

        // Perform restore of data and derived comm flags
        // Neighbouring surviving processes to recovering processes must send
        // a) backup of the recovering process's data for its recovery
        // b) replacement copy of this processes backup data
        // c) queue of flags for recovering derived communicators on recovering process
        const int nReq = 4;
        if (sendRecoveryData)
        {
            cout << "Restore: Sending " << m_dataBackup.size() << " backup items." << endl;
            cout << "Restore: Sending my backup copy to " << send_rank << endl;
            mpi::request reqs[nReq];
            reqs[0] = c.isend(send_rank, 0, m_dataBackup);
            reqs[1] = c.isend(send_rank, 1, m_data);
            reqs[2] = c.isend(send_rank, 2, m_derivedCommFlagBackup);
            reqs[3] = c.isend(send_rank, 3, m_derivedCommFlag);
            cout << "Restore: Waiting for data to be sent to " << send_rank << endl;
            mpi::wait_all(reqs, reqs + nReq);
            cout << "Restore: Complete" << endl;
            cout << "Restore: Sent " << m_data.size() << " items in queue." << endl;
        }
        // The recovering processes receive this data
        // They do not need to send anything
        else if (m_isRecovering)
        {
            mpi::request reqs[nReq];
            cout << "Restore: Receiving from " << recv_rank << endl;
            reqs[0] = c.irecv(recv_rank, 0, m_data);
            reqs[1] = c.irecv(recv_rank, 1, m_dataBackup);
            reqs[2] = c.irecv(recv_rank, 2, m_derivedCommFlag);
            reqs[3] = c.irecv(recv_rank, 3, m_derivedCommFlagBackup);
            cout << "Restore: Waiting for data from " << recv_rank << endl;
            mpi::wait_all(reqs, reqs + nReq);
            cout << "Restore: Complete" << endl;
            cout << "Restore: Received " << m_data.size() << " items in queue." << endl;
            cout << "Restore: There are " << m_derivedCommFlag.size() << " derived comms" << endl;
        }
        // Any processes not recovering, and not a neighbour, do not need to
        // participate in recovery.
    }


    // Fix derived comms
    // -----------------
    // All processes must participate in this. This includes recovering the
    // m_rowComm and m_columnComm too, as these are also just derived.
    DerivedCommFlagType vDerivedFlagRestore = m_derivedCommFlag;
    auto derivedCommIt = m_derivedComm.begin();
    cout << "About to fix derived comms." << endl;
    cout << " -> there are currently " << m_derivedComm.size() << " comms." << endl;
    for (int i = 0; i < m_derivedCommFlag.size(); ++i)
    {
        MPI_Comm tmpComm;
        int colour = vDerivedFlagRestore.front();
        vDerivedFlagRestore.pop();
cout << "Colour is " << colour << endl;
        MPI_Comm_split(m_comm, colour, 0, &tmpComm);
        if (m_isRecovering) // Recovering spare node
        {
            // Create a derived communicator from scratch
            CommMpiSharedPtr cmpi;
            if (tmpComm == MPI_COMM_NULL)
            {
                // flag == 0 => get back MPI_COMM_NULL, return a null ptr instead.
                cmpi = boost::shared_ptr<CommMpi>();
                cout << "Restored a null comm" << endl;
            }
            else
            {
                // Return a real communicator
                cmpi = boost::shared_ptr<CommMpi>(new CommMpi(tmpComm));
                cmpi->m_isRecovering = true;
                cout << "Restored a comm" << endl;
            }
            m_derivedComm.push_back(cmpi);
            cout << "New comm is of size: " << cmpi->GetSize() << endl;
            cout << "Restoring its state..." << endl << endl;
            cout << "### Restore sub-communicator ###" << endl;
            cmpi->RestoreState();
            cout << "### Finished restore of sub-communicator ###" << endl << endl;
        }
        else // surviving node
        {
            CommMpiSharedPtr cmpi = (*derivedCommIt);
            cout << "Previously comm is of size: " << cmpi->GetSize() << endl;
            cmpi->ReplaceComm(tmpComm);
            cout << "Replaced comm is of size: " <<
                    cmpi->GetSize() << endl;
            int size;
            MPI_Comm_size(cmpi->m_comm, &size);
            cout << "Underlying comm is of size: " <<  size << endl;
            cout << "Restoring its state..." << endl << endl;
            cout << "### Restore sub-communicator ###" << endl;
            cmpi->RestoreState();
            cout << "### Finished restore of sub-communicator ###" << endl << endl;
            derivedCommIt++;
        }

    }
}

static void CommMpi::HandleMpiError(MPI_Comm* pcomm, int* perr, ...)
{
    int err = *perr;

    // Get type of error and check if it is a proc failed.
    int eclass;
    MPI_Error_class(err, &eclass);
    if (MPIX_ERR_PROC_FAILED != eclass)
    {
        std::cout << "An non-proc-failed MPI error occured..." << std::endl;
        //MPI_Abort(comm, err);
    }

    int len;
    char errstr[MPI_MAX_ERROR_STRING];
    MPI_Error_string(err, errstr, &len);

    std::cout << "An MPI error occured: " << errstr << std::endl;

    throw 1; //UlfmFailureDetected(std::string(errstr));
}

void CommMpi::v_BeginTransactionLog()
{
    m_isLogging = true;
    if (m_isRecovering) {
        for (auto x : m_derivedComm)
        {
            x->m_isRecovering = true;
        }
    }
}

void CommMpi::v_EndTransactionLog()
{
    m_isLogging = false;
    m_isRecovering = false;

    for (auto x : m_derivedComm)
    {
        std::cout << "Ending transaction log on derived comm." << std::endl;
        x->EndTransactionLog();
    }

    BackupState();
}

void CommMpi::ReplaceComm(MPI_Comm commptr)
{
    m_comm = commptr;
    MPI_Comm_size(m_comm, &m_size);
    MPI_Comm_rank(m_comm, &m_rank);

    m_type = "Parallel MPI";
}

}
}
