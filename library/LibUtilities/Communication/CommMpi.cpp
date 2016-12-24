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
    if (IsRecovering())
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
    if (IsRecovering())
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
    if (IsRecovering())
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
    if (m_isRecovering)
    {
        cout << "Assuming split comm already sorted by enrolspare." << endl;
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
    m_commRow = boost::shared_ptr<Comm>(new CommMpi(newComm));

    // Split Comm into columns - all processes with same myCol are put
    // in the same communicator. The rank within this communicator is
    // the row index.
    MPI_Comm_split(m_comm, myCol, myRow, &newComm);
    m_commColumn = boost::shared_ptr<Comm>(new CommMpi(newComm));
}

/**
 * Create a new communicator if the flag is non-zero.
 */
CommSharedPtr CommMpi::v_CommCreateIf(int flag)
{
    CommSharedPtr c;
    if (m_isRecovering)
    {
        c = m_derivedComm[m_derivedRecoverIndex++];
    }
    else
    {
        MPI_Comm newComm;
        // color == MPI_UNDEF => not in the new communicator
        // key == 0 on all => use rank to order them. OpenMPI, at least,
        // implies this is faster than ordering them ourselves.
        MPI_Comm_split(m_comm, flag ? 0 : MPI_UNDEFINED, 0, &newComm);

        if (flag == 0)
        {
            // flag == 0 => get back MPI_COMM_NULL, return a null ptr instead.
            c = boost::shared_ptr<Comm>();
        }
        else
        {
            // Return a real communicator
            c = boost::shared_ptr<Comm>(new CommMpi(newComm));
        }

        m_derivedComm.push_back(c);
        m_derivedCommFlag.push_back(flag);
    }
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

    // Get number of rows and columns if split comm
    int pRows    = (m_commColumn.get() ? m_commColumn->GetSize() : 0);
    int pColumns = (m_commRow.get() ? m_commRow->GetSize() : 0);
    int pDerived = m_derivedComm.size();

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

            // Let rank 0 in the shrunk comm determine new assignments
            if ( 0 == srank )
            {
                MPI_Group cgrp, sgrp, dgrp;

                // Get the group of dead processes
                MPI_Comm_group(m_comm, &cgrp);
                MPI_Comm_group(scomm,  &sgrp);
                MPI_Group_difference(cgrp, sgrp, &dgrp);
                MPI_Group_size(dgrp, &dsize);

                for(int i = 0; i < ssize - (oldsize - dsize); i++) {
                    if( i < dsize )
                    {
                        MPI_Group_translate_ranks(dgrp, 1, &i, cgrp, &drank);
                    }
                    else
                    {
                        drank=-1; /* still a spare */
                    }
                    // send their new assignment to all spares
                    MPI_Send(&drank, 1, MPI_INT, i + oldsize - dsize, 1, scomm);
                    MPI_Send(&pRows, 1, MPI_INT, i + oldsize - dsize, 2, scomm);
                    MPI_Send(&pColumns, 1, MPI_INT, i + oldsize - dsize, 3, scomm);
                    MPI_Send(&pDerived, 1, MPI_INT, i + oldsize - dsize, 4, scomm);
                    MPI_Send(&m_derivedCommFlag[0], pDerived,
                                MPI_INT, i + oldsize - dsize, 5, scomm);
                }

                MPI_Group_free(&cgrp);
                MPI_Group_free(&sgrp);
                MPI_Group_free(&dgrp);
            }
        }
        // I am a spare waiting for my assignment
        else
        {
            MPI_Recv(&oldrank, 1, MPI_INT, 0, 1, scomm, MPI_STATUS_IGNORE);
            MPI_Recv(&pRows, 1, MPI_INT, 0, 2, scomm, MPI_STATUS_IGNORE);
            MPI_Recv(&pColumns, 1, MPI_INT, 0, 3, scomm, MPI_STATUS_IGNORE);
            MPI_Recv(&pDerived, 1, MPI_INT, 0, 4, scomm, MPI_STATUS_IGNORE);
            m_derivedCommFlag.resize(pDerived);
            MPI_Recv(&m_derivedCommFlag[0], pDerived, MPI_INT, 0, 5, scomm, MPI_STATUS_IGNORE);
            m_derivedRecoverIndex = 0;

            std::cout << "Spare received assignment: " << oldrank << std::endl;
            std::cout << "Split: " << pRows << " x " << pColumns << endl;
            std::cout << "Num of derived comms: " << pDerived << endl;

            m_isRecovering = true;
        }

        // Remove dead process and reassign spare processes to these ranks
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

    // Fix split comm
    // --------------
    // First, free old row and column communicators
    if (m_commColumn.get())
    {
        m_commColumn.reset();
    }
    if (m_commRow.get())
    {
        m_commRow.reset();
    }

    if (pRows > 0 && pColumns > 0)
    {
        ASSERTL0(pRows * pColumns == m_size,
                 "Rows/Columns do not match comm size.");

        MPI_Comm tmpComm;

        // Compute row and column in grid.
        int myCol = m_rank % pColumns;
        int myRow = (m_rank - myCol) / pColumns;

        // Split Comm into rows - all processes with same myRow are put in
        // the same communicator. The rank within this communicator is the
        // column index.
        MPI_Comm_split(m_comm, myRow, myCol, &tmpComm);
        m_commRow = boost::shared_ptr<Comm>(new CommMpi(tmpComm));

        // Split Comm into columns - all processes with same myCol are put
        // in the same communicator. The rank within this communicator is
        // the row index.
        MPI_Comm_split(m_comm, myCol, myRow, &tmpComm);
        m_commColumn = boost::shared_ptr<Comm>(new CommMpi(tmpComm));
    }

    // Fix derived comms
    // -----------------
    for (int i = 0; i < pDerived; ++i)
    {
        MPI_Comm tmpComm;
        int flag = m_derivedCommFlag[i];
        MPI_Comm_split(m_comm, flag ? 0 : MPI_UNDEFINED, 0, &tmpComm);
        if (m_isRecovering) // Recovering spare node
        {
            // Create a derived communicator from scratch
            CommSharedPtr c;
            if (flag == 0)
            {
                // flag == 0 => get back MPI_COMM_NULL, return a null ptr instead.
                c = boost::shared_ptr<Comm>();
            }
            else
            {
                // Return a real communicator
                c = boost::shared_ptr<Comm>(new CommMpi(tmpComm));
            }
            m_derivedComm.push_back(c);
        }
        else // surviving node
        {
            m_derivedComm[i]->v_ReplaceComm((void*)(tmpComm));
        }
    }

    // Perform restore
    mpi::communicator c(m_comm, mpi::comm_attach);
    int rank = c.rank();
    int size = c.size();
    int recv_rank = (rank + 1) % size;
    int send_rank = (rank + size - 1) % size;
    mpi::request reqs[2];
    cout << "Restore: Receiving from " << recv_rank << endl;
    reqs[0] = c.irecv(recv_rank, 0, m_data);
    cout << "Restore: Sending " << m_dataBackup.size() << " items." << endl;
    cout << "Restore: Sending my backup copy to " << send_rank << endl;
    reqs[1] = c.isend(send_rank, 0, m_dataBackup);
    cout << "Restore: Waiting for data from " << recv_rank << endl;
    reqs[0].wait();
    cout << "Restore: Complete" << endl;
    cout << "Restore: Received " << m_data.size() << " items in queue." << endl;

    return MPI_SUCCESS;
}

void CommMpi::v_BackupState()
{
    mpi::communicator c(m_comm, mpi::comm_attach);
    mpi::request reqs[2];
    int rank = c.rank();
    int size = c.size();
    int recv_rank = (rank + size - 1) % size;
    int send_rank = (rank + 1) % size;
    cout << "Backup: Sending " << m_data.size() << " items in queue." << endl;
    reqs[0] = c.isend(send_rank, 0, m_data);
    cout << "Backup: Sent to " << send_rank << endl;
    reqs[1] = c.irecv(recv_rank, 0, m_dataBackup);
    cout << "Backup: Waiting for data from " << recv_rank << endl;
    reqs[1].wait();
    cout << "Backup: Received " << m_dataBackup.size() << " items." << endl;
//    cout << "Backup: Sending " << m_data.size() << " items in queue." << endl;
//    c.send((rank + 1) % size, 0, m_data);
//    cout << "Backup: Sent to " << (rank + 1) % size << endl;
//    cout << "Backup: Waiting for data from " << (rank + 1) % size << endl;
//    c.recv((rank - 1) % size, 0, m_dataBackup);
//    cout << "Backup: Received " << m_dataBackup.size() << " items." << endl;

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

void CommMpi::v_ReplaceComm(void* commptr)
{
    m_comm = (MPI_Comm)(commptr);
}

}
}
