///////////////////////////////////////////////////////////////////////////////
//
// File CommMpi.h
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
// Description: CommMpi header
//
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_LIB_UTILITIES_COMMMPI_H
#define NEKTAR_LIB_UTILITIES_COMMMPI_H

#include <mpi.h>
#include <mpi-ext.h>
#include <string>
#include <queue>
#include <vector>
#include <list>


#include <LibUtilities/Communication/Comm.h>
#include <LibUtilities/Memory/NekMemoryManager.hpp>

#ifndef MPI_SYNC
#define MPISYNC 0
#else
#define MPISYNC 1
#endif

namespace Nektar
{
namespace LibUtilities
{
// Forward declarations
class CommMpi;

/// Pointer to a Communicator object.
typedef boost::shared_ptr<CommMpi> CommMpiSharedPtr;

/**
 * Exception class for Ulfm-detected failure
 */
class UlfmFailureDetected: public std::runtime_error
{
    public:
        UlfmFailureDetected(std::string && pError)
            : std::runtime_error(pError.c_str()), m_error(pError) {}

    private:
        virtual const char* what() const throw()
        {
            return m_error.c_str();
        }

        std::string m_error;
};


/// A communicator which uses MPI
class CommMpi : public Comm
{
public:
    /// Creates an instance of this class
    static CommSharedPtr create(int narg, char *arg[])
    {
        return MemoryManager<CommMpi>::AllocateSharedPtr(narg, arg);
    }

    /// Name of class
    static std::string className;

    static int nSpares;

    CommMpi(int narg, char *arg[]);
    virtual ~CommMpi();

protected:
    virtual void v_Finalise();
    virtual void* v_GetComm();
    virtual int v_GetRank();
    virtual void v_Block();
    virtual double v_Wtime();
    virtual bool v_TreatAsRankZero(void);
    virtual void v_Send(void *buf, int count, CommDataType dt, int dest);
    virtual void v_Recv(void *buf, int count, CommDataType dt, int source);
    virtual void v_SendRecv(void *sendbuf, int sendcount, CommDataType sendtype,
                            int dest, void *recvbuf, int recvcount,
                            CommDataType recvtype, int source);
    virtual void v_SendRecvReplace(void *buf, int count, CommDataType dt,
                                   int pSendProc, int pRecvProc);
    virtual void v_AllReduce(void *buf, int count, CommDataType dt,
                             enum ReduceOperator pOp);
    virtual void v_AlltoAll(void *sendbuf, int sendcount, CommDataType sendtype,
                            void *recvbuf, int recvcount,
                            CommDataType recvtype);
    virtual void v_AlltoAllv(void *sendbuf, int sendcounts[], int sensdispls[],
                             CommDataType sendtype, void *recvbuf,
                             int recvcounts[], int rdispls[],
                             CommDataType recvtype);
    virtual void v_Bcast(void *buffer, int count, CommDataType dt, int root);
    virtual void v_Exscan(Array<OneD, unsigned long long> &pData,
                          const enum ReduceOperator pOp,
                          Array<OneD, unsigned long long> &ans);

    virtual void v_Gather(void *sendbuf, int sendcount, CommDataType sendtype,
                          void *recvbuf, int recvcount, CommDataType recvtype,
                          int root);
    virtual void v_Scatter(void *sendbuf, int sendcount, CommDataType sendtype,
                           void *recvbuf, int recvcount, CommDataType recvtype,
                           int root);

    virtual void v_SplitComm(int pRows, int pColumns);
    virtual CommSharedPtr v_CommCreateIf(int flag);

private:
    typedef std::queue<std::vector<char>> StorageType;
    typedef std::list<CommMpiSharedPtr>   DerivedCommType;
    typedef std::queue<int>               DerivedCommFlagType;

    MPI_Comm m_comm;
    MPI_Comm m_agreecomm;
    int m_rank;

    bool m_isRecovering;        ///< True if we are undergoing recovery from failed process
    bool m_isLogging;           ///< True if logging MPI output
    StorageType m_data;
    StorageType m_dataBackup;
    DerivedCommType m_derivedComm; ///< Temporary derived comm list used during restore
    DerivedCommFlagType m_derivedCommFlag; ///< Log derived comm flags
    DerivedCommFlagType m_derivedCommFlagBackup; ///< Backup of neighbour flags


    static void HandleMpiError(MPI_Comm* pcomm, int* perr, ...);

    CommMpi(MPI_Comm pComm);

    virtual int v_EnrolSpare();
    virtual void v_BeginTransactionLog();
    virtual void v_EndTransactionLog();

    void BackupState();
    void RestoreState();
    void ReplaceComm(MPI_Comm commptr);

};
}
}

#endif
