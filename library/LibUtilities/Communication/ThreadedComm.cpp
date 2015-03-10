/*
 * ThreadedComm.cpp
 *
 *  Created on: 11 Mar 2013
 *      Author: simon
 */

#include <LibUtilities/Communication/ThreadedComm.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <iostream>


namespace Nektar
{
	namespace LibUtilities
	{
		using namespace Thread;
		ThreadedComm::ThreadedComm(CommSharedPtr comm, ThreadManagerSharedPtr tm) :
				m_comm(comm), m_tm(tm), m_numworkers(m_tm->GetMaxNumWorkers()),
				m_numMPI(m_comm->GetSize()), m_resDbl(m_numworkers),
				m_resInt(m_numworkers), m_ListOfThrSendDataDbl(m_numworkers),
				m_ListOfConstThrSendDataDbl(m_numworkers),
				m_ListOfThrSendDataInt(m_numworkers),
                m_ListOfThrSendDataUIntVector(m_numworkers),
				m_ListOfConstThrSendDataLon(m_numworkers),
				m_SendOffsetArr(m_numMPI),
				m_SendSizeArr(m_numMPI),
				m_RecvOffsetArr(m_numMPI), m_RecvSizeArr(m_numMPI),
				m_ListOfThrSendOffsets(m_numworkers),
				m_ListOfThrSendSizes(m_numworkers), m_ListOfThrRecvOffsets(m_numworkers),
				m_ListOfThrRecvSizes(m_numworkers),
				m_tmpSendOffsetArr(m_numworkers * m_numworkers * m_numMPI + 1),
				m_tmpRecvOffsetArr(m_numworkers * m_numworkers * m_numMPI + 1),
                m_sendToRankZeroVecInt(m_numworkers),
                m_sendToRankZeroArrInt(m_numworkers),
                m_sendToRankZeroArrDbl(m_numworkers)

		{
			ASSERTL0(m_numworkers > 1, "ThreadedComm cannot be used unless there is more than 1 thread");
            m_size = m_comm->GetSize() * m_numworkers;

            m_type = m_comm->GetType() + " + " + m_tm->GetType();

		}

		ThreadedComm::~ThreadedComm()
		{
			// empty
		}

		void ThreadedComm::v_Finalise()
		{
			if (m_tm->GetWorkerNum() == 0 || !m_tm->InThread())
			{
				m_comm->Finalise();
			}
		}

		int ThreadedComm::v_GetRank()
		{
			return m_comm->GetRank() * m_numworkers + m_tm->GetWorkerNum();
		}

        /**
         *
         */
        bool ThreadedComm::v_TreatAsRankZero(void)
        {
            return (v_GetRank() == 0);
        }

		void ThreadedComm::v_Block()
		{
			if (m_tm->GetWorkerNum() == 0 || !m_tm->InThread())
			{
				m_comm->Block();
			}
			if (m_tm->InThread())
			{
				m_tm->Hold();
			}
		}

		void ThreadedComm::v_Send(int pProc, Array<OneD, NekDouble>& pData)
		{
			ASSERTL1(m_tm->GetThrFromPartition(pProc) == 0 && m_tm->GetWorkerNum() == 0,
					"Cannot use Comm->Send() with threaded parallelism except between master threads");
			pProc = m_tm->GetRankFromPartition(pProc);
			m_comm->Send(pProc, pData);
		}

		void ThreadedComm::v_Send(int pProc, Array<OneD, int>& pData)
		{
			ASSERTL1(m_tm->GetThrFromPartition(pProc) == 0 && m_tm->GetWorkerNum() == 0,
					"Cannot use Comm->Send() with threaded parallelism except between master threads");
			pProc = m_tm->GetRankFromPartition(pProc);
			m_comm->Send(pProc, pData);
		}

        void ThreadedComm::v_Send(int pProc, std::vector<unsigned int>& pData)
		{
			ASSERTL1(m_tm->GetThrFromPartition(pProc) == 0 && m_tm->GetWorkerNum() == 0,
					"Cannot use Comm->Send() with threaded parallelism except between master threads");
			pProc = m_tm->GetRankFromPartition(pProc);
			m_comm->Send(pProc, pData);
		}

		void ThreadedComm::v_Recv(int pProc, Array<OneD, NekDouble>& pData)
		{
			ASSERTL1(m_tm->GetThrFromPartition(pProc) == 0 && m_tm->GetWorkerNum() == 0,
					"Cannot use Comm->Recv() with threaded parallelism except between master threads");
			pProc = m_tm->GetRankFromPartition(pProc);
			m_comm->Recv(pProc, pData);
		}

		void ThreadedComm::v_Recv(int pProc, Array<OneD, int>& pData)
		{
			ASSERTL1(m_tm->GetThrFromPartition(pProc) == 0 && m_tm->GetWorkerNum() == 0,
					"Cannot use Comm->Recv() with threaded parallelism except between master threads");
			pProc = m_tm->GetRankFromPartition(pProc);
			m_comm->Recv(pProc, pData);
		}

        void ThreadedComm::v_Recv(int pProc, std::vector<unsigned int>& pData)
		{
			ASSERTL1(m_tm->GetThrFromPartition(pProc) == 0 && m_tm->GetWorkerNum() == 0,
					"Cannot use Comm->Recv() with threaded parallelism except between master threads");
			pProc = m_tm->GetRankFromPartition(pProc);
			m_comm->Recv(pProc, pData);
		}

		void ThreadedComm::v_SendRecv(int pSendProc, Array<OneD, NekDouble>& pSendData,
				int pRecvProc, Array<OneD, NekDouble>& pRecvData)
		{
			ASSERTL1(m_tm->GetThrFromPartition(pSendProc) == 0 && m_tm->GetThrFromPartition(pRecvProc) == 0 && m_tm->GetWorkerNum() == 0,
					"Cannot use Comm->SendRecv() with threaded parallelism except between master threads");
			m_comm->SendRecv(pSendProc, pSendData, pRecvProc, pRecvData);
		}

		void ThreadedComm::v_SendRecv(int pSendProc, Array<OneD, int>& pSendData,
				int pRecvProc, Array<OneD, int>& pRecvData)
		{
			ASSERTL1(m_tm->GetThrFromPartition(pSendProc) == 0 && m_tm->GetThrFromPartition(pRecvProc) == 0 && m_tm->GetWorkerNum() == 0,
					"Cannot use Comm->SendRecv() with threaded parallelism except between master threads");
			m_comm->SendRecv(pSendProc, pSendData, pRecvProc, pRecvData);
		}

		void ThreadedComm::v_SendRecvReplace(int pSendProc, int pRecvProc,
				Array<OneD, NekDouble>& pSendData)
		{
			ASSERTL1(m_tm->GetThrFromPartition(pSendProc) == 0 && m_tm->GetThrFromPartition(pRecvProc) == 0 && m_tm->GetWorkerNum() == 0,
					"Cannot use Comm->SendRecvReplace() with threaded parallelism except between master threads");
			pSendProc = m_tm->GetRankFromPartition(pSendProc);
			pRecvProc = m_tm->GetRankFromPartition(pRecvProc);
			m_comm->SendRecvReplace(pSendProc, pRecvProc, pSendData);
		}

		void ThreadedComm::v_SendRecvReplace(int pSendProc, int pRecvProc,
				Array<OneD, int>& pSendData)
		{
			ASSERTL1(m_tm->GetThrFromPartition(pSendProc) == 0 && m_tm->GetThrFromPartition(pRecvProc) == 0 && m_tm->GetWorkerNum() == 0,
					"Cannot use Comm->SendRecvReplace() with threaded parallelism except between master threads");
			pSendProc = m_tm->GetRankFromPartition(pSendProc);
			pRecvProc = m_tm->GetRankFromPartition(pRecvProc);
			m_comm->SendRecvReplace(pSendProc, pRecvProc, pSendData);
		}

		void ThreadedComm::v_AllReduce(NekDouble& pData, enum ReduceOperator pOp)
		{
			GenericAllReduce(pData, pOp, m_resDbl);
		}

		void ThreadedComm::v_AllReduce(int& pData, enum ReduceOperator pOp)
		{
			GenericAllReduce(pData, pOp, m_resInt);
		}

		template <typename DataType>
		void ThreadedComm::GenericAllReduce(DataType& pData, enum ReduceOperator pOp,
				std::vector<DataType>& pRes)
		{
			unsigned int vThr = m_tm->GetWorkerNum();
			pRes[vThr] = pData;
			m_tm->Hold();
			if (vThr == 0)
			{
				switch (pOp)
				{
				case ReduceSum:
					for (unsigned int thr=1; thr < m_numworkers; ++thr)
					{
						pRes[0] += pRes[thr];
					}
					break;
				case ReduceMax:
					for (unsigned int thr=1; thr < m_numworkers; ++thr)
					{
						pRes[0] = std::max(pRes[0], pRes[thr]);
					}
					break;
				case ReduceMin:
					for (unsigned int thr=1; thr < m_numworkers; ++thr)
					{
						pRes[0] = std::min(pRes[0], pRes[thr]);
					}
					break;
				}
				m_comm->AllReduce(pRes[0], pOp);
			}
			m_tm->Hold();
			pData = pRes[0];
			m_tm->Hold();
		}

		void ThreadedComm::v_AllReduce(Array<OneD, NekDouble>& pData,
				enum ReduceOperator pOp)
		{
			GenericAllReduce(pData, pOp, m_ListOfThrSendDataDbl);
		}

		void ThreadedComm::v_AllReduce(Array<OneD, int>& pData, enum ReduceOperator pOp)
		{
			GenericAllReduce(pData, pOp, m_ListOfThrSendDataInt);
		}

		template <typename DataType>
        void ThreadedComm::GenericAllReduce(Array<OneD,DataType>& pData, enum ReduceOperator pOp,
        		std::vector<Array<OneD, DataType>*>& pRes)
        {
			unsigned int vThr = m_tm->GetWorkerNum();
			pRes[vThr] = &pData;
			m_tm->Hold();
			int vSize = pRes[0]->num_elements();
			ASSERTL1(pRes[vThr]->num_elements() == vSize, "Arrays of different length.");
			DataType* vZero = pRes[0]->data();

			int vNpp = vSize / m_numworkers;
			if (vNpp > 32) // completely pulled out of the air
			{
				// Each thread will work on a part of the array.
				int vRem = vSize - vNpp * m_numworkers;
				int vOffset = vThr * vNpp;
				if (vThr < vRem)
				{
					vOffset += vThr;
					++vNpp;
				}
				else
				{
					vOffset += vRem;
				}
				DoReduction(pOp, pRes, vOffset, vNpp);
			}
			else // job too small to bother parallelising
			{
				if (vThr == 0)
				{
					DoReduction(pOp, pRes, 0, vSize);
				}
			}

			m_tm->Hold();
			if (vThr == 0)
			{
				m_comm->AllReduce(*pRes[0], pOp);
			}
			m_tm->Hold();
			for (unsigned int i=0; i < vSize; ++i)
			{
				pData[i] = *(vZero + i);
			}
			m_tm->Hold();
        }

		template <typename DataType>
		void ThreadedComm::DoReduction(enum ReduceOperator pOp,
				std::vector<Array<OneD, DataType>*>& pRes, int pOffset, int pNpp)
		{
			DataType* vArrThr;
			DataType* vZero = pRes[0]->data() + pOffset;
			switch (pOp)
			{
			case ReduceSum:
				for (unsigned int thr=1; thr < m_numworkers; ++thr)
				{
					vArrThr = pRes[thr]->data() + pOffset;
					for (unsigned int i=0; i < pNpp; ++i)
					{
						*(vZero + i) += *(vArrThr + i);
					}
				}
				break;
			case ReduceMax:
				for (unsigned int thr=1; thr < m_numworkers; ++thr)
				{
					vArrThr = pRes[thr]->data() + pOffset;
					for (unsigned int i=0; i < pNpp; ++i)
					{
						*(vZero + i) = std::max(*(vZero + i), *(vArrThr + i));
					}
				}
				break;
			case ReduceMin:
				for (unsigned int thr=1; thr < m_numworkers; ++thr)
				{
					vArrThr = pRes[thr]->data() + pOffset;
					for (unsigned int i=0; i < pNpp; ++i)
					{
						*(vZero + i) = std::min(*(vZero + i), *(vArrThr + i));
					}
				}
				break;
			}
		}

        // Sigh.  This might be a small copy-n-paste in CommMpi but it's fairly
        // massive here.  The obvious solution is to make GenericReduce and so on
        // work on DataType* + a size argument, but that requires adding a similar
        // method to Comm itself.
        void ThreadedComm::v_AllReduce(std::vector<unsigned int>& pData,
                enum ReduceOperator pOp)

		{
			GenericAllReduceVector(pData, pOp, m_ListOfThrSendDataUIntVector);
		}

		template <typename DataType>
        void ThreadedComm::GenericAllReduceVector(std::vector<DataType>& pData,
                enum ReduceOperator pOp,
        		std::vector<std::vector<DataType>*>& pRes)
        {
			unsigned int vThr = m_tm->GetWorkerNum();
			pRes[vThr] = &pData;
			m_tm->Hold();
			int vSize = pRes[0]->size();
			ASSERTL1(pRes[vThr]->size() == vSize, "Arrays of different length.");
			DataType* vZero = &(*pRes[0])[0];

			int vNpp = vSize / m_numworkers;
			if (vNpp > 32) // completely pulled out of the air
			{
				// Each thread will work on a part of the array.
				int vRem = vSize - vNpp * m_numworkers;
				int vOffset = vThr * vNpp;
				if (vThr < vRem)
				{
					vOffset += vThr;
					++vNpp;
				}
				else
				{
					vOffset += vRem;
				}
				DoReductionVector(pOp, pRes, vOffset, vNpp);
			}
			else // job too small to bother parallelising
			{
				if (vThr == 0)
				{
					DoReductionVector(pOp, pRes, 0, vSize);
				}
			}

			m_tm->Hold();
			if (vThr == 0)
			{
				m_comm->AllReduce(*pRes[0], pOp);
			}
			m_tm->Hold();
			for (unsigned int i=0; i < vSize; ++i)
			{
				pData[i] = *(vZero + i);
			}
			m_tm->Hold();
        }

		template <typename DataType>
		void ThreadedComm::DoReductionVector(enum ReduceOperator pOp,
				std::vector<std::vector<DataType>*>& pRes, int pOffset, int pNpp)
		{
			DataType* vArrThr;
			DataType* vZero = &(*pRes[0])[0] + pOffset;
			switch (pOp)
			{
			case ReduceSum:
				for (unsigned int thr=1; thr < m_numworkers; ++thr)
				{
					vArrThr = &(*pRes[thr])[0] + pOffset;
					for (unsigned int i=0; i < pNpp; ++i)
					{
						*(vZero + i) += *(vArrThr + i);
					}
				}
				break;
			case ReduceMax:
				for (unsigned int thr=1; thr < m_numworkers; ++thr)
				{
					vArrThr = &(*pRes[thr])[0] + pOffset;
					for (unsigned int i=0; i < pNpp; ++i)
					{
						*(vZero + i) = std::max(*(vZero + i), *(vArrThr + i));
					}
				}
				break;
			case ReduceMin:
				for (unsigned int thr=1; thr < m_numworkers; ++thr)
				{
					vArrThr = &(*pRes[thr])[0] + pOffset;
					for (unsigned int i=0; i < pNpp; ++i)
					{
						*(vZero + i) = std::min(*(vZero + i), *(vArrThr + i));
					}
				}
				break;
			}
		}

		void ThreadedComm::v_AlltoAll(Array<OneD, int>& pSendData,
				Array<OneD, int>& pRecvData)
		{
			GenericAlltoAll(pSendData, pRecvData, m_tmpSendArrInt, m_tmpRecvArrInt,
					m_ListOfThrSendDataInt);
		}

		void ThreadedComm::v_AlltoAll(Array<OneD, NekDouble>& pSendData,
											Array<OneD, NekDouble>& pRecvData)
		{
			GenericAlltoAll(pSendData, pRecvData, m_tmpSendArrDbl, m_tmpRecvArrDbl,
					m_ListOfThrSendDataDbl);
		}

		template <typename DataType>
		void ThreadedComm::GenericAlltoAll(Array<OneD, DataType>& pSendData,
				Array<OneD, DataType>& pRecvData,
				Array<OneD, DataType>& pTmpSend,
				Array<OneD, DataType>& pTmpRecv,
				std::vector<Array<OneD, DataType>*>& pRes)
		{
			DataType * vTmp;
			unsigned int vThr = m_tm->GetWorkerNum();
			unsigned int vNumPerMPI = pSendData.num_elements() / m_numMPI;
			unsigned int vNumPerThr = pSendData.num_elements() / (m_numMPI * m_numworkers);
			int vRank = m_comm->GetRank();
			pRes[vThr] = &pSendData;
			if (vThr == 0)
			{
//				std::cerr << "vNumPerMPI: " << vNumPerMPI << " vNumPerThr: " << vNumPerThr << std::endl;
				pTmpSend = Array<OneD, DataType>(pSendData.num_elements() * m_numworkers *
						(m_numMPI-1) / m_numMPI);
				pTmpRecv = Array<OneD, DataType>(pSendData.num_elements() * m_numworkers *
						(m_numMPI-1) / m_numMPI);
			}
			m_tm->Hold();

			// Populate temp array to send to other MPI workers
			// with this thread's data.  Needs to be arranged carefully.
			if (m_numMPI > 1) {
				vTmp = pTmpSend.data();
				DataType* vStart = vTmp;
				vTmp += vThr * vNumPerMPI;
				DataType* vFrom = pSendData.data();
				for (int vMPI = 0; vMPI < m_numMPI; ++vMPI)
				{
					if (vThr == 0) m_SendOffsetArr[vMPI] = vTmp - vStart;
					if (vMPI == vRank)
					{
						vFrom += vNumPerMPI;
						m_SendSizeArr[vMPI] = 0;
					}
					else
					{
						for (unsigned int i=0; i < vNumPerMPI; ++i)
						{
							*vTmp = *vFrom;
							++vFrom;
							++vTmp;
						}
						vTmp += (m_numworkers - 1) * vNumPerMPI;
						m_SendSizeArr[vMPI] = m_numworkers * vNumPerMPI;
					}
				}
				m_tm->Hold();

				if (vThr == 0)
				{
//					std::cerr << "ready to send" << std::endl;
//					for (unsigned int i=0; i < pTmpSend.num_elements(); ++i)
//					{
//						std::cerr << "pTmpSend[" << i << "] = " << pTmpSend[i] << std::endl;
//					}
//
//					for (unsigned int i=0; i < m_numMPI; ++i)
//					{
//						std::cerr << "m_sizeArr[" << i << "] = " << m_SendSizeArr[i] << std::endl;
//						std::cerr << "m_offsetArr[" << i << "] = " << m_SendOffsetArr[i] << std::endl;
//					}

					m_comm->AlltoAllv(pTmpSend,
							  m_SendSizeArr,
							  m_SendOffsetArr,
							  pTmpRecv,
							  m_SendSizeArr,
							  m_SendOffsetArr);
//					std::cerr << "sent" << std::endl;

//					for (unsigned int i=0; i < pTmpRecv.num_elements(); ++i)
//					{
//						std::cerr << "pTmpRecv[" << i << "] = " << pTmpRecv[i] << std::endl;
//					}

				}
				m_tm->Hold();
			} // endif (m_numMPI > 1)

			// Now construct this thread's recv data
			// Some data from original pSendData
			vTmp = pTmpRecv.data();
			vTmp += vThr * vNumPerThr;
			DataType* vTo = pRecvData.data();
			for (int vMPI = 0; vMPI < m_numMPI; ++vMPI)
			{
				if (vMPI == vRank)
				{
					for (unsigned int thr=0; thr < m_numworkers; ++thr)
					{
						DataType* vOrigSend = pRes[thr]->data() +
								vThr * vNumPerThr +
								vMPI * vNumPerMPI;
						for (unsigned int i=0; i < vNumPerThr; ++i)
						{
							*vTo = *vOrigSend;
							++vTo;
							++vOrigSend;
						}
					}
				}
				else
				{
					for (unsigned int thr=0; thr < m_numworkers; ++thr)
					{
						for (unsigned int i=0; i < vNumPerThr; ++i)
						{
							*vTo = *vTmp;
							++vTo;
							++vTmp;
						}
						vTmp += (m_numworkers - 1) * vNumPerThr;
					}
				}
			}
			// Hold on exit because another thread might immediately re-enter and destroy common
			// data structures.
			m_tm->Hold();
		}

		void ThreadedComm::v_AlltoAllv(Array<OneD, NekDouble>& pSendData,
				Array<OneD, int>& pSendDataSizeMap,
				Array<OneD, int>& pSendDataOffsetMap, Array<OneD, NekDouble>& pRecvData,
				Array<OneD, int>& pRecvDataSizeMap,
				Array<OneD, int>& pRecvDataOffsetMap)
		{

			GenericAlltoAllv(pSendData,
					pSendDataSizeMap,
					pSendDataOffsetMap,
					pRecvData,
					pRecvDataSizeMap,
					pRecvDataOffsetMap,
					m_tmpSendArrDbl,
					m_tmpRecvArrDbl,
					m_ListOfThrSendDataDbl);
		}

		void ThreadedComm::v_AlltoAllv(Array<OneD, int>& pSendData,
				Array<OneD, int>& pSendDataSizeMap,
				Array<OneD, int>& pSendDataOffsetMap, Array<OneD, int>& pRecvData,
				Array<OneD, int>& pRecvDataSizeMap,
				Array<OneD, int>& pRecvDataOffsetMap)
		{
			GenericAlltoAllv(pSendData,
					pSendDataSizeMap,
					pSendDataOffsetMap,
					pRecvData,
					pRecvDataSizeMap,
					pRecvDataOffsetMap,
					m_tmpSendArrInt,
					m_tmpRecvArrInt,
					m_ListOfThrSendDataInt);
		}

		/**
		 * @param pTmpSend A reference to a particular class member (so that
		 * all threads can access it).  Used for actual MPI AlltoAllv.
		 * @param pTmpRecv A reference to a particular class member (so that
		 * all threads can access it).  Used for actual MPI AlltoAllv.
		 * @param pRes A list of pointers to all local threads' pSendData Arrays.
		 *
		 */
		template <typename DataType>
		void ThreadedComm::GenericAlltoAllv(Array<OneD, DataType>& pSendData,
				Array<OneD, int>& pSendDataSizeMap,
				Array<OneD, int>& pSendDataOffsetMap,
				Array<OneD, DataType>& pRecvData,
				Array<OneD, int>& pRecvDataSizeMap,
				Array<OneD, int>& pRecvDataOffsetMap,
				Array<OneD, DataType>& pTmpSend,
				Array<OneD, DataType>& pTmpRecv,
				std::vector<Array<OneD, DataType>*>& pRes)
		{
			DataType * vTmp;
			unsigned int vThr = m_tm->GetWorkerNum();
			//unsigned int vNumPerMPI = pSendData.num_elements() / m_numMPI;
			//unsigned int vNumPerThr = pSendData.num_elements() / (m_numMPI * m_numworkers);
			int vRank = m_comm->GetRank();
			pRes[vThr] = &pSendData;
			m_ListOfThrSendOffsets[vThr] = &pSendDataOffsetMap;
			m_ListOfThrSendSizes[vThr] = &pSendDataSizeMap;
			m_ListOfThrRecvOffsets[vThr] = &pRecvDataOffsetMap;
			m_ListOfThrRecvSizes[vThr] = &pRecvDataSizeMap;

			if (m_numMPI > 1)
			{
				m_tm->Hold();
				// ThreadedComm requires that there are more than one thread running.
				if (vThr == 0)
				{
				/**
				 * Construct a array of offsets for this process's threads into the temporary
				 * array we will construct later to send to MPI partners.  This is the array
				 * our MPI process will send.
				 * Also construct the MPI offsets for the MPI_AlltoAll.
				 */
					PopulateOffsets(vRank, m_tmpSendOffsetArr, m_SendOffsetArr,
							m_ListOfThrSendSizes, m_SendSizeArr, false);
					pTmpSend = Array<OneD, DataType>(
							m_tmpSendOffsetArr[m_numworkers * m_numworkers * m_numMPI]);
				}
				else if (vThr == 1)
				{
					/**
					 * This is for the array we will receive.
					 */
					PopulateOffsets(vRank, m_tmpRecvOffsetArr, m_RecvOffsetArr,
							m_ListOfThrRecvSizes, m_RecvSizeArr, true);
					pTmpRecv = Array<OneD, DataType>(
							m_tmpRecvOffsetArr[m_numworkers * m_numworkers * m_numMPI]);
				}
				m_tm->Hold();

				// Populate temp array to send to other MPI workers
				// with this thread's data.  Needs to be arranged carefully.
				vTmp = pTmpSend.data();
				DataType* vFrom = pSendData.data();
				for (int vMPI = 0; vMPI < m_numMPI; ++vMPI)
				{
					int vMPIIndex = vMPI * m_numworkers * m_numworkers +
							vThr * m_numworkers;

					if (vMPI == vRank)
					{
						//
					}
					else
					{
						for (unsigned int thr = 0; thr < m_numworkers; ++thr)
						{
							DataType* vTmpOffset = vTmp + m_tmpSendOffsetArr[vMPIIndex + thr];
							int vPartition = vMPI * m_numworkers + thr;
							DataType* vFromOffset = vFrom + pSendDataOffsetMap[vPartition];
//							std::cerr << "vMPI: " << vMPI << " vThr: " << vThr << " thr: " << thr <<
//									" vTmpOffset: " << vTmpOffset-vTmp << " vFromOffset: " <<
//									vFromOffset-vFrom << " size: " << pSendDataSizeMap[vPartition] <<
//									std::endl;
							for (int i=0; i < pSendDataSizeMap[vPartition]; ++i)
							{
								*vTmpOffset = *vFromOffset;
								++vFromOffset;
								++vTmpOffset;
							}
						}
					}
				}
				m_tm->Hold();

				/**
				 * Now make the MPI call.  Only one thread per MPI process does this.
				 */
				if (vThr == 0)
				{
//					std::cerr << "ready to send. pTmpSend.num_elements(): " << pTmpSend.num_elements() << std::endl;
//					for (unsigned int i=0; i < pTmpSend.num_elements(); ++i)
//					{
//						std::cerr << "pTmpSend[" << i << "] = " << pTmpSend[i] << std::endl;
//					}

//					for (unsigned int i=0; i < m_numMPI; ++i)
//					{
//						std::cerr << "m_sSendSizeArr[" << i << "] = " << m_SendSizeArr[i] << std::endl;
//						std::cerr << "m_SendOffsetArr[" << i << "] = " << m_SendOffsetArr[i] << std::endl;
//					}

					m_comm->AlltoAllv(pTmpSend,
							  m_SendSizeArr,
							  m_SendOffsetArr,
							  pTmpRecv,
							  m_RecvSizeArr,
							  m_RecvOffsetArr);
//					std::cerr << "sent" << std::endl;

//					for (unsigned int i=0; i < pTmpRecv.num_elements(); ++i)
//					{
//						std::cerr << "pTmpRecv[" << i << "] = " << pTmpRecv[i] << std::endl;
//					}

				}
			} // endif (m_numMPI > 1)

			m_tm->Hold();

			// Now construct this thread's recv data
			// Some data from original pSendData
			vTmp = pTmpRecv.data();
			DataType* vTo = pRecvData.data();

			for (int vMPI = 0; vMPI < m_numMPI; ++vMPI)
			{
				if (vMPI == vRank)
				{
					/**
					 * This part of the data comes from sister threads on this MPI process.
					 */
					for (unsigned int thr=0; thr < m_numworkers; ++thr)
					{
						int vRecvPartition = m_tm->GetPartitionFromRankThr(vMPI, thr);
						int vSendPartition = m_tm->GetPartitionFromRankThr(vMPI, vThr);
						Array<OneD, int>& vOffArr = *m_ListOfThrSendOffsets[thr];
						Array<OneD, int>& vSizArr = *m_ListOfThrSendSizes[thr];

						ASSERTL1(pRecvDataSizeMap[vRecvPartition] == vSizArr[vSendPartition],
								"Mismatch in sizes");
//						std::cerr << "vThr: " << vThr << " vRecvPartition: " << vRecvPartition <<
//								" vSendPartition: " << vSendPartition <<
//								" sendsiz: " << vSizArr[vSendPartition] << " recvsiz: " <<
//								pRecvDataSizeMap[vRecvPartition] << std::endl;

						DataType* vOrigSend = pRes[thr]->data() + vOffArr[vSendPartition];
						DataType* vToOffset = vTo + pRecvDataOffsetMap[vRecvPartition];
						for (unsigned int i=0; i < pRecvDataSizeMap[vRecvPartition]; ++i)
						{
							*vToOffset = *vOrigSend;
							++vToOffset;
							++vOrigSend;
						}
					}
				}
				else
				{
					/**
					 * This part of the data comes from the MPI receive array.
					 * @note We do not check to see if the sending thread (on another
					 * MPI process) has sent the same size of data as we are receiving.
					 */
					for (unsigned int thr=0; thr < m_numworkers; ++thr)
					{
						int vMPIIndex = vMPI * m_numworkers * m_numworkers +
								thr * m_numworkers;

						int vRecvPartition = vMPI * m_numworkers + thr;
						DataType* vTmpOffset = vTmp + m_tmpRecvOffsetArr[vMPIIndex + vThr];
						DataType* vToOffset = vTo + pRecvDataOffsetMap[vRecvPartition];
//						std::cerr << "vMPIIndex: " << vMPIIndex << " vTmpOffset: " << vTmpOffset - vTmp <<
//								" vToOffset: " << vToOffset - vTo << std::endl;
						for (unsigned int i=0; i < pRecvDataSizeMap[vRecvPartition]; ++i)
						{
							*vToOffset = *vTmpOffset;
							++vToOffset;
							++vTmpOffset;
						}
					}
				}
			}
			// Hold on exit because another thread might immediately re-enter and destroy common
			// data structures.
			m_tm->Hold();

		}

		void ThreadedComm::PopulateOffsets(int pRank, std::vector<int>& pTmpOffsetArr, Array<OneD, int>& pOffsetArr,
				std::vector<Array<OneD, int>* >& pSizeThrArr, Array<OneD, int>& pSizeArr, bool isRecv)
		{
			unsigned int vIndex = 0;
			pTmpOffsetArr[0] = 0;

			for (int vMPI = 0; vMPI < m_numMPI; ++vMPI)
			{
				pOffsetArr[vMPI] = pTmpOffsetArr[vIndex];
				int vSizeSum = 0;
				if (vMPI == pRank)
				{
					for (int lim=vIndex + m_numworkers * m_numworkers; vIndex < lim; )
					{
						++vIndex;
						pTmpOffsetArr[vIndex] = pTmpOffsetArr[vIndex - 1];
//						std::cerr << "pRank: " << pRank << "vIndex: " <<  vIndex <<
//								" pTmpOffsetArr: " << pTmpOffsetArr[vIndex] << std::endl;
					}
				}
				else
				{
					if (!isRecv)
					{
						for (unsigned int ithr = 0; ithr < m_numworkers; ++ithr)
						{
							Array<OneD, int>& vSizeArr = *pSizeThrArr[ithr];
							for (unsigned int thr = 0; thr < m_numworkers; ++thr)
							{
								int vPartition = vMPI * m_numworkers + thr;
								++vIndex;
								pTmpOffsetArr[vIndex] = pTmpOffsetArr[vIndex - 1] +
										vSizeArr[vPartition];
								vSizeSum += vSizeArr[vPartition];
//								std::cerr << "pRank: " << pRank << "vIndex: " <<  vIndex <<
//										" pTmpOffsetArr: " << pTmpOffsetArr[vIndex] <<
//										" " << isRecv << std::endl;
							}
						}
					}
					else
					{
						for (unsigned int ithr = 0; ithr < m_numworkers; ++ithr)
						{
							int vPartition = vMPI * m_numworkers + ithr;
							for (unsigned int thr = 0; thr < m_numworkers; ++thr)
							{
								Array<OneD, int>& vSizeArr = *pSizeThrArr[thr];
								++vIndex;
								pTmpOffsetArr[vIndex] = pTmpOffsetArr[vIndex - 1] +
										vSizeArr[vPartition];
								vSizeSum += vSizeArr[vPartition];
//								std::cerr << "pRank: " << pRank << "vIndex: " <<  vIndex <<
//										" pTmpOffsetArr: " << pTmpOffsetArr[vIndex] <<
//										" " << isRecv << std::endl;
							}
						}
					}
				}
				pSizeArr[vMPI] = vSizeSum;
			}

		}

		void ThreadedComm::v_SplitComm(int pRows, int pColumns)
		{
			if (m_tm->GetWorkerNum() == 0) {
				m_comm->SplitComm(pRows, pColumns);
				//CommSharedPtr vTmpColumn(new ThreadedComm(m_comm->GetColumnComm(), m_tm));
				CommSharedPtr vTmpColumn = m_comm->GetColumnComm();
				CommSharedPtr vTmpRow(new ThreadedComm(m_comm->GetRowComm(), m_tm));
				m_commColumn = vTmpColumn;
				m_commRow= vTmpRow;
			}
		}

		CommSharedPtr ThreadedComm::v_GetTrueComm()
		{
			return m_comm;
		}

        Gs::gs_data* ThreadedComm::v_GsInit(const Array<OneD, long> pId)
        {
			unsigned int vThr = m_tm->GetWorkerNum();
			Pack(vThr, m_ListOfConstThrSendDataLon, pId, m_tmpSendArrLon);
			m_tm->Hold();
			if (vThr == 0)
			{
				return m_comm->GsInit(m_tmpSendArrLon);
			}
			else
			{
	        	return 0;
			}
        }

        void ThreadedComm::v_GsFinalise(Gs::gs_data *pGsh)
        {
        	if (m_tm->GetWorkerNum() == 0)
        	{
        		m_comm->GsFinalise(pGsh);
        	}
        }

        void ThreadedComm::v_GsUnique(Array<OneD, long> pId)
        {
			unsigned int vThr = m_tm->GetWorkerNum();
			Pack(vThr, m_ListOfConstThrSendDataLon, pId, m_tmpSendArrLon);
			m_tm->Hold();

			if (vThr == 0)
			{
//				std::cerr << "m_tmpSendArrLon before" << std::endl;
//				for (int i=0; i < m_tmpSendArrLon.num_elements(); ++i)
//				{
//					std::cerr << "m_tmpSendArrLon[" << i << "] = " << m_tmpSendArrLon[i] << std::endl;
//				}
				m_comm->GsUnique(m_tmpSendArrLon);
//				std::cerr << "m_tmpSendArrLon before" << std::endl;
//				for (int i=0; i < m_tmpSendArrLon.num_elements(); ++i)
//				{
//					std::cerr << "m_tmpSendArrLon[" << i << "] = " << m_tmpSendArrLon[i] << std::endl;
//				}
			}
			m_tm->Hold();
			UnPack(vThr, m_ListOfConstThrSendDataLon, pId, m_tmpSendArrLon);
			m_tm->Hold();
        }

        void ThreadedComm::v_GsGather(Array<OneD, NekDouble> pU, Gs::gs_op pOp,
                Gs::gs_data *pGsh)
        {
        	unsigned int vThr = m_tm->GetWorkerNum();
        	Pack(vThr, m_ListOfConstThrSendDataDbl, pU, m_tmpSendArrDbl);
			m_tm->Hold();

			if (vThr == 0)
			{
//				std::cerr << "m_tmpSendArrDbl before" << std::endl;
//				for (int i=0; i < m_tmpSendArrDbl.num_elements(); ++i)
//				{
//					std::cerr << "m_tmpSendArrDbl[" << i << "] = " << m_tmpSendArrDbl[i] << std::endl;
//				}
				m_comm->GsGather(m_tmpSendArrDbl, pOp, pGsh);
//				std::cerr << "m_tmpSendArrDbl before" << std::endl;
//				for (int i=0; i < m_tmpSendArrDbl.num_elements(); ++i)
//				{
//					std::cerr << "m_tmpSendArrDbl[" << i << "] = " << m_tmpSendArrDbl[i] << std::endl;
//				}
			}
			m_tm->Hold();
			UnPack(vThr, m_ListOfConstThrSendDataDbl, pU, m_tmpSendArrDbl);
			m_tm->Hold();
        }

        void ThreadedComm::v_SendToRankZero(std::vector<unsigned int>&  pSendData)
        {
            int vRank = m_comm->GetRank();
            unsigned int vThr = m_tm->GetWorkerNum();
            ASSERTL0(vThr != 0 || vRank != 0, "Tried to SendToRankZero"
                " when rank zero.");
            m_sendToRankZeroVecInt[vThr] = &pSendData;
            m_tm->Hold();
            if (vThr == 0) // not rank 0 from above ASSERTL0
            {
                for (unsigned int i=0; i < m_numworkers; ++i)
                {
                    m_comm->Send(0, *m_sendToRankZeroVecInt[i]);
                }
            }
            m_tm->Hold();
        }

        void ThreadedComm::v_RecvFromAll(std::vector<std::vector<unsigned int> >&  pRecvData)
        {
            int vRank = m_comm->GetRank();
            unsigned int vThr = m_tm->GetWorkerNum();
            ASSERTL0(vThr == 0 && vRank == 0, "Tried to RecvFromAll"
                " when not rank zero.");
            m_tm->Hold();
            unsigned int j=1;
            for (unsigned int i=1; i < m_numworkers; ++i,++j)
            {
                pRecvData[j] = *m_sendToRankZeroVecInt[i];
            }
            m_tm->Hold();
            for (int r=1; r < m_numMPI; ++r)
            {
                for (unsigned int i=0; i < m_numworkers; ++i,++j)
                {
                    m_comm->Recv(r, pRecvData[j]);
                }
            }
        }

        void ThreadedComm::v_SendToRankZero(Array<OneD, int>&  pSendData)
        {
            int vRank = m_comm->GetRank();
            unsigned int vThr = m_tm->GetWorkerNum();
            ASSERTL0(vThr != 0 || vRank != 0, "Tried to SendToRankZero"
                " when rank zero.");
            m_sendToRankZeroArrInt[vThr] = &pSendData;
            m_tm->Hold();
            if (vThr == 0) // not rank 0 from above ASSERTL0
            {
                for (unsigned int i=0; i < m_numworkers; ++i)
                {
                    m_comm->Send(0, *m_sendToRankZeroArrInt[i]);
                }
            }
            m_tm->Hold();
        }

        void ThreadedComm::v_RecvFromAll(std::vector<Array<OneD, int> >&  pRecvData)
        {
            int vRank = m_comm->GetRank();
            unsigned int vThr = m_tm->GetWorkerNum();
            ASSERTL0(vThr == 0 && vRank == 0, "Tried to RecvFromAll"
                " when not rank zero.");
            m_tm->Hold();
            unsigned int j=1;
            for (unsigned int i=1; i < m_numworkers; ++i,++j)
            {
                pRecvData[j] = *m_sendToRankZeroArrInt[i];
            }
            m_tm->Hold();
            for (int r=1; r < m_numMPI; ++r)
            {
                for (unsigned int i=0; i < m_numworkers; ++i,++j)
                {
                    m_comm->Recv(r, pRecvData[j]);
                }
            }
        }

        void ThreadedComm::v_SendToRankZero(Array<OneD, NekDouble>&  pSendData)
        {
            int vRank = m_comm->GetRank();
            unsigned int vThr = m_tm->GetWorkerNum();
            ASSERTL0(vThr != 0 || vRank != 0, "Tried to SendToRankZero"
                " when rank zero.");
            m_sendToRankZeroArrDbl[vThr] = &pSendData;
            m_tm->Hold();
            if (vThr == 0) // not rank 0 from above ASSERTL0
            {
                for (unsigned int i=0; i < m_numworkers; ++i)
                {
                    m_comm->Send(0, *m_sendToRankZeroArrDbl[i]);
                }
            }
            m_tm->Hold();
        }

        void ThreadedComm::v_RecvFromAll(std::vector<Array<OneD, NekDouble> >&  pRecvData)
        {
            int vRank = m_comm->GetRank();
            unsigned int vThr = m_tm->GetWorkerNum();
            ASSERTL0(vThr == 0 && vRank == 0, "Tried to RecvFromAll"
                " when not rank zero.");
            m_tm->Hold();
            unsigned int j=1;
            for (unsigned int i=1; i < m_numworkers; ++i,++j)
            {
                pRecvData[j] = *m_sendToRankZeroArrDbl[i];
            }
            m_tm->Hold();
            for (int r=1; r < m_numMPI; ++r)
            {
                for (unsigned int i=0; i < m_numworkers; ++i,++j)
                {
                    m_comm->Recv(r, pRecvData[j]);
                }
            }
        }

        template<class DataType>
        void ThreadedComm::Pack(unsigned int pThr,
        		std::vector<Array<OneD, DataType> const *> &pRes, const Array<OneD, DataType> &pIn,
        		Array<OneD, DataType> &pPck)
        {
        	int vOff = 0;
        	const DataType* vInOff = pIn.data();
        	int vSize = 0;
        	int vNumElements = pIn.num_elements();
			pRes[pThr] = &pIn;
			m_tm->Hold();

			for (int i=0; i < m_numworkers; ++i)
			{
				int vNum = pRes[i]->num_elements();
				if (i < pThr) vOff += vNum;
				vSize += vNum;
			}

			if (pThr == 0)
			{
				pPck = Array<OneD, DataType>(vSize);
			}
			m_tm->Hold();

        	DataType* vPckOff = pPck.data() + vOff;
			for (int i=0; i < vNumElements; ++i)
			{
				*vPckOff = *vInOff;
				++vPckOff;
				++vInOff;
			}
        }

        template<class DataType>
        void ThreadedComm::UnPack(unsigned int pThr,
        		std::vector<Array<OneD, DataType> const*> &pRes, Array<OneD, DataType> &pOut,
        		Array<OneD, DataType> &pPck)
        {
        	DataType* vPckOff = pPck.data();
        	DataType* vOutOff = pOut.data();
        	int vNumElements = pOut.num_elements();
			pRes[pThr] = &pOut;
			m_tm->Hold();

			for (int i=1; i <= pThr; ++i)
			{
				vPckOff += pRes[i-1]->num_elements();
			}
			for (int i=0; i < vNumElements; ++i)
			{
				*vOutOff = *vPckOff;
				++vPckOff;
				++vOutOff;
			}
        }

	}
} /* namespace Nektar */
