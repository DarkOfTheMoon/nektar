/*
 * ThreadedComm.h
 *
 *  Created on: 11 Mar 2013
 *      Author: simon
 */

#ifndef THREADEDCOMM_H_
#define THREADEDCOMM_H_

#include <LibUtilities/Communication/Comm.h>
#include <LibUtilities/BasicUtils/Thread.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>


namespace Nektar
{
	namespace LibUtilities
	{
		using namespace Thread;
		class ThreadedComm: public Comm
		{
		public:
			ThreadedComm(CommSharedPtr comm, ThreadManagerSharedPtr tm);
			virtual ~ThreadedComm();
		private:
			CommSharedPtr 							m_comm;
			ThreadManagerSharedPtr					m_tm;
			unsigned int							m_numworkers;
			int										m_numMPI;
			std::vector<NekDouble> 					m_resDbl;
			std::vector<int>	 					m_resInt;
			std::vector<Array<OneD, NekDouble>*>	m_ListOfThrSendDataDbl;
			std::vector<Array<OneD, int>*>			m_ListOfThrSendDataInt;
			Array<OneD, NekDouble>					m_tmpSendArrDbl;
			Array<OneD, NekDouble>					m_tmpRecvArrDbl;
			Array<OneD, int>						m_tmpSendArrInt;
			Array<OneD, int>						m_tmpRecvArrInt;
			Array<OneD, int>						m_tmpArrInt;
			/**
			 * Holds offsets for MPI AlltoAllv call.
			 * These offsets are for MPI boundaries.
			 */
			Array<OneD, int>						m_SendOffsetArr;
			/**
			 * Holds sizes for MPI AlltoAllv call.
			 * These sizes are for MPI boundaries.
			 */
			Array<OneD, int>						m_SendSizeArr;
			/**
			 * Holds offsets for MPI AlltoAllv call.
			 * These offsets are for MPI boundaries.
			 */
			Array<OneD, int>						m_RecvOffsetArr;
			/**
			 * Holds sizes for MPI AlltoAllv call.
			 * These sizes are for MPI boundaries.
			 */
			Array<OneD, int>						m_RecvSizeArr;
			/**
			 * Holds Arrays of offset maps for send data for all threads.
			 */
			std::vector<Array<OneD, int>*>			m_ListOfThrSendOffsets;
			/**
			 * Holds Arrays of sizes for send data for all threads.
			 */
			std::vector<Array<OneD, int>*>			m_ListOfThrSendSizes;
			/**
			 * Holds Arrays of offset maps for recv data for all threads.
			 */
			std::vector<Array<OneD, int>*>			m_ListOfThrRecvOffsets;
			/**
			 * Holds Arrays of sizes for recv data for all threads.
			 */
			std::vector<Array<OneD, int>*>			m_ListOfThrRecvSizes;
			/**
			 * Holds offsets for temporary send Array to help with
			 * its construction by master thread to send over MPI.
			 * Offsets are for thread boundaries.
			 */
			std::vector<int>						m_tmpSendOffsetArr;
			/**
			 * Holds sizes for temporary send Array to help with
			 * its construction by master thread to send over MPI.
			 * Sizes are for thread boundaries.
			 */
			std::vector<int>						m_tmpSendSizeArr;
			/**
			 * Holds offsets for temporary recv Array to help with
			 * its deconstruction by threads after it's been sent over MPI.
			 * Offsets are for thread boundaries.
			 *
			 */
			std::vector<int>						m_tmpRecvOffsetArr;

        protected:
            virtual void v_Finalise();
            virtual int  v_GetRank();
            virtual void v_Block();
            virtual void v_Send(int pProc, Array<OneD, NekDouble>& pData);
            virtual void v_Send(int pProc, Array<OneD, int>& pData);
            virtual void v_Recv(int pProc, Array<OneD, NekDouble>& pData);
            virtual void v_Recv(int pProc, Array<OneD, int>& pData);
            virtual void v_SendRecv(int pSendProc,
                                    Array<OneD, NekDouble>& pSendData,
                                    int pRecvProc,
                                    Array<OneD, NekDouble>& pRecvData);
            virtual void v_SendRecv(int pSendProc,
                                    Array<OneD, int>& pSendData,
                                    int pRecvProc,
                                    Array<OneD, int>& pRecvData);
			virtual void v_SendRecvReplace(int pSendProc,
										   int pRecvProc,
										   Array<OneD, NekDouble>& pSendData);
			virtual void v_SendRecvReplace(int pSendProc,
										   int pRecvProc,
										   Array<OneD, int>& pSendData);
            virtual void v_AllReduce(NekDouble& pData,
                                     enum ReduceOperator pOp);
            virtual void v_AllReduce(int& pData,
                                     enum ReduceOperator pOp);
            virtual void v_AllReduce(Array<OneD, NekDouble>& pData,
                                     enum ReduceOperator pOp);
            virtual void v_AllReduce(Array<OneD, int      >& pData,
                                     enum ReduceOperator pOp);
			virtual void v_AlltoAll(Array<OneD, NekDouble>& pSendData,
									Array<OneD, NekDouble>& pRecvData);
            virtual void v_AlltoAll(Array<OneD, int>& pSendData,
									Array<OneD, int>& pRecvData);
			virtual void v_AlltoAllv(Array<OneD, NekDouble>& pSendData,
									Array<OneD, int>& pSendDataSizeMap,
									Array<OneD, int>& pSendDataOffsetMap,
									Array<OneD, NekDouble>& pRecvData,
									Array<OneD, int>& pRecvDataSizeMap,
									Array<OneD, int>& pRecvDataOffsetMap);
			virtual void v_AlltoAllv(Array<OneD, int>& pSendData,
									Array<OneD, int>& pSendDataSizeMap,
									Array<OneD, int>& pSendDataOffsetMap,
									Array<OneD, int>& pRecvData,
									Array<OneD, int>& pRecvDataSizeMap,
									Array<OneD, int>& pRecvDataOffsetMap);
            virtual void v_SplitComm(int pRows, int pColumns);
            virtual CommSharedPtr v_GetTrueComm();
            virtual Gs::gs_data* v_GsInit(const Array<OneD, long> pId);
            virtual void v_GsFinalise(Gs::gs_data *pGsh);
            virtual void v_GsUnique(const Array<OneD, long> pId);
            virtual void v_GsGather(Array<OneD, NekDouble> pU, Gs::gs_op pOp,
                    Gs::gs_data *pGsh, Array<OneD, NekDouble> pBuffer
                                                     = NullNekDouble1DArray);

        private:
            template <typename DataType>
            void GenericAllReduce(DataType& pData, enum ReduceOperator pOp,
            		std::vector<DataType>& pRes);
            template <typename DataType>
            void GenericAllReduce(Array<OneD,DataType>& pData, enum ReduceOperator pOp,
            		std::vector<Array<OneD, DataType>*>& pRes);
            template <typename DataType>
            void DoReduction(enum ReduceOperator pOp, std::vector<Array<OneD, DataType>*>& pRes, int pOffset,
            		int pNpp);
    		template <typename DataType>
    		void GenericAlltoAll(Array<OneD, DataType>& pSendData,
    				Array<OneD, DataType>& pRecvData,
    				Array<OneD, DataType>& pTmpSend,
    				Array<OneD, DataType>& pTmpRecv,
    				std::vector<Array<OneD, DataType>*>& pRes);
    		template <typename DataType>
    		void GenericAlltoAllv(Array<OneD, DataType>& pSendData,
					Array<OneD, int>& pSendDataSizeMap,
					Array<OneD, int>& pSendDataOffsetMap,
    				Array<OneD, DataType>& pRecvData,
					Array<OneD, int>& pRecvDataSizeMap,
					Array<OneD, int>& pRecvDataOffsetMap,
    				Array<OneD, DataType>& pTmpSend,
    				Array<OneD, DataType>& pTmpRecv,
    				std::vector<Array<OneD, DataType>*>& pRes);
    		void PopulateOffsets(int pRank, std::vector<int>& pTmpOffsetArr, Array<OneD, int>& pOffsetArr,
    				std::vector<Array<OneD, int>* >& pSizeThrArr, Array<OneD, int>& pSizeArr, bool isRecv);

		};

	} /* namespace Thread */
} /* namespace Nektar */
#endif /* THREADEDCOMM_H_ */
