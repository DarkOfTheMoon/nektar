///////////////////////////////////////////////////////////////////////////////
//
// File Comm.cpp
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
// Description: Base communication class
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Communication/Comm.h>
#include <loki/Singleton.h>             // for CreateUsingNew, NoDestroy, etc

namespace Nektar
{
    namespace LibUtilities
    {
        Comm::Comm(int narg, char* arg[])
        {

        }

        Comm::Comm()
        {

        }

        Comm::~Comm()
        {

        }

        CommSharedPtr Comm::v_GetTrueComm()
        {
        	return shared_from_this();
        }

        // Method to enforce that all existing files are removed
        bool Comm::v_RemoveExistingFiles(void)
        {
            return true;
        }

        CommFactory& GetCommFactory()
        {
            typedef Loki::SingletonHolder<CommFactory,
                Loki::CreateUsingNew,
                Loki::NoDestroy,
                Loki::SingleThreaded> Type;
            return Type::Instance();
        }

        void Comm::v_SendToRankZero(std::vector<unsigned int>&  pSendData)
        {
            ASSERTL0(GetRank() != 0,
                "Tried to SendToRankZero from rank zero.");
            v_Send(0, pSendData);
        }

        void Comm::v_SendToRankZero(Array<OneD, int>&  pSendData)
        {
            ASSERTL0(GetRank() != 0,
                "Tried to SendToRankZero from rank zero.");
            v_Send(0, pSendData);
        }

        void Comm::v_SendToRankZero(Array<OneD, NekDouble>&  pSendData)
        {
            ASSERTL0(GetRank() != 0,
                "Tried to SendToRankZero from rank zero.");
            v_Send(0, pSendData);
        }

        void Comm::v_RecvFromAll(std::vector<std::vector<unsigned int> >&  pRecvData)
        {
            ASSERTL0(GetRank() == 0,
                "Tried to RecvFromAll from not rank zero.");
            ASSERTL0(pRecvData.size() >= GetSize(),
                "Recv vector too small in RecvFromAll.");
            for (unsigned int i=1; i<GetSize(); ++i)
            {
                v_Recv(i, pRecvData[i]);
            }
        }

        void Comm::v_RecvFromAll(std::vector<Array<OneD, int> >&  pRecvData)
        {
            ASSERTL0(GetRank() == 0,
                "Tried to RecvFromAll from not rank zero.");
            ASSERTL0(pRecvData.size() >= GetSize(),
                "Recv vector too small in RecvFromAll.");
            for (unsigned int i=1; i<GetSize(); ++i)
            {
                v_Recv(i, pRecvData[i]);
            }
        }

        void Comm::v_RecvFromAll(std::vector<Array<OneD, NekDouble> >&  pRecvData)
        {
            ASSERTL0(GetRank() == 0,
                "Tried to RecvFromAll from not rank zero.");
            ASSERTL0(pRecvData.size() >= GetSize(),
                "Recv vector too small in RecvFromAll.");
            for (unsigned int i=1; i<GetSize(); ++i)
            {
                v_Recv(i, pRecvData[i]);
            }
        }
    }
}
