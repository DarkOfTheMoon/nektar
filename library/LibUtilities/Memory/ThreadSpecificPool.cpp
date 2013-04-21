////////////////////////////////////////////////////////////////////////////////
//
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
////////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Memory/ThreadSpecificPool.hpp>
#include <LibUtilities/BasicUtils/Thread.h>
#include <boost/thread/shared_mutex.hpp>
#include <boost/thread/locks.hpp>

namespace Nektar
{

    MemPool& GetMemoryPool()
    {
//        typedef Loki::SingletonHolder<MemPool ,
//                Loki::CreateUsingNew,
//                Loki::NoDestroy > Type;
//        return Type::Instance();

    	/*
    	 * Now have a single MemPool per thread (because MemPool is *not* thread safe).
    	 * Should not need locking here because std::map is supposed to be thread safe
    	 * when used this way.  Locking here slows things down, too.
    	 *
    	 * HOWEVER!  It seems std::map isn't behaving (or I'm misreading the docs).
    	 * Occasionally get crashes without the locks.
    	 */
    	static boost::shared_mutex mutex;
    	static std::map<unsigned int, MemPool *> s_threadPools;
    	static Nektar::Thread::ThreadManagerSharedPtr sThrMan; // stored to prevent ThreadManager from
    														   // going out of scope while it's still needed.
    	Nektar::Thread::ThreadManagerSharedPtr vThrMan = Nektar::Thread::ThreadManager::GetInstance();
    	unsigned int vThr = vThrMan ? vThrMan->GetWorkerNum() : 0;
		if (vThrMan)
		{
			sThrMan = vThrMan;
		}
		boost::shared_lock<boost::shared_mutex> RLock(mutex);
    	if (s_threadPools.count(vThr) == 0)
    	{
    		RLock.unlock();
    		boost::unique_lock<boost::shared_mutex> WLock(mutex);
    		s_threadPools[vThr] = new MemPool();
        	return *(s_threadPools[vThr]);
    	}
    	return *(s_threadPools[vThr]);
    }
}
