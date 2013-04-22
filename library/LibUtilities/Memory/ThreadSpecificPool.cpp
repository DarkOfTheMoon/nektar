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
//	namespace detail
//	{
	static Nektar::Thread::ThreadManagerSharedPtr sThrMan;
	typedef boost::shared_ptr<std::vector<MemPool *> > MemoryPoolPool;
	static bool s_threadPoolsEnabled = false;
//	}

	boost::shared_ptr<std::vector<MemPool *> >& GetMemoryPoolPool()
	{
		typedef Loki::SingletonHolder<MemoryPoolPool ,
				Loki::CreateUsingNew,
				Loki::NoDestroy > Type;
		MemoryPoolPool &p = Type::Instance();
		if (!p)
		{
			p = boost::shared_ptr<std::vector<MemPool *> >(new std::vector<MemPool *>(1,static_cast<MemPool*>(0)));
		}
		return p;
	}
    MemPool& GetMemoryPool()
    {
//        typedef Loki::SingletonHolder<MemPool ,
//                Loki::CreateUsingNew,
//                Loki::NoDestroy > Type;
//        return Type::Instance();

    	/*
    	 * Now have a single MemPool per thread (because MemPool is *not* thread safe).
    	 * We want to avoid locking in here because this function is used heavily.
    	 * Even shared locking causes a large slowdown.
    	 */
    	Nektar::Thread::ThreadManagerSharedPtr vThrMan = Nektar::Thread::ThreadManager::GetInstance();
    	unsigned int vThr = vThrMan ? vThrMan->GetWorkerNum() : 0;
		if (!sThrMan && vThrMan)
		{
			sThrMan = vThrMan;
		}
    	static MemoryPoolPool &p = GetMemoryPoolPool();
    	if (!s_threadPoolsEnabled)
    	{
    		ASSERTL1(vThr == 0, "Using threaded memory pool before it's inited");
        	ASSERTL1(p, "Done goofed");
    		if ((*p)[0] == 0)
    		{
    			(*p)[0] = new MemPool();
    		}
    	}
    	return *((*p)[vThr]);
    }

    /**
     * @brief Sets up memory pools for all threads.
     *
     * @note This must be called by a single thread.
     */
    void InitMemoryPools(unsigned int pNumThr, bool pEnabled)
    {
    	Nektar::Thread::ThreadManagerSharedPtr vThrMan = Nektar::Thread::ThreadManager::GetInstance();
    	MemoryPoolPool &p = GetMemoryPoolPool();

    	p->resize(pNumThr);
    	for (unsigned int i=0; i < pNumThr; ++i)
    	{
    		if ((*p)[i] == 0) // for thread 0
    		{
    			(*p)[i] = new MemPool();
    		}
    	}
    	s_threadPoolsEnabled = pEnabled;
    }
}
