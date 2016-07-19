#ifndef LIBUTILITIES_FOUNDATIONS_INTERP
#define LIBUTILITIES_FOUNDATIONS_INTERP

#include <LibUtilities/Foundations/Points.hpp>

namespace Nektar
{
namespace LibUtilities
{
namespace Foundations
{

template<typename TData, typename TPts1, typename TPts2>
class Interp
{
    public:
        Interp(PointsKey p1, PointsKey p2) : t1(p1), t2(p2) {}

        boost::shared_ptr<NekMatrix<TData> > GetI()
        {
            Array<OneD, TData> interp
                = t1.GetInterpMatrix(t2.GetZ(), t2.GetNumPoints());
            TData* t = interp.data();
            unsigned int np1 = t1.GetNumPoints();
            unsigned int np2 = t2.GetNumPoints();
            boost::shared_ptr< NekMatrix<NekDouble> > returnval(MemoryManager<NekMatrix<NekDouble> >::AllocateSharedPtr(np1,np2,t));
        }

    private:
        TPts1 t1;
        TPts2 t2;
};

}
}
}

#endif
