#include <iostream>
#include <type_traits>
using namespace std;

#include <LibUtilities/BasicUtils/NekManager.hpp>
#include <LibUtilities/LinearAlgebra/NekTypeDefs.hpp>
#include <LibUtilities/Foundations/Foundations.hpp>
#include <LibUtilities/Foundations/Points.hpp>
#include <LibUtilities/Foundations/Basis.hpp>
#include <LibUtilities/Foundations/Basis/BasisTemplate.hpp>
#include <LibUtilities/Foundations/Basis/BasisGauss.hpp>
#include <LibUtilities/Foundations/Basis/BasisBernstein.hpp>
#include <LibUtilities/Foundations/Interp.hpp>
#include <LibUtilities/Foundations/Points/PointsTraits.hpp>

using namespace Nektar;
using namespace Nektar::LibUtilities::Foundations;

// VecMat
struct SingleVector {};
struct MultipleVector {};

// Impl
struct SumFac {};
struct LocalMatrix {};
struct IterPerExp {};
struct StdMat {};

// Op
struct BwdTransOp {};

// Multiplicity policies
class MultiplicityPolicySingle {
    protected:
        void dgemv(DNekMatSharedPtr pMat, Array<OneD, NekDouble>& pIn, Array<OneD, NekDouble>& pOut)
        {
            cout << "Act on single vector" << endl;
        }
};

class MultiplicityPolicyMultiple {
    protected:
        void dgemv(DNekMatSharedPtr pMat, Array<OneD, NekDouble>& pIn, Array<OneD, NekDouble>& pOut)
        {
            cout << "Act on multiple vectors" << endl;
//            Blas::Dgemm('N', 'N', pMat->GetRows(), pIn.num_elements() / pMat->GetColumns(),
//                    pMat->GetColumns(), 1.0, pMat->GetRawPtr(),
//                    pMat->GetRows(), pIn.get(), pMat->GetColumns(),
//                    0.0, pOut.get(), pMat->GetRows());
        }
};

// Operator policies
// These can be defined in any way as long as there is the createMatrix function
class OperatorPolicyBwdTrans {
    public:
        DNekMatSharedPtr createMatrix(BasisBase<NekDouble>& pBasis) {
            const PointsBase<NekDouble>& pts = pBasis.GetPoints();
            return MemoryManager<DNekMat>::AllocateSharedPtr(
                    pBasis.GetNumModes(),
                    pts.GetNumPoints(),
                    pBasis.GetBdata().get());
        }
};

class OperatorPolicyIProductWRTBase {
    public:
        DNekMatSharedPtr createMatrix(BasisBase<NekDouble>& pBasis) {

        }
};

template<typename, typename T>
struct has_createMatrix {
    static_assert(
        std::integral_constant<T, false>::value,
        "Second template parameter needs to be of function type.");
};

template<typename C, typename Ret, typename... Args>
struct has_createMatrix<C, Ret(Args...)> {
private:
    template<typename T>
    static constexpr auto check(T*)
    -> typename
        std::is_same<
            decltype( std::declval<T>().createMatrix( std::declval<Args>()... ) ),
            Ret
        >::type;  // attempt to call it and see if the return type is correct

    template<typename>
    static constexpr std::false_type check(...);

    typedef decltype(check<C>(0)) type;

public:
    static constexpr bool value = type::value;
};


template<typename TData>
class OperatorBase {

    public:
        OperatorBase(BasisBase<TData>& pBasis) : m_basis(pBasis) {}
        virtual ~OperatorBase() {}

        virtual void operator()(Array<OneD, TData>& in, Array<OneD, TData>& out) = 0;

    protected:
        BasisBase<TData>& m_basis;
};


/**
 * @brief Primary template for Operator class
 */
template<typename TData, typename TShape, typename TImpl, typename OpPolicy, typename MultPolicy>
class Operator;


/**
 * @brief LocalMatrix implementation (only valid on single elements)
 */
template<typename TData, typename TShape, typename OpPolicy>
class Operator<TData, TShape, LocalMatrix, OpPolicy, MultiplicityPolicySingle>
        : public OperatorBase<TData>, private OpPolicy, private MultiplicityPolicySingle
{
        using OperatorBase<TData>::m_basis;

        static_assert(has_createMatrix<OpPolicy, DNekMatSharedPtr(BasisBase<NekDouble>&)>::value, "No createMatrix");

    public:
        Operator(BasisBase<TData>& pBasis) : OperatorBase<TData>(pBasis) {}
        virtual ~Operator() {}

        virtual void operator()(Array<OneD, TData>& in, Array<OneD, TData>& out) {
            DNekMatSharedPtr mat = this->createMatrix(m_basis);
            this->dgemv(mat, in, out);
        }
};


/**
 * @brief SumFac implementation
 */
template<typename TData, typename OpPolicy, typename MultPolicy>
class Operator<TData, Quadrilateral, SumFac, OpPolicy, MultPolicy>
        : public OperatorBase<TData>, private OpPolicy, private MultPolicy
{
        using OperatorBase<TData>::m_basis;

        static_assert(has_createMatrix<OpPolicy, DNekMatSharedPtr(BasisBase<NekDouble>&)>::value, "No createMatrix");

    public:
        Operator(BasisBase<TData>& pBasis) : OperatorBase<TData>(pBasis) {}
        virtual ~Operator() {}

        virtual void operator()(Array<OneD, TData>& in, Array<OneD, TData>& out) {
            for (int i = 0; i < m_basis.GetNumConstituentBases(); ++i)
            {
                BasisBase<TData> b = m_basis.GetConstitutentBasis(i);
            }
//            DNekMatSharedPtr mat = this->createMatrix(m_basis);
//            this->dgemv(mat, in, out);
        }
};


/**
 * @brief IterPerExp implementation
 */
template<typename TData, typename OpPolicy, typename MultPolicy>
class Operator<TData, Quadrilateral, IterPerExp, OpPolicy, MultPolicy>
        : public OperatorBase<TData>, private OpPolicy, private MultPolicy
{
        using OperatorBase<TData>::m_basis;

        static_assert(has_createMatrix<OpPolicy, DNekMatSharedPtr(BasisBase<NekDouble>&)>::value, "No createMatrix");

    public:
        Operator(BasisBase<TData>& pBasis) : OperatorBase<TData>(pBasis) {}
        virtual ~Operator() {}

        virtual void operator()(Array<OneD, TData>& in, Array<OneD, TData>& out) {
//            DNekMatSharedPtr mat = this->createMatrix(m_basis);
//            this->dgemv(mat, in, out);
        }
};


template<typename TData>
class TestQuadExp {
    public:
        TestQuadExp(BasisBase<NekDouble>& pB0) {
            m_BwdTransOpImpl = new Operator<TData, Quadrilateral, LocalMatrix,
                                            OperatorPolicyBwdTrans,
                                            MultiplicityPolicySingle>(pB0);
            m_IProdWRTBaseOpImpl = new Operator<TData, Quadrilateral, LocalMatrix,
                                            OperatorPolicyIProductWRTBase,
                                            MultiplicityPolicySingle>(pB0);
        }

        void BwdTrans(Array<OneD, NekDouble>& in, Array<OneD, NekDouble>& out)
        {
            (*m_BwdTransOpImpl)(in, out);
        }

    private:
        OperatorBase<NekDouble>* m_BwdTransOpImpl;
        OperatorBase<NekDouble>* m_IProdWRTBaseOpImpl;
};



int main () {
    PointsKey key;
    key.m_numpoints[0] = 5;
    key.m_numpoints[1] = 6;
    key.m_numpoints[2] = 7;

    BasisKey bkey;
    bkey.m_id = "MOD,GGL,Segment";
    bkey.m_ptsKey = key;
    bkey.m_nummodes[0] = 4;
    bkey.m_nummodes[1] = 5;
    bkey.m_nummodes[2] = 6;

    // Test instantiation of a Basis
    Basis<double, Segment, GaussGaussLegendre, ModifiedLegendre> B0(bkey);

    // Create test element and perform operations
    cout << "Create test element" << endl;
    Array<OneD, NekDouble> testarray;
    TestQuadExp<NekDouble> quadexp(B0);
    cout << "Perform bwdtrans" << endl;
    quadexp.BwdTrans(testarray, testarray);

}
