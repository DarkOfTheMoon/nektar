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





BasisSharedPtr<double> CreateBasisObject(const BasisKey& bkey)
{
    return GetBasisFactory().CreateInstance(bkey.m_id, bkey);
}


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

    // Test some points keys
    Points<double, Segment,       GaussGaussLegendre> P(key);
    Points<double, Quadrilateral, std::tuple<GaussGaussLegendre, GaussGaussLegendre>> Pquad(key);
    Points<double, Hexahedron,    std::tuple<GaussGaussLegendre, GaussGaussLegendre, GaussGaussLegendre>> Phex(key);

    // These should fail
    //Points<double, Segment,       Fekete> Pseg_fek(key);
    //Points<double, Triangle,      GaussGaussLegendre> Ptri_gauss(key);

    // Test extraction of constituent points distributions
    Array<OneD, double> p1 = Pquad.GetZ(0);
    Array<OneD, double> p2 = Pquad.GetZ(1);
    cout << p1.num_elements() << endl;
    cout << p2.num_elements() << endl;

    Array<OneD, double> h1 = Phex.GetZ(0);
    Array<OneD, double> h2 = Phex.GetZ(1);
    Array<OneD, double> h3 = Phex.GetZ(2);
    cout << h1.num_elements() << endl;
    cout << h2.num_elements() << endl;
    cout << h3.num_elements() << endl;

    // Test access through base class pointer
    PointsBase<double>* ptr1 = new Points<double, Segment,       GaussGaussLegendre>(key);
    PointsBase<double>* ptr2 = new Points<double, Quadrilateral, std::tuple<GaussGaussLegendre, GaussGaussLegendre>>(key);
    cout << ptr1->GetZ(0).num_elements() << endl;
    cout << ptr2->GetZ(0).num_elements() << endl;
    cout << ptr2->GetZ(1).num_elements() << endl;
    try {
        cout << ptr2->GetZ(2).num_elements() << endl;
    }
    catch (...) {
        cout << "Could not access third points distribution in 2D tensor product" << endl;
    }

    // Test instantiation of a Basis
    Basis<double, Segment, GaussGaussLegendre, ModifiedLegendre> B0(bkey);
    Basis<double, Quadrilateral, std::tuple<GaussGaussLegendre, GaussGaussLegendre>, std::tuple<ModifiedLegendre, ModifiedLegendre>> Bquad(bkey);

    // This should fail
    //Basis<double, Triangle, GaussGaussLegendre, BernsteinTriangle> Bbtri(bkey);

//    cout << sizeof(Array<OneD, double>) << endl;
    cout << "Size of points key: " << sizeof(PointsKey) << endl;
    cout << "Size of GGL points: " << sizeof(Points<double, Segment, GaussGaussLegendre>) << endl;
    cout << "Size of GGL quad pts: " << sizeof(Points<double, Quadrilateral, std::tuple<GaussGaussLegendre, GaussGaussLegendre>>) << endl;
    cout << "Size of Modified basis: " << sizeof(Basis<double, Segment, GaussGaussLegendre, ModifiedLegendre>) << endl;

    // Create a points object using the factory
    cout << "Points factory creation." << endl;
    GetPointsFactory().PrintAvailableClasses(std::cout);
    PointsSharedPtr<double> facptr = GetPointsFactory().CreateInstance("GGL,Segment", key);
    cout << facptr->GetNumPoints() << endl;

    cout << "Basis creation" << endl;
    GetBasisFactory().PrintAvailableClasses(std::cout);
    BasisSharedPtr<double> bfacptr = GetBasisFactory().CreateInstance("MOD,GGL,Segment", bkey);
    BasisSharedPtr<double> bfacquadptr = GetBasisFactory().CreateInstance("MOD_MOD,GGL_GGL,Quadrilateral", bkey);

    cout << "Create using manager" << endl;
    LibUtilities::NekManager<BasisKey, BasisBase<double>, BasisKey::opLess> bmanager;
    bmanager.RegisterGlobalCreator(boost::bind(CreateBasisObject, _1));
    BasisSharedPtr<double> bmanptr = bmanager[bkey];

    // Create test element and perform operations
    cout << "Create test element" << endl;
    Array<OneD, NekDouble> testarray;
    TestQuadExp<NekDouble> quadexp(B0);
    cout << "Perform bwdtrans" << endl;
    quadexp.BwdTrans(testarray, testarray);

}
