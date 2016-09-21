///////////////////////////////////////////////////////////////////////////////
//
// File: BasisRegistration.cpp
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
// Description: Populate the Basis factory
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Foundations/Basis.hpp>
#include <LibUtilities/Foundations/Basis/BasisGauss.hpp>
#include <LibUtilities/Foundations/Basis/BasisBernstein.hpp>

namespace Nektar
{
namespace LibUtilities
{
namespace Foundations
{

/**
 * @brief Get the singleton Basis factory.
 */
BasisFactory<NekDouble>& GetBasisFactory()
{
    typedef Loki::SingletonHolder<BasisFactory<NekDouble>,
                                  Loki::CreateUsingNew,
                                  Loki::NoDestroy> Type;
    // Putting Loki::ClassLevelLockable causes an assertion!
    return Type::Instance();
}

/**
 * Register available basis classes with factory
 */
static std::string bases[] = {
        GetBasisFactory().RegisterCreatorFunction("MOD,GGL,Segment", Basis<NekDouble, Segment, std::tuple<GaussGaussLegendre>, std::tuple<ModifiedLegendre>>::create,
                                                  "Modified basis with Gauss-Gauss-Legendre on Segment"),
        GetBasisFactory().RegisterCreatorFunction("MOD,GRM,Segment", Basis<NekDouble, Segment, std::tuple<GaussRadauMLegendre>, std::tuple<ModifiedLegendre>>::create,
                                                  "Modified basis with Gauss-RadauM-Legendre on Segment"),
        GetBasisFactory().RegisterCreatorFunction("MOD_MOD,GGL_GGL,Quadrilateral", Basis<NekDouble, Quadrilateral, std::tuple<GaussGaussLegendre, GaussGaussLegendre>, std::tuple<ModifiedLegendre, ModifiedLegendre>>::create,
                                                  "Modified basis with Gauss-Gauss-Legendre on Quadrilateral")
};

/**
 * @brief Get the singleton Basis factory.
 */
BasisManager<NekDouble>& GetBasisManager()
{
    typedef Loki::SingletonHolder<BasisManager<NekDouble>,
                                  Loki::CreateUsingNew,
                                  Loki::NoDestroy> Type;
    // Putting Loki::ClassLevelLockable causes an assertion!
    return Type::Instance();
}

BasisSharedPtr<NekDouble> BasisCreator(const BasisKey& pKey)
{
    return GetBasisFactory().CreateInstance(pKey.m_id, pKey);
}

static bool basismaninit = GetBasisManager().RegisterGlobalCreator(boost::bind(BasisCreator, _1));

}
}
}
