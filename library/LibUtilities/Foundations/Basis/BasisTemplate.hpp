///////////////////////////////////////////////////////////////////////////////
//
// File: BasisTemplate.cpp
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
// Description: Primary template for Basis classes
//
///////////////////////////////////////////////////////////////////////////////

#ifndef LIBRARY_LIBUTILITIES_FOUNDATIONS_BASIS_BASISTEMPLATE_HPP_
#define LIBRARY_LIBUTILITIES_FOUNDATIONS_BASIS_BASISTEMPLATE_HPP_

#include <LibUtilities/Foundations/Basis.hpp>

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/LibUtilitiesDeclspec.h>
#include <LibUtilities/Foundations/Foundations.hpp>
#include <LibUtilities/Foundations/Shape.hpp>
#include <LibUtilities/Foundations/ShapeTypes.hpp>
#include <LibUtilities/Foundations/Points/PointsTypes.hpp>
#include <LibUtilities/Foundations/Basis/BasisTypes.hpp>
#include <LibUtilities/Foundations/Points.hpp>

#include <LibUtilities/Foundations/Points/PointsTemplate.hpp>
#include <LibUtilities/Foundations/Points/PointsGauss.hpp>
#include <LibUtilities/Foundations/Points/PointsFekete.hpp>

namespace Nektar
{
namespace LibUtilities
{
namespace Foundations
{

template<typename TData>
class BasisBase;

/// A shared pointer to a BasisBase object.
template<typename TData>
using BasisSharedPtr = boost::shared_ptr<BasisBase<TData>>;

/**
 * @brief Primary template for composite Points classes.
 *
 * No implementation is provided, since this non-specialised template should
 * never be instantiated.
 */
template<typename TData, typename TShape, typename... TTuple>
class Basis;

/**
 * @brief Specialisation for composite basis.
 *
 * @details A composite basis is a multi-dimensional basis made up of two or
 * more basis types. For each constituent basis type, the corresponding points
 * type of matching dimension must also be provided.
 *
 * At least one points type \TPts and one basis type \TBasis must be given. The
 * number of points types must match the number of basis types. There may be at
 * most three basis and points types specified.
 */
template<typename TData, typename TShape, typename... TPts, typename... TBasis>
class Basis<TData, TShape, std::tuple<TPts...>, std::tuple<TBasis...>>
            : public BasisBase<TData>
{
        /// The type of this class
        typedef Basis<TData,
                      TShape,
                      std::tuple<TPts...>,
                      std::tuple<TBasis...>
                     > ThisType;
        /// The type of the base class
        typedef BasisBase<TData> BaseType;
        /// The type of the std::tuple which represents the list of constituent
        /// basis types which make up this composite basis.
        typedef std::tuple<
                Basis<TData,
                      typename traits::basis_traits<TBasis>::native_shape,
                      TPts,
                      TBasis
                     >...> TupleType;

        /// At least one basis type must be given.
        static_assert(sizeof...(TBasis) > 0,
                      "No basis type given.");
        /// No more than three basis types can be given
        static_assert(sizeof...(TBasis) < 4,
                      "Too many basis types given. Must be 3 or less.");
        /// The same number of points types and basis types must be given
        static_assert(sizeof...(TPts) == sizeof...(TBasis),
                      "Basis and Points types do not match.");
        /// The dimension of the composite basis type must match the dimension
        /// of the shape type
        static_assert(traits::basis_traits<TBasis...>::dimension ==
                      traits::shape_traits<TShape>::dimension,
                      "Basis dimension and shape dimension do not agree.");
        /// The dimension of the composite points type must match the dimension
        /// of the shape type
        static_assert(traits::points_traits<TPts...>::dimension ==
                      traits::shape_traits<TShape>::dimension,
                      "Points dimension does not match shape dimension,");

    public:
        /**
         * @brief Create an instance of this composite basis type.
         * @param pKey The key describing the parameters for the basis.
         */
        static BasisSharedPtr<TData> create(const BasisKey& pKey)
        {
            BasisSharedPtr<TData> b = MemoryManager<ThisType>::AllocateSharedPtr(pKey);
            return b;
        }


        /**
         * @brief Constructs a new instance of this basis using given key.
         * @param pKey The key describing the parameters for the basis.
         */
        Basis(const BasisKey& pKey) : BasisBase<TData>(pKey) {
            BasisBase<TData>::m_typehash = {typeid(TBasis).hash_code()...};
            // Compute full basis data here...
        }

    protected:
        /// A composite Points object on which the basis is constructed
        Points<TData, TShape, std::tuple<TPts...>> m_points;
        /// A tuple with entries corresponding to the constituent basis types
        TupleType x;

    private:
        virtual int v_GetShapeDimension() const
        {
            return traits::shape_traits<TShape>::dimension;
        }
        virtual int v_GetShapeNumBoundaryElements() const
        {
            return traits::shape_traits<TShape>::num_bnd_elmts;
        }
        virtual bool v_IsBoundaryInterior() const
        {
            return traits::basis_traits<TBasis...>::is_boundary_interior;
        }
        virtual bool v_IsCollocation() const
        {
            return traits::basis_traits<TBasis...>::is_collocation;
        }

        /**
         * @copydoc BasisBase::GetNumConstitutentBases()
         */
        virtual int v_GetNumConstituentBases() const
        {
            return traits::basis_traits<TBasis...>::num_constituent_bases;
        }

        /**
         * @copydoc BasisBase::GetConstituentBasis(int)
         * @details This function returns a reference to a base class object
         *          for the ith entry in the basis tuple.
         */
        virtual const BasisBase<TData>& v_GetConstituentBasis(int i) const
        {
            switch (i)
            {
                case 0: return GetTupleEntry<0>();
                case 1: return GetTupleEntry<1>();
                case 2: return GetTupleEntry<2>();
                default: throw i;
            }
        }

        /**
         * @copydoc BasisBase::GetPoints()
         */
        virtual const PointsBase<TData>& v_GetPoints() const
        {
            return m_points;
        }

        /**
         * @copydoc BasisBase::GetShapeType()
         */
        virtual LibUtilities::Foundations::ShapeType v_GetShapeType() const
        {
            return traits::shape_traits<TShape>::type;
        }

        /**
         * @copydoc BasisBase::GetShapeName()
         */
        virtual std::string v_GetShapeName() const
        {
            return std::string(traits::shape_traits<TShape>::name);
        }

        /**
         * @brief Retrieve a given entry of the tuple by induction.
         * @tparam i The index to retrieve
         * @returns A reference to a BasisBase object which is the ith entry.
         */
        template<size_t i>
        inline typename std::enable_if<i < sizeof...(TBasis), const BasisBase<TData>&>::type
        GetTupleEntry() const
        {
            return std::get<i>(x);
        }

        /**
         * @brief Retrieve a given entry of the tuple by induction.
         * Terminating case.
         * @tparam i The index to retrieve
         * @returns A reference to a BasisBase object which is the ith entry.
         *
         */
        template<size_t i>
        inline typename std::enable_if<i >= sizeof...(TBasis), const BasisBase<TData>&>::type
        GetTupleEntry() const
        {
            throw std::logic_error("Out of range.");
        }
};

/**
 * @brief This defines the ThisType and BaseType data types for the fundamental
 * basis class as well as performing validation on the provided shape and
 * points distributions.
 */
#define BASIS_DEFINES_VALIDATION(x) \
public: \
    typedef Basis<TData, TShape, TPts, x> ThisType; \
    typedef BasisBase<TData> BaseType; \
    static_assert(is_not_tuple<TPts>::value, "Is a tuple."); \
    static_assert(traits::points_traits<TPts>::dimension == \
                  traits::shape_traits<TShape>::dimension, \
                  "Points dimension does not match shape dimension,"); \
    static_assert(traits::basis_traits<ModifiedLegendre>::dimension == \
                  traits::shape_traits<TShape>::dimension, \
                  "Basis dimension and shape dimension do not agree."); \
    static_assert(traits::expansion_traits<TShape, x>::is_valid, \
                  "Not a valid combination of shape and basis type."); \
    static BasisSharedPtr<TData> create(const BasisKey& pKey) \
    { \
        BasisSharedPtr<TData> p \
                = MemoryManager<ThisType>::AllocateSharedPtr(pKey); \
        return p; \
    }

/**
 * @brief This defines the core virtual functions which follow the same pattern
 * for all basis classes.
 */
#define BASIS_CORE_FUNCTIONS(x) \
protected: \
    virtual int v_GetShapeDimension() const \
    { \
        return traits::shape_traits<TShape>::dimension; \
    } \
    virtual int v_GetShapeNumBoundaryElements() const \
    { \
        return traits::shape_traits<TShape>::num_bnd_elmts; \
    } \
    virtual bool v_IsBoundaryInterior() const \
    { \
        return traits::basis_traits<x>::is_boundary_interior; \
    } \
    virtual bool v_IsCollocation() const \
    { \
        return traits::basis_traits<x>::is_collocation; \
    } \
    virtual int v_GetNumConstituentBases() const \
    { \
        return traits::basis_traits<x>::num_constituent_bases; \
    } \
    virtual const PointsBase<TData>& v_GetPoints() const \
    { \
        return m_points; \
    } \
    virtual LibUtilities::Foundations::ShapeType v_GetShapeType() const \
    { \
        return traits::shape_traits<TShape>::type; \
    } \
    virtual std::string v_GetShapeName() const \
    { \
        return std::string(traits::shape_traits<TShape>::name); \
    } \

/**
 * @brief This defines the Points object within the fundamental Basis class.
 */
#define BASIS_CORE_DATA \
private: \
    Points<TData, \
           typename traits::points_traits<TPts>::native_shape, \
           TPts \
          > m_points;

}
}
}
#endif
