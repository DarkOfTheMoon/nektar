///////////////////////////////////////////////////////////////////////////////
//
// File: Basis.hpp
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
// Description: Base class for Basis types
//
///////////////////////////////////////////////////////////////////////////////

#ifndef LIBUTILITIES_FOUNDATIONS_BASIS
#define LIBUTILITIES_FOUNDATIONS_BASIS

#include <iostream>
#include <type_traits>
#include <typeinfo>

#include <loki/Singleton.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/NekManager.hpp>
#include <LibUtilities/LibUtilitiesDeclspec.h>
#include <LibUtilities/Foundations/Basis/BasisTypes.hpp>
#include <LibUtilities/Foundations/Basis/BasisKey.h>
#include <LibUtilities/Foundations/Basis/BasisTemplate.hpp>
#include <LibUtilities/Foundations/Points.hpp>
#include <LibUtilities/Foundations/Shape.hpp>

namespace Nektar
{
namespace LibUtilities
{
namespace Foundations
{

/**
 * @brief   Basis base class defining interface and storage for basis data.
 * @details This class provides a base class for all bases. Bases may be
 *          a single fundamental basis type, or a composition of multiple basis
 *          types. For example, a 2D basis may be a fundamentally 2D basis of
 *          2D modes, or a composition of two (possibly different) 1D bases.
 *
 *          The implementation for different fundamental bases is provided in
 *          template specialisations of the Basis class. Each basis has an
 *          associated Points instance which implements the point distribution
 *          on which quadrature is performed.
 *
 * @tparam  TData The data type used for defining the basis.
 */
template<typename TData>
class BasisBase
{
    public:
        /**
         * @brief Destructor.
         */
        virtual ~BasisBase()
        {
        };

        /**
         * @brief Get the dimension of the native shape for this basis.
         */
        inline int GetShapeDimension() const
        {
            return v_GetShapeDimension();
        }

        /**
         * @brief Get the number of boundary elements for the native shape.
         */
        inline int GetShapeNumBoundaryElements() const
        {
            return v_GetShapeNumBoundaryElements();
        }

        /**
         * @brief Get whether the basis is of boundary-interior type.
         */
        inline bool IsBoundaryInterior() const
        {
            return v_IsBoundaryInterior();
        }

        /**
         * @brief Get whether the basis is of collocation type
         */
        inline bool IsCollocation() const
        {
            return v_IsCollocation();
        }

        /**
         * @brief Get number of constituent bases.
         */
        inline int GetNumConstituentBases() const
        {
            return v_GetNumConstituentBases();
        }

        /**
         * @brief Get a constituent basis in a basis composition.
         */
        inline const BasisBase<TData>& GetConstituentBasis(int i) const
        {
            return v_GetConstituentBasis(i);
        }

        /**
         * @copydoc BasisBasis::GetConstituentBasis(int)
         */
        inline const BasisBase<TData>& operator[](const int i) const
        {
            return GetConstituentBasis(i);
        }

        /**
         * Compares the typeid hashcode of this polymorphic
         */
        template<typename... TBasisRhs>
        inline const bool IsOfType() const
        {
            return (m_typehash == std::vector<size_t>({typeid(TBasisRhs).hash_code()...}));
        }

        /**
         * @brief Get the points object associated with this basis.
         */
        inline const PointsBase<TData>& GetPoints() const
        {
            return v_GetPoints();
        }

        /**
         * @brief Get the basis key
         */
        inline const BasisKey& GetKey() const
        {
            return m_key;
        }

        /**
         * @brief Get the shape type associated with this basis.
         */
        inline LibUtilities::Foundations::ShapeType GetShapeType() const
        {
            return v_GetShapeType();
        }

        /**
         * @brief Get the shape name as a string.
         */
        inline std::string GetShapeName() const
        {
            return v_GetShapeName();
        }

        /**
         * @brief Return total number of modes from the basis specification.
         */
        inline int GetNumModes() const
        {
            return m_key.m_nummodes[0];
        }

        /**
         * @brief Returns the maximum number of modes of all constituent
         *        bases.
         */
        inline int GetNumModesMax() const
        {
            unsigned int& a = m_key.m_nummodes[0];
            unsigned int& b = m_key.m_nummodes[1];
            unsigned int& c = m_key.m_nummodes[2];
            return (a > b ? (a > c ? a : (b > c ? b : c)) : (b > c ? b : (a > c ? a : c)));
        }

        /**
         * @brief Get basis matrix data \f$ \phi_j(\xi_i) \f$.
         */
        const Array<OneD, const TData>& GetBdata() const
        {
            return m_bdata;
        }

        /**
         * @brief Get derivative of basis matrix data.
         */
        const Array<OneD, const TData>& GetDbdata() const
        {
            return m_dbdata;
        }

    protected:
        /// Holds the basis key for this basis.
        BasisKey                m_key;

        /// Holds the evaluation of the basis polynomials at quadrature points.
        Array<OneD, TData>      m_bdata;

        /// Holds the derivation matrix for the basis.
        Array<OneD, TData>      m_dbdata;

        /// Type index for comparison
        std::vector<std::size_t> m_typehash;

        /**
         * @brief Constructor. This base class can not be instantiated directly.
         */
        BasisBase(const BasisKey& pKey)
                : m_key(pKey), m_typehash({0})
        {
        }

        /**
         * @brief Default constructor.
         */
        BasisBase() {}

        /**
         * @brief For a multi-constituent basis, compute total number of modes
         * by induction.
         */
        template<typename TBasis1, typename TBasis2, typename... TBasisOther>
        unsigned int GetNumberOfModes(int i = 0)
        {
            return traits::basis_traits<TBasis1>::get_total_modes(m_key.m_nummodes[i])
                    * GetNumberOfModes<TBasis2, TBasisOther...>(i + 1);
        }

        /**
         * @brief For a multi-constituent basis, compute total number of modes
         * by induction. Terminating case.
         */
        template<typename TBasis>
        unsigned int GetNumberOfModes(int i = 0)
        {
            return traits::basis_traits<TBasis>::get_total_modes(m_key.m_nummodes[i]);
        }


        /**
         * @copydoc BasisBase::GetShapeDimension
         */
        virtual int v_GetShapeDimension() const = 0;

        /**
         * @copydoc BasisBase::GetShapeNumBoundaryElements
         */
        virtual int v_GetShapeNumBoundaryElements() const = 0;

        /**
         * @copydoc BasisBase::IsBoundaryInterior()
         */
        virtual bool v_IsBoundaryInterior() const = 0;

        /**
         * @copydoc BasisBase::IsCollocation()
         */
        virtual bool v_IsCollocation() const = 0;

        /**
         * @copydoc BasisBase::GetNumConstitutentBases()
         * Defaults to 1.
         */
        virtual int v_GetNumConstituentBases() const
        {
            return 1;
        }

        /**
         * @copydoc BasisBase::GetConstituentBasis(int)
         */
        virtual const BasisBase<TData>& v_GetConstituentBasis(int i) const
        {
            if (i > 0)
            {
                throw i;
            }
            return *this;
        }

        /**
         * @copydoc BasisBase::GetPoints()
         */
        virtual const PointsBase<TData>& v_GetPoints() const = 0;

        /**
         * @copydoc BasisBase::GetShapeType()
         */
        virtual LibUtilities::Foundations::ShapeType v_GetShapeType() const = 0;

        /**
         * @copydoc BasisBase::GetShapeName()
         */
        virtual std::string v_GetShapeName() const = 0;
};

/// A shared pointer to a BasisBase object.
template<typename TData>
using BasisSharedPtr = boost::shared_ptr<BasisBase<TData>>;

/// A factory pattern for creating Basis objects.
template<typename TData>
using BasisFactory = LibUtilities::NekFactory<
        std::string, BasisBase<TData>, const BasisKey&>;

//template<typename TData>
LIB_UTILITIES_EXPORT BasisFactory<NekDouble>& GetBasisFactory();

/// A manager for Basis objects
template<typename TData>
using BasisManager = LibUtilities::NekManager<
        BasisKey, BasisBase<TData>, BasisKey::opLess>;

LIB_UTILITIES_EXPORT BasisManager<NekDouble>& GetBasisManager();

BasisSharedPtr<NekDouble> BasisCreator(const BasisKey& pKey);

}
}
}

#endif
