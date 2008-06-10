///////////////////////////////////////////////////////////////////////////////
//
// File: MatrixStorageType.h
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
// Description: Defines the classes used in the NekMatrix template parameter list
// to define different storage types for matrices.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_MATRIX_STORAGE_TYPE_H
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_MATRIX_STORAGE_TYPE_H

namespace Nektar
{
    enum MatrixStorage 
    {
        eFULL, 
        eDIAGONAL,
        eUPPER_TRIANGULAR,
        eLOWER_TRIANGULAR,
        eSYMMETRIC,
        eBANDED,
        eSYMMETRIC_BANDED,
        eUPPER_TRIANGULAR_BANDED,
        eLOWER_TRIANGULAR_BANDED
    };
      

    /// \brief Tag for matrices using full storage.
    class FullMatrixTag {};

    /// \brief Tag for matrices which only store data on the diagonal.
    class DiagonalMatrixTag {};

    /// \brief Tag for matrices which store data in the upper triangular portion.
    class UpperTriangularMatrixTag {};

    /// \brief Tag for matrices which store data in the lower triangular portion.
    class LowerTriangularMatrixTag {};

    /// \brief Tag for matrices which are symmetric.
    class SymmetricMatrixTag {};

    /// \brief Tag for matrices which are banded.
    class BandedMatrixTag {};

    class SymmetricBandedMatrixTag{};
    class UpperTriangularBandedMatrixTag {};
    class LowerTriangularBandedMatrixTag {};

    template<typename T>
    class ConvertToMatrixStorageEnum;
    
    template<>
    class ConvertToMatrixStorageEnum<FullMatrixTag>
    {
        public:
            enum { Value = eFULL };
    };
        
    template<>
    class ConvertToMatrixStorageEnum<DiagonalMatrixTag>
    {
        public:
            enum { Value = eDIAGONAL };
    };

    template<>
    class ConvertToMatrixStorageEnum<UpperTriangularMatrixTag>
    {
        public:
            enum { Value = eUPPER_TRIANGULAR };
    };
    
    template<>
    class ConvertToMatrixStorageEnum<LowerTriangularMatrixTag>
    {
        public:
            enum { Value = eLOWER_TRIANGULAR };
    };

    template<>
    class ConvertToMatrixStorageEnum<SymmetricMatrixTag>
    {
        public:
            enum { Value = eSYMMETRIC };
    };

    template<>
    class ConvertToMatrixStorageEnum<BandedMatrixTag>
    {
        public:
            enum { Value = eBANDED };
    };

    template<>
    class ConvertToMatrixStorageEnum<SymmetricBandedMatrixTag>
    {
        public:
            enum { Value = eSYMMETRIC_BANDED };
    };

    template<>
    class ConvertToMatrixStorageEnum<UpperTriangularBandedMatrixTag>
    {
        public:
            enum { Value = eUPPER_TRIANGULAR_BANDED };
    };

    template<>
    class ConvertToMatrixStorageEnum<LowerTriangularBandedMatrixTag>
    {
        public:
            enum { Value = eLOWER_TRIANGULAR_BANDED };
    };
        
 }   
    
#endif //NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_MATRIX_STORAGE_TYPE_H
