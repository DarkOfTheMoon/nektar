///////////////////////////////////////////////////////////////////////////////
//
// File PetscEnv.cpp
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
// Description: Wrap up the PETSc environment
//
///////////////////////////////////////////////////////////////////////////////
#include <iostream>

#include "PetscEnv.h"
#ifdef NEKTAR_USING_PETSC
#include "petscsys.h"

namespace Nektar
{
    namespace Utilities
    {

        PetscEnv::PetscEnv(int* argc, char*** argv) :
                doesOwnPetsc(false)
        {
            if (!Initialized())
            {
                PetscErrorCode ierr = PetscInitialize(argc, argv, NULL, NULL);
                if (ierr)
                    std::cerr << "Error in PetscInitialize. Code = " << ierr << std::endl;
                doesOwnPetsc = true;
            }
        }

        PetscEnv::~PetscEnv()
        {
            if (doesOwnPetsc)
            {
                PetscFinalize();
            }
        }

        bool PetscEnv::Initialized()
        {
            PetscBool flag = PETSC_FALSE;
            PetscInitialized(&flag);
            return flag == PETSC_TRUE;
        }

        bool PetscEnv::Finalized()
        {
            PetscBool flag = PETSC_FALSE;
            PetscFinalized(&flag);
            return flag == PETSC_TRUE;
        }
    }
}
#else
namespace Nektar
{
    namespace Utilities
    {

        PetscEnv::PetscEnv(int* argc, char*** argv) :
                doesOwnPetsc(false)
        {
        }

        PetscEnv::~PetscEnv()
        {
        }

        bool PetscEnv::Initialized()
        {
            return false;
        }

        bool PetscEnv::Finalized()
        {
            return false;
        }
    }
}
#endif
