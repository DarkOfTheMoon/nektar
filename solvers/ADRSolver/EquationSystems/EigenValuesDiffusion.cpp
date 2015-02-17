///////////////////////////////////////////////////////////////////////////////
//
// File EigenValuesDiffusion.cpp
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
// Description:
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>

#include <ADRSolver/EquationSystems/EigenValuesDiffusion.h>

namespace Nektar
{
    string EigenValuesDiffusion::className = GetEquationSystemFactory().
        RegisterCreatorFunction(
            "EigenValuesDiffusion", EigenValuesDiffusion::create,
            "Eigenvalues of the weak advection operator.");

    EigenValuesDiffusion::EigenValuesDiffusion(
        const LibUtilities::SessionReaderSharedPtr& pSession)
        : UnsteadyDiffusion(pSession)
    {
    }

    void EigenValuesDiffusion::v_InitObject()
    {
        UnsteadyDiffusion::v_InitObject();
    }

    EigenValuesDiffusion::~EigenValuesDiffusion()
    {

    }

    void EigenValuesDiffusion::v_DoSolve()
    {
        const int nVariables = 1;
        const int npoints = GetNpoints();

        Array<OneD, Array<OneD, NekDouble> > inarray (nVariables);
        Array<OneD, Array<OneD, NekDouble> > tmp     (nVariables);
        Array<OneD, Array<OneD, NekDouble> > outarray(nVariables);
        Array<OneD, NekDouble> weakMatrix(npoints * npoints,0.0);

        inarray [0] = Array<OneD, NekDouble>(npoints, 0.0);
        outarray[0] = Array<OneD, NekDouble>(npoints, 0.0);
        tmp     [0] = Array<OneD, NekDouble>(npoints, 0.0);

        for (int j = 0; j < npoints; j++)
        {
            inarray[0][j] = 1.0;

            /// Feeding the weak Diffusion oprator with a vector (inarray)
            /// Looping on inarray and changing the position of the only
            /// non-zero entry we simulate the multiplication by the identity
            /// matrix.  The results stored in outarray is one of the columns of
            /// the weak advection oprators which are then stored in MATRIX for
            /// the futher eigenvalues calculation.
            m_diffusion->Diffuse(nVariables, m_fields, inarray, outarray);

            Vmath::Vcopy(npoints, &outarray[0][0], 1, &weakMatrix[j*npoints], 1);

            /// Set the j-th entry of inarray back to zero
            inarray[0][j] = 0.0;
        }

        // Calulating the eigenvalues of the weak advection operator stored in
        // (MATRIX) using Lapack routines

        char jobvl = 'N';
        char jobvr = 'N';
        int info = 0, lwork = 3*npoints;
        NekDouble dum;

        Array<OneD, NekDouble> eigReal(npoints);
        Array<OneD, NekDouble> eigImag(npoints);

        Array<OneD, NekDouble> work(lwork);

        Lapack::Dgeev(jobvl,jobvr,npoints,weakMatrix.get(),npoints,eigReal.get(),eigImag.get(),&dum,1,&dum,1,&work[0],lwork,info);

        ofstream eigFile("eig-diffusion.txt");

        for(int j = 0; j < npoints; j++)
        {
            eigFile << scientific << eigReal[j] << " " << eigImag[j] << endl;
        }

        eigFile.close();
    }
}
