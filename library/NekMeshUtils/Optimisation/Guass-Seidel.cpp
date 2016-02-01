////////////////////////////////////////////////////////////////////////////////
//
//  File: SurfaceMeshing.cpp
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
//  Description: surfacemeshing object methods.
//
////////////////////////////////////////////////////////////////////////////////

#include <NekMeshUtils/Optimisation/Guass-Seidel.h>

using namespace std;
namespace Nektar
{
namespace NekMeshUtils
{

    Array<OneD, NekDouble> gsOptimise(NekDouble alpha, Array<OneD, NekDouble> x, DNekMat H, DNekMat J)
    {
        //Array<OneD, NekDouble> eig_r(x.num_elements());
        //Array<OneD, NekDouble> eig_i(x.num_elements());

        //H.EigenSolve(eig_r, eig_i);

        //cout << eig_r[0] << " " << eig_r[1] << " " << eig_r[2] << endl;

        Array<OneD, NekDouble> dx(x.num_elements(),0.0);
        NekDouble Diff;
        int ct = 0;
        do
        {
            Array<OneD, NekDouble> dxn(x.num_elements(),0.0);

            for(int i = 0; i < x.num_elements(); i++)
            {
                NekDouble presum = 0.0;
                for(int j = 0; j < i; j++)
                {
                    presum += H(i,j)*dxn[j];
                }
                NekDouble postsum = 0.0;
                for(int j = i + 1; j < x.num_elements(); j++)
                {
                    postsum += H(i,j)*dx[j];
                }

                dxn[i] = 1.0/H(i,i) * (-J(i,0) - presum - postsum);
            }

            Diff = 0.0;
            for(int i = 0; i < x.num_elements(); i++)
            {
                Diff += (dxn[i] - dx[i]) * (dxn[i] - dx[i]);
            }
            Diff = sqrt(Diff);

            for(int i = 0; i < x.num_elements(); i++)
            {
                dx[i] = dxn[i];
            }

            if(ct>1000000)
            {
                cout << endl << endl << J << endl << endl << H << endl;
                cout << "Failed to converge" << endl;
                abort();
            }
            ct++;

        }
        while(Diff > 1E-6);

        return dx;
    }

}
}
