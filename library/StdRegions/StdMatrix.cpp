///////////////////////////////////////////////////////////////////////////////
//
// File StdMatrix.cpp
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
// Description: Matrix definitions 
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>

#include <LibUtilities/Polylib/Polylib.h>

#include <StdRegions/StdRegions.hpp>
#include <StdRegions/StdMatrix.h>
#include <StdRegions/StdExpansion.h>

#ifdef max
#undef max
#endif

namespace Nektar
{
    namespace StdRegions 
    {

        void StdMatContainer::FillPackedMatrix()
        {
            int i,j,cnt;

            ASSERTL1(m_lda, "lda is not set");

            // check to see if memory is declared and if not setup
            if(!m_packed_matrix)
            {
                SetMemPackedMatrix();
            }

            switch(m_matform)
            {
            case eSymmetric_Positive: 
            case eSymmetric:
                // store matrix in symmetric form 
                cnt = 0;
                for(i = 0; i < m_lda; ++i)
                {
                    for(j = i; j < m_lda; ++j)
                    {
                        m_packed_matrix[cnt++] = m_matrix[j*m_lda+i];
                    }
                }
                break;
            case eSymmetric_Positive_Banded:
                // store matrix in symmetric banded 
                cnt = 0;
                for(i = 0; i < m_lda; ++i)
                {
                    for(j  = i; j < i+m_bwidth; ++j)
                    {
                        m_packed_matrix[cnt++] = m_matrix[i*m_lda+j];
                    }
                }
                break;
            case eGeneral_Banded:
                for(i = 0; i < m_lda; ++i)
                {
                    for(j  = std::max(i-m_ldiag,0); j < i+m_bwidth; ++j)
                    {
                        m_packed_matrix[j*(2*m_ldiag+m_bwidth)+m_ldiag+m_bwidth-1+(i-j)]
                        = m_matrix[i*m_lda+j];
                    }
                }
                break;
            case eGeneral_Full:
                // pack full matrix transposing matrix 
                for(i = 0; i < m_lda; ++i)
                {
                    Blas::Dcopy(m_lda,m_matrix+i*m_lda,1,m_packed_matrix+i,m_lda);
                }
                break;
            }
        }

        void StdMatContainer::Mxv(const double *in, double *out)
        {
            BstShrDArray wsp; 
            double *tmp;

            ASSERTL1(m_lda, "_lda is not set");

            if(in == out)
            {
                wsp = GetDoubleTmpSpace(m_lda);
                tmp = wsp.get();
                Blas::Dcopy(m_lda,in,1,tmp,1);
            }
            else
            {
                tmp = (double *) in;
            }

            // assume row major matrix 
            Blas::Dgemv('T',m_lda,m_lda,1.0,m_matrix,m_lda,tmp,1,0.0,out,1);      
        }

        void StdMatContainer::Factor()
        {
            int info;

            if(m_factored)
            {
                return;
            }

            ASSERTL1(m_lda, "m_lda is not set");

            if(!m_packed_matrix)
            {
                FillPackedMatrix();
            }

            switch(m_matform)
            {
            case eSymmetric:
                m_ipiv = new int[m_lda];
                Lapack::Dsptrf('L',m_lda,m_packed_matrix,m_ipiv,info);
                ASSERTL0(info==0, "matrix did not factor");
                break;
            case eSymmetric_Positive:
                Lapack::Dpptrf('L', m_lda, m_packed_matrix, info);
                ASSERTL0(info==0, "matrix did not factor");
                break;
            case eSymmetric_Positive_Banded:
                Lapack::Dpbtrf('L',m_lda,m_bwidth-1,m_packed_matrix,m_bwidth,info);
                ASSERTL0(info==0, "matrix did not factor");
                break;
            case eGeneral_Banded:
                m_ipiv = new int[m_lda];
                Lapack::Dgbtrf(m_lda,m_lda,m_ldiag,m_bwidth-1,m_packed_matrix,
                    2*m_ldiag+m_bwidth,m_ipiv,info);
                ASSERTL0(info==0, "matrix did not factor");
                break;
            case eGeneral_Full:
                m_ipiv = new int[m_lda];
                Lapack::Dgetrf(m_lda,m_lda,m_packed_matrix,m_lda,m_ipiv,info);     
                ASSERTL0(info==0, "matrix did not factor");
                break;
            }

            m_factored = true;    
        }

        void StdMatContainer::Solve(double *u, const int nrhs)
        {
            int info;

            if(!m_factored) 
            {
                Factor();
            }

            switch(m_matform)
            {
            case eSymmetric:
                Lapack::Dsptrs('L',m_lda,nrhs,m_packed_matrix,m_ipiv,u,m_lda,info);
                ASSERTL0(info==0, "matrix did not solve");
                break;
            case eSymmetric_Positive:
                Lapack::Dpptrs('L', m_lda,nrhs,m_packed_matrix,u,m_lda,info);
                ASSERTL0(info==0, "matrix did not solve");
                break;
            case eSymmetric_Positive_Banded:
                Lapack::Dpbtrs('L',m_lda,m_bwidth-1,nrhs,m_packed_matrix,m_bwidth,u,
                    m_lda,info);
                ASSERTL0(info==0, "matrix did not solve");
                break;
            case eGeneral_Banded:
                Lapack::Dgbtrs('N',m_lda,m_ldiag,m_bwidth-1,nrhs,m_packed_matrix,
                    2*m_ldiag+m_bwidth,m_ipiv,u,m_lda,info);
                ASSERTL0(info==0, "matrix did not solve");
                break;
            case eGeneral_Full:
                Lapack::Dgetrs('N',m_lda,nrhs,m_packed_matrix,m_lda,m_ipiv,u,
                    m_lda,info);     
                ASSERTL0(info==0, "matrix did not solve");
                break;
            }
        }

        void StdMatContainer::SolveTranspose(double *u, const int nrhs)
        {
            int info;

            if(!m_factored) 
            {
                Factor();
            }

            switch(m_matform)
            {
            case eSymmetric:
                Lapack::Dsptrs('L',m_lda,nrhs,m_packed_matrix,m_ipiv,u,m_lda,info);
                ASSERTL0(info==0, "matrix did not solve");
                break;
            case eSymmetric_Positive:
                Lapack::Dpptrs('L', m_lda,nrhs,m_packed_matrix,u,m_lda,info);
                ASSERTL0(info==0, "matrix did not solve");
                break;
            case eSymmetric_Positive_Banded:
                Lapack::Dpbtrs('L',m_lda,m_bwidth-1,nrhs,m_packed_matrix,m_bwidth,u,
                    m_lda,info);
                ASSERTL0(info==0, "matrix did not solve");
                break;
            case eGeneral_Banded:
                Lapack::Dgbtrs('T',m_lda,m_ldiag,m_bwidth-1,nrhs,m_packed_matrix,
                    2*m_ldiag+m_bwidth,m_ipiv,u,m_lda,info);
                ASSERTL0(info==0, "matrix did not solve");
                break;
            case eGeneral_Full:
                Lapack::Dgetrs('T',m_lda,nrhs,m_packed_matrix,m_lda,m_ipiv,u,
                    m_lda,info);     
                ASSERTL0(info==0, "matrix did not solve");
                break;
            }
        }

        // declare the memory of _packed_matrix depending upon its definition
        void StdMatContainer::SetMemPackedMatrix()
        {

            if(m_packed_matrix) 
            {
                return;
            }

            ASSERTL0((m_lda > 0), "m_lda not defined ");

            switch(m_matform)
            {
            case eSymmetric:
            case eSymmetric_Positive: 
                m_packed_matrix = new double [m_lda*(m_lda+1)/2];
                Vmath::Zero(m_lda*(m_lda+1)/2,m_packed_matrix,1);
                break;
            case eSymmetric_Positive_Banded:
                ASSERTL0((m_bwidth > 0 ), "m_bwidth  not set");
                m_packed_matrix = new double [m_lda*m_bwidth];
                Vmath::Zero(m_lda*m_bwidth,m_packed_matrix,1);
                break;
            case eGeneral_Banded:
                ASSERTL0((m_bwidth > 0 )&&(m_ldiag > 0),
                    "m_bwidth or m_ldiag is  not set");
                m_packed_matrix = new double [m_lda*2*(m_ldiag+m_bwidth)];
                Vmath::Zero(m_lda*2*(m_ldiag+m_bwidth),m_packed_matrix,1);
                break;
            case eGeneral_Full:
                m_packed_matrix = new double [m_lda*m_lda];
                Vmath::Zero(m_lda*m_lda,m_packed_matrix,1);
                break;
            } 
        }

        double StdMatContainer::L2ConditionNo()
        {

            BstShrDArray wsp  = GetDoubleTmpSpace(m_lda);
            double *er = wsp.get();
            BstShrDArray wsp1 = GetDoubleTmpSpace(m_lda);
            double *ei = wsp1.get();
            double max,min;

            EigenValues(er,ei,(double *)NULL);

            Vmath::Vmul (m_lda,er,1,er,1,er,1);
            Vmath::Vmul (m_lda,ei,1,ei,1,ei,1);
            Vmath::Vadd (m_lda,er,1,ei,1,er,1);

            max = sqrt(er[Vmath::Imax(m_lda,er,1)]);
            min = sqrt(er[Vmath::Imin(m_lda,er,1)]);

            if(min < 1e-12){ // if min < 1e-12 find second smallest ev
                fprintf(stderr,"Min ev < 1e-12 using second ev\n");
                er[Vmath::Imin(m_lda,er,1)] += max;
                min = sqrt(er[Vmath::Imin(m_lda,er,1)]);
            }

            return max/min;
        }

        double StdMatContainer::MaxEigenValue()
        {
            BstShrDArray wsp  = GetDoubleTmpSpace(m_lda);
            double *er = wsp.get();
            BstShrDArray wsp1 = GetDoubleTmpSpace(m_lda);
            double *ei = wsp1.get();
            double max;

            EigenValues(er,ei,(double *)NULL);

            Vmath::Vmul(m_lda,er,1,er,1,er,1);
            Vmath::Vmul(m_lda,ei,1,ei,1,ei,1);
            Vmath::Vadd(m_lda,er,1,ei,1,er,1);

            max = sqrt(er[Vmath::Imax(m_lda,er,1)]);

            return max;
        }

        int StdMatContainer::NullSpaceDim(const double tol)
        {
            BstShrDArray wsp  = GetDoubleTmpSpace(m_lda);
            double *er = wsp.get();
            BstShrDArray wsp1 = GetDoubleTmpSpace(m_lda);
            double *ei = wsp1.get();
            int i,ndim; 

            EigenValues(er,ei,(double *)NULL);

            Vmath::Vmul (m_lda,er,1,er,1,er,1);
            Vmath::Vmul (m_lda,ei,1,ei,1,ei,1);
            Vmath::Vadd (m_lda,er,1,ei,1,er,1);
            Vmath::Vsqrt(m_lda,er,1,er,1);

            ndim = 0;
            for(i = 0; i < m_lda; ++i)
            {
                if(er[i] < tol)
                {
                    ++ndim;
                }
            }

            return ndim;
        }

        void StdMatContainer::EigenValues(const char file[])
        {
            int i;
            BstShrDArray wsp  = GetDoubleTmpSpace(m_lda);
            double *er = wsp.get();
            BstShrDArray wsp1 = GetDoubleTmpSpace(m_lda);
            double *ei = wsp1.get();
            FILE *fp;

            fp = fopen(file,"w");

            EigenValues(er,ei,(double *)NULL);

            fprintf(fp,"# Real Imag Magnitude\n");

            for(i = 0; i < m_lda; ++i)
            {
                fprintf(fp,"%lg %lg %lg \n",er[i],ei[i],sqrt(er[i]*er[i]+ei[i]*ei[i]));
            }

            fclose(fp);
        }

        void StdMatContainer::EigenValues(double *er, double *ei, double *evecs)
        {
            double dum;
            int    info;

            switch(m_matform){
      case eSymmetric_Positive:
      case eSymmetric:
          {
              BstShrDArray work = GetDoubleTmpSpace(3*m_lda);
              Vmath::Zero(m_lda,ei,1);
              Lapack::Dspev('N','L',m_lda,m_packed_matrix,er,&dum,1,work.get(),info);
              ASSERTL0(info==0, "info is not zero");
              break;
          }
      case eSymmetric_Positive_Banded:
          {
              BstShrDArray work = GetDoubleTmpSpace(3*m_lda);
              Vmath::Zero(m_lda,ei,1);
              Lapack::Dsbev('N','L',m_lda,m_bwidth-1,m_packed_matrix,m_bwidth,er,
                  &dum,1,work.get(),info);
              ASSERTL0(info==0, "info is not zero");
              break;
          }
      case eGeneral_Banded:
          ASSERTL0(false, "Eigenvalue evaluation for genaral baneded matrix needs coding");
          break;
      case eGeneral_Full:
          {
              BstShrDArray work = GetDoubleTmpSpace(4*m_lda);
              if(evecs)
              {
                  Lapack::Dgeev('N','V',m_lda,m_packed_matrix,m_lda,er,ei,&dum,1,
                      evecs,m_lda,work.get(),4*m_lda,info);
              }
              else
              {
                  double dum1;
                  Lapack::Dgeev('N','N',m_lda,m_packed_matrix,m_lda,er,ei,&dum,1,
                      &dum1,1,work.get(),4*m_lda,info);
              }	  
              ASSERTL0(info==0,"info is not zero");
              break;
          }
            }
        }

        StdMatContainer * StdMatrix::GetLocalMass(StdExpansion * E)
        {
            std::vector<StdMatContainer*>::iterator def;
            StdMatContainer *M;       

            def = find(m_local_mass.begin(),m_local_mass.end(),m_local_mass_cur,*E);

            if(def != m_local_mass.end())
            {
                M = def[0];
                m_local_mass_cur = def;
            }
            else
            {
                double * outarray = new double[E->GetNcoeffs()*E->GetNcoeffs()];
                E->GenMassMatrix(outarray);

                StdMatContainer * tmp = new StdMatContainer(E,outarray);

                // set up how matrix should be inverted 
		E->SetInvInfo(tmp,eMassMatrix);

                m_local_mass.push_back(tmp);
                m_local_mass_cur = (--m_local_mass.end());
		
                M = tmp; 
            }
	    
            return M;
        }

        StdMatContainer * StdMatrix::GetLocalLap(StdExpansion * E)
        {
            std::vector<StdMatContainer*>::iterator def;
            StdMatContainer *M;

            def = find(m_local_lap.begin(),m_local_lap.end(),m_local_lap_cur,*E);

            if(def != m_local_lap.end())
            {
                M = def[0];
                m_local_lap_cur = def;
            }
            else
            {
                double * outarray = new double[E->GetNcoeffs()*E->GetNcoeffs()];
                E->GenLapMatrix(outarray);

                StdMatContainer * tmp = new StdMatContainer(E,outarray);

                // set up how matrix should be inverted 
		E->SetInvInfo(tmp,eLapMatrix);

                m_local_lap.push_back(tmp);
                m_local_lap_cur = (--m_local_lap.end());
                M = tmp;
            }
            return M;
        }

        StdMatContainer * StdMatrix::GetNBasisTrans(StdExpansion * E)
        {
            std::vector<StdMatContainer*>::iterator def;
            StdMatContainer *M;       
	    
            def = find(m_nbasis_trans.begin(),m_nbasis_trans.end(),
		       m_nbasis_trans_cur,*E);

            if(def != m_nbasis_trans.end())
            {
                M = def[0];
                m_nbasis_trans_cur = def;
            }
            else
            {
                double * outarray = new double[E->GetNcoeffs()*E->GetNcoeffs()];
                E->GenNBasisTransMatrix(outarray);

                StdMatContainer * tmp = new StdMatContainer(E,outarray);
                // set up how matrix should be inverted 
		E->SetInvInfo(tmp,eNBasisTrans);

                m_nbasis_trans.push_back(tmp);
                m_nbasis_trans_cur = (--m_nbasis_trans.end());

                M = tmp; 
            }
            return M;
        }

        //-----------------------------------
        //  I/O Routines 
        //-----------------------------------

        void StdMatContainer::DumpMatrix(FILE *out)
        {
            int i,j;

            ASSERTL0(m_lda, "m_lda not defined");
            ASSERTL0(m_matrix, "m_matrix not defined");

            for(i = 0; i < m_lda; ++i)
            {
                for(j = 0; j < m_lda; ++j)
                {
                    fprintf(out, "%9.6lf ",m_matrix[i*m_lda + j]);
                }
                fputc('\n',out);
            }
        }
        void StdMatContainer::ShowMatrixStructure(FILE *out)
        {
            int i,j;
            ASSERTL0(m_lda, "m_lda not defined");
            ASSERTL0(m_matrix, "m_matrix not defined");

            fprintf(out,"\nMatrix entries are defined as:\n");
            fprintf(out,"\t+ is a positive entry,\n"
                "\t* is a negative entry,\n"
                "\t- is a zero entry (< 1e-12)\n\n");

            for(i = 0; i < m_lda; ++i)
            {
                for(j = 0; j < m_lda; ++j)
                {
                    if(fabs(m_matrix[i*m_lda + j]) > 1e-12)
                        if(m_matrix[i*m_lda + j] > 0)
                            fputc('+',out);
                        else
                            fputc('*',out);
                    else
                        fputc('-',out);
                } 
                fputc('\n',out);
            }
            fputc('\n',out);
        }


        bool operator  == (const StdExpansion &x, const StdMatContainer *y)
        {
            bool var = true;

            if(x.m_numbases != (*y).m_numbases)
            {
                var = false;
            }
            else
            {
                var = true;
            }

            for(int i=0;i<x.m_numbases;i++)
            {
                var = var && (x.m_base[i]->GetBasisOrder() == 
                    (*y).m_base[i]->GetBasisOrder());
                var = var && (x.m_base[i] == (*y).m_base[i]);
            }

            return var;
        }

        bool operator  != (const StdExpansion &x, const StdMatContainer *y)
        {
            bool var = true;

            if(x.m_numbases != (*y).m_numbases)
            {
                var = false;
            }
            else
            {
                var = true;
            }

            for(int i=0;i<x.m_numbases;i++)
            {
                var = var && (x.m_base[i]->GetBasisOrder() == 
                    (*y).m_base[i]->GetBasisOrder());
                var = var && (x.m_base[i] == (*y).m_base[i]);
            }

            return !var;
        }


        bool operator  == (const StdMatContainer *x, const StdExpansion &y)
        {
            bool var = true;

            if((*x).m_numbases != y.m_numbases)
            {
                var = false;
            }
            else
            {
                var = true;
            }

            for(int i=0;i<y.m_numbases;i++)
            {
                var = var && ((*x).m_base[i]->GetBasisOrder() == 
                    y.m_base[i]->GetBasisOrder());
                var = var && ((*x).m_base[i] == y.m_base[i]);
            }

            return var;
        }


        bool operator  != (const StdMatContainer *x, const StdExpansion &y)
        {
            bool var = true;

            if((*x).m_numbases != y.m_numbases)
            {
                var = false;
            }
            else
            {
                var = true;
            }

            for(int i=0;i<y.m_numbases;i++)
            {
                var = var && ((*x).m_base[i]->GetBasisOrder() == 
                    y.m_base[i]->GetBasisOrder());
                var = var && ((*x).m_base[i] == y.m_base[i]);
            }

            return !var;
        }

    } // end of namespace stdregion    
} // end of namespace 


/** 
* $Log: StdMatrix.cpp,v $
* Revision 1.9  2007/01/17 16:05:40  pvos
* updated doxygen documentation
*
* Revision 1.8  2007/01/10 20:17:31  sherwin
* Accidentally left in a ShowMatrix in Factor
*
* Revision 1.7  2007/01/10 20:13:25  sherwin
* Fixed bug in NodalPointReorder2D
*
* Revision 1.6  2006/12/10 19:00:54  sherwin
* Modifications to handle nodal expansions
*
* Revision 1.5  2006/07/02 17:16:18  sherwin
*
* Modifications to make MultiRegions work for a connected domain in 2D (Tris)
*
* Revision 1.4  2006/06/06 15:25:21  jfrazier
* Removed unreferenced variables and replaced ASSERTL0(false, ....) with
* NEKERROR.
*
* Revision 1.3  2006/06/01 15:03:59  sherwin
* Updated header file to include StdRegions.hpp and removed NekMemoryManager.hpp
*
* Revision 1.2  2006/06/01 14:46:16  kirby
* *** empty log message ***
*
* Revision 1.1  2006/05/04 18:58:31  kirby
* *** empty log message ***
*
* Revision 1.29  2006/04/25 20:23:33  jfrazier
* Various fixes to correct bugs, calls to ASSERT, etc.
*
* Revision 1.28  2006/04/01 21:59:27  sherwin
* Sorted new definition of ASSERT
*
* Revision 1.27  2006/03/21 09:21:32  sherwin
* Introduced NekMemoryManager
*
* Revision 1.26  2006/03/01 22:59:12  sherwin
*
* First working version of Project1D
*
* Revision 1.25  2006/03/01 17:07:33  sherwin
*
* Added new location of polylib in header definitions
*
* Revision 1.24  2006/03/01 08:25:04  sherwin
*
* First compiling version of StdRegions
*
* Revision 1.23  2006/02/26 23:37:30  sherwin
*
* Updates and compiling checks upto StdExpansions1D
*
**/ 






