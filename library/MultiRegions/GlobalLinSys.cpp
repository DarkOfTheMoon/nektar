///////////////////////////////////////////////////////////////////////////////
//
// File GlobalLinSys.cpp
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
// Description: GlobalLinSys definition
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/GlobalLinSys.h>
#include <MultiRegions/Preconditioner.h>
#include <LocalRegions/MatrixKey.h>
#include <LocalRegions/Expansion.h>
#include <LibUtilities/BasicUtils/SessionReader.h>

#include <LibUtilities/LinearAlgebra/NekTypeDefs.hpp>
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>

namespace Nektar
{
    namespace MultiRegions
    {
        std::string GlobalLinSys::lookupIds[12] = {
            LibUtilities::SessionReader::RegisterEnumValue(
                "GlobalSysSoln", "DirectFull",
                MultiRegions::eDirectFullMatrix),
            LibUtilities::SessionReader::RegisterEnumValue(
                "GlobalSysSoln", "DirectStaticCond",
                MultiRegions::eDirectStaticCond),
            LibUtilities::SessionReader::RegisterEnumValue(
                "GlobalSysSoln", "DirectMultiLevelStaticCond",
                MultiRegions::eDirectMultiLevelStaticCond),
            LibUtilities::SessionReader::RegisterEnumValue(
                "GlobalSysSoln", "IterativeFull",
                MultiRegions::eIterativeFull),
            LibUtilities::SessionReader::RegisterEnumValue(
                "GlobalSysSoln", "IterativeStaticCond",
                MultiRegions::eIterativeStaticCond),
            LibUtilities::SessionReader::RegisterEnumValue(
                "GlobalSysSoln", "IterativeMultiLevelStaticCond",
                MultiRegions::eIterativeMultiLevelStaticCond),
            LibUtilities::SessionReader::RegisterEnumValue(
                "GlobalSysSoln", "XxtFull",
                MultiRegions::eXxtFullMatrix),
            LibUtilities::SessionReader::RegisterEnumValue(
                "GlobalSysSoln", "XxtStaticCond",
                MultiRegions::eXxtStaticCond),
            LibUtilities::SessionReader::RegisterEnumValue(
                "GlobalSysSoln", "XxtMultiLevelStaticCond",
                MultiRegions::eXxtMultiLevelStaticCond),
            LibUtilities::SessionReader::RegisterEnumValue(
                "GlobalSysSoln", "PETScFull",
                MultiRegions::ePETScFullMatrix),
            LibUtilities::SessionReader::RegisterEnumValue(
                "GlobalSysSoln", "PETScStaticCond",
                MultiRegions::ePETScStaticCond),
            LibUtilities::SessionReader::RegisterEnumValue(
                "GlobalSysSoln", "PETScMultiLevelStaticCond",
                MultiRegions::ePETScMultiLevelStaticCond)
        };

        std::string GlobalLinSys::def = LibUtilities::SessionReader::
            RegisterDefaultSolverInfo("GlobalSysSoln",
                                      "DirectMultiLevelStaticCond");

        /**
         * @class GlobalLinSys
         *
         * Consider the linear system
         * \f$\boldsymbol{M\hat{u}}_g=\boldsymbol{\hat{f}}\f$.
         * Distinguishing between the boundary and interior components of
         * \f$\boldsymbol{\hat{u}}_g\f$ and \f$\boldsymbol{\hat{f}}\f$ using
         * \f$\boldsymbol{\hat{u}}_b\f$,\f$\boldsymbol{\hat{u}}_i\f$ and
         * \f$\boldsymbol{\hat{f}}_b\f$,\f$\boldsymbol{\hat{f}}_i\f$
         * respectively, this system can be split into its constituent parts as
         * \f[\left[\begin{array}{cc}
         * \boldsymbol{M}_b&\boldsymbol{M}_{c1}\\
         * \boldsymbol{M}_{c2}&\boldsymbol{M}_i\\
         * \end{array}\right]
         * \left[\begin{array}{c}
         * \boldsymbol{\hat{u}_b}\\
         * \boldsymbol{\hat{u}_i}\\
         * \end{array}\right]=
         * \left[\begin{array}{c}
         * \boldsymbol{\hat{f}_b}\\
         * \boldsymbol{\hat{f}_i}\\
         * \end{array}\right]\f]
         * where \f$\boldsymbol{M}_b\f$ represents the components of
         * \f$\boldsymbol{M}\f$ resulting from boundary-boundary mode
         * interactions,
         * \f$\boldsymbol{M}_{c1}\f$ and \f$\boldsymbol{M}_{c2}\f$ represent the
         * components resulting from coupling between the boundary-interior
         * modes, and \f$\boldsymbol{M}_i\f$ represents the components of
         * \f$\boldsymbol{M}\f$ resulting from interior-interior mode
         * interactions.
         *
         * The solution of the linear system can now be determined in two steps:
         * \f{eqnarray*}
         * \mathrm{step 1:}&\quad&(\boldsymbol{M}_b-\boldsymbol{M}_{c1}
         * \boldsymbol{M}_i^{-1}\boldsymbol{M}_{c2}) \boldsymbol{\hat{u}_b} =
         * \boldsymbol{\hat{f}}_b - \boldsymbol{M}_{c1}\boldsymbol{M}_i^{-1}
         * \boldsymbol{\hat{f}}_i,\nonumber \\
         * \mathrm{step 2:}&\quad&\boldsymbol{\hat{u}_i}=\boldsymbol{M}_i^{-1}
         * \left( \boldsymbol{\hat{f}}_i
         *      - \boldsymbol{M}_{c2}\boldsymbol{\hat{u}_b}
         * \right). \nonumber \\ \f}
         * As the inverse of \f$\boldsymbol{M}_i^{-1}\f$ is
         * \f[ \boldsymbol{M}_i^{-1} = \left [\underline{\boldsymbol{M}^e_i}
         * \right ]^{-1} = \underline{[\boldsymbol{M}^e_i]}^{-1} \f]
         * and the following operations can be evaluated as,
         * \f{eqnarray*}
         * \boldsymbol{M}_{c1}\boldsymbol{M}_i^{-1}\boldsymbol{\hat{f}}_i &
         * =& \boldsymbol{\mathcal{A}}_b^T \underline{\boldsymbol{M}^e_{c1}}
         * \underline{[\boldsymbol{M}^e_i]}^{-1} \boldsymbol{\hat{f}}_i \\
         * \boldsymbol{M}_{c2} \boldsymbol{\hat{u}_b} &=&
         * \underline{\boldsymbol{M}^e_{c2}} \boldsymbol{\mathcal{A}}_b
         * \boldsymbol{\hat{u}_b}.\f}
         * where \f$\boldsymbol{\mathcal{A}}_b \f$ is the permutation matrix
         * which scatters from global to local degrees of freedom, only the
         * following four matrices should be constructed:
         * - \f$\underline{[\boldsymbol{M}^e_i]}^{-1}\f$
         * - \f$\underline{\boldsymbol{M}^e_{c1}}
         *                          \underline{[\boldsymbol{M}^e_i]}^{-1}\f$
         * - \f$\underline{\boldsymbol{M}^e_{c2}}\f$
         * - The Schur complement: \f$\boldsymbol{M}_{\mathrm{Schur}}=
         *   \quad\boldsymbol{M}_b-\boldsymbol{M}_{c1}\boldsymbol{M}_i^{-1}
         *   \boldsymbol{M}_{c2}\f$
         *
         * The first three matrices are just a concatenation of the
         * corresponding local matrices and they can be created as such. They
         * also allow for an elemental evaluation of the operations concerned.
         *
         * The global Schur complement however should be assembled from the
         * concatenation of the local elemental Schur complements, that is,
         * \f[ \boldsymbol{M}_{\mathrm{Schur}}=\boldsymbol{M}_b
         *          - \boldsymbol{M}_{c1}
         * \boldsymbol{M}_i^{-1} \boldsymbol{M}_{c2} =
         * \boldsymbol{\mathcal{A}}_b^T \left [\underline{\boldsymbol{M}^e_b -
         * \boldsymbol{M}^e_{c1} [\boldsymbol{M}^e_i]^{-1}
         * (\boldsymbol{M}^e_{c2})} \right ] \boldsymbol{\mathcal{A}}_b \f]
         * and it is the only matrix operation that need to be evaluated on a
         * global level when using static condensation.
         * However, due to the size and sparsity of the matrix
         * \f$\boldsymbol{\mathcal{A}}_b\f$, it is more efficient to assemble
         * the global Schur matrix using the mapping array bmap\f$[e][i]\f$
         * contained in the input argument \a locToGloMap. The global Schur
         * complement is then constructed as:
         * \f[\boldsymbol{M}_{\mathrm{Schur}}\left[\mathrm{\a bmap}[e][i]\right]
         * \left[\mathrm{\a bmap}[e][j]\right]=\mathrm{\a bsign}[e][i]\cdot
         * \mathrm{\a bsign}[e][j]
         * \cdot\boldsymbol{M}^e_{\mathrm{Schur}}[i][j]\f]
         * All four matrices are stored in the \a GlobalLinSys returned by this
         * function.
         */

        /**
         * Given a block matrix, construct a global matrix system according to
         * a local to global mapping. #m_linSys is constructed by
         * AssembleFullMatrix().
         * @param   pkey        Associated linear system key.
         * @param   locToGloMap Local to global mapping.
         */
        GlobalLinSys::GlobalLinSys(const GlobalLinSysKey &pKey,
                const boost::weak_ptr<ExpList> &pExpList,
                const boost::shared_ptr<AssemblyMap>
                                   &pLocToGloMap):
            m_linSysKey(pKey),
            m_expList(pExpList),
            m_robinBCInfo(m_expList.lock()->GetRobinBCInfo()),
            m_weakDirichletBCInfo(m_expList.lock()->GetWeakDirichletBCInfo())
        {
        }

        /**
         *
         */
        GlobalLinSysFactory& GetGlobalLinSysFactory()
        {
            typedef Loki::SingletonHolder<GlobalLinSysFactory,
                Loki::CreateUsingNew,
                Loki::NoDestroy,
                Loki::SingleThreaded> Type;
            return Type::Instance();
        }

        /**
         * @brief Create a preconditioner object from the parameters defined in
         * the supplied assembly map.
         *
         * @param asmMap  Assembly map used to construct the global system.
         */
        PreconditionerSharedPtr GlobalLinSys::CreatePrecon(AssemblyMapSharedPtr asmMap)
        {
            PreconditionerType pType = asmMap->GetPreconType();
            std::string PreconType = MultiRegions::PreconditionerTypeMap[pType];
            return GetPreconFactory().CreateInstance(
                PreconType, GetSharedThisPtr(), asmMap);
        }


        /**
         * @brief Get the number of blocks in this system. 
         *
         * At the top level this corresponds to the number of elements in the
         * expansion list.
         */
        int GlobalLinSys::v_GetNumBlocks()
        {
            return m_expList.lock()->GetExpSize();
        }

        /**
         * @brief Retrieves the block matrix from n-th expansion using the
         * matrix key provided by the #m_linSysKey.
         *
         * @param   n           Number of the expansion.
         * @return              Block matrix for the specified expansion.
         */
        DNekScalMatSharedPtr GlobalLinSys::v_GetBlock(unsigned int n)
        {
            boost::shared_ptr<MultiRegions::ExpList> expList = m_expList.lock();
            int cnt = 0;
            DNekScalMatSharedPtr loc_mat;

            LocalRegions::ExpansionSharedPtr vExp = 
                boost::dynamic_pointer_cast<LocalRegions::Expansion>(
                    expList->GetExp(n));

            // need to be initialised with zero size for non variable
            // coefficient case
            StdRegions::VarCoeffMap vVarCoeffMap;

            // retrieve variable coefficients
            if(m_linSysKey.GetNVarCoeffs() > 0)
            {
                StdRegions::VarCoeffMap::const_iterator x;
                cnt = expList->GetPhys_Offset(n);
                
                for (x = m_linSysKey.GetVarCoeffs().begin(); 
                     x != m_linSysKey.GetVarCoeffs().end(); ++x)
                {
                    vVarCoeffMap[x->first] = x->second + cnt;
                }
            }

            LocalRegions::MatrixKey matkey(m_linSysKey.GetMatrixType(),
                                           vExp->DetShapeType(),
                                           *vExp, m_linSysKey.GetConstFactors(),
                                           vVarCoeffMap);
            loc_mat = vExp->GetLocMatrix(matkey);

            // apply robin boundary conditions to the matrix.
            if(m_robinBCInfo.count(n) != 0) // add robin mass matrix
            {
                RobinBCInfoSharedPtr rBC;

                // declare local matrix from scaled matrix.
                int rows = loc_mat->GetRows();
                int cols = loc_mat->GetColumns();
                const NekDouble *dat = loc_mat->GetRawPtr();
                DNekMatSharedPtr new_mat = MemoryManager<DNekMat>::
                    AllocateSharedPtr(rows,cols,dat);
                Blas::Dscal(rows*cols,loc_mat->Scale(),new_mat->GetRawPtr(),1);

                // add local matrix contribution
                for(rBC = m_robinBCInfo.find(n)->second;rBC; rBC = rBC->next)
                {
                    vExp->AddRobinMassMatrix(
                        rBC->m_robinID, rBC->m_robinPrimitiveCoeffs, new_mat);
                }

                // redeclare loc_mat to point to new_mat plus the scalar.
                loc_mat = MemoryManager<DNekScalMat>::AllocateSharedPtr(
                    1.0, new_mat);
            }

            if(m_weakDirichletBCInfo.count(n) != 0) // Add matrices for weak Dirichlet BCs
            {
                WeakDirichletBCInfoSharedPtr wDBC;

                // declare local matrix from scaled matrix.
                int rows = loc_mat->GetRows();
                int cols = loc_mat->GetColumns();

                /*
                std::cout << "***" << std::endl;
                std::cout << "matrix rank in bc = " << rows << " x " << cols << std::endl;
                std::cout << "***" << std::endl;
                */

                const NekDouble *dat = loc_mat->GetRawPtr();
                DNekMatSharedPtr new_mat = MemoryManager<DNekMat>::
                    AllocateSharedPtr(rows,cols,dat);
                Blas::Dscal(rows*cols,loc_mat->Scale(),new_mat->GetRawPtr(),1);

                int numBdryFacets = 0;
                for(wDBC = m_weakDirichletBCInfo.find(n)->second; wDBC; wDBC = wDBC->next)
                {
                    numBdryFacets++;
                }

                Array<OneD, int> edgeids(numBdryFacets);

                numBdryFacets = 0;
                for(wDBC = m_weakDirichletBCInfo.find(n)->second; wDBC; wDBC = wDBC->next)
                {
                    edgeids[numBdryFacets] = wDBC->m_weakDirichletID;
                    numBdryFacets++;
                }

                vExp->AddWeakDirichletElementContribution(edgeids, *new_mat);

                // redeclare loc_mat to point to new_mat plus the scalar.
                loc_mat = MemoryManager<DNekScalMat>::AllocateSharedPtr(
                    1.0, new_mat);
            }

            // finally return the matrix.
            return loc_mat;
        }

        /**
         * @brief Retrieves a the static condensation block matrices from n-th
         * expansion using the matrix key provided by the #m_linSysKey.
         * 
         * @param   n           Number of the expansion
         * @return              2x2 Block matrix holding the static condensation
         *                      matrices for the n-th expansion.
         */
        DNekScalBlkMatSharedPtr GlobalLinSys::v_GetStaticCondBlock(
            unsigned int n)
        {
            boost::shared_ptr<MultiRegions::ExpList> expList = m_expList.lock();
            int cnt = 0;
            DNekScalBlkMatSharedPtr loc_mat;
            DNekScalMatSharedPtr    tmp_mat;

            StdRegions::StdExpansionSharedPtr vExp = expList->GetExp(n);

            // need to be initialised with zero size for non variable
            // coefficient case
            StdRegions::VarCoeffMap vVarCoeffMap;

            // retrieve variable coefficients
            if(m_linSysKey.GetNVarCoeffs() > 0)
            {
                StdRegions::VarCoeffMap::const_iterator x;
                cnt = expList->GetPhys_Offset(n);
                for (x  = m_linSysKey.GetVarCoeffs().begin(); 
                     x != m_linSysKey.GetVarCoeffs().end  (); ++x)
                {
                    vVarCoeffMap[x->first] = x->second + cnt;
                }
            }

            LocalRegions::MatrixKey matkey(m_linSysKey.GetMatrixType(),
                                           vExp->DetShapeType(),
                                           *vExp,
                                           m_linSysKey.GetConstFactors(),
                                           vVarCoeffMap);

            loc_mat = vExp->GetLocStaticCondMatrix(matkey);

            if(m_robinBCInfo.count(n) != 0) // add robin mass matrix
            {
                RobinBCInfoSharedPtr rBC;

                tmp_mat = loc_mat->GetBlock(0,0);

                // declare local matrix from scaled matrix.
                int rows = tmp_mat->GetRows();
                int cols = tmp_mat->GetColumns();
                const NekDouble *dat = tmp_mat->GetRawPtr();
                DNekMatSharedPtr new_mat = MemoryManager<DNekMat>::
                    AllocateSharedPtr(rows, cols, dat);
                Blas::Dscal(rows*cols,tmp_mat->Scale(),new_mat->GetRawPtr(),1);

                // add local matrix contribution
                for(rBC = m_robinBCInfo.find(n)->second;rBC; rBC = rBC->next)
                {
                    vExp->AddRobinMassMatrix(
                        rBC->m_robinID, rBC->m_robinPrimitiveCoeffs, new_mat);
                }

                // redeclare loc_mat to point to new_mat plus the scalar.
                tmp_mat = MemoryManager<DNekScalMat>::AllocateSharedPtr(
                    1.0, new_mat);
                DNekScalBlkMatSharedPtr new_loc_mat;
                unsigned int exp_size[] = {tmp_mat->GetRows(), loc_mat->GetBlock(1,1)->GetRows()};
                unsigned int nblks = 2;
                new_loc_mat = MemoryManager<DNekScalBlkMat>::AllocateSharedPtr(nblks, nblks, exp_size, exp_size);


                new_loc_mat->SetBlock(0,0,tmp_mat);
                new_loc_mat->SetBlock(0,1,loc_mat->GetBlock(0,1));
                new_loc_mat->SetBlock(1,0,loc_mat->GetBlock(1,0));
                new_loc_mat->SetBlock(1,1,loc_mat->GetBlock(1,1));
                loc_mat = new_loc_mat;
            }

            if(m_weakDirichletBCInfo.count(n) != 0) // add weak Dirichlet matrix
            {
               WeakDirichletBCInfoSharedPtr wDBC;

               /*
               DNekScalMatSharedPtr loc_fullmat = vExp->GetLocMatrix(matkey);
               */

               const DNekScalMatSharedPtr loc_fullmat =
                  vExp->as<LocalRegions::Expansion>()->GetLocMatrix(matkey);

               const int rows = loc_fullmat->GetRows();
               const int cols = loc_fullmat->GetColumns();

               const NekDouble *dat = loc_fullmat->GetRawPtr();
               DNekMatSharedPtr new_mat = MemoryManager<DNekMat>::
                   AllocateSharedPtr(rows,cols,dat);
               Blas::Dscal(rows*cols,loc_fullmat->Scale(),new_mat->GetRawPtr(),1);

               // add local matrix contribution
               int numBdryFacets = 0;
               for(wDBC = m_weakDirichletBCInfo.find(n)->second; wDBC; wDBC = wDBC->next)
               {
                   numBdryFacets++;
               }

               Array<OneD, int> edgeids(numBdryFacets);

               for(wDBC = m_weakDirichletBCInfo.find(n)->second; wDBC; wDBC = wDBC->next)
               {
                   edgeids[numBdryFacets] = wDBC->m_weakDirichletID;
                   numBdryFacets++;
               }

                vExp->AddWeakDirichletElementContribution(edgeids, *new_mat);

                // set up block matrix system
                unsigned int nbdry = vExp->NumBndryCoeffs();
                unsigned int nint = (unsigned int)(vExp->GetNcoeffs() - nbdry);
                unsigned int exp_size[] = {nbdry, nint};
                unsigned int nblks=2;
                loc_mat = MemoryManager<DNekScalBlkMat>::AllocateSharedPtr(nblks, nblks, exp_size, exp_size);
                {
                    int i,j;
                    const NekDouble factor = 1.0;
                    NekDouble            invfactor = 1.0/factor;
                    NekDouble            one = 1.0;
                    DNekMat &mat = *new_mat;
                    DNekMatSharedPtr A = MemoryManager<DNekMat>::AllocateSharedPtr(nbdry,nbdry);
                    DNekMatSharedPtr B = MemoryManager<DNekMat>::AllocateSharedPtr(nbdry,nint);
                    DNekMatSharedPtr C = MemoryManager<DNekMat>::AllocateSharedPtr(nint,nbdry);
                    DNekMatSharedPtr D = MemoryManager<DNekMat>::AllocateSharedPtr(nint,nint);

                    Array<OneD,unsigned int> bmap(nbdry);
                    Array<OneD,unsigned int> imap(nint);
                    vExp->GetBoundaryMap(bmap);
                    vExp->GetInteriorMap(imap);

                    for(i = 0; i < nbdry; ++i)
                    {
                        for(j = 0; j < nbdry; ++j)
                        {
                            (*A)(i,j) = mat(bmap[i],bmap[j]);
                        }

                        for(j = 0; j < nint; ++j)
                        {
                            (*B)(i,j) = mat(bmap[i],imap[j]);
                        }
                    }

                    for(i = 0; i < nint; ++i)
                    {
                        for(j = 0; j < nbdry; ++j)
                        {
                            (*C)(i,j) = mat(imap[i],bmap[j]);
                        }

                        for(j = 0; j < nint; ++j)
                        {
                            (*D)(i,j) = mat(imap[i],imap[j]);
                        }
                    }

                    // Calculate static condensed system
                    if(nint)
                    {
                        D->Invert();
                        (*B) = (*B)*(*D);
                        (*A) = (*A) - (*B)*(*C);
                    }

                    DNekScalMatSharedPtr     Atmp;

                    loc_mat->SetBlock(0,0,Atmp = MemoryManager<DNekScalMat>::AllocateSharedPtr(factor,A));
                    loc_mat->SetBlock(0,1,Atmp = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,B));
                    loc_mat->SetBlock(1,0,Atmp = MemoryManager<DNekScalMat>::AllocateSharedPtr(factor,C));
                    loc_mat->SetBlock(1,1,Atmp = MemoryManager<DNekScalMat>::AllocateSharedPtr(invfactor,D));

                }
            }

            return loc_mat;
        }

        /**
         * @brief Releases the static condensation block matrices from NekManager
         * of n-th expansion using the matrix key provided by the #m_linSysKey.
         * 
         * @param   n           Number of the expansion
         */
        void GlobalLinSys::v_DropStaticCondBlock(unsigned int n)
        {
            boost::shared_ptr<MultiRegions::ExpList> expList = m_expList.lock();

            StdRegions::StdExpansionSharedPtr vExp = expList->GetExp(n);

            // need to be initialised with zero size for non variable
            // coefficient case
            StdRegions::VarCoeffMap vVarCoeffMap;

            // retrieve variable coefficients
            if(m_linSysKey.GetNVarCoeffs() > 0)
            {
                StdRegions::VarCoeffMap::const_iterator x;
                int cnt = expList->GetPhys_Offset(n);
                for (x  = m_linSysKey.GetVarCoeffs().begin(); 
                     x != m_linSysKey.GetVarCoeffs().end  (); ++x)
                {
                    vVarCoeffMap[x->first] = x->second + cnt;
                }
            }

            LocalRegions::MatrixKey matkey(m_linSysKey.GetMatrixType(),
                                           vExp->DetShapeType(),
                                           *vExp,
                                           m_linSysKey.GetConstFactors(),
                                           vVarCoeffMap);

            vExp->DropLocStaticCondMatrix(matkey);
        }

        void GlobalLinSys::v_InitObject()
        {
            NEKERROR(ErrorUtil::efatal, "Method does not exist" );
	}

        void GlobalLinSys::v_Initialise(
            const boost::shared_ptr<AssemblyMap>& pLocToGloMap)
        {
            NEKERROR(ErrorUtil::efatal, "Method does not exist" );
	}
    } //end of namespace
} //end of namespace

