///////////////////////////////////////////////////////////////////////////////
//
// File CoupledLinearNS.cpp
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
// Description: Coupled  Solver for the Linearised Incompressible
// Navier Stokes equations
///////////////////////////////////////////////////////////////////////////////

#include <ViscoelasticSolver/EquationSystems/CoupledLinearNS.h>
#include <LibUtilities/BasicUtils/Timer.h>
#include <LocalRegions/MatrixKey.h>
#include <MultiRegions/GlobalLinSysDirectStaticCond.h>

namespace Nektar
{
    
    string CoupledLinearNS::className = SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction("CoupledLinearisedNS", CoupledLinearNS::create);
    
    
    /**
     *  @class CoupledLinearNS 
     *
     * Set up expansion field for velocity and pressure, the local to
     * global mapping arrays and the basic memory definitions for
     * coupled matrix system
     */ 
    CoupledLinearNS::CoupledLinearNS(const LibUtilities::SessionReaderSharedPtr &pSession):
        m_singleMode(false),
        m_zeroMode(false),
        LogConformation(pSession)
    {
    }

    void CoupledLinearNS::v_InitObject()
    {
        EquationSystem::v_InitObject();
        LogConformation::v_InitObject();

        int  n,i,j,k,eid;
        int  expdim = m_graph->GetMeshDimension();
        int  n_exp  = m_fields[m_velocity[0]]->GetNumElmts();
        int  nvel   = m_velocity.num_elements();

        // Get Expansion list for orthogonal expansion at p-2
        const SpatialDomains::ExpansionMap &pressure_exp = GenPressureExp(m_graph->GetExpansions("u"));

        if(m_equationType==eViscoelasticDEVSSG || m_equationType==eViscoelasticLogc || m_equationType==eViscoelasticDEVSSGALE)
        {
        	m_nConvectiveFields= m_velocity.num_elements() + m_stress.num_elements();
        }
        else
        {
        	m_nConvectiveFields = m_velocity.num_elements();
        }
        
        if(NoCaseStringCompare(m_boundaryConditions->GetVariable(m_nConvectiveFields-1),"p") == 0)
        {
            ASSERTL0(false,"Last field is defined as pressure but this is not suitable for this solver, please remove this field as it is implicitly defined");
        }
        // Decide how to declare explist for pressure. 
        if(expdim == 2)
        {
            int nz; 

            if(m_HomogeneousType == eHomogeneous1D)
            {

                ASSERTL0(nvel > 2,"Expect to have three at least three components of velocity variables");
                LibUtilities::BasisKey Homo1DKey = m_fields[0]->GetHomogeneousBasis()->GetBasisKey();
                
                m_pressure = MemoryManager<MultiRegions::ExpList3DHomogeneous1D>::AllocateSharedPtr(m_session, Homo1DKey, m_LhomZ, m_useFFT,m_dealiasing, pressure_exp);

                ASSERTL1(m_npointsZ%2==0,"Non binary number of planes have been specified");
                nz = m_npointsZ/2;                
            }
            else
            {
                m_pressure = MemoryManager<MultiRegions::ExpList2D>::AllocateSharedPtr(m_session, pressure_exp);
                nz = 1;
            }
            
            Array<OneD, MultiRegions::ExpListSharedPtr> velocity(m_velocity.num_elements());
            for(i =0 ; i < m_velocity.num_elements(); ++i)
            {
                velocity[i] = m_fields[m_velocity[i]];
            }

            // Set up Array of mappings
            m_locToGloMap = Array<OneD, CoupledLocalToGlobalC0ContMapSharedPtr> (nz);
            
            if(m_session->DefinesSolverInfo("SingleMode"))
            {

                ASSERTL0(nz <=2 ,"For single mode calculation can only have  nz <= 2");
                m_singleMode = true;
                if(m_session->DefinesSolverInfo("BetaZero"))
                {
                    m_zeroMode = true;
                }
                int nz_loc = 2;
                m_locToGloMap[0] = MemoryManager<CoupledLocalToGlobalC0ContMap>::AllocateSharedPtr(m_session,m_graph,m_boundaryConditions,velocity,m_pressure,nz_loc,m_zeroMode);
            }
            else 
            {
                // base mode
                int nz_loc = 1;
                m_locToGloMap[0] = MemoryManager<CoupledLocalToGlobalC0ContMap>::AllocateSharedPtr(m_session,m_graph,m_boundaryConditions,velocity,m_pressure,nz_loc);
                
                if(nz > 1)
                {
                    nz_loc = 2;
                    // Assume all higher modes have the same boundary conditions and re-use mapping
                    m_locToGloMap[1] = MemoryManager<CoupledLocalToGlobalC0ContMap>::AllocateSharedPtr(m_session,m_graph,m_boundaryConditions,velocity,m_pressure->GetPlane(2),nz_loc,false);
                    // Note high order modes cannot be singular
                    for(i = 2; i < nz; ++i)
                    {
                        m_locToGloMap[i] = m_locToGloMap[1];
                    }
                }
            }
        }
        else if (expdim == 3)
        {
            //m_pressure = MemoryManager<MultiRegions::ExpList3D>::AllocateSharedPtr(pressure_exp);
            ASSERTL0(false,"Setup mapping aray");
        }
        else
        {
            ASSERTL0(false,"Exp dimension not recognised");
        }
        
		switch(m_equationType)
		{
		 case eSteadyStokes:
		 case eSteadyOseen:
		 case eUnsteadyNavierStokes:
		 case eUnsteadyStokes:
			 // for dimensionless case
			 m_laplacefac = 1.0;
			 break;
		 case eUnsteadyViscoelastic:
		 case eViscoelasticLogc:
			 // for dimensionless case
			 m_laplacefac = m_beta;
			 break;
		 case eViscoelasticDEVSSG:
		 case eViscoelasticDEVSSGALE:
			 // for dimensionless case
			 m_laplacefac = m_beta + m_theta;
			 break;
		 case eNoEquationType:
		 default:
			 ASSERTL0(false,"Unknown or undefined equation type");
			 break;
		}
    }

    /**
     * Set up a coupled linearised Naviers-Stokes solve. Primarily a
     * wrapper fuction around the full mode by mode version of
     * SetUpCoupledMatrix.
     * 
     */
    void CoupledLinearNS::SetUpCoupledMatrix(const NekDouble lambda,  const Array< OneD, Array< OneD, NekDouble > > &Advfield, bool IsLinearNSEquation, bool AddGraduT)
   {
        
        int nz;
        if(m_singleMode)
	   {
            
            NekDouble lambda_imag; 
            
            // load imaginary componont of any potential shift
            // Probably should be called from DriverArpack but not yet
            // clear how to do this
            m_session->LoadParameter("imagShift",lambda_imag,NekConstants::kNekUnsetDouble);
            
            nz = 1;
		   m_mat  = Array<OneD, CoupledSolverMatrices> (nz);

		   ASSERTL1(m_npointsZ <=2,"Only expected a maxmimum of two planes in single mode linear NS solver");

		   if(m_zeroMode)
		   {
			   SetUpCoupledMatrix(lambda,Advfield,IsLinearNSEquation,AddGraduT,0,m_mat[0],m_locToGloMap[0],lambda_imag);
            }
		   else
		   {
			   NekDouble beta = 2*M_PI/m_LhomZ;
			   NekDouble lam = lambda + m_kinvis*beta*beta;

			   SetUpCoupledMatrix(lam,Advfield,IsLinearNSEquation,AddGraduT,1,m_mat[0],m_locToGloMap[0],lambda_imag);
            }
	   }
	   else
	   {
		   int n;
		   if(m_npointsZ > 1)
		   {
			   nz = m_npointsZ/2;
		   }
		   else
		   {
			   nz =  1;
		   }

		   m_mat  = Array<OneD, CoupledSolverMatrices> (nz);

		   // mean mode or 2D mode.
		   SetUpCoupledMatrix(lambda,Advfield,IsLinearNSEquation,AddGraduT,0,m_mat[0],m_locToGloMap[0]);

		   for(n = 1; n < nz; ++n)
		   {
			   NekDouble beta = 2*M_PI*n/m_LhomZ;

			   NekDouble lam = lambda + m_kinvis*beta*beta;

			   SetUpCoupledMatrix(lam,Advfield,IsLinearNSEquation,AddGraduT,n,m_mat[n],m_locToGloMap[n]);
		   }
	   }

   }

    /**
     * Set up a coupled linearised Naviers-Stokes solve in the
     * following manner: 
     * 
     * Consider the linearised Navier-Stokes element matrix denoted as
     * 
     * \f[ \left [ \begin{array}{ccc} A
     * & -D_{bnd}^T & B\\ -D_{bnd}& 0 & -D_{int}\\ \tilde{B}^T &
     * -D_{int}^T & C \end{array} \right ] \left [ \begin{array}{c} {\bf
     * v}_{bnd} \\ p\\ {\bf v}_{int} \end{array} \right ] = \left [
     * \begin{array}{c} {\bf f}_{bnd} \\ 0\\ {\bf f}_{int} \end{array}
     * \right ] \f] 
     * 
     * where \f${\bf v}_{bnd}, {\bf f}_{bnd}\f$ denote the degrees of
     * freedom of the elemental velocities on the boundary of the
     * element, \f${\bf v}_{int}, {\bf f}_{int}\f$ denote the degrees
     * of freedom of the elemental velocities on the interior of the
     * element and \f$p\f$ is the piecewise continuous pressure.  The
     * matrices have the interpretation
     * 
     * \f[  A[n,m]  = (\nabla \phi^b_n, \nu \nabla
     * \phi^b_m) + (\phi^b_n,{\bf U \cdot \nabla} \phi^b_m) +
     * (\phi^b_n \nabla^T {\bf U} \phi^b_m), \f] 
     * 
     * \f[ B[n,m] = (\nabla \phi^b_n, \nu \nabla \phi^i_m) +
     * (\phi^b_n,{\bf U \cdot \nabla} \phi^i_m) + (\phi^b_n \nabla^T
     * {\bf U} \phi^i_m),\f]
     * 
     * \f[\tilde{B}^T[n,m] = (\nabla \phi^i_n, \nu \nabla \phi^b_m) +
     * (\phi^i_n,{\bf U \cdot \nabla} \phi^b_m) + (\phi^i_n \nabla^T
     * {\bf U} \phi^b_m) \f]
     * 
     * \f[ C[n,m] = (\nabla \phi^i_n, \nu \nabla
     * \phi^i_m) + (\phi^i_n,{\bf U \cdot \nabla} \phi^i_m) +
     * (\phi^i_n \nabla^T {\bf U} \phi^i_m),\f]
     * 
     * \f[ D_{bnd}[n,m] = (\psi_m,\nabla \phi^b),\f]
     * 
     * \f[ D_{int}[n,m] = (\psi_m,\nabla \phi^i) \f]
     * 
     * where \f$\psi\f$ is the space of pressures typically at order
     * \f$P-2\f$ and \f$\phi\f$ is the velocity vector space of
     * polynomial order \f$P\f$. (Note the above definitions for the
     * transpose terms shoudl be interpreted with care since we have
     * used a tensor differential for the \f$\nabla^T \f$ terms)
     * 
     * Note \f$B = \tilde{B}^T\f$ if just considering the stokes
     * operator and then \f$C\f$ is also symmetric. Also note that
     * \f$A,B\f$ and \f$C\f$ are block diagonal in the Oseen equation
     * when \f$\nabla^T {\bf U}\f$ are zero.
     * 
     * Since \f$C\f$ is invertible we can premultiply the governing
     * elemental equation by the following matrix: 
     * 
     * \f[  \left [ \begin{array}{ccc} I & 0 &
     * -BC^{-1}\\ 0& I & D_{int}C^{-1}\\ 0 & 0 & I \end{array}
     * \right ] \left \{ \left [ \begin{array}{ccc} A & -D_{bnd}^T &
     * B\\ -D_{bnd}& 0 & -D_{int}\\ \tilde{B}^T & -D_{int}^T & C
     * \end{array} \right ] \left [ \begin{array}{c} {\bf v}_{bnd} \\
     * p\\ {\bf v_{int}} \end{array} \right ] = \left [
     * \begin{array}{c} {\bf f}_{bnd} \\ 0\\ {\bf f_{int}} \end{array}
     * \right ] \right \} \f]
     * 
     *  which if we multiply out the matrix equation we obtain
     * 
     * \f[ \left [ \begin{array}{ccc} A - B
     * C^{-1}\tilde{B}^T & -D_{bnd}^T +B C^{-1} D_{int}^T& 0\\
     * -D_{bnd}+D_{int} C^{-1} \tilde{B}^T & -D_{int} C^{-1}
     * D_{int}^T & 0\\ \tilde{B}^T & -D_{int}^T & C \end{array} \right ]
     * \left [ \begin{array}{c} {\bf v}_{bnd} \\ p\\ {\bf v_{int}}
     * \end{array} \right ] = \left [ \begin{array}{c} {\bf f}_{bnd}
     * -B C^{-1} {\bf f}_{int}\\ f_p = D_{int} C^{-1} {\bf
     * f}_{int}\\ {\bf f_{int}} \end{array} \right ] \f]
     *
     * In the above equation the \f${\bf v}_{int}\f$ degrees of freedom
     * are decoupled and so we need to solve for the \f${\bf v}_{bnd},
     * p\f$ degrees of freedom.  The final step is to perform a second
     * level of static condensation but where we will lump the mean
     * pressure mode (or a pressure degree of freedom containing a
     * mean component) with the velocity boundary degrees of
     * freedom. To do we define \f${\bf b} = [{\bf v}_{bnd}, p_0]\f$ where
     * \f$p_0\f$ is the mean pressure mode and \f$\hat{p}\f$ to be the
     * remainder of the pressure space.  We now repartition the top \f$2
     * \times 2\f$ block of matrices of previous matrix equation  as
     * 
     * \f[ \left [ \begin{array}{cc} \hat{A} & \hat{B}\\ \hat{C} &
     * \hat{D} \end{array} \right ] \left [ \begin{array}{c} {\bf b}
     * \\ \hat{p} \end{array} \right ] = \left [ \begin{array}{c}
     * \hat{\bf f}_{bnd} \\ \hat{f}_p \end{array} \right ]
     * \label{eqn.linNS_fac2} \f]
     * 
     * where 
     * 
     * \f[ \hat{A}_{ij} = \left [ \begin{array}{cc} A - B
     * C^{-1}\tilde{B}^T & [-D_{bnd}^T +B C^{-1} D_{int}^T]_{i0}\\
     * {[}-D_{bnd}+D_{int} C^{-1} \tilde{B}^T]_{0j} & -[D_{int}
     * C^{-1} D_{int}^T ]_{00} \end{array} \right ] \f] 
     * 
     * \f[ \hat{B}_{ij} = \left [ \begin{array}{c} [-D_{bnd}^T +B
     * C^{-1} D_{int}^T]_{i,j+1} \\ {[} -D_{int} C^{-1} D^T_{int}
     * ]_{0j}\end{array} \right ] \f]
     * 
     * \f[ \hat{C}_{ij}  =  \left [\begin{array}{cc} -D_{bnd} + 
     * D_{int} C^{-1} \tilde{B}^T, & {[} -D_{int} C^{-1} D^T_{int}
     * ]_{i+1,0}\end{array} \right ] \f]
     * 
     * \f[ \hat{D}_{ij}  =  \left [\begin{array}{c} {[} -D_{int}
     * C^{-1} D^T_{int} ]_{i+1,j+1}\end{array} \right ] \f]
     *
     * and 
     *
     * \f[ fh\_{bnd} = \left [ \begin{array}{c} {\bf
     * f}_{bnd} -B C^{-1} {\bf f}_{int}\\ {[}D_{int} C^{-1} {\bf
     * f}_{int}]_0 \end{array}\right ] \hspace{1cm} [fh\_p_{i} =
     * \left [ \begin{array}{c} {[}D_{int} C^{-1} {\bf f}_{iint}]_{i+1}
     * \end{array}\right ] \f]
     * 
     * Since the \f$\hat{D}\f$ is decoupled and invertible we can now
     * statically condense the previous matrix equationto decouple
     * \f${\bf b}\f$ from \f$\hat{p}\f$ by solving the following system
     * 
     * \f[ \left [ \begin{array}{cc} \hat{A} - \hat{B} \hat{D}^{-1}
     * \hat{C} & 0 \\ \hat{C} & \hat{D} \end{array} \right ] \left
     * [ \begin{array}{c} {\bf b} \\ \hat{p} \end{array} \right ] =
     * \left [ \begin{array}{c} \hat{\bf f}_{bnd} - \hat{B}
     * \hat{D}^{-1} \hat{f}_p\\ \hat{f}_p \end{array} \right ] \f] 
     * 
     * The matrix \f$\hat{A} - \hat{B} \hat{D}^{-1} \hat{C}\f$ has
     * to be globally assembled and solved iteratively or
     * directly. One we obtain the solution to \f${\bf b}\f$ we can use
     * the second row of equation fourth matrix equation to solve for
     * \f$\hat{p}\f$ and finally the last row of equation second
     * matrix equation to solve for \f${\bf v}_{int}\f$.
     *
     */

    void CoupledLinearNS::SetUpCoupledMatrix(const NekDouble lambda,  const Array< OneD, Array< OneD, NekDouble > > &Advfield, bool IsLinearNSEquation
      		 ,bool AddGraduT,const int HomogeneousMode, CoupledSolverMatrices &mat, CoupledLocalToGlobalC0ContMapSharedPtr &locToGloMap, const NekDouble lambda_imag)
  	 {
  		 int  n,i,j,k,eid;
  		 int  expdim = m_graph->GetMeshDimension();
  		 int  nel  = m_fields[m_velocity[0]]->GetNumElmts();
  		 int  nvel   = m_velocity.num_elements();
  		 static int nofcalls=1;

  		 // if Advfield is defined can assume it is an Oseen or LinearNS equation
  		 bool AddAdvectionTerms = (Advfield ==  NullNekDoubleArrayofArray)? false: true;
  		 //bool AddAdvectionTerms = true; // Temporary debugging trip

  		 // call velocity space Static condensation and fill block
  		 // matrices.  Need to set this up differently for Oseen and
  		 // Lin NS.  Ideally should make block diagonal for Stokes and
  		 // Oseen problems.
  		 DNekScalMatSharedPtr loc_mat;
  		 StdRegions::StdExpansionSharedPtr locExp;
  		 NekDouble one = 1.0;
  		 int nint,nbndry;
  		 int rows, cols;
  		 NekDouble zero = 0.0;
  		 Array<OneD, unsigned int> bmap,imap;

  		 Array<OneD,unsigned int> nsize_bndry   (nel);
  		 Array<OneD,unsigned int> nsize_bndry_p1(nel);
  		 Array<OneD,unsigned int> nsize_int     (nel);
  		 Array<OneD,unsigned int> nsize_p       (nel);
  		 Array<OneD,unsigned int> nsize_p_m1    (nel);

  		 int nz_loc;

  		 if(HomogeneousMode) // Homogeneous mode flag
  		 {
  			 nz_loc = 2;
  		 }
  		 else
  		 {
  			 if(m_singleMode)
  			 {
  				 nz_loc = 2;
  			 }
  			 else
  			 {
  				 nz_loc = 1;
  			 }
  		 }

  		 // Set up block matrix sizes -
  		 for(n = 0; n < nel; ++n)
  		 {
  			 eid = m_fields[m_velocity[0]]->GetOffset_Elmt_Id(n);
  			 nsize_bndry[n] = nvel*m_fields[m_velocity[0]]->GetExp(eid)->NumBndryCoeffs()*nz_loc;
  			 nsize_bndry_p1[n] = nsize_bndry[n]+nz_loc;
  			 nsize_int  [n] = (nvel*m_fields[m_velocity[0]]->GetExp(eid)->GetNcoeffs()*nz_loc - nsize_bndry[n]);
  			 nsize_p[n] = m_pressure->GetExp(eid)->GetNcoeffs()*nz_loc;
  			 nsize_p_m1[n] = nsize_p[n]-nz_loc;
  		 }

  		 MatrixStorage blkmatStorage = eDIAGONAL;
  		 if(nofcalls==1)
  		 {

  		 mat.m_BCinv = MemoryManager<DNekScalBlkMat>
  			 ::AllocateSharedPtr(nsize_bndry,nsize_int,blkmatStorage);
  		 mat.m_Btilde = MemoryManager<DNekScalBlkMat>
  			 ::AllocateSharedPtr(nsize_bndry,nsize_int,blkmatStorage);
  		 mat.m_Cinv = MemoryManager<DNekScalBlkMat>
  			 ::AllocateSharedPtr(nsize_int,nsize_int,blkmatStorage);

  		 mat.m_D_bnd = MemoryManager<DNekScalBlkMat>
  			 ::AllocateSharedPtr(nsize_p,nsize_bndry,blkmatStorage);

  		 mat.m_D_int = MemoryManager<DNekScalBlkMat>
  			 ::AllocateSharedPtr(nsize_p,nsize_int,blkmatStorage);
  		 }

  		 DNekScalBlkMatSharedPtr pAh = MemoryManager<DNekScalBlkMat>
  			 ::AllocateSharedPtr(nsize_bndry_p1,nsize_bndry_p1,blkmatStorage);
  		 // Final level static condensation matrices.
  		 DNekScalBlkMatSharedPtr pBh = MemoryManager<DNekScalBlkMat>
  			 ::AllocateSharedPtr(nsize_bndry_p1,nsize_p_m1,blkmatStorage);
  		 DNekScalBlkMatSharedPtr pCh = MemoryManager<DNekScalBlkMat>
  			 ::AllocateSharedPtr(nsize_p_m1,nsize_bndry_p1,blkmatStorage);
  		 DNekScalBlkMatSharedPtr pDh = MemoryManager<DNekScalBlkMat>
  			 ::AllocateSharedPtr(nsize_p_m1,nsize_p_m1,blkmatStorage);


  		 for(n = 0; n < nel; ++n)
  		 {
  			 eid = m_fields[m_velocity[0]]->GetOffset_Elmt_Id(n);
  			 nbndry = nsize_bndry[n];
  			 nint   = nsize_int[n];
  			 k = nsize_bndry_p1[n];

  			 DNekMatSharedPtr Ah = MemoryManager<DNekMat>::AllocateSharedPtr(k,k,zero);
  			 DNekMatSharedPtr B  = MemoryManager<DNekMat>::AllocateSharedPtr(nbndry,nint,zero);
  			 DNekMatSharedPtr C  = MemoryManager<DNekMat>::AllocateSharedPtr(nbndry,nint,zero);
  			 DNekMatSharedPtr D  = MemoryManager<DNekMat>::AllocateSharedPtr(nint, nint,zero);

  			 DNekMatSharedPtr Dbnd = MemoryManager<DNekMat>::AllocateSharedPtr(nsize_p[n],nsize_bndry[n],zero);
  			 DNekMatSharedPtr Dint = MemoryManager<DNekMat>::AllocateSharedPtr(nsize_p[n],nsize_int[n],zero);

  			 locExp = m_fields[m_velocity[0]]->GetExp(eid);
  			 locExp->GetBoundaryMap(bmap);
  			 locExp->GetInteriorMap(imap);
             StdRegions::ConstFactorMap factors;
             factors[StdRegions::eFactorLambda] = (lambda*m_reynolds)/(m_laplacefac);
  			 LocalRegions::MatrixKey helmkey(StdRegions::eHelmholtz,
  					 	 	 	 	 	 	 locExp->DetShapeType(),
  											 *locExp,
  											 factors);


  			 int ncoeffs = m_fields[m_velocity[0]]->GetExp(eid)->GetNcoeffs();
  			 int nbmap = bmap.num_elements();
  			 int nimap = imap.num_elements();

  			 Array<OneD, NekDouble> coeffs  = m_fields[m_velocity[0]]->GetExp(eid)->UpdateCoeffs();
  			 Array<OneD, NekDouble> phys    = m_fields[m_velocity[0]]->GetExp(eid)->UpdatePhys();
  			 int psize   = m_pressure->GetExp(eid)->GetNcoeffs();

  			 Array<OneD, NekDouble> deriv   = m_pressure->GetExp(eid)->UpdatePhys();
  			 Array<OneD, NekDouble> pcoeffs = m_pressure->GetExp(eid)->UpdateCoeffs();

  			 if(AddAdvectionTerms == false && AddGraduT == false) // use static condensed managed matrices
  			 {
  				 // construct velocity matrices using statically
  				 // condensed elemental matrices and then construct
  				 // pressure matrix systems
  				 DNekScalBlkMatSharedPtr CondMat;
  				 CondMat = locExp->GetLocStaticCondMatrix(helmkey);

  				 for(k = 0; k < nvel*nz_loc; ++k)
  				 {
  					 DNekScalMat &SubBlock = *CondMat->GetBlock(0,0);
  					 rows = SubBlock.GetRows();
  					 cols = SubBlock.GetColumns();
  					 for(i = 0; i < rows; ++i)
  					 {
  						 for(j = 0; j < cols; ++j)
  						 {
  							 (*Ah)(i+k*rows,j+k*cols) = m_laplacefac*SubBlock(i,j);
  						 }
  					 }
  				 }

  				 for(k = 0; k < nvel*nz_loc; ++k)
  				 {
  					 DNekScalMat &SubBlock  = *CondMat->GetBlock(0,1);
  					 DNekScalMat &SubBlock1 = *CondMat->GetBlock(1,0);
  					 rows = SubBlock.GetRows();
  					 cols = SubBlock.GetColumns();
  					 for(i = 0; i < rows; ++i)
  					 {
  						 for(j = 0; j < cols; ++j)
  						 {
  							 (*B)(i+k*rows,j+k*cols) = SubBlock(i,j);
  													(*C)(i+k*rows,j+k*cols) = m_laplacefac*SubBlock1(j,i);
    }
  					 }
  				 }

  				 for(k = 0; k < nvel*nz_loc; ++k)
  				 {
  					 DNekScalMat &SubBlock = *CondMat->GetBlock(1,1);
  					 NekDouble inv_kinvis = 1.0/m_laplacefac;
  					 rows = SubBlock.GetRows();
  					 cols = SubBlock.GetColumns();
  					 for(i = 0; i < rows; ++i)
  					 {
  						 for(j = 0; j < cols; ++j)
  						 {
  							 (*D)(i+k*rows,j+k*cols) = inv_kinvis*SubBlock(i,j);
  						 }
  					 }
  				 }


  				 // Loop over pressure space and construct boundary block matrices.
  				 for(i = 0; i < bmap.num_elements(); ++i)
  				 {
  					 // Fill element with mode
  					 Vmath::Zero(ncoeffs,coeffs,1);
  					 coeffs[bmap[i]] = 1.0;
  					 m_fields[m_velocity[0]]->GetExp(eid)->BwdTrans(coeffs,phys);

  					 // Differentiation & Inner product wrt base.
  					 for(j = 0; j < nvel; ++j)
  					 {
  						 if( (nz_loc == 2)&&(j == 2)) // handle d/dz derivative
  						 {
  							 NekDouble beta =  2*M_PI*HomogeneousMode/m_LhomZ;

  							 Vmath::Smul(m_fields[m_velocity[0]]->GetExp(eid)->GetTotPoints(), beta, phys,1,deriv,1);

  							 m_pressure->GetExp(eid)->IProductWRTBase(deriv,pcoeffs);

  							 Blas::Dcopy(psize,&(pcoeffs)[0],1,
  										 Dbnd->GetRawPtr() +
  										 ((nz_loc*j+1)*bmap.num_elements()+i)*nsize_p[n],1);

  							 Vmath::Neg(psize,pcoeffs,1);
  							 Blas::Dcopy(psize,&(pcoeffs)[0],1,
  										 Dbnd->GetRawPtr() +
  										 ((nz_loc*j)*bmap.num_elements()+i)*nsize_p[n]+psize,1);

  						 }
  						 else
  						 {
  							 if(j < 2) // required for mean mode of homogeneous expansion
  							 {
  				locExp->PhysDeriv(MultiRegions::DirCartesianMap[j],phys,deriv);
  								 m_pressure->GetExp(eid)->IProductWRTBase(deriv,pcoeffs);
  								 // copy into column major storage.
  								 for(k = 0; k < nz_loc; ++k)
  								 {
  									 Blas::Dcopy(psize,&(pcoeffs)[0],1,
  												 Dbnd->GetRawPtr() +
  												 ((nz_loc*j+k)*bmap.num_elements()+i)*nsize_p[n]+ k*psize,1);
  								 }
  							 }
  						 }
  					 }
  				 }

  				 for(i = 0; i < imap.num_elements(); ++i)
  				 {
  					 // Fill element with mode
  					 Vmath::Zero(ncoeffs,coeffs,1);
  					 coeffs[imap[i]] = 1.0;
  					 m_fields[m_velocity[0]]->GetExp(eid)->BwdTrans(coeffs,phys);

  					 // Differentiation & Inner product wrt base.
  					 for(j = 0; j < nvel; ++j)
  					 {
  						 if( (nz_loc == 2)&&(j == 2)) // handle d/dz derivative
  						 {

  							 NekDouble beta = 2*M_PI*HomogeneousMode/m_LhomZ;

  							 Vmath::Smul(m_fields[m_velocity[0]]->GetExp(eid)->GetTotPoints(), beta, phys,1,deriv,1);

  							 m_pressure->GetExp(eid)->IProductWRTBase(deriv,pcoeffs);

  							 Blas::Dcopy(psize,&(pcoeffs)[0],1,
  										 Dint->GetRawPtr() +
  										 ((nz_loc*j+1)*imap.num_elements()+i)*nsize_p[n],1);

  							 Vmath::Neg(psize,pcoeffs,1);
  							 Blas::Dcopy(psize,&(pcoeffs)[0],1,
  										 Dint->GetRawPtr() +
  										 ((nz_loc*j)*imap.num_elements()+i)*nsize_p[n]+psize,1);

  						 }
  						 else
  						 {
  							 if(j < 2) // required for mean mode of homogeneous expansion
  							 {
  								 //m_fields[m_velocity[0]]->GetExp(eid)->PhysDeriv(j,phys, deriv);
  								 locExp->PhysDeriv(MultiRegions::DirCartesianMap[j],phys,deriv);

  								 m_pressure->GetExp(eid)->IProductWRTBase(deriv,pcoeffs);

  								 // copy into column major storage.
  								 for(k = 0; k < nz_loc; ++k)
  								 {
  									 Blas::Dcopy(psize,&(pcoeffs)[0],1,
  												 Dint->GetRawPtr() +
  												 ((nz_loc*j+k)*imap.num_elements()+i)*nsize_p[n]+ k*psize,1);
  								 }
  							 }
  						 }
  					 }
  				 }

   #if 0
  			 if(nz_loc  == 2)
  			 {
  				 //cout << "Ah" << endl;
  				 //cout << *Ah << endl;
  				 cout << "B" << endl;
  				 cout << *B << endl;
  				 cout << "C" << endl;
  				 cout << *C << endl;
  				 cout << "D" << endl;
  				 cout << *D << endl;
  				 cout << "Dbnd" << endl;
  				 cout << *Dbnd << endl;
  				 cout << "Dint" << endl;
  				 cout << *Dint << endl;
  				 exit(1);
  			 }
   #endif

  			 }
  			 else //include Advection and GraduT
  			 {
  				 // construct velocity matrices and pressure systems at
  				 // the same time resusing differential of velocity
  				 // space

  				 DNekScalMat &HelmMat = *(boost::dynamic_pointer_cast<LocalRegions::Expansion>(locExp)->GetLocMatrix(helmkey));

  				 //if(AddGraduT)
  				 //{
  					 LocalRegions::MatrixKey lap00(StdRegions::eLaplacian00,
  																locExp->DetShapeType(),
  																*locExp);

  					 LocalRegions::MatrixKey lap11(StdRegions::eLaplacian11,
  													locExp->DetShapeType(),
  													*locExp);

  					 DNekScalMat &Lap00 = *(boost::dynamic_pointer_cast<LocalRegions::Expansion>(locExp)->GetLocMatrix(lap00));
  					 DNekScalMat &Lap11 = *(boost::dynamic_pointer_cast<LocalRegions::Expansion>(locExp)->GetLocMatrix(lap11));
  				 //}

  				 Array<OneD, NekDouble> Advtmp;
  				 Array<OneD, Array<OneD, NekDouble> > AdvDeriv(nvel*nvel);
  				 // Use ExpList phys array for temporaary storage
  				 Array<OneD, NekDouble> tmpphys = m_fields[0]->UpdatePhys();
  				 int phys_offset = m_fields[m_velocity[0]]->GetPhys_Offset(eid);
  				 int nv;
  				 int npoints = locExp->GetTotPoints();

  				 // Calculate derivative of base flow
  				 if(IsLinearNSEquation)
  				 {
  					 int nv1;
  					 int cnt = 0;
  					 AdvDeriv[0] = Array<OneD, NekDouble>(nvel*nvel*npoints);
  					 for(nv = 0; nv < nvel; ++nv)
  					 {
  						 for(nv1 = 0; nv1 < nvel; ++nv1)
  						 {
  							 if(cnt < nvel*nvel-1)
  							 {
  								 AdvDeriv[cnt+1] = AdvDeriv[cnt] + npoints;
  								 ++cnt;
  							 }

  							 if((nv1 == 2)&&(m_HomogeneousType == eHomogeneous1D))
  							 {
  								 Vmath::Zero(npoints,AdvDeriv[nv*nvel+nv1],1); // dU/dz = 0
  							 }
  							 else
  							 {
  								 locExp->PhysDeriv(MultiRegions::DirCartesianMap[nv1],Advfield[nv] + phys_offset, AdvDeriv[nv*nvel+nv1]);
  							 }
  						 }
  					 }
  				 }


  				 for(i = 0; i < nbmap; ++i)
  				 {

  					 // Fill element with mode
  					 Vmath::Zero(ncoeffs,coeffs,1);
  					 coeffs[bmap[i]] = 1.0;
  					 locExp->BwdTrans(coeffs,phys);

  					 for(k = 0; k < nvel*nz_loc; ++k)
  					 {
  						 for(j = 0; j < nbmap; ++j)
  						 {
  							 (*Ah)(i+k*nbmap,j+k*nbmap) += m_laplacefac*HelmMat(bmap[i],bmap[j]);

  							 if(AddGraduT)
  							 {
  							   if(k == 0)
  							   {
  								   (*Ah)(i+k*nbmap,j+k*nbmap) += m_laplacefac*Lap00(bmap[i],bmap[j]);
  							   }
  							   else if(k == 1)
  							   {
  								   (*Ah)(i+k*nbmap,j+k*nbmap) += m_laplacefac*Lap11(bmap[i],bmap[j]);
  							   }
  							 }
  						 }

  						 for(j = 0; j < nimap; ++j)
  						 {
  							 (*B)(i+k*nbmap,j+k*nimap) += m_laplacefac*HelmMat(bmap[i],imap[j]);

  							 if(AddGraduT)
  							 {
  							   if(k == 0)
  							   {
  								   (*B)(i+k*nbmap,j+k*nimap) += m_laplacefac*Lap00(bmap[i],imap[j]);
  							   }
  							   else if(k == 1)
  							   {
  								   (*B)(i+k*nbmap,j+k*nimap) += m_laplacefac*Lap11(bmap[i],imap[j]);
  							   }
  							 }
  						 }
  					 }

  					 for(k = 0; k < nvel; ++k)
  					 {
  						 if((nz_loc == 2)&&(k == 2)) // handle d/dz derivative
  						 {
  							 NekDouble beta = 2*M_PI*HomogeneousMode/m_LhomZ;

  							 // Real Component
  							 Vmath::Smul(npoints,beta,phys,1,deriv,1);

  							 m_pressure->GetExp(eid)->IProductWRTBase(deriv,pcoeffs);
  							 Blas::Dcopy(psize,&(pcoeffs)[0],1,
  										 Dbnd->GetRawPtr() +
  										 ((nz_loc*k+1)*bmap.num_elements()+i)*nsize_p[n],1);

  							 // Imaginary Component
  							 Vmath::Neg(psize,pcoeffs,1);
  							 Blas::Dcopy(psize,&(pcoeffs)[0],1,
  										 Dbnd->GetRawPtr() +
  										 ((nz_loc*k)*bmap.num_elements()+i)*nsize_p[n]+psize,1);

  							 // now do advection terms
  							 Vmath::Vmul(npoints,
  										 Advtmp = Advfield[k] + phys_offset,
  										 1,deriv,1,tmpphys,1);

  							 locExp->IProductWRTBase(tmpphys,coeffs);


  							 // real contribution
  							 for(nv = 0; nv < nvel; ++nv)
  							 {
  								 for(j = 0; j < nbmap; ++j)
  								 {
  									 (*Ah)(j+2*nv*nbmap,i+(2*nv+1)*nbmap) +=
  										 m_reynolds*coeffs[bmap[j]];
  								 }

  								 for(j = 0; j < nimap; ++j)
  								 {
  									 (*C)(i+(2*nv+1)*nbmap,j+2*nv*nimap) +=
  										m_reynolds*coeffs[imap[j]];
  								 }
  							 }

  							 Vmath::Neg(ncoeffs,coeffs,1);
  							 // imaginary contribution
  							 for(nv = 0; nv < nvel; ++nv)
  							 {
  								 for(j = 0; j < nbmap; ++j)
  								 {
  									 (*Ah)(j+(2*nv+1)*nbmap,i+2*nv*nbmap) +=
  										m_reynolds*coeffs[bmap[j]];
  								 }

  								 for(j = 0; j < nimap; ++j)
  								 {
  									 (*C)(i+2*nv*nbmap,j+(2*nv+1)*nimap) +=
  									   m_reynolds*coeffs[imap[j]];
  								 }
  							 }
  						 }
  						 else
  						 {
  							 if(k < 2)
  							 {
  								 //locExp->PhysDeriv(k,phys, deriv);
  								 locExp->PhysDeriv(MultiRegions::DirCartesianMap[k],phys, deriv);
  								 Vmath::Vmul(npoints,
  											 Advtmp = Advfield[k] + phys_offset,
  											 1,deriv,1,tmpphys,1);
  								 locExp->IProductWRTBase(tmpphys,coeffs);


  								 for(nv = 0; nv < nvel*nz_loc; ++nv)
  								 {
  									 for(j = 0; j < nbmap; ++j)
  									 {
  										 (*Ah)(j+nv*nbmap,i+nv*nbmap) +=
  												 m_reynolds*coeffs[bmap[j]];
  									 }

  									 for(j = 0; j < nimap; ++j)
  									 {
  										 (*C)(i+nv*nbmap,j+nv*nimap) +=
  												 m_reynolds*coeffs[imap[j]];
  									 }
  								 }
  								 // copy into column major storage.
  								 m_pressure->GetExp(eid)->IProductWRTBase(deriv,pcoeffs);
  								 for(j = 0; j < nz_loc; ++j)
  								 {
  									 Blas::Dcopy(psize,&(pcoeffs)[0],1,
  												 Dbnd->GetRawPtr() +
  												 ((nz_loc*k+j)*bmap.num_elements() + i)*nsize_p[n]+ j*psize,1);
  								 }
  							 }
  						 }

  						 if(IsLinearNSEquation)
  						 {
  							 for(nv = 0; nv < nvel; ++nv)
  							 {
  								 // u' . Grad U terms
  								 Vmath::Vmul(npoints,phys,1,
  											 AdvDeriv[k*nvel+nv],
  											 1,tmpphys,1);
  								 locExp->IProductWRTBase(tmpphys,coeffs);

  								 for(int n1 = 0; n1 < nz_loc; ++n1)
  								 {
  									 for(j = 0; j < nbmap; ++j)
  									 {
  										 (*Ah)(j+(k*nz_loc+n1)*nbmap,
  											   i+(nv*nz_loc+n1)*nbmap) +=
  													   m_reynolds*coeffs[bmap[j]];
  									 }

  									 for(j = 0; j < nimap; ++j)
  									 {
  										 (*C)(i+(nv*nz_loc+n1)*nbmap,
  											  j+(k*nz_loc+n1)*nimap) +=
  													  m_reynolds*coeffs[imap[j]];
  									 }
  								 }
  							 }
  						 }
  					 }
  				 }


  				 for(i = 0; i < nimap; ++i)
  				 {
  					 // Fill element with mode
  					 Vmath::Zero(ncoeffs,coeffs,1);
  					 coeffs[imap[i]] = 1.0;
  					 locExp->BwdTrans(coeffs,phys);

  					 for(k = 0; k < nvel*nz_loc; ++k)
  					 {
  						 for(j = 0; j < nbmap; ++j) // C set up as transpose
  						 {
  							 (*C)(j+k*nbmap,i+k*nimap) += m_laplacefac*HelmMat(imap[i],bmap[j]);
  							 if(AddGraduT)
  							 {
  							   if(k == 0)
  							   {
  								   (*C)(j+k*nbmap,i+k*nimap) += m_laplacefac*Lap00(imap[i],bmap[j]);
  							   }
  							   else if(k == 1)
  							   {
  								   (*C)(j+k*nbmap,i+k*nimap) += m_laplacefac*Lap11(imap[i],bmap[j]);
  							   }
  							 }
  						 }

  						 for(j = 0; j < nimap; ++j)
  						 {
  							 (*D)(i+k*nimap,j+k*nimap) += m_laplacefac*HelmMat(imap[i],imap[j]);
  							 if(AddGraduT)
  							 {
  							   if(k == 0)
  							   {
  								   (*D)(i+k*nimap,j+k*nimap) += m_laplacefac*Lap00(imap[i],imap[j]);
  							   }
  							   else if(k == 1)
  							   {
  								   (*D)(i+k*nimap,j+k*nimap) += m_laplacefac*Lap11(imap[i],imap[j]);
  							   }
  							 }
  						 }
  					 }


  					 for(k = 0; k < nvel; ++k)
  					 {
  						 if((nz_loc == 2)&&(k == 2)) // handle d/dz derivative
  						 {
  							 NekDouble beta = 2*M_PI*HomogeneousMode/m_LhomZ;

  							 // Real Component
  							 Vmath::Smul(npoints,beta,phys,1,deriv,1);
  							 m_pressure->GetExp(eid)->IProductWRTBase(deriv,pcoeffs);
  							 Blas::Dcopy(psize,&(pcoeffs)[0],1,
  										 Dint->GetRawPtr() +
  										 ((nz_loc*k+1)*imap.num_elements()+i)*nsize_p[n],1);
  							 // Imaginary Component
  							 Vmath::Neg(psize,pcoeffs,1);
  							 Blas::Dcopy(psize,&(pcoeffs)[0],1,
  										 Dint->GetRawPtr() +
  										 ((nz_loc*k)*imap.num_elements()+i)*nsize_p[n]+psize,1);

  							 // Advfield[k] *d/dx_k to all velocity
  							 // components on diagonal
  							 Vmath::Vmul(npoints, Advtmp = Advfield[k] + phys_offset,1,deriv,1,tmpphys,1);
  							 locExp->IProductWRTBase(tmpphys,coeffs);

  							 // Real Components
  							 for(nv = 0; nv < nvel; ++nv)
  							 {
  								 for(j = 0; j < nbmap; ++j)
  								 {
  									 (*B)(j+2*nv*nbmap,i+(2*nv+1)*nimap) +=
  											 m_reynolds*coeffs[bmap[j]];
  								 }

  								 for(j = 0; j < nimap; ++j)
  								 {
  									 (*D)(j+2*nv*nimap,i+(2*nv+1)*nimap) +=
  											 m_reynolds*coeffs[imap[j]];
  								 }
  							 }
  							 Vmath::Neg(ncoeffs,coeffs,1);
  							 // Imaginary
  							 for(nv = 0; nv < nvel; ++nv)
  							 {
  								 for(j = 0; j < nbmap; ++j)
  								 {
  									 (*B)(j+(2*nv+1)*nbmap,i+2*nv*nimap) +=
  											 m_reynolds*coeffs[bmap[j]];
  								 }

  								 for(j = 0; j < nimap; ++j)
  								 {
  									 (*D)(j+(2*nv+1)*nimap,i+2*nv*nimap) +=
  											 m_reynolds*coeffs[imap[j]];
  								 }
  							 }

  						 }
  						 else
  						 {
  							 if(k < 2)
  							 {
  								 // Differentiation & Inner product wrt base.
  								 //locExp->PhysDeriv(k,phys, deriv);
  								 locExp->PhysDeriv(MultiRegions::DirCartesianMap[k],phys,deriv);
  								 Vmath::Vmul(npoints,
  											 Advtmp = Advfield[k] + phys_offset,
  											 1,deriv,1,tmpphys,1);
  								 locExp->IProductWRTBase(tmpphys,coeffs);


  								 for(nv = 0; nv < nvel*nz_loc; ++nv)
  								 {
  									 for(j = 0; j < nbmap; ++j)
  									 {
  										 (*B)(j+nv*nbmap,i+nv*nimap) +=
  												 m_reynolds*coeffs[bmap[j]];
  									 }

  									 for(j = 0; j < nimap; ++j)
  									 {
  										 (*D)(j+nv*nimap,i+nv*nimap) +=
  												 m_reynolds*coeffs[imap[j]];
  									 }
  								 }
  								 // copy into column major storage.
  								 m_pressure->GetExp(eid)->IProductWRTBase(deriv,pcoeffs);
  								 for(j = 0; j < nz_loc; ++j)
  								 {
  									 Blas::Dcopy(psize,&(pcoeffs)[0],1,
  												 Dint->GetRawPtr() +
  												 ((nz_loc*k+j)*imap.num_elements() + i)*nsize_p[n]+j*psize,1);
  								 }
  							 }
  						 }

  						 if(IsLinearNSEquation)
  						 {
  							 int n1;
  							 for(nv = 0; nv < nvel; ++nv)
  							 {
  								 // u'.Grad U terms
  								 Vmath::Vmul(npoints,phys,1,
  											 AdvDeriv[k*nvel+nv],
  											 1,tmpphys,1);
  								 locExp->IProductWRTBase(tmpphys,coeffs);

  								 for(n1 = 0; n1 < nz_loc; ++n1)
  								 {
  									 for(j = 0; j < nbmap; ++j)
  									 {
  										 (*B)(j+(k*nz_loc+n1)*nbmap,
  											  i+(nv*nz_loc+n1)*nimap) +=
  											 m_reynolds*coeffs[bmap[j]];
  									 }

  									 for(j = 0; j < nimap; ++j)
  									 {
  										 (*D)(j+(k*nz_loc+n1)*nimap,
  											  i+(nv*nz_loc+n1)*nimap) +=
  											 m_reynolds*coeffs[imap[j]];
  									 }
  								 }
  							 }
  						 }
  					 }
  				 }
   #if 0
  			 if(nz_loc  == 1)
  			 {
  				 //cout << "Ah" << endl;
  				 //cout << *Ah << endl;
  				 cout << "B" << endl;
  				 cout >> *B << endl;
  				 cout << "C" << endl;
  				 cout >> *C << endl;
  				 cout << "D" << endl;
  				 cout >> *D << endl;
  				 cout << *D << endl;
  				 //cout << "Dbnd" << endl;
  				 //cout >> *Dbnd << endl;
  				 //cout << "Dint" << endl;
  				 //cout >> *Dint << endl;
  				 exit(1);
  			 }
   #endif

  				 D->Invert();
  				 (*B) = (*B)*(*D);


  				 // perform (*Ah) = (*Ah) - (*B)*(*C) but since size of
  				 // Ah is larger than (*B)*(*C) easier to call blas
  				 // directly
  				 Blas::Dgemm('N','T', B->GetRows(), C->GetRows(),
  							 B->GetColumns(), -1.0, B->GetRawPtr(),
  							 B->GetRows(), C->GetRawPtr(),
  							 C->GetRows(), 1.0,
  							 Ah->GetRawPtr(), Ah->GetRows());
  			 }

  			 mat.m_BCinv->SetBlock(n,n,loc_mat = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,B));
  			 mat.m_Btilde->SetBlock(n,n,loc_mat = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,C));
  			 mat.m_Cinv->SetBlock(n,n,loc_mat = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,D));
  			 mat.m_D_bnd->SetBlock(n,n,loc_mat = MemoryManager<DNekScalMat>::AllocateSharedPtr(one, Dbnd));
  			 mat.m_D_int->SetBlock(n,n,loc_mat = MemoryManager<DNekScalMat>::AllocateSharedPtr(one, Dint));


  			 // Do matrix manipulations and get final set of block matries
  			 // reset boundary to put mean mode into boundary system.

  			 DNekMatSharedPtr Cinv,BCinv,Btilde;
  			 DNekMat  DintCinvDTint, BCinvDTint_m_DTbnd, DintCinvBTtilde_m_Dbnd;

  			 Cinv   = D;
  			 BCinv  = B;
  			 Btilde = C;


  			 DintCinvDTint = (*Dint)*(*Cinv)*Transpose(*Dint);

  			 BCinvDTint_m_DTbnd =  (*BCinv)*Transpose(*Dint) - Transpose(*Dbnd);

  			 // This could be transpose of BCinvDint in some cases
  			 DintCinvBTtilde_m_Dbnd =(*Dint)*(*Cinv)*Transpose(*Btilde) -  (*Dbnd);

  			 // Set up final set of matrices.
  			 DNekMatSharedPtr Bh = MemoryManager<DNekMat>::AllocateSharedPtr(nsize_bndry_p1[n],nsize_p_m1[n],zero);
  			 DNekMatSharedPtr Ch = MemoryManager<DNekMat>::AllocateSharedPtr(nsize_p_m1[n],nsize_bndry_p1[n],zero);
  			 DNekMatSharedPtr Dh = MemoryManager<DNekMat>::AllocateSharedPtr(nsize_p_m1[n], nsize_p_m1[n],zero);


  			 // Copy matrices into final structures.
  			 int n1,n2;
  			 for(n1 = 0; n1 < nz_loc; ++n1)
  			 {
  				 for(i = 0; i < psize-1; ++i)
  				 {
  					 for(n2 = 0; n2 < nz_loc; ++n2)
  					 {
  						 for(j = 0; j < psize-1; ++j)
  						 {
  							 (*Dh)(i+n1*(psize-1),j+n2*(psize-1)) =
  								 -DintCinvDTint(i+1+n1*psize,j+1+n2*psize);
  						 }
  					 }
  				 }
  			 }

  			 for(n1 = 0; n1 < nz_loc; ++n1)
  			 {
  				 for(i = 0; i < nsize_bndry_p1[n]-nz_loc; ++i)
  				 {
  					 (*Ah)(i,nsize_bndry_p1[n]-nz_loc+n1) = BCinvDTint_m_DTbnd(i,n1*psize);
  					 (*Ah)(nsize_bndry_p1[n]-nz_loc+n1,i) = DintCinvBTtilde_m_Dbnd(n1*psize,i);
  				 }
  			 }

  			 for(n1 = 0; n1 < nz_loc; ++n1)
  			 {
  				 for(n2 = 0; n2 < nz_loc; ++n2)
  				 {
  					 (*Ah)(nsize_bndry_p1[n]-nz_loc+n1,nsize_bndry_p1[n]-nz_loc+n2) =
  						 -DintCinvDTint(n1*psize,n2*psize);
  				 }
  			 }

  			 for(n1 = 0; n1 < nz_loc; ++n1)
  			 {
  				 for(j = 0; j < psize-1; ++j)
  				 {
  					 for(i = 0; i < nsize_bndry_p1[n]-nz_loc; ++i)
  					 {
  						 (*Bh)(i,j+n1*(psize-1)) = BCinvDTint_m_DTbnd(i,j+1+n1*psize);
  						 (*Ch)(j+n1*(psize-1),i) = DintCinvBTtilde_m_Dbnd(j+1+n1*psize,i);
  					 }

  					 (*Bh)(nsize_bndry_p1[n]-nz_loc+n1,j+n1*(psize-1)) = -DintCinvDTint(n1*psize,j+1+n1*psize);
  					 (*Ch)(j+n1*(psize-1),nsize_bndry_p1[n]-nz_loc+n1) = -DintCinvDTint(j+1+n1*psize,n1*psize);
  				 }
  			 }

   #if 0
  			 if(nz_loc  == 2)
  			 {
  				 //cout << "Cinv" << endl;
  				 //cout << *Cinv << endl;
  				 //cout << "Dbnd" << endl;
  				 //cout << *Dbnd << endl;
  				 //cout << "Dint" << endl;
  				 //cout << *Dint << endl;
  				 cout << "Ah" << endl;
  				 cout << *Ah << endl;
  				 cout << "Bh" << endl;
  				 cout << *Bh << endl;
  				 cout << "Ch" << endl;
  				 cout << *Ch << endl;
  				 cout << "Dh" << endl;
  				 cout << *Dh << endl;
  				 exit(1);
  			 }
   #endif
  			 // Do static condensation
  			 Dh->Invert();
  			 (*Bh) = (*Bh)*(*Dh);
  			 (*Ah) = (*Ah) - (*Bh)*(*Ch);

  			 // Set matrices for later inversion. Probably do not need to be
  			 // attached to class
  			 pAh->SetBlock(n,n,loc_mat = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,Ah));
  			 pBh->SetBlock(n,n,loc_mat = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,Bh));
  			 pCh->SetBlock(n,n,loc_mat = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,Ch));
  			 pDh->SetBlock(n,n,loc_mat = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,Dh));
  		 }

  		 // Set up global coupled boundary solver.
  		 // This is a key to define the solution matrix type
  		 // currently we are giving it a argument of eLInearAdvectionReaction
  		 // since this then makes the matrix storage of type eFull
  		 MultiRegions::GlobalLinSysKey key(StdRegions::eLinearAdvectionReaction,locToGloMap);
  		 mat.m_CoupledBndSys.reset();
  		 mat.m_CoupledBndSys = MemoryManager<MultiRegions::GlobalLinSysDirectStaticCond>::AllocateSharedPtr(key,m_fields[0],pAh,pBh,pCh,pDh,locToGloMap);
  	 }
    
    void CoupledLinearNS::v_PrintSummary(std::ostream &out)
    {
        cout <<  "\tSolver Type     : Coupled Linearised NS" <<endl;
    }

    void CoupledLinearNS::v_DoInitialise(void)
    {
        switch(m_equationType)
        {
        case eUnsteadyStokes:
        case eUnsteadyNavierStokes:
            {

                LibUtilities::TimeIntegrationMethod intMethod;
                std::string TimeIntStr = m_session->GetSolverInfo("TIMEINTEGRATIONMETHOD");
                int i;
                for(i = 0; i < (int) LibUtilities::SIZE_TimeIntegrationMethod; ++i)
                {
                    if(NoCaseStringCompare(LibUtilities::TimeIntegrationMethodMap[i],TimeIntStr) == 0 )
                    {
                        intMethod = (LibUtilities::TimeIntegrationMethod)i; 
                        break;
                    }
                }
                
                ASSERTL0(i != (int) LibUtilities::SIZE_TimeIntegrationMethod, "Invalid time integration type.");
                
                switch(intMethod)
                {
                case LibUtilities::eIMEXOrder1: 
                    {
                        m_intSteps = 1;
                        m_integrationScheme = Array<OneD, LibUtilities::TimeIntegrationSchemeSharedPtr> (m_intSteps);
                        LibUtilities::TimeIntegrationSchemeKey       IntKey0(intMethod);
                        m_integrationScheme[0] = LibUtilities::TimeIntegrationSchemeManager()[IntKey0];
                    }
                    break;
                case LibUtilities::eIMEXOrder2: 
                    {
                        m_intSteps = 2;
                        m_integrationScheme = Array<OneD, LibUtilities::TimeIntegrationSchemeSharedPtr> (m_intSteps);
                        LibUtilities::TimeIntegrationSchemeKey       IntKey0(LibUtilities::eIMEXOrder1);
                        m_integrationScheme[0] = LibUtilities::TimeIntegrationSchemeManager()[IntKey0];
                        LibUtilities::TimeIntegrationSchemeKey       IntKey1(intMethod);
                        m_integrationScheme[1] = LibUtilities::TimeIntegrationSchemeManager()[IntKey1];
                    }
                    break;
                default:
                    ASSERTL0(0,"Integration method not setup: Options include ImexOrder1, ImexOrder2");
                    break;
                }
                
                // Could defind this from IncNavierStokes class? 
                m_integrationOps.DefineOdeRhs(&CoupledLinearNS::EvaluateAdvection, this);

                m_integrationOps.DefineImplicitSolve(&CoupledLinearNS::SolveUnsteadyStokesSystem,this);

                // Set initial condition using time t=0
                SetInitialConditions(0.0);
            }
        case eSteadyStokes:
            SetUpCoupledMatrix(0.0);
            break;
        case eUnsteadyNavierStokesALE:
        {
        	SetInitialConditions(0.0);
            ImportFldPressure();
            int phystot=m_fields[0]->GetTotPoints();
			cout << endl;
			for(int i = 0; i < m_fields.num_elements(); ++i)
			{
			cout << GetVariable(i) << ": Max :" << Vmath::Vmax(phystot,m_fields[i]->GetPhys(),1);
			cout << " Min: " << Vmath::Vmin(phystot,m_fields[i]->GetPhys(),1)<< "\t";
			}
			cout << endl;

			cout << "p" << ": Max :" << Vmath::Vmax(phystot,m_pressure->GetPhys(),1);
			cout << " Min: " << Vmath::Vmin(phystot,m_pressure->GetPhys(),1)<< endl;
			SetUpCoupledMatrix(0.0);
        }
		break;
        case eSteadyOseen:
            {
                Array<OneD, Array<OneD, NekDouble> > AdvField(m_velocity.num_elements());
                for(int i = 0; i < m_velocity.num_elements(); ++i)
                {
                    AdvField[i] = Array<OneD, NekDouble> (m_fields[m_velocity[i]]->GetTotPoints(),0.0);
                }

                ASSERTL0(m_session->DefinesFunction("AdvectionVelocity"),
                            "Advection Velocity section must be defined in "
                            "session file.");

                std::vector<std::string> fieldStr;
                for(int i = 0; i < m_velocity.num_elements(); ++i)
                {
                    fieldStr.push_back(m_boundaryConditions->GetVariable(m_velocity[i]));
                }
                EvaluateFunction(fieldStr,AdvField,"AdvectionVelocity");
                
                SetUpCoupledMatrix(0.0,AdvField,false);
            }
            break;
        case eSteadyLinearisedNS:
            {                
                
                Array<OneD, Array<OneD, NekDouble> > AdvField(m_velocity.num_elements());
                for(int i = 0; i < m_velocity.num_elements(); ++i)
                {
                    AdvField[i] = Array<OneD, NekDouble> (m_fields[m_velocity[i]]->GetTotPoints(),0.0);
                }
                
                ASSERTL0(m_session->DefinesFunction("AdvectionVelocity"),
                            "Advection Velocity section must be defined in "
                            "session file.");

                std::vector<std::string> fieldStr;
                for(int i = 0; i < m_velocity.num_elements(); ++i)
                {
                    fieldStr.push_back(m_boundaryConditions->GetVariable(m_velocity[i]));
                }
                EvaluateFunction(fieldStr,AdvField,"AdvectionVelocity");

                SetUpCoupledMatrix(m_lambda,AdvField,true);
            }
            break;
            case eViscoelasticDEVSSG:
        case eUnsteadyViscoelastic:
        case eViscoelasticLogc:
        {
        	 LibUtilities::TimeIntegrationMethod intMethod;
        	                std::string TimeIntStr = m_session->GetSolverInfo("TIMEINTEGRATIONMETHOD");
        	                int i;
        	                for(i = 0; i < (int) LibUtilities::SIZE_TimeIntegrationMethod; ++i)
        	                {
        	                    if(NoCaseStringCompare(LibUtilities::TimeIntegrationMethodMap[i],TimeIntStr) == 0 )
        	                    {
        	                        intMethod = (LibUtilities::TimeIntegrationMethod)i;
        	                        break;
        	                    }
        	                }

        	                ASSERTL0(i != (int) LibUtilities::SIZE_TimeIntegrationMethod, "Invalid time integration type.");

        	                switch(intMethod)
        	                {
        	                case LibUtilities::eIMEXOrder1:
        	                    {
        	                        m_intSteps = 1;
        	                        m_integrationScheme = Array<OneD, LibUtilities::TimeIntegrationSchemeSharedPtr> (m_intSteps);
        	                        LibUtilities::TimeIntegrationSchemeKey       IntKey0(intMethod);
        	                        m_integrationScheme[0] = LibUtilities::TimeIntegrationSchemeManager()[IntKey0];
        	                    }
        	                    break;
        	                case LibUtilities::eIMEXOrder2:
        	                    {
        	                        m_intSteps = 2;
        	                        m_integrationScheme = Array<OneD, LibUtilities::TimeIntegrationSchemeSharedPtr> (m_intSteps);
        	                        LibUtilities::TimeIntegrationSchemeKey       IntKey0(LibUtilities::eIMEXOrder1);
        	                        m_integrationScheme[0] = LibUtilities::TimeIntegrationSchemeManager()[IntKey0];
        	                        LibUtilities::TimeIntegrationSchemeKey       IntKey1(intMethod);
        	                        m_integrationScheme[1] = LibUtilities::TimeIntegrationSchemeManager()[IntKey1];
        	                    }
        	                    break;
        	                default:
        	                    ASSERTL0(0,"Integration method not setup: Options include ImexOrder1, ImexOrder2");
        	                    break;
        	                }

        	                // Could defind this from IncNavierStokes class?
        	                m_integrationOps.DefineOdeRhs(&CoupledLinearNS::EvaluateAdvection, this);

        	                m_integrationOps.DefineImplicitSolve(&CoupledLinearNS::SolveUnsteadyStokesSystem,this);

        	                // Set initial condition using time t=0
        	                //SetInitialConditions(0.0);

            // Storage array for weak stress divergence of current and past time levels
               m_stressdivergenceweaku = Array<OneD, Array<OneD, NekDouble> > (m_intSteps);
               m_stressdivergenceweakv = Array<OneD, Array<OneD, NekDouble> > (m_intSteps);

               int npoints = m_fields[0]->GetNcoeffs();
              // int npoints = m_fields[0]->GetTotPoints();

               for(int n = 0; n < m_intSteps; ++n)
               {
               	m_stressdivergenceweaku[n] = Array<OneD, NekDouble>(npoints);
               	m_stressdivergenceweakv[n] = Array<OneD, NekDouble>(npoints);
               }

       	 // Set initial condition using time t=0
			SetInitialConditions(0.0);
		//	SetFixedCylinderVelocityField();
			SetInitialConditionsStress(0.0);
			SetUpCoupledMatrix(0.0);
		//	ImportFldPressure();

			int phystot=m_fields[0]->GetTotPoints();

             cout << endl;
             for(i = 0; i < m_fields.num_elements(); ++i)
             {
            cout << GetVariable(i) << ": Max :" << Vmath::Vmax(phystot,m_fields[i]->GetPhys(),1);
            cout << " Min: " << Vmath::Vmin(phystot,m_fields[i]->GetPhys(),1)<< "\t";
             }
             cout << endl;

             cout << "p" << ": Max :" << Vmath::Vmax(phystot,m_pressure->GetPhys(),1);
             cout << " Min: " << Vmath::Vmin(phystot,m_pressure->GetPhys(),1)<< "\t";

        }
        break;
        case eViscoelasticDEVSSGALE:
        {
       	    LibUtilities::TimeIntegrationMethod intMethod;
			std::string TimeIntStr = m_session->GetSolverInfo("TIMEINTEGRATIONMETHOD");
			int i;
			for(i = 0; i < (int) LibUtilities::SIZE_TimeIntegrationMethod; ++i)
			{
				if(NoCaseStringCompare(LibUtilities::TimeIntegrationMethodMap[i],TimeIntStr) == 0 )
				{
					intMethod = (LibUtilities::TimeIntegrationMethod)i;
					break;
				}
			}

            ASSERTL0(i != (int) LibUtilities::SIZE_TimeIntegrationMethod, "Invalid time integration type.");

            switch(intMethod)
            {
				case LibUtilities::eIMEXOrder1:
					{
						m_intSteps = 1;
					}
					break;
				case LibUtilities::eIMEXOrder2:
					{
						m_intSteps = 2;
					}
					break;
				default:
					ASSERTL0(0,"Integration method not setup: Options include ImexOrder1, ImexOrder2");
					break;
            }

            m_stressdivergenceweaku = Array<OneD, Array<OneD, NekDouble> > (m_intSteps);
            m_stressdivergenceweakv = Array<OneD, Array<OneD, NekDouble> > (m_intSteps);

            int npoints = m_fields[0]->GetNcoeffs();
           // int npoints = m_fields[0]->GetTotPoints();

            for(int n = 0; n < m_intSteps; ++n)
            {
            	m_stressdivergenceweaku[n] = Array<OneD, NekDouble>(npoints);
            	m_stressdivergenceweakv[n] = Array<OneD, NekDouble>(npoints);
            }

        	SetInitialConditions(0.0);
            ImportFldPressure();
            SetInitialConditionsStress(0.0);

            int phystot=m_fields[0]->GetTotPoints();
			cout << endl;
			for(int i = 0; i < m_fields.num_elements(); ++i)
			{
			cout << GetVariable(i) << ": Max :" << Vmath::Vmax(phystot,m_fields[i]->GetPhys(),1);
			cout << " Min: " << Vmath::Vmin(phystot,m_fields[i]->GetPhys(),1)<< "\t";
			}
			cout << endl;

			cout << "p" << ": Max :" << Vmath::Vmax(phystot,m_pressure->GetPhys(),1);
			cout << " Min: " << Vmath::Vmin(phystot,m_pressure->GetPhys(),1)<< endl;
			SetUpCoupledMatrix(0.0);
        }
		break;
        case eNoEquationType:
        default:
            ASSERTL0(false,"Unknown or undefined equation type for CoupledLinearNS");
        }
    }

    void CoupledLinearNS::EvaluateAdvection(const Array<OneD, const Array<OneD, NekDouble> > &inarray, 
                                            Array<OneD, Array<OneD, NekDouble> > &outarray,
                                            const NekDouble time)
    {
        // evaluate convection terms
        if(m_reynolds>0)
    	{
    		EvaluateAdvectionTerms(inarray,outarray);
    	}
    	else
    	{
    		int physTot = m_fields[0]->GetTotPoints();
    		Vmath::Zero(physTot,outarray[m_velocity[0]],1);
    		Vmath::Zero(physTot,outarray[m_velocity[1]],1);
    	}
        int nqtot  =m_fields[0]->GetTotPoints(); 
        //add the force
        if(m_session->DefinesFunction("BodyForce"))
        {
            for(int i = 0; i < m_nConvectiveFields; ++i)
            {
                Vmath::Vadd(nqtot,outarray[i],1,(m_forces[i]->GetPhys()),1,outarray[i],1);
            }        
        }
        
        if(m_equationType==eUnsteadyViscoelastic || m_equationType==eViscoelasticDEVSSG)
        {
        	EvaluateRhsViscoelasticEquation(inarray,outarray);
        }
        else if(m_equationType==eViscoelasticLogc)
		{
        	 if(m_reynolds>0)
        	 {
        		 AddStressDivergenceLogc(inarray,outarray);
        	 }

			 EvaluateRhsViscoelasticEquationLogc(inarray,outarray);
		}
    }

    void CoupledLinearNS::SolveUnsteadyStokesSystem(const Array<OneD, const Array<OneD, NekDouble> > &inarray, 
                                                    Array<OneD, Array<OneD, NekDouble> > &outarray, 
                                                    const NekDouble time, 
                                                    const NekDouble aii_Dt)
    {
        int i,n;
        Array<OneD, Array< OneD, NekDouble> > F(m_nConvectiveFields);
        NekDouble  lambda = 1.0/aii_Dt; 
        static NekDouble lambda_store; 
        Array <OneD, Array<OneD, NekDouble> > forcing(m_velocity.num_elements());
        static int count=0;
        static int cntoutput = 1;
        // Matrix solution

     //  if(count==0)
     //   {
        //TODO: undo this
        if(m_reynolds>0)
        {
        // Matrix solution
        if(fabs(lambda_store - lambda) > 1e-10)
        {
            cout << "Setting up Stokes matrix problem [.";
            fflush(stdout);
            SetUpCoupledMatrix(lambda);
            cout << "]" << endl;
            lambda_store = lambda;
        }
        
        SetBoundaryConditions(time);

        // Forcing for advection solve 
        for(int i = 0; i < m_velocity.num_elements(); ++i)
        {
            m_fields[m_velocity[i]]->IProductWRTBase(inarray[i],m_fields[m_velocity[i]]->UpdateCoeffs());
            Vmath::Smul(m_fields[m_velocity[i]]->GetNcoeffs(),lambda*m_reynolds,m_fields[m_velocity[i]]->GetCoeffs(), 1,m_fields[m_velocity[i]]->UpdateCoeffs(),1);
            forcing[i] = m_fields[m_velocity[i]]->GetCoeffs();
        }
       }
       else
       {
     	    for(int i = 0; i < m_velocity.num_elements(); ++i)
            {
                forcing[i] = Array<OneD, NekDouble> (m_fields[m_velocity[0]]->GetNcoeffs(),0.0);
            }
       }
        if(m_equationType==eUnsteadyViscoelastic || m_equationType==eViscoelasticDEVSSG)
        {
        	AddStressDivergenceWeakFormExplicit(inarray,forcing);
        }
        if(m_equationType==eViscoelasticLogc && m_reynolds==0.0)
        {
    	   AddStressDivergenceWeakFormLogc(inarray,forcing);
        }
        
        SolveLinearNS(forcing);
        
        for(int i = 0; i < m_velocity.num_elements(); ++i)
        {
            m_fields[m_velocity[i]]->BwdTrans(m_fields[m_velocity[i]]->GetCoeffs(),outarray[i]);
            m_fields[m_velocity[i]]->BwdTrans(m_fields[m_velocity[i]]->GetCoeffs(),m_fields[m_velocity[i]]->UpdatePhys());
        }

        //static int cntoutput = 1;
        //Checkpoint_TecplotOutput(cntoutput++, true);
        
        if(m_equationType==eViscoelasticDEVSSG)
        {
        	CalcGL2Projection(outarray);
        }
       // }

        // TODO: remove this part again just a test case !!!
     /*   for(int i = 0; i < m_velocity.num_elements(); ++i)
		{
			Vmath::Vcopy(m_fields[0]->GetTotPoints(),m_fields[m_velocity[i]]->GetPhys(),1,outarray[i],1);
		} */

        if(m_equationType==eUnsteadyViscoelastic || m_equationType==eViscoelasticDEVSSG)
        {
			if(m_alpha)
			{
				CalcSourceTermImplicitGiesekus(inarray,outarray,aii_Dt);
			}
			else
			{
				CalcSourceTermImplicitOldB(inarray,outarray,aii_Dt);
			}

			for(int i = 0; i < m_stress.num_elements(); ++i)
			{
			   if(m_projectionType == MultiRegions::eGalerkin)
			   {
				   int VelDim     = m_velocity.num_elements();
				   const Array<OneD,  Array<OneD,const NekDouble> > velocity(VelDim);
				   Array<OneD,  NekDouble> stressin;
				   for(int j = 0; j < VelDim; ++j)
					{
						//velocity[j] = outarray[m_velocity[j]];
					}

				  // m_fields[m_stress[i]]->FwdTransSUPG(velocity,outarray[m_stress[i]],m_fields[m_stress[i]]->UpdateCoeffs(),false);
					// m_fields[m_stress[i]]->IProductWRTSUPG(velocity,outarray[m_stress[i]],m_fields[m_stress[i]]->UpdateCoeffs());
					// m_fields[m_stress[i]]->MultiplyByInvMassMatrix(m_fields[m_stress[i]]->GetCoeffs(), m_fields[m_stress[i]]->UpdateCoeffs());
				   m_fields[m_stress[i]]->FwdTrans(outarray[m_stress[i]],m_fields[m_stress[i]]->UpdateCoeffs());
			   }
			   else
			   {
				   m_fields[m_stress[i]]->FwdTrans(outarray[m_stress[i]],m_fields[m_stress[i]]->UpdateCoeffs());
			   }

				 m_fields[m_stress[i]]->BwdTrans(m_fields[m_stress[i]]->GetCoeffs(),outarray[m_stress[i]]);

				 m_stressold[i] = outarray[m_stress[i]];
				 // Update physical space
				 m_fields[m_stress[i]]->BwdTrans(m_fields[m_stress[i]]->GetCoeffs(),m_fields[m_stress[i]]->UpdatePhys());
				 //m_fields[m_stress[i]]->SetPhys(outarray[m_stress[i]]);
				 //m_fields[m_stress[i]]->SetPhysState(true);
			}
		}
        else if(m_equationType==eViscoelasticLogc)
        {
        	for(i = 0; i < m_logc.num_elements(); i++)
			{
                 if(m_projectionType == MultiRegions::eGalerkin)
                 {
				 /* For calculations from the steady state solution this transformation introduces an error */
				 m_fields[m_logc[i]]->FwdTrans(inarray[m_logc[i]],m_fields[m_logc[i]]->UpdateCoeffs());

				 // This causes an error at times for channel flow of 10^{-7}
				 m_fields[m_logc[i]]->BwdTrans(m_fields[m_logc[i]]->GetCoeffs(),outarray[m_logc[i]]);
				// m_fields[m_logc[i]]->BwdTrans(m_fields[m_logc[i]]->GetCoeffs(),m_fields[m_logc[i]]->UpdatePhys());
                 }
                else
                 {
                	 outarray[m_logc[i]] = inarray[m_logc[i]];
                 }
			}

        	/* Compute the elastic stresses at the new time level */
			  LogcToElasticL2Projection(inarray, outarray);
        	//TODO: dont know if this is necessary
        /*	for(i = 0; i < m_stress.num_elements(); ++i)
			{//TODO: activate this again
        		m_fields[m_stress[i]]->FwdTrans(m_fields[m_stress[i]]->GetPhys(),m_fields[m_stress[i]]->UpdateCoeffs());
        		 m_fields[m_stress[i]]->BwdTrans(m_fields[m_stress[i]]->GetCoeffs(),m_fields[m_stress[i]]->UpdatePhys());

        		//m_fields[m_stress[i]]->FwdTrans(outarray[m_stress[i]],m_fields[m_stress[i]]->UpdateCoeffs());
				// m_fields[m_stress[i]]->BwdTrans(m_fields[m_stress[i]]->GetCoeffs(),outarray[m_stress[i]]);
			} */

        }

        // project pressure field to velocity space
		//m_pressure->BwdTrans(m_pressure->GetCoeffs(), m_pressure->UpdatePhys());


		 if(!((count+1)%m_checksteps))
		 {
			  Checkpoint_TecplotOutput(cntoutput, true);
			  CheckPointOutput(cntoutput);
			  cntoutput++;
		 }

    	if(m_equationType==eUnsteadyViscoelastic|| m_equationType==eViscoelasticDEVSSG)
    	{


    		//cout << "Max" << m_fields[m_]
            if(count==0)
            {
            //Checkpoint_MatlabOutput(count, false);
            }

    		if(!((count+1)%m_infosteps))
			{
    	      /*  for(int i = 0; i < m_velocity.num_elements(); ++i)
    	        {
    	            m_fields[m_velocity[i]]->BwdTrans(m_fields[m_velocity[i]]->GetCoeffs(),m_fields[m_velocity[i]]->UpdatePhys());
    	        } */
    		            Viscoelastic::WriteDragToFile(count,outarray);

    		//Checkpoint_TecplotOutput(count, false);
    		          //  Checkpoint_MatlabOutput(count, false);
			}

    	}
    	count++;
    }
    
    void CoupledLinearNS::SolveUnsteadyStokesSystemALE(const Array<OneD, const Array<OneD, NekDouble> > &inarray,
    												   Array<OneD, Array<OneD, NekDouble> > &outarray,
    												   const NekDouble time,
    												   const NekDouble aii_Dt)
	{
	   int i,n;
	   int phystot = m_fields[0]->GetTotPoints();
	   int ncoeffs= m_fields[m_velocity[0]]->GetNcoeffs();
	   int VelDim = m_velocity.num_elements();
	   Array<OneD, Array< OneD, NekDouble> > F(m_nConvectiveFields);
	   NekDouble  lambda = 1.0/aii_Dt;
	   static NekDouble lambda_store;
	   NekDouble gstar = 9.81*m_reynolds;
	   Array<OneD, NekDouble> gx(phystot,gstar);
	   Array<OneD, NekDouble> gxweak(ncoeffs,0.0);
	   static int cntoutput = 0;
	   static int nsteps = 0;
	   NekDouble maximum=1.0;
	   int it=1;
	   Array<OneD, NekDouble> max(m_velocity.num_elements());
	   for(i = 0; i<m_velocity.num_elements(); i++)
	   {
		   max[i]=1.0;
	   }

	   Array <OneD, Array<OneD, NekDouble> > forcing(m_velocity.num_elements());
	   for(i = 0; i < m_velocity.num_elements(); ++i)
	   {
		   forcing[i] = Array<OneD, NekDouble> (ncoeffs,0.0);
	   }
	   Array <OneD, Array<OneD, NekDouble> > fieldsold(m_nConvectiveFields);
	   Array <OneD, Array<OneD, NekDouble> > fieldsn(m_nConvectiveFields);
	   Array <OneD, Array<OneD, NekDouble> > rhs(m_nConvectiveFields);

	   for(i = 0; i < m_nConvectiveFields; ++i)
	   {
		   fieldsold[i] = Array<OneD, NekDouble> (phystot,0.0);
		   fieldsn[i] = Array<OneD, NekDouble> (phystot,0.0);
		   rhs[i] = Array<OneD, NekDouble> (phystot,0.0);
	   }

	   for(i = 0; i < m_nConvectiveFields; ++i)
	   {
		 fieldsn[i] = Array<OneD, NekDouble> (phystot,0.0);
		 Vmath::Vcopy(phystot,inarray[i],1,fieldsn[i],1);
	   }

	   Array<OneD, Array<OneD, NekDouble> > Advfield(VelDim);
	   Array<OneD, NekDouble > Deriv;

	   while(maximum>1e-5)
	   {

		for(i = 0; i < m_nConvectiveFields; ++i)
		{
			Vmath::Vcopy(phystot,m_fields[i]->GetPhys(),1,fieldsold[i],1);
		}

	   for(i = 0; i < VelDim; ++i)
	   {
		   Advfield[i] = Array<OneD, NekDouble>(phystot);
		   Vmath::Vcopy(phystot,m_fields[m_velocity[i]]->GetPhys(),1,Advfield[i],1);
	   }

	   if(m_ALE)
	   {
//		  ComputeVelocityMinusMeshVelocity(Advfield);
	   }

	   // Matrix solution
	   SetUpCoupledMatrix(lambda, Advfield, false, true);

	   SetBoundaryConditions(time);

	   // Forcing for advection solve
	   for(i = 0; i < m_velocity.num_elements(); ++i)
	   {
		   m_fields[m_velocity[i]]->IProductWRTBase(inarray[i],forcing[i]);
		   Vmath::Smul(ncoeffs,lambda*m_reynolds,forcing[i], 1,forcing[i],1);
	   }

       if(m_equationType==eViscoelasticDEVSSGALE)
       {
       	  AddStressDivergenceWeakFormExplicit(inarray,forcing);
       }

	   SubGraduTRHS(forcing);

	   SolveLinearNS(forcing);

	   for(i = 0; i < m_velocity.num_elements(); ++i)
	   {
		   m_fields[m_velocity[i]]->BwdTrans(m_fields[m_velocity[i]]->GetCoeffs(),m_fields[m_velocity[i]]->UpdatePhys());
		   m_fields[m_velocity[i]]->BwdTrans(m_fields[m_velocity[i]]->GetCoeffs(),outarray[i]);
	   }

	   m_pressure->BwdTrans(m_pressure->GetCoeffs(),m_pressure->UpdatePhys());

       if(m_equationType==eViscoelasticDEVSSGALE)
       {
       	  CalcGL2Projection(outarray);
       	  SolveOlroydBPseudoSpectral(outarray,fieldsn,aii_Dt);
       }

	   for(i = 0; i < m_velocity.num_elements(); ++i)
		{
		 Vmath::Vsub(phystot,m_fields[m_velocity[i]]->GetPhys(),1,fieldsold[i],1,fieldsold[i],1);
		 Vmath::Vabs(phystot,fieldsold[i],1,fieldsold[i],1);
		 max[i] = Vmath::Vmax(phystot,fieldsold[i],1);
		}
	   maximum = (max[0]<max[1]) ? max[1] : max[0];
	   if(!(((nsteps+1)%m_infosteps))||!(it%100))
		{
			cout << "it:" << it << " max:" << maximum << endl;
		}

	   it++;
	   }

	   for(int i = 0; i < m_velocity.num_elements(); ++i)
	   {
		   m_fields[m_velocity[i]]->BwdTrans(m_fields[m_velocity[i]]->GetCoeffs(),outarray[i]);
	   }

	   nsteps++;
	}

	
	void CoupledLinearNS::v_TransCoeffToPhys(void)
    {
		int nfields = m_fields.num_elements();
		for (int k=0 ; k < nfields; ++k)
		{
			//Backward Transformation in physical space for time evolution
			m_forces[k]->BwdTrans_IterPerExp(m_forces[k]->GetCoeffs(),
										     m_forces[k]->UpdatePhys());			
		}
		
	}
	
	void CoupledLinearNS::v_TransPhysToCoeff(void)
    {
		int nfields = m_fields.num_elements();
		for (int k=0 ; k < nfields; ++k)
		{
			//Forward Transformation in physical space for time evolution
			m_forces[k]->FwdTrans_IterPerExp(m_forces[k]->GetPhys(),
										     m_forces[k]->UpdateCoeffs());
			
		}
	}
	
    void CoupledLinearNS::v_DoSolve(void)
    {
        switch(m_equationType)
        {
        case eUnsteadyViscoelastic:
        case eViscoelasticLogc:
        case eViscoelasticDEVSSG:
			AdvanceInTime(m_steps);
			break;
        case eUnsteadyStokes:
        case eUnsteadyNavierStokes:
            AdvanceInTime(m_steps);
            break;
        case eUnsteadyNavierStokesALE:
        case eViscoelasticDEVSSGALE:
            AdvanceInTimeALE(m_steps);
            break;
        case eSteadyStokes:
        case eSteadyOseen:
        case eSteadyLinearisedNS:
            Solve();
            break;
        case eNoEquationType:
        default:
            ASSERTL0(false,"Unknown or undefined equation type for CoupledLinearNS");
        }
    }

    void CoupledLinearNS::AdvanceInTimeALE(int nsteps)
	{
	   int i,n;
	   int it;
	   int phystot = m_fields[0]->GetTotPoints();
	   static int nchk = 0;

	   int n_fields = m_fields.num_elements();
	   Array<OneD, NekDouble> max(m_nConvectiveFields);
	   NekDouble maximum=1.0;

	   // Set up wrapper to fields data storage.
	   Array<OneD, Array<OneD, NekDouble> >   fields(m_nConvectiveFields);
	   for(i = 0; i < m_nConvectiveFields; ++i)
	   {
		   fields[i]  = m_fields[i]->UpdatePhys();
	   }

	   Array<OneD, Array<OneD, NekDouble> >   fieldsold(m_nConvectiveFields);
	   for(i = 0; i < m_nConvectiveFields; ++i)
	   {
		   fieldsold[i]  = Array<OneD, NekDouble> (m_fields[0]->GetTotPoints(),0.0);
		   Vmath::Vcopy(phystot,fields[i],1,fieldsold[i],1);
	   }

	   //Time advance
	   for(n = 0; n < nsteps; ++n)
	   {
		   // Advance velocity fields
		   if(m_ALE)
		   {
		/*	if(n){
				DoMovingMesh(m_timestep, m_time);
				UpdateFields();
				UpdatePressure();
				if(m_equationType == eViscoelasticDEVSSGALE)
				{
					UpdateG();
				}
			}
			else
			{
				DoMovingMeshInitialise(m_timestep);
				UpdateFields();
				UpdatePressure();
				if(m_equationType == eViscoelasticDEVSSGALE)
				{
					UpdateG();
				}
			} */
		   }

			if(n==0)
			{
				maximum =1;
				it =1;
				while(maximum>1e-5)
				{
					SolveUnsteadyStokesSystemALE(fieldsold,fields, m_time, m_timestep);

					for(i = 0; i < m_nConvectiveFields; ++i)
					{
						Vmath::Vsub(m_fields[0]->GetTotPoints(),fields[i],1,fieldsold[i],1,fieldsold[i],1);
						max[i] = Vmath::Vmax(m_fields[0]->GetTotPoints(),fieldsold[i],1);
					}

					maximum = (max[0] > max[1] ? max[0] : max[1]);
					//if(!((it+1)%m_infosteps))
					if(m_infosteps && !((it+1)%m_infosteps))
					{
						for(i = 0; i < m_nConvectiveFields; ++i)
						{
							cout << "outer it: "<< it << ": Max( " << m_session->GetVariable(i) << ")= " << max[i] << endl;
						}
					}

					for(i = 0; i < m_nConvectiveFields; ++i)
					{
						Vmath::Vcopy(phystot,fields[i],1,fieldsold[i],1);
					}
					it++;
				}
			}

			SolveUnsteadyStokesSystemALE(fieldsold,fields, m_time, m_timestep);
			for(i = 0; i < m_nConvectiveFields; ++i)
			{
				Vmath::Vcopy(phystot,fields[i],1,fieldsold[i],1);
			}

			// Write out current time step
			if(m_infosteps && !((n+1)%m_infosteps))
			{
				cout << "Step: " << n+1 << "  Time: " << m_time << endl;
				for(i = 0; i < m_fields.num_elements(); ++i)
				{
				   cout << GetVariable(i) << ": Max :" << Vmath::Vmax(phystot,m_fields[i]->GetPhys(),1);
				   cout << " Min: " << Vmath::Vmin(phystot,m_fields[i]->GetPhys(),1) << "\t";
				   if(i==1) cout << endl;
				   if(isnan(Vmath::Vmax(phystot,m_fields[i]->GetPhys(),1)))
				   {
					  ASSERTL0(false,"Number too big stop computing");
				   }
				}
				cout << endl;
				cout << "p: Max :" << Vmath::Vmax(phystot,m_pressure->GetPhys(),1);
				cout << " Min: " << Vmath::Vmin(phystot,m_pressure->GetPhys(),1)<< "\t";
				cout << endl;
			}

			// dump data in m_fields->m_coeffs to file.
			if(m_checksteps && n&&(!((n+1)%m_checksteps)))
			{
				for(i = 0; i < m_nConvectiveFields; ++i)
				{
					m_fields[i]->SetPhys(fields[i]);
					m_fields[i]->SetPhysState(true);
				}
				Checkpoint_Output(nchk);
				Checkpoint_TecplotOutput(nchk, true);
		//		WriteXmlFile(nchk);
				nchk++;
			}

			m_time += m_timestep;

	   }

	//Array<OneD, NekDouble>   test(256,0.0);
	   for(i = 0; i < m_nConvectiveFields; ++i)
	   {
		   m_fields[i]->SetPhys(fields[i]);
		   m_fields[i]->SetPhysState(true);
	   }
	}

    void CoupledLinearNS::Solve(void)
    {
        Array <OneD, Array<OneD, NekDouble> > forcing(m_velocity.num_elements());
        
        if(m_session->DefinesFunction("BodyForce"))
        {
            for(int i = 0; i < m_velocity.num_elements(); ++i)
            {
                m_forces[i]->IProductWRTBase(m_forces[i]->GetPhys(),
                                             m_forces[i]->UpdateCoeffs());
                if(m_HomogeneousType == eHomogeneous1D)
                {
                    if(!m_singleMode) // Do not want to Fourier transoform in case of single mode stability analysis. 
                    {
                        m_forces[i]->HomogeneousFwdTrans(m_forces[i]->GetCoeffs(),m_forces[i]->UpdateCoeffs());
                    }
                }
                forcing[i] = m_forces[i]->GetCoeffs();
            }
        }
        else
        {
            // Should put in read forcing in here. 
            for(int i = 0; i < m_velocity.num_elements(); ++i)
            {
                forcing[i] = Array<OneD, NekDouble> (m_fields[m_velocity[0]]->GetNcoeffs(),0.0);
            }
        }
        
        SolveLinearNS(forcing);

    }

    const SpatialDomains::ExpansionMap &CoupledLinearNS::GenPressureExp(const SpatialDomains::ExpansionMap &VelExp)
    {
        int i;
        SpatialDomains::ExpansionMapShPtr returnval;

        returnval = MemoryManager<SpatialDomains::ExpansionMap>::AllocateSharedPtr();
        
        SpatialDomains::ExpansionMap::const_iterator  expMapIter;
        int nummodes;

        for (expMapIter = VelExp.begin(); expMapIter != VelExp.end(); ++expMapIter)
        {
            LibUtilities::BasisKeyVector BasisVec;

            for(i = 0; i <  expMapIter->second->m_basisKeyVector.size(); ++i)
            {
                LibUtilities::BasisKey B = expMapIter->second->m_basisKeyVector[i];
                nummodes = B.GetNumModes();
                ASSERTL0(nummodes > 3,"Velocity polynomial space not sufficiently high (>= 4)");
                // Should probably set to be an orthogonal basis. 
                LibUtilities::BasisKey newB(B.GetBasisType(),nummodes-2,B.GetPointsKey());
                BasisVec.push_back(newB);
            }

            // Put new expansion into list. 
            SpatialDomains::ExpansionShPtr expansionElementShPtr =
                MemoryManager<SpatialDomains::Expansion>::AllocateSharedPtr(expMapIter->second->m_geomShPtr, BasisVec);
            (*returnval)[expMapIter->first] = expansionElementShPtr;
        }
        
        // Save expansion into graph. 
        m_graph->SetExpansions("p",returnval);
            
        return *returnval;
    }
    
    /** 
     *  @param forcing A list of forcing functions for each velocity
     *  component
     *   
     *  The routine involves two levels of static
     *  condensations. Initially we require a statically condensed
     *  forcing function which requires the following manipulation
     *
     * \f[ {F\_bnd} = {\bf f}_{bnd} -m\_B \,m\_Cinv\, {\bf f}_{int},
     * \hspace{1cm} F\_p = m\_D\_{int}\, m\_Cinv\, {\bf f}_{int} \f]
     * 
     *  Where \f${\bf f}_{bnd}\f$ denote the forcing degrees of
     *  freedom of the elemental velocities on the boundary of the
     *  element, \f${\bf f}_{int}\f$ denote the forcing degrees of
     *  freedom of the elemental velocities on the interior of the
     *  element. (see detailed description for more details).  
     * 
     * This vector is further manipulated into 
     *
     * \f[ Fh\_{bnd} = \left [ \begin{array}{c} f\_{bnd} -m\_B \,
     * m\_Cinv\, {\bf f}_{int}\\ \left [m\_D\_{int} \, m\_Cinv \,{\bf
     * f}_{int} \right]_0 \end{array}\right ] \hspace{1cm} [Fh\_p]_{i} =
     * \begin{array}{c} [m\_D\_{int} \, m\_Cinv \, {\bf
     * f}_{int}]_{i+1} \end{array} \f]
     *
     * where \f$-{[}m\_D\_{int}^T\, m\_Cinv \,{\bf f}_{int}]_0\f$
     * which is corresponds to the mean mode of the pressure degrees
     * of freedom and is now added to the boundary system and the
     * remainder of the block becomes the interior forcing for the
     * inner static condensation (see detailed description for more
     * details) which is set up in a GlobalLinSysDirectStaticCond
     * class.
     * 
     * Finally we perform the final maniplation of the forcing to 
     * using hte 
     * \f[ Fh\_{bnd} = Fh\_{bnd} - m\_Bh \,m\_Chinv \, Fh\_p \f]
     *
     * We can now call the solver to the global coupled boundary
     * system (through the call to #m_CoupledBndSys->Solve) to obtain
     * the velocity boundary solution as the mean pressure solution,
     * i.e.
     * 
     * \f[ {\cal A}^T(\hat{A} - \hat{C}^T \hat{D}^{-1} \hat{B} ){\cal
     * A} \, Bnd =  Fh\_{bnd} \f]
     * 
     * Once we know the solution to the above the rest of the pressure
     * modes are recoverable thorugh
     * 
     *  \f[ Ph = m\_Dhinv\, (Bnd  - m\_Ch^T \, Fh_{bnd}) \f]
     * 
     * We can now unpack \f$ Fh\_{bnd} \f$ (last elemental mode) and
     * \f$ Ph \f$ into #m_pressure and \f$ F_p\f$ and \f$ Fh\_{bnd}\f$
     * into a closed pack list of boundary velocoity degrees of
     * freedom stored in \f$ F\_bnd\f$.
     * 
     * Finally the interior velocity degrees of freedom are then
     * obtained through the relationship
     * 
     *   \f[ F\_{int} = m\_Cinv\ ( F\_{int} + m\_D\_{int}^T\,
     *   F\_p - m\_Btilde^T\, Bnd) \f]
     * 
     * We then unpack the solution back to the MultiRegion structures
     * #m_velocity and #m_pressure
     */
    void CoupledLinearNS::SolveLinearNS(const Array<OneD, Array<OneD, NekDouble> > &forcing)
    {
        int i,n;
        Array<OneD,  MultiRegions::ExpListSharedPtr> vel_fields(m_velocity.num_elements());
        Array<OneD, Array<OneD, NekDouble> > force(m_velocity.num_elements());

        if(m_HomogeneousType == eHomogeneous1D)
        {
            int ncoeffsplane = m_fields[m_velocity[0]]->GetPlane(0)->GetNcoeffs();
            for(n = 0; n < m_npointsZ/2; ++n)
            {
                // Get the a Fourier mode of velocity and forcing components. 
                for(i = 0; i < m_velocity.num_elements(); ++i)
                {
                    vel_fields[i] = m_fields[m_velocity[i]]->GetPlane(2*n);
                    // Note this needs to correlate with how we pass forcing
                    force[i] = forcing[i] + 2*n*ncoeffsplane;
                }
                
                SolveLinearNS(force,vel_fields,m_pressure->GetPlane(2*n),n);
            }
            for(i = 0; i < m_velocity.num_elements(); ++i)
            {
                m_fields[m_velocity[i]]->SetPhysState(false);
            }
            m_pressure->SetPhysState(false);
        }
        else
        {
            for(i = 0; i < m_velocity.num_elements(); ++i)
            {
                vel_fields[i] = m_fields[m_velocity[i]];
                // Note this needs to correlate with how we pass forcing
                force[i] = forcing[i];
            }
            
            SolveLinearNS(force,vel_fields,m_pressure);
        }
    }

    void CoupledLinearNS::SolveLinearNS(const Array<OneD, Array<OneD, NekDouble> > &forcing,  Array<OneD, MultiRegions::ExpListSharedPtr> &fields, MultiRegions::ExpListSharedPtr &pressure,  const int mode)
    {
        int i,j,k,n,eid,cnt,cnt1;
        int nbnd,nint,offset;
        int nvel = m_velocity.num_elements();
        int nel  = fields[0]->GetNumElmts();
        Array<OneD, unsigned int> bmap, imap; 

        Array<OneD, NekDouble > f_bnd(m_mat[mode].m_BCinv->GetRows());
        NekVector< NekDouble  > F_bnd(f_bnd.num_elements(), f_bnd, eWrapper);
        Array<OneD, NekDouble > f_int(m_mat[mode].m_BCinv->GetColumns());
        NekVector< NekDouble  > F_int(f_int.num_elements(),f_int, eWrapper);

        int nz_loc;
        int  nplanecoeffs = fields[m_velocity[0]]->GetNcoeffs();// this is fine since we pass the nplane coeff data. 
 
        if(mode) // Homogeneous mode flag
        {
            nz_loc = 2;
         }
        else
        {
            if(m_singleMode)
            {
                nz_loc = 2;
            }
            else
            {
                nz_loc = 1;
                if(m_HomogeneousType == eHomogeneous1D)
                {
                    // Zero fields to set complex mode to zero;
                    for(i = 0; i < fields.num_elements(); ++i)
                    {
                        Vmath::Zero(2*fields[i]->GetNcoeffs(),fields[i]->UpdateCoeffs(),1);
                    }
                    Vmath::Zero(2*pressure->GetNcoeffs(),pressure->UpdateCoeffs(),1);
                }
            }
        }
        
        // Assemble f_bnd and f_int
        cnt = cnt1 = 0;
        for(i = 0; i < nel; ++i) // loop over elements
        {
            eid = fields[m_velocity[0]]->GetOffset_Elmt_Id(i);
            fields[m_velocity[0]]->GetExp(eid)->GetBoundaryMap(bmap);
            fields[m_velocity[0]]->GetExp(eid)->GetInteriorMap(imap);
            nbnd   = bmap.num_elements();
            nint   = imap.num_elements();
            offset = fields[m_velocity[0]]->GetCoeff_Offset(eid);
            
            for(j = 0; j < nvel; ++j) // loop over velocity fields 
            {
                for(n = 0; n < nz_loc; ++n)
                {
                    for(k = 0; k < nbnd; ++k)
                    {
                        f_bnd[cnt+k] = forcing[j][n*nplanecoeffs + 
                                                  offset+bmap[k]];
                    }
                    for(k = 0; k < nint; ++k)
                    {
                        f_int[cnt1+k] = forcing[j][n*nplanecoeffs + 
                                                   offset+imap[k]];
                    }
                    cnt  += nbnd;
                    cnt1 += nint;
                }
            }
        }

        Array<OneD, NekDouble > f_p(m_mat[mode].m_D_int->GetRows());
        NekVector<  NekDouble > F_p(f_p.num_elements(),f_p,eWrapper);
        NekVector<  NekDouble > F_p_tmp(m_mat[mode].m_Cinv->GetRows());

        // fbnd does not currently hold the pressure mean
        F_bnd = F_bnd - (*m_mat[mode].m_BCinv)*F_int;
        F_p_tmp = (*m_mat[mode].m_Cinv)*F_int;
        F_p = (*m_mat[mode].m_D_int) * F_p_tmp;
        
        // construct inner forcing 
        Array<OneD, NekDouble > bnd   (m_locToGloMap[mode]->GetNumGlobalCoeffs(),0.0);
        Array<OneD, NekDouble > fh_bnd(m_locToGloMap[mode]->GetNumGlobalCoeffs(),0.0);
        
        const Array<OneD,const int>& loctoglomap
            = m_locToGloMap[mode]->GetLocalToGlobalMap();
        const Array<OneD,const NekDouble>& loctoglosign
            = m_locToGloMap[mode]->GetLocalToGlobalSign();

        
        offset = cnt = 0; 
        for(i = 0; i < nel; ++i)
        {
            eid  = fields[0]->GetOffset_Elmt_Id(i);
            nbnd = nz_loc*fields[0]->GetExp(eid)->NumBndryCoeffs(); 
            
            for(j = 0; j < nvel; ++j)
            {
                for(k = 0; k < nbnd; ++k)
                {
                    fh_bnd[loctoglomap[offset+j*nbnd+k]] += 
                        loctoglosign[offset+j*nbnd+k]*f_bnd[cnt+k];
                }
                cnt += nbnd;
            }

            nint    = pressure->GetExp(eid)->GetNcoeffs();
            offset += nvel*nbnd + nint*nz_loc; 
        }

        offset = cnt1 = 0; 
        for(i = 0; i <  nel; ++i)
        {
            eid  = fields[0]->GetOffset_Elmt_Id(i);
            nbnd = nz_loc*fields[0]->GetExp(eid)->NumBndryCoeffs(); 
            nint = pressure->GetExp(eid)->GetNcoeffs(); 

            for(n = 0; n < nz_loc; ++n)
            {
                for(j = 0; j < nint; ++j)
                {
                    fh_bnd[loctoglomap[offset + nvel*nbnd + n*nint+j]] = f_p[cnt1+j];
                }
                cnt1   += nint;
            }
            offset += nvel*nbnd + nz_loc*nint; 
        }

        //  Set Weak BC into f_bnd and Dirichlet Dofs in bnd
        const Array<OneD,const int>& bndmap
            = m_locToGloMap[mode]->GetBndCondCoeffsToGlobalCoeffsMap();
        
        // Forcing function with weak boundary conditions and
        // Dirichlet conditions
        int bndcnt=0;
        
        for(k = 0; k < nvel; ++k)
        {
            const Array<OneD, SpatialDomains::BoundaryConditionShPtr> bndConds = fields[k]->GetBndConditions();
            Array<OneD, const MultiRegions::ExpListSharedPtr> bndCondExp;
            if(m_HomogeneousType == eHomogeneous1D) 
            {
                bndCondExp = m_fields[k]->GetPlane(2*mode)->GetBndCondExpansions();
            }
            else
            {
                bndCondExp = m_fields[k]->GetBndCondExpansions();
            }

            for(i = 0; i < bndCondExp.num_elements(); ++i)
            {
                const Array<OneD, const NekDouble > bndCondCoeffs = bndCondExp[i]->GetCoeffs();
                cnt = 0;
                for(n = 0; n < nz_loc; ++n)
                {
                    if(bndConds[i]->GetBoundaryConditionType() 
                       == SpatialDomains::eDirichlet)
                    {
                        for(j = 0; j < (bndCondExp[i])->GetNcoeffs(); j++)
                        {
                            bnd[bndmap[bndcnt++]]
                                = bndCondCoeffs[cnt++];
                        }
                    }
                    else
                    {                    
                        for(j = 0; j < (bndCondExp[i])->GetNcoeffs(); j++)
                        {
                            fh_bnd[bndmap[bndcnt++]]
                                += bndCondCoeffs[cnt++];
                        }
                    }
                }
            }
        }
        
        m_mat[mode].m_CoupledBndSys->Solve(fh_bnd,bnd,m_locToGloMap[mode]);
        
        // unpack pressure and velocity boundary systems. 
        offset = cnt = 0; 
        int totpcoeffs = pressure->GetNcoeffs();
        Array<OneD, NekDouble> p_coeffs = pressure->UpdateCoeffs();
        for(i = 0; i <  nel; ++i)
        {
            eid  = fields[0]->GetOffset_Elmt_Id(i);
            nbnd = nz_loc*fields[0]->GetExp(eid)->NumBndryCoeffs(); 
            nint = pressure->GetExp(eid)->GetNcoeffs(); 
            
            for(j = 0; j < nvel; ++j)
            {
                for(k = 0; k < nbnd; ++k)
                {
                    f_bnd[cnt+k] = loctoglosign[offset+j*nbnd+k]*bnd[loctoglomap[offset + j*nbnd + k]];
                }
                cnt += nbnd;
            }
            offset += nvel*nbnd + nint*nz_loc;
        }

        pressure->SetPhysState(false);

        offset = cnt = cnt1 = 0;
        for(i = 0; i < nel; ++i)
        {
            eid  = fields[0]->GetOffset_Elmt_Id(i);
            nint = pressure->GetExp(eid)->GetNcoeffs(); 
            nbnd = fields[0]->GetExp(eid)->NumBndryCoeffs(); 
            cnt1 = pressure->GetCoeff_Offset(eid);
            
            for(n = 0; n < nz_loc; ++n)
            {
                for(j = 0; j < nint; ++j)
                {
                    p_coeffs[n*totpcoeffs + cnt1+j] = 
                        f_p[cnt+j] = bnd[loctoglomap[offset + 
                                                     (nvel*nz_loc)*nbnd + 
                                                     n*nint +j]];
                }
                cnt += nint;
            }
            offset += (nvel*nbnd + nint)*nz_loc;
        }

        // Back solve first level of static condensation for interior
        // velocity space and store in F_int
        F_int = F_int + Transpose(*m_mat[mode].m_D_int)*F_p
                                       - Transpose(*m_mat[mode].m_Btilde)*F_bnd;
        F_int = (*m_mat[mode].m_Cinv)*F_int;
    
        // Unpack solution from Bnd and F_int to v_coeffs 
        cnt = cnt1 = 0;
        for(i = 0; i < nel; ++i) // loop over elements
        {
            eid  = fields[m_velocity[0]]->GetOffset_Elmt_Id(i);
            fields[0]->GetExp(eid)->GetBoundaryMap(bmap);
            fields[0]->GetExp(eid)->GetInteriorMap(imap);
            nbnd   = bmap.num_elements();
            nint   = imap.num_elements();
            offset = fields[0]->GetCoeff_Offset(eid);

            for(j = 0; j < nvel; ++j) // loop over velocity fields 
            {
                for(n = 0; n < nz_loc; ++n)
                {
                    for(k = 0; k < nbnd; ++k)
                    {
                        fields[j]->SetCoeff(n*nplanecoeffs + 
                                            offset+bmap[k],
                                            f_bnd[cnt+k]);
                    }
                    
                    for(k = 0; k < nint; ++k)
                    {
                        fields[j]->SetCoeff(n*nplanecoeffs + 
                                            offset+imap[k],
                                            f_int[cnt1+k]);
                    }
                    cnt  += nbnd;
                    cnt1 += nint;
                }
            }
        }
        
        for(j = 0; j < nvel; ++j) 
        {
            fields[j]->SetPhysState(false);
        }
    }

    void CoupledLinearNS::Checkpoint_TecplotOutput(const int n, bool IsInPhysicalSpace)
    {
    	 int nq = m_fields[0]->GetTotPoints();

    	    std::string var = "";
    	    for(int j = 0; j < m_fields.num_elements(); ++j)
    	      {
    		var = var + ", " + m_boundaryConditions->GetVariable(j);
    	      }

    	    var = var + ", " + "p";

    	    char chkout[16] = "";
    	    sprintf(chkout, "%d", n);
    	    std::string outname = m_sessionName + "_" + chkout + ".dat";
    	    ofstream outfile(outname.c_str());

    	    // put inarray in m_phys
    	    if (IsInPhysicalSpace == false)
    	      {
    		for(int i = 0; i < m_fields.num_elements(); ++i)
    		  {
    		    m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(),m_fields[i]->UpdatePhys());
    		  }
    	      }


           // m_fields[0]->FwdTrans_IterPerExp(m_pressure->GetPhys(),fieldcoeffs[i]);

    	    m_fields[0]->WriteTecplotHeader(outfile,var);

    	    for(int i = 0; i < m_fields[0]->GetExpSize(); ++i)
    	      {
    		m_fields[0]->WriteTecplotZone(outfile,i);
    		for(int j = 0; j < m_fields.num_elements(); ++j)
    	    	  {
    	    	    m_fields[j]->WriteTecplotField(outfile,i);
    	    	  }
    		m_pressure->WriteTecplotField(outfile,i);
    	      }
    }

    void CoupledLinearNS::CheckPointOutput(const int n)
    {
        Array<OneD, Array<OneD, NekDouble > > fieldcoeffs(m_fields.num_elements()+1);
        Array<OneD, std::string> variables(m_fields.num_elements()+1);
        int i;

        for(i = 0; i < m_fields.num_elements(); ++i)
        {
            fieldcoeffs[i] = m_fields[i]->UpdateCoeffs();
            variables[i]   = m_boundaryConditions->GetVariable(i);
        }

        fieldcoeffs[i] = Array<OneD, NekDouble>(m_fields[0]->GetNcoeffs());

        // project pressure field to velocity space

        m_pressure->BwdTrans(m_pressure->GetCoeffs(), m_pressure->UpdatePhys());

        m_fields[0]->FwdTrans_IterPerExp(m_pressure->GetPhys(),fieldcoeffs[i]);
        variables[i] = "p";

        char chkout[16] = "";
	    char procout[16] = "";

        sprintf(chkout, "%d", n);
        std::string outname = m_sessionName + "_" + chkout + ".chk";

       // cout << "Adding p!!!!!!!!!!!" << endl;
        WriteFld(outname,m_fields[0],fieldcoeffs,variables);

    }




    void CoupledLinearNS::v_Output(void)
    {    
        Array<OneD, Array<OneD, NekDouble > > fieldcoeffs(m_fields.num_elements()+1);
        Array<OneD, std::string> variables(m_fields.num_elements()+1);
        int i;
        
        for(i = 0; i < m_fields.num_elements(); ++i)
        {
            fieldcoeffs[i] = m_fields[i]->UpdateCoeffs();
            variables[i]   = m_boundaryConditions->GetVariable(i);
        }
        
        fieldcoeffs[i] = Array<OneD, NekDouble>(m_fields[0]->GetNcoeffs());

        // project pressure field to velocity space
      
        m_pressure->BwdTrans(m_pressure->GetCoeffs(), m_pressure->UpdatePhys());
        
        m_fields[0]->FwdTrans_IterPerExp(m_pressure->GetPhys(),fieldcoeffs[i]);
        variables[i] = "p"; 
  
        std::string outname = m_sessionName + ".fld";

        WriteFld(outname,m_fields[0],fieldcoeffs,variables);
  
    }

    int CoupledLinearNS::v_GetForceDimension()
    {
           return m_session->GetVariables().size();
    }

    void CoupledLinearNS::ImportFldPressure()
    {

    if (m_session->DefinesFunction("InitialConditions"))
        {
            // Check for restart file.
            if (m_session->GetFunctionType("InitialConditions",0)
                    == LibUtilities::eFunctionTypeFile)
            {
                std::string restartfile
                            = m_session->GetFunctionFilename("InitialConditions",0);
                cout << "\tRestart file: "<< restartfile << endl;

                std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef;
                std::vector<std::vector<NekDouble> > FieldData;
 
                LibUtilities::Import(restartfile,FieldDef,FieldData);

                int j = m_fields.num_elements();
                // last FieldDate is pressure (written last to Fld file)
                // Copy FieldData into m_fields

				Vmath::Zero(m_pressure->GetNcoeffs(),m_pressure->UpdateCoeffs(),1);

				for(int i = 0; i < FieldDef.size(); ++i)
				{
					m_pressure->ExtractDataToCoeffs(FieldDef[i], FieldData[i],
													FieldDef[i]->m_fields[j],
													m_pressure->UpdateCoeffs());
				}
				m_pressure->BwdTrans(m_pressure->GetCoeffs(),
									m_pressure->UpdatePhys());

            }
        }
    }


    void CoupledLinearNS::UpdatePressure()
    {
    	// project pressure field to velocity space
    	m_pressure->BwdTrans(m_pressure->GetCoeffs(), m_pressure->UpdatePhys());
    	Array<OneD, NekDouble> pressurephys = m_pressure->GetPhys();

        // Get Expansion list for orthogonal expansion at p-2
          const SpatialDomains::ExpansionMap &pressure_exp = GenPressureExp(m_graph->GetExpansions("u"));

     	 MultiRegions::ExpList2DSharedPtr field;

   	     if(!(field = boost::dynamic_pointer_cast<
   	        			   MultiRegions::ExpList2D>(m_pressure)))
		  {
			  ASSERTL0(false,"Dynamics cast failed");
		  }

//     	 field->UpdateExpList2D(m_session,pressure_exp);

         // m_pressure = MemoryManager<MultiRegions::ExpList2D>::AllocateSharedPtr(m_comm, pressure_exp);
          //SetUp2DExpansionC0ContMap();

          m_pressure->SetPhys(pressurephys);
			// Fill Coeffs
			m_pressure->FwdTrans(m_pressure->GetPhys(),m_pressure->UpdateCoeffs());
			m_pressure->BwdTrans(m_pressure->GetCoeffs(),m_pressure->UpdatePhys());
    }
  
}

/**
* $Log: CoupledLinearNS.cpp,v $
**/
