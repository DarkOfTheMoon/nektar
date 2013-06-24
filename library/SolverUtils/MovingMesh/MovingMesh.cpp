///////////////////////////////////////////////////////////////////////////////
//
// File AdvectionTerm.cpp
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
// Description: Base class for Navier-Stokes advection term
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/MovingMesh/MovingMesh.h>
#include <cstdio>
#include <cstdlib>

#include <cmath>

#include <string>
namespace Nektar
{

	MovingMeshFactory& GetMovingMeshFactory()
	{
		typedef Loki::SingletonHolder<MovingMeshFactory,
			Loki::CreateUsingNew > Type;
		return Type::Instance();
	}

	string MovingMesh::className  = GetMovingMeshFactory().RegisterCreatorFunction("Elliptic", MovingMesh::create);

    /**
     * Constructor. Creates ...
     *
     * \param
     * \param
     */

    MovingMesh::MovingMesh(const LibUtilities::SessionReaderSharedPtr &pSession):
		    UnsteadySystem(pSession),
            m_session(pSession)
	{

	}

    void MovingMesh::v_InitObject()
    {
    	int i,j;
              m_spacedim =  m_graph->GetSpaceDimension();
              int numfields = m_fields.num_elements();
              int ncoords =3;
              std::string velids[] = {"wx","wy"};

              // Set up Velocity field to point to the first m_expdim of m_fields;
              m_meshvelocity = Array<OneD,int>(m_spacedim);

              for(i = 0; i < m_spacedim; ++i)
              {
                  for(j = 0; j < numfields; ++j)
                  {
                      std::string var = m_boundaryConditions->GetVariable(j);
                      if(NoCaseStringCompare(velids[i],var) == 0)
                      {
                          m_meshvelocity[i] = j;
                          break;
                      }

                      if(j == numfields)
                      {
                          std::string error = "Failed to find field: " + var;
                          ASSERTL0(false,error.c_str());
                      }
                  }
              }

              std::string velids2[] = {"u","v"};

              // Set up Velocity field to point to the first m_expdim of m_fields;
              m_velocity = Array<OneD,int>(m_spacedim);

              for(i = 0; i < m_spacedim; ++i)
              {
                  for(j = 0; j < numfields; ++j)
                  {
                      std::string var = m_session->GetVariable(j);
                      if(NoCaseStringCompare(velids2[i],var) == 0)
                      {
                          m_velocity[i] = j;
                          break;
                      }

                      if(j == numfields)
                      {
                          std::string error = "Failed to find field: " + var;
                          ASSERTL0(false,error.c_str());
                      }
                  }
              }

        	//  m_meshvelocity   = Array<OneD, MultiRegions::ExpListSharedPtr>(m_spacedim);
        	  m_meshcoords= Array <OneD, Array<OneD, NekDouble> >(ncoords);
        	  m_meshcoordsold= Array <OneD, Array<OneD, NekDouble> >(ncoords);

    		  if(!(m_mesh2D = boost::dynamic_pointer_cast<
    			   SpatialDomains::MeshGraph2D>(m_graph)))
    		  {
    			  ASSERTL0(false,"Dynamics cast failed");
    		  }

            m_phystot = m_fields[0]->GetTotPoints();
            m_meshcoords[0] = Array<OneD, NekDouble> (m_phystot*ncoords);
            m_meshcoordsold[0] = Array<OneD, NekDouble> (m_phystot*ncoords);
            for(i = 1; i < ncoords; ++i)
    	   {
            	m_meshcoords[i] = m_meshcoords[i-1] + m_phystot;
            	m_meshcoordsold[i] = m_meshcoordsold[i-1] + m_phystot;
    	   }


            m_fields[0]->GetCoords(m_meshcoords[0],m_meshcoords[1],m_meshcoords[2]);
            m_fields[0]->GetCoords(m_meshcoordsold[0],m_meshcoordsold[1],m_meshcoordsold[2]);

           // m_DisField = MemoryManager<MultiRegions::DisContField2D>::AllocateSharedPtr(m_comm,*m_mesh2D,*m_boundaryConditions,m_meshvelocity[0]);
            const std::string var = "wx";
            m_ExpField = MemoryManager<MultiRegions::ExpList2D>::AllocateSharedPtr(m_session,
                           m_graph,
                            true,
                            m_session->GetVariable(m_meshvelocity[0]));

            int intSteps = 3;
            m_meshvelwx =  Array<OneD, Array<OneD, NekDouble> >  (intSteps);
            m_meshvelwy =  Array<OneD, Array<OneD, NekDouble> >  (intSteps);
            for(int i=0;i<intSteps;i++)
            {
            	 m_meshvelwx[i] = Array<OneD, NekDouble>(m_phystot);
            	 m_meshvelwy[i] = Array<OneD, NekDouble>(m_phystot);
            }

            int BCintSteps = 1;
            int n,cnt;
    		Array<OneD, MultiRegions::ExpListSharedPtr>  BndExp     = m_fields[m_meshvelocity[0]]->GetBndCondExpansions();
    		m_FBCwx =  Array<OneD, Array<OneD, NekDouble> >  (BCintSteps);
    		m_FBCwy =  Array<OneD, Array<OneD, NekDouble> >  (BCintSteps);

    		for(cnt = n = 0; n < m_fields[0]->GetBndConditions().num_elements(); ++n)
    		{
    		  string type = m_fields[0]->GetBndConditions()[n]->GetUserDefined().GetEquation();
    		  if(type == "FreeSurface")
    		  {
    			 cnt += BndExp[n]->GetTotPoints();;
    		  }
    		}
    		for(int i=0;i<BCintSteps;i++)
    		{
    			 m_FBCwx[i] = Array<OneD, NekDouble>(cnt);
    			 m_FBCwy[i] = Array<OneD, NekDouble>(cnt);
    		}

            //Setting up free surface spline
            SetUpFreeSurfaceSplines();

    //        SetInitialMeshSloshing();
      /*      SetInitialMeshDieSwell();
            UpdateMesh();
            UpdateFields(); */
    }


	MovingMesh::~MovingMesh()
	{
	}

	void MovingMesh::SetUpFreeSurfaceSplines()
	{
		int cnt,n;
		Array<OneD, MultiRegions::ExpListSharedPtr>  BndExp     = m_fields[m_meshvelocity[0]]->GetBndCondExpansions();
		Array<OneD, const SpatialDomains::BoundaryConditionShPtr > BndConds = m_fields[m_meshvelocity[0]]->GetBndConditions();

		// Count number of free surface regions
        for(cnt = n = 0; n < BndConds.num_elements(); ++n)
         {
             string type = BndConds[n]->GetUserDefined().GetEquation();
     //for the ones marked as Drag calculate drag along boundary edge
             if(type == "FreeSurface")
             {
            	 cnt++;
             }
         }

        m_freesurfacesplines = Array<OneD, LibUtilities::CubicSplineSharedPtr>(cnt);
        LibUtilities::SplineBoundaryType bcleft,bcright;
		bcleft = LibUtilities::eClampedQuadraticLagrange;
		bcright = LibUtilities::eNotaKnot;
		int bcl=0,bcr=0;

        for(cnt = n = 0; n < BndConds.num_elements(); ++n)
	   {
		   string type = BndConds[n]->GetUserDefined().GetEquation();
		   //for the ones marked as Drag calculate drag along boundary edge
		   if(type == "FreeSurface")
		   {
			   m_freesurfacesplines[cnt] = MemoryManager<LibUtilities::CubicSpline>::AllocateSharedPtr();
			   m_freesurfacesplines[cnt] = m_freesurfacesplines[cnt]->SetUpSpline(BndExp[n],bcleft,bcright,bcl,bcr);
			 //  m_freesurfacesplines[cnt] = m_freesurfacesplines[cnt]->SetUpSplineCornerPoints(BndExp[n]);
			   cnt++;
		   }
	   }
	}

	// Set up dirichlet bc at all boundaries and userdefined free surface
/*	void ALE::GenerateBoundaryConditionsMeshVelocity(Array<OneD, MultiRegions::ExpListSharedPtr > m_fields)
	{
		Array<OneD, const SpatialDomains::BoundaryConditionShPtr > BndConds= m_fields[0]->GetBndConditions();

		// iterate over boundary regions
		for(cnt = n = 0; n < m_fields[0]->GetBndConditions().num_elements(); ++n)
	   {


	   }


	  	std::string equation, userDefined;
	        	for(i = 0; i < in.m_boundaryRegions.size();i++)
	        	{
	        	    BoundaryConditionMapShPtr newboundaryConditions = MemoryManager<BoundaryConditionMap>::AllocateSharedPtr();
	        		for (Variable::iterator varIter = in.m_variables.begin();
	        		                                varIter != in.m_variables.end(); ++varIter)
	        		{
	        			BoundaryConditionType type = (*in.m_boundaryConditions[i])[*varIter]->GetBoundaryConditionType();
	        			LibUtilities::Equation userdefined = (*in.m_boundaryConditions[i])[*varIter]->GetUserDefined();
	        			string userdefinedeqn = userdefined.GetEquation();

	        			// Problem might be that eqntype is not build up clearly
	        			(*newboundaryConditions)[*varIter] = (*in.m_boundaryConditions[i])[*varIter];
	        		}

	        		(*returnval).m_boundaryConditions[i] = newboundaryConditions;

	        	}
	} */

	 void MovingMesh::DoMovingMeshInitialise(const NekDouble timestep)
	 {
	    	int outputcount = 0;

	    	//WriteValuesAlongFreeSurfaceToFile("afterHelmholtz",outputcount);

	    	// Reshuffle wx time vector
			int  nint    = 1;
			int  nlevels = m_meshvelwx.num_elements();
			m_meshvelwx[0] = m_fields[m_meshvelocity[0]]->GetPhys();
			m_meshvelwy[0] = m_fields[m_meshvelocity[1]]->GetPhys();

	        // Extrapolate to n+1
	        Vmath::Smul(m_phystot,kFreeSurfaceBCsExtrapolation[nint-1][nint-1],
	        		m_meshvelwx[nint-1],1,m_meshvelwx[nlevels-1],1);
	        Vmath::Smul(m_phystot,kFreeSurfaceBCsExtrapolation[nint-1][nint-1],
	                		m_meshvelwy[nint-1],1,m_meshvelwy[nlevels-1],1);

	        DeformMesh2(timestep);
	        DeformMeshBoundary(timestep);

	        // Updates Mesh, Fields and Coordinates of Mesh
	        UpdateMesh();
	        // Compute Mesh Velocity from Coordinates of Mesh
	        //ComputeMeshVelocity(aii_Dt);

	        //WriteValuesAlongFreeSurfaceToFile("ComputedValues",outputcount);
	}


	void MovingMesh::DoMovingMesh(NekDouble timestep,NekDouble time)
    {
    	static int outputcount = 1;

        //CreateFreeSurfaceFunction(aii_Dt);
    	CreateFreeSurfaceSpline(timestep, time);
    	//SetOutflowBCALE();
    	//SetPlateBCALE();

    	//WriteValuesAlongFreeSurfaceToFile("beforeHelmholtz",outputcount);

        //int phystot = m_fields[m_meshvelocity[0]]->GetTotPoints();
        Array<OneD,NekDouble> forcing(m_phystot,0.0);

        // At this place we could solve a poisson equation with RHS=0
    	for(int i=0;i<m_spacedim;i++)
		{
    		// Set Meshvelocity to fluid velocity for caber pulling phase
			//Vmath::Vcopy(m_fields[m_velocity[i]]->GetTotPoints(),m_fields[m_velocity[i]]->GetPhys(),1,m_fields[m_meshvelocity[i]]->UpdatePhys(),1);
			//m_fields[m_meshvelocity[i]]->FwdTrans(m_fields[m_meshvelocity[i]]->GetPhys(),m_fields[m_meshvelocity[i]]->UpdateCoeffs());
    		m_fields[m_meshvelocity[i]]->HelmSolve(forcing,m_fields[m_meshvelocity[i]]->UpdateCoeffs(),0.0);
			m_fields[m_meshvelocity[i]]->BwdTrans(m_fields[m_meshvelocity[i]]->GetCoeffs(),m_fields[m_meshvelocity[i]]->UpdatePhys());
		}

    	//WriteValuesAlongFreeSurfaceToFile("afterHelmholtz",outputcount);

    	// Reshuffle wx time vector
    	static int ncalls = 1;
		Array<OneD, NekDouble> tmpx,tmpy;
		int  nlevels = m_meshvelwx.num_elements();
		int  nint    = min(ncalls++,nlevels);
		// Reshuffle Bc Storage vector
		tmpx = m_meshvelwx[nlevels-1];
		tmpy = m_meshvelwy[nlevels-1];
		for(int n = nlevels-1; n > 0; --n)
		{
			m_meshvelwx[n] = m_meshvelwx[n-1];
			m_meshvelwy[n] = m_meshvelwy[n-1];
		}
		m_meshvelwx[0] = m_fields[m_meshvelocity[0]]->GetPhys();
		m_meshvelwy[0] = m_fields[m_meshvelocity[1]]->GetPhys();

        // Extrapolate to n+1
        Vmath::Smul(m_phystot,kFreeSurfaceBCsExtrapolation[nint-1][nint-1],
        		m_meshvelwx[nint-1],1,m_meshvelwx[nlevels-1],1);
        Vmath::Smul(m_phystot,kFreeSurfaceBCsExtrapolation[nint-1][nint-1],
                		m_meshvelwy[nint-1],1,m_meshvelwy[nlevels-1],1);
        for(int n = 0; n < nint-1; ++n)
        {
            Vmath::Svtvp(m_phystot,kFreeSurfaceBCsExtrapolation[nint-1][n],
            		m_meshvelwx[n],1,m_meshvelwx[nlevels-1],1,
            		m_meshvelwx[nlevels-1],1);
            Vmath::Svtvp(m_phystot,kFreeSurfaceBCsExtrapolation[nint-1][n],
                      		m_meshvelwy[n],1,m_meshvelwy[nlevels-1],1,
                      		m_meshvelwy[nlevels-1],1);
        }

        DeformMesh2(timestep);
        UpdateMeshFreeSurfaceBoundary(timestep);
       // DeformMeshBoundary(aii_Dt);

      //  DeformMeshBoundary(aii_Dt);

        // Updates Mesh, Fields and Coordinates of Mesh
        UpdateMesh();
        // Compute Mesh Velocity from Coordinates of Mesh
        ComputeMeshVelocity(timestep);

        //WriteValuesAlongFreeSurfaceToFile("ComputedValues",outputcount);
    	//UpdateMeshandFields();
        //WARNING: This Produced a Zero Meshvelocity at free surface boundary!!!!!!!!!!
        // But might be necessary to avoid jumps between elements
    /*    for(int i=0;i<m_spacedim;i++)
		{
        	m_fields[m_meshvelocity[i]]->FwdTrans(m_fields[m_meshvelocity[i]]->GetPhys(),m_fields[m_meshvelocity[i]]->UpdateCoeffs());
        	m_fields[m_meshvelocity[i]]->BwdTrans(m_fields[m_meshvelocity[i]]->GetCoeffs(),m_fields[m_meshvelocity[i]]->UpdatePhys());
		} */
        outputcount++;
    }

	void MovingMesh::DoALEAxiSymmetric(const NekDouble timestep,
			const NekDouble time)
	{
		static int outputcount = 1;

		//CreateFreeSurfaceFunction(aii_Dt);
		CreateFreeSurfaceSpline(timestep, time);
		//SetOutflowBCALE();
		//SetPlateBCALE();

		//WriteValuesAlongFreeSurfaceToFile("beforeHelmholtz",outputcount);

		//int phystot = m_fields[m_meshvelocity[0]]->GetTotPoints();
		Array<OneD,NekDouble> forcing(m_phystot,0.0);
		Array<OneD,NekDouble> r(m_phystot,0.0);
		Array<OneD,NekDouble> z(m_phystot,0.0);
		Array<OneD,NekDouble> dummy(m_phystot,0.0);
		m_fields[m_meshvelocity[0]]->GetCoords(z,r,dummy);

		// At this place we could solve a poisson equation with RHS=0
		for(int i=0;i<m_spacedim;i++)
		{
			// Set Meshvelocity to fluid velocity for caber pulling phase
			//Vmath::Vcopy(m_fields[m_velocity[i]]->GetTotPoints(),m_fields[m_velocity[i]]->GetPhys(),1,m_fields[m_meshvelocity[i]]->UpdatePhys(),1);
			//m_fields[m_meshvelocity[i]]->FwdTrans(m_fields[m_meshvelocity[i]]->GetPhys(),m_fields[m_meshvelocity[i]]->UpdateCoeffs());
			//m_fields[m_meshvelocity[i]]->HelmSolve(forcing,m_fields[m_meshvelocity[i]]->UpdateCoeffs(),0.0);
			Vmath::Vmul(m_phystot,r,1,m_fields[m_meshvelocity[i]]->GetPhys(),1,m_fields[m_meshvelocity[i]]->UpdatePhys(),1);
			m_fields[m_meshvelocity[i]]->FwdTrans(m_fields[m_meshvelocity[i]]->GetPhys(),m_fields[m_meshvelocity[i]]->UpdateCoeffs());
			m_fields[m_meshvelocity[i]]->BwdTrans(m_fields[m_meshvelocity[i]]->GetCoeffs(),m_fields[m_meshvelocity[i]]->UpdatePhys());
		}

		//WriteValuesAlongFreeSurfaceToFile("afterHelmholtz",outputcount);

		// Reshuffle wx time vector
		static int ncalls = 1;
		Array<OneD, NekDouble> tmpx,tmpy;
		int  nlevels = m_meshvelwx.num_elements();
		int  nint    = min(ncalls++,nlevels);
		// Reshuffle Bc Storage vector
		tmpx = m_meshvelwx[nlevels-1];
		tmpy = m_meshvelwy[nlevels-1];
		for(int n = nlevels-1; n > 0; --n)
		{
			m_meshvelwx[n] = m_meshvelwx[n-1];
			m_meshvelwy[n] = m_meshvelwy[n-1];
		}
		m_meshvelwx[0] = m_fields[m_meshvelocity[0]]->GetPhys();
		m_meshvelwy[0] = m_fields[m_meshvelocity[1]]->GetPhys();

		// Extrapolate to n+1
		Vmath::Smul(m_phystot,kFreeSurfaceBCsExtrapolation[nint-1][nint-1],
				m_meshvelwx[nint-1],1,m_meshvelwx[nlevels-1],1);
		Vmath::Smul(m_phystot,kFreeSurfaceBCsExtrapolation[nint-1][nint-1],
						m_meshvelwy[nint-1],1,m_meshvelwy[nlevels-1],1);
		for(int n = 0; n < nint-1; ++n)
		{
			Vmath::Svtvp(m_phystot,kFreeSurfaceBCsExtrapolation[nint-1][n],
					m_meshvelwx[n],1,m_meshvelwx[nlevels-1],1,
					m_meshvelwx[nlevels-1],1);
			Vmath::Svtvp(m_phystot,kFreeSurfaceBCsExtrapolation[nint-1][n],
							m_meshvelwy[n],1,m_meshvelwy[nlevels-1],1,
							m_meshvelwy[nlevels-1],1);
		}

        DeformMesh2(timestep);
        UpdateMeshFreeSurfaceBoundary(timestep);
       // DeformMeshBoundary(aii_Dt);

      //  DeformMeshBoundary(aii_Dt);

        // Updates Mesh, Fields and Coordinates of Mesh
        UpdateMesh();
        // Compute Mesh Velocity from Coordinates of Mesh
        ComputeMeshVelocity(timestep);

		//WriteValuesAlongFreeSurfaceToFile("ComputedValues",outputcount);
		//UpdateMeshandFields();
		//WARNING: This Produced a Zero Meshvelocity at free surface boundary!!!!!!!!!!
		// But might be necessary to avoid jumps between elements
	/*    for(int i=0;i<m_spacedim;i++)
		{
			m_fields[m_meshvelocity[i]]->FwdTrans(m_fields[m_meshvelocity[i]]->GetPhys(),m_fields[m_meshvelocity[i]]->UpdateCoeffs());
			m_fields[m_meshvelocity[i]]->BwdTrans(m_fields[m_meshvelocity[i]]->GetCoeffs(),m_fields[m_meshvelocity[i]]->UpdatePhys());
		} */
		outputcount++;
	}

/*
 * Updates m_mesh2D with new mesh information
 */
    void MovingMesh::SetMeshVelocityatBoundary(const NekDouble timestep)
    {
        Array<OneD, const SpatialDomains::BoundaryConditionShPtr > BndConds;
        Array<OneD, MultiRegions::ExpListSharedPtr>  BndExp, BndExpy;

        // Grab boundary regions and expansions of velocity u
        //BndConds   = m_fields[0]->GetBndConditions();
        //BndExp     = m_fields[0]->GetBndCondExpansions();
        BndConds   = m_fields[m_meshvelocity[0]]->GetBndConditions();
        BndExp     = m_fields[m_meshvelocity[0]]->GetBndCondExpansions();
        BndExpy     = m_fields[m_meshvelocity[1]]->GetBndCondExpansions();

        StdRegions::StdExpansion1DSharedPtr EdgeExp;
        SpatialDomains::SegGeomSharedPtr edge;
        StdRegions::StdExpansionSharedPtr ElExp;
        Array<OneD, Array<OneD, NekDouble> > normals, locnormals,weaknormals;
		int vertID1,vertID2;
		SpatialDomains::VertexComponentSharedPtr vertex1,vertex2;
		SpatialDomains::MeshGraph2DSharedPtr mesh2D;

		MultiRegions::ExpList1DSharedPtr  FreeSurfaceFct;

        Array<OneD, int> ElmtID,EdgeID;
        int cnt,n,el,edgeID,elmID;
        Array<OneD, const NekDouble> Uphyselement, Vphyselement, wxel,wyel;

        //Fills ElmtID and EdgeID with global ids of elements id number and edges id number of the boundary expansion
        m_fields[m_meshvelocity[0]]->GetBoundaryToElmtMap(ElmtID,EdgeID);

        int cnt2=0;
        for(cnt = n = 0; n < m_fields[0]->GetBndConditions().num_elements(); ++n)
         {
             string type = m_fields[0]->GetBndConditions()[n]->GetUserDefined().GetEquation();
     //for the ones marked as Drag calculate drag along boundary edge
             if(type == "FreeSurface")
             {
            	 int nq= BndExp[n]->GetTotPoints();
            	 normals = Array<OneD, Array<OneD, NekDouble> > (2);
				 for (int k = 0; k < 2; ++k)
				 {
					normals[k] = Array<OneD, NekDouble>(nq,0.0);
				 }

			   for(el = 0; el < BndExp[n]->GetExpSize(); ++el,cnt2++)
				{
				   int offset = BndExp[n]->GetPhys_Offset(el);
				   int nqel = BndExp[n]->GetExp(el)->GetTotPoints();
				   EdgeExp =  boost::dynamic_pointer_cast<StdRegions::StdExpansion1D> (BndExp[n]->GetExp(el));
				   ElExp = m_fields[m_meshvelocity[0]]->GetExp(ElmtID[cnt2]);
				   //EdgeExp->SetUpPhysNormals(ElExp,EdgeID[cnt2]);
				   //locnormals	=  EdgeExp->GetMetricInfo()->GetNormal();
				   ElExp->ComputeEdgeNormal(EdgeID[cnt2]);
				   locnormals = ElExp->GetEdgeNormal(EdgeID[cnt2]);
				   for (int k = 0; k < 2; ++k)
					{
						Vmath::Vcopy(nqel, &(locnormals[k][0]), 1,
												&(normals[k][offset]), 1);
					}
					int nquad_e = EdgeExp->GetNumPoints(0);
		            Array<OneD,NekDouble> xedge(nquad_e,0.0);
					Array<OneD,NekDouble> yedge(nquad_e,0.0);
					Array<OneD,NekDouble> zedge(nquad_e,0.0);

					EdgeExp->GetCoords(xedge,yedge,zedge);

					/*cout<< "Before smoothing" << endl;
				   					for(int i=0;i<nquad_e;i++)
				   					{
				   						//if(xedge[i]==11.0 && yedge[i]<=0.0)
				   														{
				   						cout << i << ": (" << xedge[i] << "," << yedge[i] << ")"  << "\t";
				   						cout << "(n1,n2) = (" << locnormals[0][i] << ", " << locnormals[1][i] << ")" << endl;
				   														}
				   					}*/

				}

			/*	if(!(FreeSurfaceFct =  boost::dynamic_pointer_cast<MultiRegions::ExpList1D>(BndExp[n])))
					{
						ASSERTL0(false,"Dynamics cast failed");
					}
					 MultiRegions::ContField1DSharedPtr contfield =
				    	                  	 				MemoryManager<MultiRegions::ContField1D>
				    	                  	 				::AllocateSharedPtr(m_comm,*FreeSurfaceFct,m_solnType);
					 contfield->FwdTrans(normals[0],contfield->UpdateCoeffs());
					 contfield->BwdTrans(contfield->GetCoeffs(),normals[0]);
					 contfield->FwdTrans(normals[1],contfield->UpdateCoeffs());
					 contfield->BwdTrans(contfield->GetCoeffs(),normals[1]); */

            	 //loop over all elements along the boundary region
            	bool first=true;
				for(el = 0; el < BndExp[n]->GetExpSize(); ++el,cnt++)
				{
					Uphyselement = (m_fields[0]->GetPhys())+ m_fields[0]->GetPhys_Offset(ElmtID[cnt]);
					Vphyselement = (m_fields[1]->GetPhys())+ m_fields[1]->GetPhys_Offset(ElmtID[cnt]);
					wxel = (m_fields[m_meshvelocity[0]]->GetPhys())+ m_fields[m_meshvelocity[0]]->GetPhys_Offset(ElmtID[cnt]);
					wyel = (m_fields[m_meshvelocity[1]]->GetPhys())+ m_fields[m_meshvelocity[1]]->GetPhys_Offset(ElmtID[cnt]);

					ElExp = m_fields[m_meshvelocity[0]]->GetExp(ElmtID[cnt]);
				/*	int nq = locExp->GetTotPoints();
					Array<OneD,NekDouble> xel(nq,0.0);
					Array<OneD,NekDouble> yel(nq,0.0);
					Array<OneD,NekDouble> zel(nq,0.0);
					ElExp->GetCoords(xel,yel,zel);
					for(int i=0;i<nq;i++)
					{
						cout << i << ": (" << xel[i] << "," << yel[i] << ")=";
						cout << "(u,v)=(" << Uphyselement[i] <<","<< Vphyselement[i] << ")"<< endl;
					} */

					EdgeExp =  boost::dynamic_pointer_cast<StdRegions::StdExpansion1D> (BndExp[n]->GetExp(el));
					int nquad_e = EdgeExp->GetNumPoints(0);
		            Array<OneD,NekDouble> xedge(nquad_e,0.0);
					Array<OneD,NekDouble> yedge(nquad_e,0.0);
					Array<OneD,NekDouble> zedge(nquad_e,0.0);
					Array<OneD,NekDouble> Uphysedge(nquad_e,0.0);
					Array<OneD,NekDouble> Vphysedge(nquad_e,0.0);
					Array<OneD,NekDouble> wxedge(nquad_e,0.0);
					Array<OneD,NekDouble> wyedge(nquad_e,0.0);

					// EdgeID in within element
					edgeID = EdgeID[cnt];
					ElExp->GetEdgePhysVals(edgeID,EdgeExp,Uphyselement,Uphysedge);
					ElExp->GetEdgePhysVals(edgeID,EdgeExp,Vphyselement,Vphysedge);
					ElExp->GetEdgePhysVals(edgeID,EdgeExp,wxel,wxedge);
					ElExp->GetEdgePhysVals(edgeID,EdgeExp,wyel,wyedge);

					EdgeExp->GetCoords(xedge,yedge,zedge);

					//EdgeExp->SetUpPhysNormals(ElExp,edgeID);

					//locnormals	=  EdgeExp->GetMetricInfo()->GetNormal();
					  ElExp->ComputeEdgeNormal(EdgeID[cnt2]);
					  normals = ElExp->GetEdgeNormal(EdgeID[cnt2]);
					int offset = BndExp[n]->GetPhys_Offset(el);
					 int nqel = BndExp[n]->GetExp(el)->GetTotPoints();
					 for (int k = 0; k < 2; ++k)
					{
						Vmath::Vcopy(nqel, &(normals[k][offset]), 1,
								&(locnormals[k][0]), 1);
					}
					if(ElExp->GetEorient(edgeID) == StdRegions::eBackwards)
					{
						Vmath::Neg(nquad_e,locnormals[0],1);
						Vmath::Neg(nquad_e,locnormals[1],1);
					}

				/*	cout<< "After smoothing" << endl;
					for(int i=0;i<nquad_e;i++)
					{
						//if(xedge[i]==11.0 && yedge[i]<=0.0)
														{
						cout << i << ": (" << xedge[i] << "," << yedge[i] << ")"  << "\t";
						cout << "(n1,n2) = (" << locnormals[0][i] << ", " << locnormals[1][i] << ")" << endl;
														}
					} */



					// compute wx/wy on the edge
					for(int i=0;i<nquad_e;i++)
					{
						//wxedge[i]= (Uphysedge[i]*locnormals[0][i]+Vphysedge[i]*locnormals[1][i])*locnormals[0][i];
						wxedge[i]=0.0;
						wyedge[i]= (Uphysedge[i]*locnormals[0][i]+Vphysedge[i]*locnormals[1][i])*locnormals[1][i];
					}

					// Copy result in boundary expansion
					 int id1  = BndExp[n]->GetPhys_Offset(el);


					Vmath::Vcopy(nquad_e,&wxedge[0], 1,&(BndExp[n]->UpdatePhys())[id1],1);
					Vmath::Vcopy(nquad_e,&wyedge[0], 1,&(BndExpy[n]->UpdatePhys())[id1],1);
					// MOVE edge, i.e. calculate X^{n+1} at free surface boundary
					// TODO: Implement time integration instead

				/*	if(first)
					{ */

					/*	for(int i=0;i<nquad_e;i++)
						{
							//if(xedge[i]==11.0 && yedge[i]<=0.0)
												    		{
							cout << i << ": (" << xedge[i] << "," << yedge[i] << ")"  << "\t";
							cout << "(n1,n2) = (" << locnormals[0][i] << ", " << locnormals[1][i] << ")" << endl;
							cout << "Edge: (u,v)=(" << Uphysedge[i] <<","<< Vphysedge[i] << ")"<< "\t";
							cout << "(wx,wy)=(" << wxedge[i] <<","<< wyedge[i] << ")"<< endl;
												    		}
						} */
				//	}


				}
            	BndExp[n]->FwdTrans_BndConstrained(BndExp[n]->GetPhys(),BndExp[n]->UpdateCoeffs());
            	BndExpy[n]->FwdTrans_BndConstrained(BndExpy[n]->GetPhys(),BndExpy[n]->UpdateCoeffs());


             }
             else if(type == "" || type == "TimeDependent")  // for all other boundary conditi
			 {
				 cnt += BndExp[n]->GetExpSize();
				 cnt2 += BndExp[n]->GetExpSize();
			 }
         }
    }

    /// Update Positions of inner vertices to new position
    void MovingMesh::DeformMesh(const NekDouble timestep)
    {
    	// number of elements
    	int  nlevels = m_meshvelwx.num_elements();
    	int eid, id, diff,el,j,nqel;
    	int vertId, edgeId;
    	StdRegions::StdExpansionSharedPtr ElExp;
    	int globaledgeID;
        SpatialDomains::SegGeomSharedPtr edge;
        StdRegions::StdExpansion1DSharedPtr EdgeExp;
		Array<OneD, Array<OneD, NekDouble> > normals, locnormals,weaknormals;
		Array<OneD, const NekDouble> wxel,wyel;
		int vertID1,vertID2;
		 SpatialDomains::VertexComponentSharedPtr vertex1,vertex2;
		 SpatialDomains::Geometry1DSharedPtr ElmtSegGeom;

		 // local expansion vector
		 const StdRegions::StdExpansionVector &locExpVector = *(m_DisField->GetExp());

		 int nel = locExpVector.size();

        // Loop over elemente and collect forward expansion
		int nexp = m_fields[m_meshvelocity[0]]->GetExpSize();
		int nquad_e,n,e,offset,phys_offset;
		Array<OneD,NekDouble> e_tmp;
		//Array<OneD, Array<OneD, StdRegions::StdExpansion1DSharedPtr> > elmtToTrace = m_traceMap->GetElmtToTrace();

    	// Loop over all elements
        for(el = 0; el < nel; ++el)
        {
        	// Get element id
            eid = m_fields[m_meshvelocity[0]]->GetOffset_Elmt_Id(el);
            //int nqel =

            //ElExp = locExpVector[eid]; element standard expansion
            for(j = 0; j < locExpVector[eid]->GetNedges(); ++j)
            {
                vertId = (locExpVector[eid]->GetGeom2D())->GetVid(j);
                // global ID
                edgeId = (locExpVector[eid]->GetGeom2D())->GetEid(j);
                nquad_e = locExpVector[eid]->GetEdgeNumPoints(j);
                nqel  = locExpVector[eid]->GetTotPoints();

                Array<OneD,NekDouble> xel(nqel,0.0);
				Array<OneD,NekDouble> yel(nqel,0.0);
				Array<OneD,NekDouble> zel(nqel,0.0);

				locExpVector[eid]->GetCoords(xel,yel,zel);

				LocalRegions::Expansion2DSharedPtr exp2d = boost::dynamic_pointer_cast<LocalRegions::Expansion2D>(m_DisField->GetExp(eid));

               // wxel = (m_fields[m_meshvelocity[0]]->GetPhys())+ m_fields[m_meshvelocity[0]]->GetPhys_Offset(eid);
               // wyel = (m_fields[m_meshvelocity[1]]->GetPhys())+ m_fields[m_meshvelocity[1]]->GetPhys_Offset(eid);

				wxel = m_meshvelwx[nlevels-1] + m_fields[m_meshvelocity[0]]->GetPhys_Offset(eid);
				wyel = m_meshvelwy[nlevels-1] + m_fields[m_meshvelocity[1]]->GetPhys_Offset(eid);

                EdgeExp = exp2d->GetEdgeExp(j);
               // ElmtSegGeom  = (locExpVector[eid]->GetGeom2D())->GetEdge(j);

                Array<OneD,NekDouble> Uedge(nquad_e,0.0);
                Array<OneD,NekDouble> xedge(nquad_e,0.0);
                Array<OneD,NekDouble> yedge(nquad_e,0.0);
                Array<OneD,NekDouble> zedge(nquad_e,0.0);

                EdgeExp->GetCoords(xedge,yedge,zedge);
                //locExpVectorwx[eid]->GetEdgePhysVals(j,xel,xedge);
                //locExpVectorwx[eid]->GetEdgePhysVals(j,xel,yedge);
                //locExpVectorwx[eid]->GetEdgePhysVals(edgeId,U,Uedge);

				Array<OneD,NekDouble> wxedge(nquad_e,0.0);
				Array<OneD,NekDouble> wyedge(nquad_e,0.0);

				locExpVector[eid]->GetEdgePhysVals(j,EdgeExp,wxel,wxedge);
				locExpVector[eid]->GetEdgePhysVals(j,EdgeExp,wyel,wyedge);

				for(int i=0;i<nquad_e;i++)
				{
					xedge[i]= xedge[i] + (timestep)*wxedge[i];
					yedge[i]= yedge[i] + (timestep)*wyedge[i];
				}

				edge = m_graph->GetEdge(edgeId);

				vertID1 = edge->GetVid(0);
				//vertID2 = edge->GetVid(1);
				//TODO: move every vertex just once
				vertex1 = m_graph->GetVertex(vertID1);
				//vertex2 = m_graph->GetVertex(vertID1);
				//cout << "vertID: " << vertID1 << endl;
				vertex1->UpdatePosition(xedge[0],yedge[0],zedge[0]);
				//vertex2->UpdatePosition(xedge[nquad_e-1],yedge[nquad_e-1],zedge[nquad_e-1]);

				/*for(int k=0;k<nquad_e;k++)
							{
								cout << "(x,y)=(" << xedge[k] << "," << yedge[k] << ")" << "\t";
								cout << "(wx,wy)=(" << wxedge[k] << "," << wyedge[k] << ")" << "\t";
							}
							cout << endl; */
            }
        }
    }

    /// Move GLL nodes along the free surface
    void MovingMesh::DeformMesh2(const NekDouble timestep)
     {
  	// number of elements
  	int  nlevels = m_meshvelwx.num_elements();
  	NekDouble x,y,z;
  	int eid, id, diff,el,j,nqel;
  	int vertId, edgeId;
  	StdRegions::StdExpansionSharedPtr ElExp;
  	int globaledgeID;
  	   SpatialDomains::SegGeomSharedPtr edge;
  	   StdRegions::StdExpansion1DSharedPtr EdgeExp;
  	Array<OneD, Array<OneD, NekDouble> > normals, locnormals,weaknormals;
  	Array<OneD, const NekDouble> wxel,wyel;
  	int vertID1,vertID2;
  	 SpatialDomains::VertexComponentSharedPtr vertex1,vertex2;
  	 SpatialDomains::Geometry1DSharedPtr ElmtSegGeom;

  	 // local expansion vector
  	 const StdRegions::StdExpansionVector &locExpVector = *(m_ExpField->GetExp());

  	 int nel = locExpVector.size();

  	   // Loop over elemente and collect forward expansion
  	int nexp = m_fields[m_meshvelocity[0]]->GetExpSize();
  	int nquad_e,n,e,offset,phys_offset;
  	Array<OneD,NekDouble> e_tmp;
  	map<int,int> VertexDone;

  	// Loop over all elements
  	   for(el = 0; el < nel; ++el)
  	   {
  		// Get element id
  		   eid = m_fields[m_meshvelocity[0]]->GetOffset_Elmt_Id(el);
  		   //int nqel =

  		   //ElExp = locExpVector[eid]; element standard expansion
  		   for(j = 0; j < locExpVector[eid]->GetNedges(); ++j)
  		   {
  			   //vertId = (locExpVector[eid]->GetGeom2D())->GetVid(j);
  			   // global ID
  			   edgeId = (locExpVector[eid]->GetGeom2D())->GetEid(j);
  			   edge = m_graph->GetEdge(edgeId);
  			   vertID1 = edge->GetVid(0);
  			   // TODO: introduce check if vertex ID has already been treated or not!!
  			   if(VertexDone.count(vertID1)==0)
  			   {
  				   VertexDone[vertID1] = el;
  				   nquad_e = locExpVector[eid]->GetEdgeNumPoints(j);
  				   //nqel  = locExpVector[eid]->GetTotPoints();

  				   wxel = m_meshvelwx[nlevels-1] + m_fields[m_meshvelocity[0]]->GetPhys_Offset(eid);
  				   wyel = m_meshvelwy[nlevels-1] + m_fields[m_meshvelocity[1]]->GetPhys_Offset(eid);


  				   vertex1 = m_graph->GetVertex(vertID1);
  				   vertex1->GetCoords(x, y, z);

  				   Array<OneD,NekDouble> wxedge(nquad_e,0.0);
  				   Array<OneD,NekDouble> wyedge(nquad_e,0.0);

  				   locExpVector[eid]->GetEdgePhysVals(j,wxel,wxedge);
  				   locExpVector[eid]->GetEdgePhysVals(j,wyel,wyedge);

  				   x = x + (timestep)*wxedge[0];
  				   y = y + (timestep)*wyedge[0];

  				   vertex1->UpdatePosition(x,y,z);
  			   }

  		   }
  	   }
     }

      void MovingMesh::DeformMeshBoundary(const NekDouble aii_Dt)
      {
      	int  nlevels = m_meshvelwx.num_elements();

      	Array<OneD, const SpatialDomains::BoundaryConditionShPtr > BndConds;
      	        Array<OneD, MultiRegions::ExpListSharedPtr>  BndExp, BndExpy;

      	        // Grab boundary regions and expansions of velocity u
      	        //BndConds   = m_fields[0]->GetBndConditions();
      	        //BndExp     = m_fields[0]->GetBndCondExpansions();
      	        BndConds   = m_fields[m_meshvelocity[0]]->GetBndConditions();
      	        BndExp     = m_fields[m_meshvelocity[0]]->GetBndCondExpansions();
      	        BndExpy     = m_fields[m_meshvelocity[1]]->GetBndCondExpansions();

      	        StdRegions::StdExpansion1DSharedPtr EdgeExp;
      	        SpatialDomains::SegGeomSharedPtr edge;
      	        StdRegions::StdExpansionSharedPtr ElExp;
      	        Array<OneD, Array<OneD, NekDouble> > normals, locnormals,weaknormals;
      			int vertID1,vertID2;
      			 SpatialDomains::VertexComponentSharedPtr vertex1,vertex2;
      			 SpatialDomains::MeshGraph2DSharedPtr mesh2D;

      	        Array<OneD, int> ElmtID,EdgeID;
      	        int cnt,n,el,edgeID,elmID;
      	        Array<OneD, const NekDouble> Uphyselement, Vphyselement, wxel,wyel;

      	        //Fills ElmtID and EdgeID with global ids of elements id number and edges id number of the boundary expansion
      	        m_fields[m_meshvelocity[0]]->GetBoundaryToElmtMap(ElmtID,EdgeID);

      	        for(cnt = n = 0; n < m_fields[m_meshvelocity[0]]->GetBndConditions().num_elements(); ++n)
      	              {
      	                  string type = m_fields[m_meshvelocity[0]]->GetBndConditions()[n]->GetUserDefined().GetEquation();
      	          //for the ones marked as Drag calculate drag along boundary edge
      	                  if(type == "FreeSurface")
      	                  {
      	                 	 //loop over all elements along the boundary region
      	                 	bool first=true;
      	     				for(el = 0; el < BndExp[n]->GetExpSize(); ++el,cnt++)
      	     				{
      	     					Uphyselement = (m_fields[0]->GetPhys())+ m_fields[0]->GetPhys_Offset(ElmtID[cnt]);
      	     					Vphyselement = (m_fields[1]->GetPhys())+ m_fields[1]->GetPhys_Offset(ElmtID[cnt]);
      	     					//wxel = (m_fields[m_meshvelocity[0]]->GetPhys())+ m_fields[m_meshvelocity[0]]->GetPhys_Offset(ElmtID[cnt]);
      	     					//wyel = (m_fields[m_meshvelocity[1]]->GetPhys())+ m_fields[m_meshvelocity[1]]->GetPhys_Offset(ElmtID[cnt]);

      	     					wxel = m_meshvelwx[nlevels-1] + m_fields[m_meshvelocity[0]]->GetPhys_Offset(ElmtID[cnt]);
      	     					wyel = m_meshvelwy[nlevels-1] + m_fields[m_meshvelocity[1]]->GetPhys_Offset(ElmtID[cnt]);

      	     					ElExp = m_fields[m_meshvelocity[0]]->GetExp(ElmtID[cnt]);
      	     					EdgeExp =  boost::dynamic_pointer_cast<StdRegions::StdExpansion1D> (BndExp[n]->GetExp(el));

      	     					int nquad_e = EdgeExp->GetNumPoints(0);
      	     		            Array<OneD,NekDouble> xedge(nquad_e,0.0);
      	     					Array<OneD,NekDouble> yedge(nquad_e,0.0);
      	     					Array<OneD,NekDouble> zedge(nquad_e,0.0);
      	     					Array<OneD,NekDouble> wxedge(nquad_e,0.0);
      	     					Array<OneD,NekDouble> wyedge(nquad_e,0.0);

      	     					// EdgeID in within element
      	     					edgeID = EdgeID[cnt];
      	     					ElExp->GetEdgePhysVals(edgeID,EdgeExp,wxel,wxedge);
      	     					ElExp->GetEdgePhysVals(edgeID,EdgeExp,wyel,wyedge);

      	     					EdgeExp->GetCoords(xedge,yedge,zedge);

      	     					 int id1  = BndExp[n]->GetPhys_Offset(el);

      	     						for(int i=0;i<nquad_e;i++)
      	     						{
      	     							xedge[i]= xedge[i] + (aii_Dt)*wxedge[i];
      	     							yedge[i]= yedge[i] + (aii_Dt)*wyedge[i];
      	     						}


      	     						edgeID= m_mesh2D->GetEidFromElmt(ElExp->DetExpansionType(),
      	     						                edgeID, ElmtID[cnt]);

      	     						edge = m_mesh2D->GetEdge(edgeID);


      	     						vertID1=edge->GetVid(0);
      	     						vertID2=edge->GetVid(1);

      	     						// Make changes to m_graph

      	     						vertex1 = m_graph->GetVertex(vertID1);
      	     						vertex2 = m_graph->GetVertex(vertID2);

      	     						vertex1->UpdatePosition(xedge[0],yedge[0],zedge[0]);
      	     						vertex2->UpdatePosition(xedge[nquad_e-1],yedge[nquad_e-1],zedge[nquad_e-1]);

      	     						// mesh->UpdateVertex(vertID1,xedge[0],yedge[0],0);
      	     						// mesh->UpdateVertex(vertID2,xedge[nquad_e-1],yedge[nquad_e-1],0);

      	     						 m_graph->CreateCurvedEdge(edgeID,xedge,yedge,zedge);
      	     				}

      	                  }
      	                  else   // for all other boundary conditi
      	     			 {
      	     				 cnt += BndExp[n]->GetExpSize();
      	     			 }
      	              }

  		        // Reset Values into Pressure BCs
  				/*       for(cnt = n = 0; n < PBndConds.num_elements(); ++n)
  				        {
  				            // High order boundary condition;
  				            if(PBndConds[n]->GetUserDefined().GetEquation() == "H")
  				            {
  				                int nq = PBndExp[n]->GetNcoeffs();
  				                Vmath::Vcopy(nq,&(m_pressureHBCs[nlevels-1])[cnt],1,&(PBndExp[n]->UpdateCoeffs()[0]),1);
  				                cnt += nq;
  				            }
  				        } */
      }

    void MovingMesh::ComputeNewCoordsatFreeSurfaceBoundary(Array<OneD, Array<OneD, NekDouble> > &newcoords, const NekDouble timestep)
    {
       	int  nlevels = m_meshvelwx.num_elements();

       	Array<OneD, const SpatialDomains::BoundaryConditionShPtr > BndConds;
		Array<OneD, MultiRegions::ExpListSharedPtr>  BndExp, BndExpy;

		// Grab boundary regions and expansions of velocity u
		//BndConds   = m_fields[0]->GetBndConditions();
		//BndExp     = m_fields[0]->GetBndCondExpansions();
		BndConds   = m_fields[m_meshvelocity[0]]->GetBndConditions();
		BndExp     = m_fields[m_meshvelocity[0]]->GetBndCondExpansions();
		BndExpy     = m_fields[m_meshvelocity[1]]->GetBndCondExpansions();

		StdRegions::StdExpansion1DSharedPtr EdgeExp;
		SpatialDomains::SegGeomSharedPtr edge;
		StdRegions::StdExpansionSharedPtr ElExp;
		Array<OneD, Array<OneD, NekDouble> > normals, locnormals,weaknormals;
		int vertID1,vertID2;
		SpatialDomains::VertexComponentSharedPtr vertex1,vertex2;
		SpatialDomains::MeshGraph2DSharedPtr mesh2D;

		Array<OneD, int> ElmtID,EdgeID;
		int cnt,n,el,edgeID,elmID;
		Array<OneD, const NekDouble> Uphyselement, Vphyselement, wxel,wyel;

		//Fills ElmtID and EdgeID with global ids of elements id number and edges id number of the boundary expansion
		m_fields[m_meshvelocity[0]]->GetBoundaryToElmtMap(ElmtID,EdgeID);

		for(cnt = n = 0; n < m_fields[m_meshvelocity[0]]->GetBndConditions().num_elements(); ++n)
		  {
			  string type = m_fields[m_meshvelocity[0]]->GetBndConditions()[n]->GetUserDefined().GetEquation();
	  //for the ones marked as Drag calculate drag along boundary edge
			  if(type == "FreeSurface")
			  {
				 //loop over all elements along the boundary region
				bool first=true;
				for(el = 0; el < BndExp[n]->GetExpSize(); ++el,cnt++)
				{
					Uphyselement = (m_fields[0]->GetPhys())+ m_fields[0]->GetPhys_Offset(ElmtID[cnt]);
					Vphyselement = (m_fields[1]->GetPhys())+ m_fields[1]->GetPhys_Offset(ElmtID[cnt]);
					//wxel = (m_fields[m_meshvelocity[0]]->GetPhys())+ m_fields[m_meshvelocity[0]]->GetPhys_Offset(ElmtID[cnt]);
					//wyel = (m_fields[m_meshvelocity[1]]->GetPhys())+ m_fields[m_meshvelocity[1]]->GetPhys_Offset(ElmtID[cnt]);

					wxel = m_meshvelwx[nlevels-1] + m_fields[m_meshvelocity[0]]->GetPhys_Offset(ElmtID[cnt]);
					wyel = m_meshvelwy[nlevels-1] + m_fields[m_meshvelocity[1]]->GetPhys_Offset(ElmtID[cnt]);

					ElExp = m_fields[m_meshvelocity[0]]->GetExp(ElmtID[cnt]);
					EdgeExp =  boost::dynamic_pointer_cast<StdRegions::StdExpansion1D> (BndExp[n]->GetExp(el));

					int nquad_e = EdgeExp->GetNumPoints(0);
					Array<OneD,NekDouble> xedge(nquad_e,0.0);
					Array<OneD,NekDouble> yedge(nquad_e,0.0);
					Array<OneD,NekDouble> zedge(nquad_e,0.0);
					Array<OneD,NekDouble> wxedge(nquad_e,0.0);
					Array<OneD,NekDouble> wyedge(nquad_e,0.0);

					// EdgeID in within element
					edgeID = EdgeID[cnt];
					ElExp->GetEdgePhysVals(edgeID,EdgeExp,wxel,wxedge);
					ElExp->GetEdgePhysVals(edgeID,EdgeExp,wyel,wyedge);

					EdgeExp->GetCoords(xedge,yedge,zedge);

					int id1  = BndExp[n]->GetPhys_Offset(el);

					for(int i=0;i<nquad_e;i++)
					{
						newcoords[0][id1+i]= xedge[i] + (timestep)*wxedge[i];
						newcoords[1][id1+i]= yedge[i] + (timestep)*wyedge[i];
					}

					/*cout << "New Coords= " << endl;
					for(int i=0;i<nquad_e;i++)
					{
						cout << "(x,y)=(" << newcoords[0][id1+i] << "," << newcoords[1][id1+i] << ")" << endl;
						cout << "(wx,wy)=(" << wxedge[i] << "," << wyedge[i] << ")" << endl;
					} */
				}

			  }
			  else   // for all other boundary conditi
			 {
				 cnt += BndExp[n]->GetExpSize();
			 }
		  }
       }

    void MovingMesh::UpdateMeshFreeSurfaceBoundary(const NekDouble timestep)
	{
    	int cnt,n,i;
    	int el,edgeID,elmID;
		Array<OneD, MultiRegions::ExpListSharedPtr>  BndExp     = m_fields[m_meshvelocity[0]]->GetBndCondExpansions();
		Array<OneD, const SpatialDomains::BoundaryConditionShPtr > BndConds = m_fields[m_meshvelocity[0]]->GetBndConditions();
		Array<OneD, Array<OneD, NekDouble> > newcoords(3);
		LibUtilities::CubicSplineSharedPtr CubicSpline;
		LibUtilities::SplineBoundaryType bcleft,bcright;
		bcleft = LibUtilities::eClampedQuadraticLagrange;
		bcright = LibUtilities::eNotaKnot;
		int bcl=0,bcr=0;

		StdRegions::StdExpansion1DSharedPtr EdgeExp;
		StdRegions::StdExpansionSharedPtr ElExp;

		SpatialDomains::SegGeomSharedPtr edge;
		int vertID1,vertID2;
		SpatialDomains::VertexComponentSharedPtr vertex1,vertex2;
		SpatialDomains::MeshGraph2DSharedPtr mesh2D;

		Array<OneD, int> ElmtID,EdgeID;

		m_fields[m_meshvelocity[0]]->GetBoundaryToElmtMap(ElmtID,EdgeID);

		// For every free surface region compute gll points along new free surface boundary
		for(cnt = n = 0; n < BndConds.num_elements(); ++n)
		 {
			 string type = BndConds[n]->GetUserDefined().GetEquation();
	 //for the ones marked as Drag calculate drag along boundary edge
			 if(type == "FreeSurface")
			 {
				 int nq = BndExp[n]->GetTotPoints();
				 for(i=0;i<newcoords.num_elements();i++)
				 {
					 newcoords[i] = Array<OneD, NekDouble>(nq,0.0);
				 }
				 /* Compute new coordinates of free surface boundary */
				 ComputeNewCoordsatFreeSurfaceBoundary(newcoords,timestep);
				 /* Construct cubic spline through new coordinates of free surface */
				 CubicSpline = MemoryManager<LibUtilities::CubicSpline>::AllocateSharedPtr(newcoords[0],newcoords[1],bcleft,bcright,bcl,bcr);
				 /* Compute GLL Points along each edge */
				 for(el = 0; el < BndExp[n]->GetExpSize(); ++el,cnt++)
				 {
						ElExp = m_fields[m_meshvelocity[0]]->GetExp(ElmtID[cnt]);
						EdgeExp =  boost::dynamic_pointer_cast<StdRegions::StdExpansion1D> (BndExp[n]->GetExp(el));
						edgeID = EdgeID[cnt];

						int nquad_e = EdgeExp->GetNumPoints(0);
						Array<OneD,NekDouble> xedge(nquad_e,0.0);
						Array<OneD,NekDouble> yedge(nquad_e,0.0);
						Array<OneD,NekDouble> xtmp(nquad_e,0.0);
						Array<OneD,NekDouble> ytmp(nquad_e,0.0);
						Array<OneD,NekDouble> zedge(nquad_e,0.0);
						Array<OneD,const NekDouble> e_tmp;

						// EdgeID in within element
						int id1  = BndExp[n]->GetPhys_Offset(el);

						/* Copy values of new coordinates into xedge, yedge */
						for(int i=0;i<nquad_e;i++)
						{
							xedge[i] = newcoords[0][id1+i];
							yedge[i] = newcoords[1][id1+i];
						}

					/*	cout << "Coords before= " << endl;
						for(int i=0;i<nquad_e;i++)
						{
							cout << "(x,y)=(" << xedge[i] << "," << yedge[i] << ")" << endl;
						} */
						/* Compute GLL positions of xedge,yedge and save them in xedge, yedge*/
						/* For Backward edges result must be copied in reverse order to xtmp to compute GLL points */
						if(ElExp->GetEorient(edgeID) == StdRegions::eForwards)
						{
							Vmath::Vcopy(nquad_e,xedge,1,xtmp,1);
							Vmath::Vcopy(nquad_e,yedge,1,ytmp,1);
						}
						else
						{
							Vmath::Vcopy(nquad_e,e_tmp = xedge+(nquad_e-1),-1,xtmp,1);
							Vmath::Vcopy(nquad_e,e_tmp = yedge+(nquad_e-1),-1,ytmp,1);
						}

						CubicSpline->ComputeGLLPointsAlongEdge(xtmp,ytmp,nquad_e);
						if(ElExp->GetEorient(edgeID) == StdRegions::eForwards)
						{
							Vmath::Vcopy(nquad_e,xtmp,1,xedge,1);
							Vmath::Vcopy(nquad_e,ytmp,1,yedge,1);
						}
						else
						{
							Vmath::Vcopy(nquad_e,e_tmp = xtmp+(nquad_e-1),-1,xedge,1);
							Vmath::Vcopy(nquad_e,e_tmp = ytmp+(nquad_e-1),-1,yedge,1);
						}
						/*cout << "Coords after=" << endl;
						for(int i=0;i<nquad_e;i++)
						{
							cout << "(x,y)=(" << xedge[i] << "," << yedge[i] << ")" << endl;
						} */

						/* Move the according edges and vertices to new coordinates */
						edgeID= m_mesh2D->GetEidFromElmt(ElExp->DetExpansionType(),
						    	     						                edgeID, ElmtID[cnt]);

						edge = m_mesh2D->GetEdge(edgeID);

						vertID1=edge->GetVid(0);
						vertID2=edge->GetVid(1);

						// Make changes to m_graph
						 vertex1 = m_graph->GetVertex(vertID1);
						 vertex2 = m_graph->GetVertex(vertID2);

						vertex1->UpdatePosition(xedge[0],yedge[0],zedge[0]);
						vertex2->UpdatePosition(xedge[nquad_e-1],yedge[nquad_e-1],zedge[nquad_e-1]);

						// mesh->UpdateVertex(vertID1,xedge[0],yedge[0],0);
						// mesh->UpdateVertex(vertID2,xedge[nquad_e-1],yedge[nquad_e-1],0);

						 m_graph->CreateCurvedEdge(edgeID,xedge,yedge,zedge);
				 }
			 }
			 else   // for all other boundary conditi
			 {
				 cnt += BndExp[n]->GetExpSize();
			 }
		 }

	}

    void MovingMesh::UpdateFields()
    {
    	// velocities plus pressure
    	int i;

        for(i = 0 ; i < m_fields.num_elements(); i++)
		{
        	MultiRegions::ContField2DSharedPtr ContField;
        	MultiRegions::DisContField2DSharedPtr DisContField;

        	if(ContField = boost::dynamic_pointer_cast<
        			MultiRegions::ContField2D>(m_fields[i]))
  			{
        		ContField->UpdateContField2D(m_session,m_graph,*m_boundaryConditions,
        	    			 m_session->GetVariable(i));
  			}
        	else if(DisContField = boost::dynamic_pointer_cast<
        			MultiRegions::DisContField2D>(m_fields[i]))
  			{
        		DisContField->UpdateDisContField2D(m_session,m_graph,*m_boundaryConditions,
        	    			 m_session->GetVariable(i),true);
  			}
        	else
        	{
        		cout << "Warning Field " << i << "has not been updated" << endl;
        	}
		}

    }


    // Update m_graph and m_boundarycondtions objects
    void MovingMesh::UpdateMesh()
    {
    	m_mesh2D = m_mesh2D->UpdateMesh(m_graph);

	    if(!(m_graph = boost::dynamic_pointer_cast<
							   SpatialDomains::MeshGraph>(m_mesh2D)))
	    {
		  ASSERTL0(false,"Dynamics cast failed");
	    }

	    m_boundaryConditions = m_boundaryConditions->UpdateBoundaryConditions(m_graph,m_boundaryConditions);
		 // Dummy to deal with mesh properties
		/*  m_DisField->UpdateDisContField2D(m_comm,*m_mesh2D,*m_boundaryConditions,
						 m_boundaryConditions->GetVariable(0),
										 m_solnType,
										 true); */
	    m_ExpField->UpdateExpList2D(m_session,m_mesh2D,m_session->GetVariable(m_meshvelocity[0]));

		for(int i = 0; i < 3; ++i)
		{
			Vmath::Vcopy(m_phystot,m_meshcoords[i],1,m_meshcoordsold[i],1);
		}

		// new coordinates
		// TODO: we need a field that saves the new coords
		//m_DisField->GetCoords(m_meshcoords[0],m_meshcoords[1],m_meshcoords[2]);
		m_ExpField->GetCoords(m_meshcoords[0],m_meshcoords[1],m_meshcoords[2]);
    }

    void MovingMesh::ComputeMeshVelocity(const NekDouble timestep)
    {
    	int VelDim = m_spacedim;
    	for(int i=0;i<VelDim;i++)
    	{
        	for(int j=0;j<m_phystot;j++)
        	{
        		m_fields[m_meshvelocity[i]]->UpdatePhys()[j]=(m_meshcoords[i][j]-m_meshcoordsold[i][j])/timestep;
        	}
    	}
    }

    void MovingMesh::ComputeVelocityMinusMeshVelocity(
			Array<OneD, Array<OneD, NekDouble> > &velocity)
    {
    	int VelDim=velocity.num_elements();

    	for(int i=0;i<VelDim;i++)
    	{
    		Vmath::Vsub(m_phystot,velocity[i],1,m_fields[m_meshvelocity[i]]->GetPhys(),1,velocity[i],1);

    	}
    	//Vmath::Vsub(m_phystot,velocity[0],1,m_meshvelwx[1],1,velocity[0],1);
    	//Vmath::Vsub(m_phystot,velocity[1],1,m_meshvelwy[1],1,velocity[1],1);

    	//static int cntwhatever = 1;
    	//WriteArrayAlongFreeSurfaceToFile(m_sessionName,cntwhatever++,velocity);
    }

    void MovingMesh::CreateFreeSurfaceFunction(const NekDouble timestep)
    {
    	 Array<OneD, const SpatialDomains::BoundaryConditionShPtr > BndConds;
		Array<OneD, MultiRegions::ExpListSharedPtr>  BndExp, BndExpy;
		MultiRegions::ExpList1DSharedPtr  FreeSurfaceFct;
		Array<OneD, Array<OneD, NekDouble> > normals, locnormals;

		 BndConds   = m_fields[m_meshvelocity[0]]->GetBndConditions();
		BndExp     = m_fields[m_meshvelocity[0]]->GetBndCondExpansions();
		BndExpy     = m_fields[m_meshvelocity[1]]->GetBndCondExpansions();

		StdRegions::StdExpansionSharedPtr ElExp;
		StdRegions::StdExpansion1DSharedPtr EdgeExp;
		Array<OneD, int> ElmtID,EdgeID;
		 int cnt,n,el,edgeID,elmID;
		 Array<OneD, const NekDouble> Uphyselement, Vphyselement;


		 //Fills ElmtID and EdgeID with global ids of elements id number and edges id number of the boundary expansion
		 m_fields[m_meshvelocity[0]]->GetBoundaryToElmtMap(ElmtID,EdgeID);

    	        /*  for(int i=0;i<EdgeID.num_elements();i++)
    	          {
    	          	cout<< "elementIDs" << ElmtID[i] << endl;
    	          	cout<< "edgeIDs" << EdgeID[i] << endl;
    	          } */
    	          // ensure continuity in normals
    	          map<int, int> elements;
    	          int k=0;
    	          for(cnt = n = 0; n < m_fields[m_meshvelocity[0]]->GetBndConditions().num_elements(); ++n)
    	                 {
    	                     string type = m_fields[m_meshvelocity[0]]->GetBndConditions()[n]->GetUserDefined().GetEquation();
    	             //for the ones marked as Drag calculate drag along boundary edge
    	                     if(type == "FreeSurface")
    	                     {
    	                  	  // BndExp[n]->SetUpTangents();
    	                  	   //int size=BndExp[n]->GetExpSize();
    	                  	   int phystot = BndExp[n]->GetTotPoints();
    	                  	//  MultiRegions::ExpList1DSharedPtr BndExp1D;
    	                  	/*   if(!(FreeSurfaceFct =  boost::dynamic_pointer_cast<MultiRegions::ExpList1D>(BndExp[n])))
    	  						{
    	  							ASSERTL0(false,"Dynamics cast failed");
    	  						}
    	                  	 MultiRegions::ContField1DSharedPtr dfdxfct =
    	                  	 				MemoryManager<MultiRegions::ContField1D>
    	                  	 				::AllocateSharedPtr(m_comm,*FreeSurfaceFct,m_solnType); */


    	                  	   int nq= BndExp[n]->GetTotPoints();
							   Array<OneD,NekDouble> x(nq,0.0);
								Array<OneD,NekDouble> y(nq,0.0);
								Array<OneD,NekDouble> z(nq,0.0);
								Array<OneD,NekDouble> dfdx(nq,0.0);
								BndExp[n]->GetCoords(x,y,z);

								// f(x,t) = y
								FreeSurfaceFct->SetPhys(y);
								FreeSurfaceFct->PhysDeriv(0,FreeSurfaceFct->GetPhys(),dfdx);
							//	dfdxfct->SetPhys(dfdx);

    	                  /*	   cout << "Setting up free surface function:" << endl;
    	  					for(int i=0;i<nq;i++)
							{
								cout << i << ": (" << x[i] << "," << y[i] << ")=";
								//cout << "(n1,n2)=(" << normals[0][i] <<","<< normals[1][i] << ")"<< endl;
								cout << "f(x,t)=" << y[i] << ", df/dx=" << dfdx[i] << endl;
							} */

							//	dfdxfct->FwdTrans(dfdx,dfdxfct->UpdateCoeffs());
							//	dfdxfct->BwdTrans(dfdxfct->GetCoeffs(),dfdx);
						  /* 	cout << "Smoothing of dfdx (C0):" << endl;
							for(int i=0;i<nq;i++)
							{
								cout << i << ": (" << x[i] << "," << y[i] << ")=";
								//cout << "(n1,n2)=(" << normals[0][i] <<","<< normals[1][i] << ")"<< endl;
								cout << "f(x,t)=" << y[i] << ", df/dx=" << dfdx[i] << endl;
							} */


    	  			     	bool first=true;
    	  							for(el = 0; el < BndExp[n]->GetExpSize(); ++el,cnt++)
    	  							{
    	  								Uphyselement = (m_fields[0]->GetPhys())+ m_fields[0]->GetPhys_Offset(ElmtID[cnt]);
    	  								Vphyselement = (m_fields[1]->GetPhys())+ m_fields[1]->GetPhys_Offset(ElmtID[cnt]);

    	  								ElExp = m_fields[m_meshvelocity[0]]->GetExp(ElmtID[cnt]);
    	  							/*	int nq = locExp->GetTotPoints();
    	  								Array<OneD,NekDouble> xel(nq,0.0);
    	  								Array<OneD,NekDouble> yel(nq,0.0);
    	  								Array<OneD,NekDouble> zel(nq,0.0);
    	  								ElExp->GetCoords(xel,yel,zel);
    	  								for(int i=0;i<nq;i++)
    	  								{
    	  									cout << i << ": (" << xel[i] << "," << yel[i] << ")=";
    	  									cout << "(u,v)=(" << Uphyselement[i] <<","<< Vphyselement[i] << ")"<< endl;
    	  								} */

    	  								EdgeExp =  boost::dynamic_pointer_cast<StdRegions::StdExpansion1D> (BndExp[n]->GetExp(el));
    	  								int nquad_e = EdgeExp->GetNumPoints(0);
    	  					            Array<OneD,NekDouble> xedge(nquad_e,0.0);
    	  								Array<OneD,NekDouble> yedge(nquad_e,0.0);
    	  								Array<OneD,NekDouble> zedge(nquad_e,0.0);
    	  								Array<OneD,NekDouble> Uphysedge(nquad_e,0.0);
    	  								Array<OneD,NekDouble> Vphysedge(nquad_e,0.0);
    	  								Array<OneD,NekDouble> wxedge(nquad_e,0.0);
    	  								Array<OneD,NekDouble> wyedge(nquad_e,0.0);

    	  								// EdgeID in within element
    	  								edgeID = EdgeID[cnt];
    	  								ElExp->GetEdgePhysVals(edgeID,EdgeExp,Uphyselement,Uphysedge);
    	  								ElExp->GetEdgePhysVals(edgeID,EdgeExp,Vphyselement,Vphysedge);

    	  								EdgeExp->GetCoords(xedge,yedge,zedge);
    	  								// Copy result in boundary expansion
										int id1  = BndExp[n]->GetPhys_Offset(el);
										int id  = FreeSurfaceFct->GetPhys_Offset(el);

    	  								// compute wx/wy on the edge
    	  								for(int i=0;i<nquad_e;i++)
    	  								{
    	  									//wyedge[i]= (1.0/(1.0+(dfdx[id+i]*dfdx[id+i])))*(Vphysedge[i]-Uphysedge[i]*dfdx[id+i]);
    	  									//wxedge[i]= -1.0*dfdx[id+i]*wyedge[i];
    	  									// wx = 0
    	  									wyedge[i] = Vphysedge[i]-Uphysedge[i]*dfdx[id+i];
    	  									wxedge[i]=0.0;
    	  								}

    	  								Vmath::Vcopy(nquad_e,&wxedge[0], 1,&(BndExp[n]->UpdatePhys())[id1],1);
    	  								Vmath::Vcopy(nquad_e,&wyedge[0], 1,&(BndExpy[n]->UpdatePhys())[id1],1);
    	  								// MOVE edge, i.e. calculate X^{n+1} at free surface boundary
    	  								// TODO: Implement time integration instead

    	  							/*	if(first)
    	  								{ */

    	  							/*		for(int i=0;i<nquad_e;i++)
    	  									{
    	  										//if(xedge[i]==11.0 && yedge[i]<=0.0)
    	  															    		{
    	  										cout << i << ": (" << xedge[i] << "," << yedge[i] << ")"  << "\t";
    	  										cout << "Edge: (u,v)=(" << Uphysedge[i] <<","<< Vphysedge[i] << ")"<< "\t";
    	  										cout << "dfdx=(" << dfdx[id+i] << "\t";
    	  										cout << "(wx,wy)=(" << wxedge[i] <<","<< wyedge[i] << ")"<< endl;
    	  															    		}
    	  									} */
    	  							//	}


    	  							}
    	  			            	BndExp[n]->FwdTrans_BndConstrained(BndExp[n]->GetPhys(),BndExp[n]->UpdateCoeffs());
    	  			            	BndExpy[n]->FwdTrans_BndConstrained(BndExpy[n]->GetPhys(),BndExpy[n]->UpdateCoeffs());
    	                     }
    	                     else  // for all other boundary conditi
    	  					 {
    	  						 cnt += BndExp[n]->GetExpSize();
    	  					 }
    	                 }
    }

    void MovingMesh::CreateFreeSurfaceSpline(const NekDouble timestep,

    										const NekDouble time)
    {
     	Array<OneD, const SpatialDomains::BoundaryConditionShPtr > BndConds;
 		Array<OneD, MultiRegions::ExpListSharedPtr>  BndExp, BndExpy;
 		MultiRegions::ExpList1DSharedPtr  FreeSurfaceFct;
 		Array<OneD, Array<OneD, NekDouble> > normals, locnormals;
 		static int outputcnt = 1;
 		NekDouble Uplate;
 		if(time < 10)
 		{
 			Uplate =1.0;
 		}
 		else
 		{
 			Uplate = 0.0;
 		}

 		BndConds   = m_fields[m_meshvelocity[0]]->GetBndConditions();
 		BndExp     = m_fields[m_meshvelocity[0]]->GetBndCondExpansions();
 		BndExpy     = m_fields[m_meshvelocity[1]]->GetBndCondExpansions();

 		StdRegions::StdExpansionSharedPtr ElExp;
 		StdRegions::StdExpansion1DSharedPtr EdgeExp;
 		Array<OneD, int> ElmtID,EdgeID;
 		int cnt, n,el,edgeID,elmID;
 		Array<OneD, const NekDouble> Uphyselement, Vphyselement;

 		/*stringstream convert;
 		convert << "FSBC" << outputcnt << ".txt";
 		string filename = convert.str();
 		ofstream outfile(filename.c_str());
 		outputcnt++; */

 		//Fills ElmtID and EdgeID with global ids of elements id number and edges id number of the boundary expansion
 		m_fields[m_meshvelocity[0]]->GetBoundaryToElmtMap(ElmtID,EdgeID);

 		NekDouble xmax = 0.0;
 		int j,k,l;
 		l =0;
 		// Determine xmax because its outflow value that needs to be set to value from before
		for(cnt = n = 0; n < m_fields[m_meshvelocity[0]]->GetBndConditions().num_elements(); ++n)
		{
			 string type = m_fields[m_meshvelocity[0]]->GetBndConditions()[n]->GetUserDefined().GetEquation();
			 if(type == "FreeSurface")
			 {
				 int nquadbc = BndExp[n]->GetTotPoints();
				 Array<OneD,NekDouble> xbc(nquadbc,0.0);
				 Array<OneD,NekDouble> ybc(nquadbc,0.0);
				 Array<OneD,NekDouble> zbc(nquadbc,0.0);

				 BndExp[n]->GetCoords(xbc,ybc,zbc);

				 xmax = Vmath::Vmax(nquadbc,xbc,1);
			 }
			 else // for all other boundary conditi
			 {
				 cnt += BndExp[n]->GetExpSize();
			 }
		}

		int ispline=0;
	    for(cnt = n = 0; n < m_fields[m_meshvelocity[0]]->GetBndConditions().num_elements(); ++n)
	    {
			 string type = m_fields[m_meshvelocity[0]]->GetBndConditions()[n]->GetUserDefined().GetEquation();
	 //for the ones marked as Drag calculate drag along boundary edge
			 if(type == "FreeSurface")
			 {
				m_freesurfacesplines[ispline]->UpdateSpline(BndExp[n]);
				// m_freesurfacesplines[ispline]->UpdateSplineCornerPoints(BndExp[n]);
				//m_freesurfacesplines[ispline]->WriteMatlabFiles(m_sessionName, outputcnt++);

				for(el = 0; el < BndExp[n]->GetExpSize(); ++el,cnt++)
				{
					Uphyselement = (m_fields[0]->GetPhys())+ m_fields[0]->GetPhys_Offset(ElmtID[cnt]);
					Vphyselement = (m_fields[1]->GetPhys())+ m_fields[1]->GetPhys_Offset(ElmtID[cnt]);

					ElExp = m_fields[m_meshvelocity[0]]->GetExp(ElmtID[cnt]);

					EdgeExp =  boost::dynamic_pointer_cast<StdRegions::StdExpansion1D> (BndExp[n]->GetExp(el));
					int nquad_e = EdgeExp->GetNumPoints(0);
					Array<OneD,NekDouble> xedge(nquad_e,0.0);
					Array<OneD,NekDouble> yedge(nquad_e,0.0);
					Array<OneD,NekDouble> zedge(nquad_e,0.0);
					Array<OneD,NekDouble> Uphysedge(nquad_e,0.0);
					Array<OneD,NekDouble> Vphysedge(nquad_e,0.0);
					Array<OneD,NekDouble> wxedge(nquad_e,0.0);
					Array<OneD,NekDouble> wyedge(nquad_e,0.0);

					EdgeExp->GetCoords(xedge,yedge,zedge);
					edgeID = EdgeID[cnt];
					ElExp->GetEdgePhysVals(edgeID,EdgeExp,Uphyselement,Uphysedge);
					ElExp->GetEdgePhysVals(edgeID,EdgeExp,Vphyselement,Vphysedge);

					EdgeExp->GetCoords(xedge,yedge,zedge);
					NekDouble dfdx,nx,ny;

					// compute wx/wy on the edge
					for(int i=0;i<nquad_e;i++)
					{
						// wx = 0
						m_freesurfacesplines[ispline]->GetSplineNormal(xedge[i],nx,ny);
						dfdx = -nx/ny;
						wyedge[i] = Vphysedge[i]+Uphysedge[i]*nx/ny;
						wxedge[i]=0.0;
						//wyedge[i]= (1.0/(1.0+(dfdx*dfdx)))*(Vphysedge[i]-Uphysedge[i]*dfdx);
						//wxedge[i]= -1.0*dfdx*wyedge[i];
						//wyedge[i] = Vphysedge[i];
						//wxedge[i]= Uphysedge[i];

					}

				/*	for(int i=0;i<nquad_e;i++)
					{

						// Compute wy wx for outflow corner!!!
						if(xedge[i] == xmax)
						{
							StdRegions::EdgeOrientation edgedir = ElExp->GetEorient(edgeID);
							//if(edgedir == StdRegions::eForwards)
							{
								wxedge[i]= Uplate;
								wyedge[i]= 0.0;
								//cout << "Set wxedge to one"  << endl;
								//wyedge[i]= wyedge[i+1];
								//wxedge[i]= wxedge[i+1];
							}
							//else
							{
								//wyedge[i]= wyedge[i-1];
								//wxedge[i]= wxedge[i-1];
							}
						}
					} */




				/*	for(int i=0;i<nquad_e;i++)
					{
						outfile << xedge[i] << "\t"<< yedge[i] <<  "\t";
						outfile <<  Uphysedge[i] <<"\t"<< Vphysedge[i] <<  "\t";
						outfile << wxedge[i] << "\t" << wyedge[i] << endl;
					} */

				/*	for(int i=0;i<nquad_e;i++)
					{
						cout<< xedge[i] << "\t"<< yedge[i] <<  "\t";
						cout <<  Uphysedge[i] <<"\t"<< Vphysedge[i] <<  "\t";
						cout << wxedge[i] << "\t" << wyedge[i] << endl;
					} */

					// Copy result in boundary expansion
					int id1  = BndExp[n]->GetPhys_Offset(el);
					Vmath::Vcopy(nquad_e,&wxedge[0], 1,&(BndExp[n]->UpdatePhys())[id1],1);
					Vmath::Vcopy(nquad_e,&wyedge[0], 1,&(BndExpy[n]->UpdatePhys())[id1],1);
				}
				BndExp[n]->FwdTrans_BndConstrained(BndExp[n]->GetPhys(),BndExp[n]->UpdateCoeffs());
				BndExpy[n]->FwdTrans_BndConstrained(BndExpy[n]->GetPhys(),BndExpy[n]->UpdateCoeffs());
				ispline++;
			 }
			 else // for all other boundary conditi
			 {
				 cnt += BndExp[n]->GetExpSize();
			 }
		 }
     }

    void MovingMesh::SmoothNormalsAlongFreeSurface(const NekDouble timestep)
    {
		Array<OneD, const SpatialDomains::BoundaryConditionShPtr > BndConds;
		Array<OneD, MultiRegions::ExpListSharedPtr>  BndExp, BndExpy;
		MultiRegions::ExpList1DSharedPtr  FreeSurfaceFct;
		Array<OneD, Array<OneD, NekDouble> > normals, locnormals;

		BndConds   = m_fields[m_meshvelocity[0]]->GetBndConditions();
		BndExp     = m_fields[m_meshvelocity[0]]->GetBndCondExpansions();
		BndExpy     = m_fields[m_meshvelocity[1]]->GetBndCondExpansions();

		StdRegions::StdExpansionSharedPtr ElExp;
		StdRegions::StdExpansion1DSharedPtr EdgeExp;
		Array<OneD, int> ElmtID,EdgeID;
		int cnt,n,el,edgeID,elmID;
		Array<OneD, const NekDouble> Uphyselement, Vphyselement;

		//Fills ElmtID and EdgeID with global ids of elements id number and edges id number of the boundary expansion
		m_fields[m_meshvelocity[0]]->GetBoundaryToElmtMap(ElmtID,EdgeID);

		// Loop over Boundary Regions
		for(cnt = n = 0; n < m_fields[m_meshvelocity[0]]->GetBndConditions().num_elements(); ++n)
		{
			string type = m_fields[m_meshvelocity[0]]->GetBndConditions()[n]->GetUserDefined().GetEquation();
			//for the ones marked as Drag calculate drag along boundary edge
			if(type == "FreeSurface")
			{
				 //  if(!(m_freesurfacenormals[0] =  boost::dynamic_pointer_cast<MultiRegions::ContField1D>(BndExp[n])))
					{
						ASSERTL0(false,"Dynamics cast failed");
					}
			}
		}
     }


    void MovingMesh::AddFreeSurfaceTension(Array<OneD, Array<OneD, NekDouble> > &Weakoutarray, const NekDouble tension)
      {
			Array<OneD, const SpatialDomains::BoundaryConditionShPtr > BndConds;
			Array<OneD, MultiRegions::ExpListSharedPtr>  BndExp, BndExpv;
			MultiRegions::ExpList1DSharedPtr  FreeSurfaceFct;
			Array<OneD, Array<OneD, NekDouble> > normals, locnormals;
			Array<OneD, NekDouble> e_outarrayu, e_outarrayv;
			Array<OneD, NekDouble> curv;
			bool NegateNormals;

			BndConds   = m_fields[0]->GetBndConditions();
			BndExp     = m_fields[0]->GetBndCondExpansions();
			BndExpv    = m_fields[1]->GetBndCondExpansions();

			StdRegions::StdExpansionSharedPtr ElExp;
			StdRegions::StdExpansion1DSharedPtr EdgeExp,EdgeExpv;
			Array<OneD, int> ElmtID,EdgeID;
			int cnt,n,el,edgeID,elmtid;
			Array<OneD, const NekDouble> U,V;
			Array<OneD, NekDouble> Uvals,Vvals;


			//Fills ElmtID and EdgeID with global ids of elements id number and edges id number of the boundary expansion
			m_fields[0]->GetBoundaryToElmtMap(ElmtID,EdgeID);

			for(cnt = n = 0; n < m_fields[0]->GetBndConditions().num_elements(); ++n)
			{
				 string type = m_fields[0]->GetBndConditions()[n]->GetUserDefined().GetEquation();
		 //for the ones marked as Drag calculate drag along boundary edge
				 if(type == "FreeSurface")
				 {
				   int phystot = BndExp[n]->GetTotPoints();
				   if(!(FreeSurfaceFct =  boost::dynamic_pointer_cast<MultiRegions::ExpList1D>(BndExp[n])))
					{
						ASSERTL0(false,"Dynamics cast failed");
					}

				   int nq= BndExp[n]->GetTotPoints();
				   Array<OneD,NekDouble> x(nq,0.0);
				   Array<OneD,NekDouble> y(nq,0.0);
				   Array<OneD,NekDouble> temp(nq,0.0);
				   Array<OneD,NekDouble> dfdx(nq,0.0);
				   Array<OneD,NekDouble> dfdx2(nq,0.0);
				   Array<OneD,NekDouble> curvature(nq,0.0);
				   BndExp[n]->GetCoords(x,y,temp);

					// f(x,t) = y
					FreeSurfaceFct->SetPhys(y);
					FreeSurfaceFct->PhysDeriv(0,FreeSurfaceFct->GetPhys(),dfdx);
					FreeSurfaceFct->PhysDeriv(0,dfdx,dfdx2);

					Vmath::Vabs(nq,dfdx2,1,dfdx2,1);
					Vmath::Vmul(nq,dfdx,1,dfdx,1,temp,1);
					Vmath::Sadd(nq,1.0,temp,1,temp,1);
					for(int i=0;i<nq;i++)
					{
						temp[i] = pow(temp[i],1.5);
					}
					Vmath::Vdiv(nq,dfdx2,1,temp,1,curvature,1);

			  	   /*cout << "Calc curvature:" << endl;
				for(int i=0;i<nq;i++)
				{
					cout << i << ": (" << x[i] << "," << y[i] << ")=";
					//cout << "(n1,n2)=(" << normals[0][i] <<","<< normals[1][i] << ")"<< endl;
					cout << "cur=" << curvature[i] << endl;
				} */

				//NekDouble tension = m_reynolds/m_weber;
					Vmath::Smul(nq,tension,curvature,1,curvature,1);

						for(el = 0; el < BndExp[n]->GetExpSize(); ++el,cnt++)
						{
							EdgeExp =  boost::dynamic_pointer_cast<StdRegions::StdExpansion1D> (BndExp[n]->GetExp(el));
							EdgeExpv =  boost::dynamic_pointer_cast<StdRegions::StdExpansion1D> (BndExpv[n]->GetExp(el));
							int nquad_e = EdgeExp->GetNumPoints(0);
							int id1  = BndExp[n]->GetPhys_Offset(el);
							edgeID = EdgeID[cnt];
							ElExp = m_fields[0]->GetExp(ElmtID[cnt]);
							int nqel = ElExp->GetTotPoints();

							Array<OneD,NekDouble> dUdx(nqel,0.0);
							Array<OneD,NekDouble> dUdy(nqel,0.0);
							Array<OneD,NekDouble> dVdx(nqel,0.0);
							Array<OneD,NekDouble> dVdy(nqel,0.0);


						    U = (m_fields[0]->GetPhys())+ m_fields[0]->GetPhys_Offset(ElmtID[cnt]);
						    V = (m_fields[1]->GetPhys())+ m_fields[1]->GetPhys_Offset(ElmtID[cnt]);


						    // Calculating deformation tensor
						    ElExp->PhysDeriv(0,U,dUdx);
						    ElExp->PhysDeriv(1,U,dUdy);
						    ElExp->PhysDeriv(0,V,dVdx);
						    ElExp->PhysDeriv(1,V,dVdy);

							e_outarrayu = Weakoutarray[0] + m_fields[0]->GetCoeff_Offset(ElmtID[cnt]);
							e_outarrayv = Weakoutarray[1] + m_fields[1]->GetCoeff_Offset(ElmtID[cnt]);
							// TODO:test that right curvature values are taken for comparison
							Array<OneD,NekDouble> curvfromnormal(nquad_e,0.0);
							Array<OneD,NekDouble> curvfromnormaltimesnormalx(nquad_e,0.0);
							Array<OneD,NekDouble> curvfromnormaltimesnormaly(nquad_e,0.0);
							Array<OneD,NekDouble> dnxdx(nquad_e,0.0);
							Array<OneD,NekDouble> dnydy(nquad_e,0.0);
							Array<OneD,NekDouble> dUdxedge(nquad_e,0.0);
							Array<OneD,NekDouble> dUdyedge(nquad_e,0.0);
							Array<OneD,NekDouble> dVdxedge(nquad_e,0.0);
							Array<OneD,NekDouble> dVdyedge(nquad_e,0.0);


							//normals	=  EdgeExp->GetMetricInfo()->GetNormal();
							normals = ElExp->GetEdgeNormal(edgeID);

							if(ElExp->GetEorient(edgeID) == StdRegions::eBackwards)
							{
								 Vmath::Neg(nquad_e,normals[0],1);
								 Vmath::Neg(nquad_e,normals[1],1);
							}

							ElExp->GetEdgePhysVals(edgeID,EdgeExp,dUdx,dUdxedge);
							ElExp->GetEdgePhysVals(edgeID,EdgeExp,dUdy,dUdyedge);
							ElExp->GetEdgePhysVals(edgeID,EdgeExp,dVdx,dVdxedge);
							ElExp->GetEdgePhysVals(edgeID,EdgeExp,dVdy,dVdyedge);

							EdgeExp->PhysDeriv(0,normals[0],dnxdx);
							EdgeExp->PhysDeriv(1,normals[1],dnydy);

							Vmath::Vadd(nquad_e,dnxdx,1,dnydy,1,curvfromnormal,1);
							Vmath::Neg(nquad_e,curvfromnormal,1);
							Vmath::Smul(nquad_e,tension,curvfromnormal,1,curvfromnormal,1);


							//Array<OneD,NekDouble> curvfromnormaledge(nquad_e,0.0);
							//ElExp->GetEdgePhysVals(edgeID,EdgeExp,curvfromnormal,curvfromnormaledge);
							curv = curvature + BndExp[n]->GetPhys_Offset(el);

							/*for(int j=0;j<nquad_e;j++)
							{
								cout << "curv=" << curv[j] << ", curvfromnormal=" << curvfromnormal[j] << endl;
							} */

						/*	Vmath::Vmul(nquad_e,normals[0],1,curvfromnormal,1,curvfromnormaltimesnormalx,1);
							Vmath::Vmul(nquad_e,normals[1],1,curvfromnormal,1,curvfromnormaltimesnormaly,1);



							Vmath::Vsub(nquad_e,curvfromnormal,1,dUdxedge,1,dUdxedge,1);
							//Vmath::Vsub(nquad_e,curvfromnormal,1,dUdyedge,1,dUdyedge,1);
							//Vmath::Vsub(nquad_e,curvfromnormal,1,dVdxedge,1,dVdxedge,1);
							Vmath::Neg(nquad_e,dUdyedge,1);
							Vmath::Neg(nquad_e,dVdxedge,1);
							Vmath::Vsub(nquad_e,curvfromnormal,1,dVdyedge,1,dVdyedge,1); */

							Vmath::Vmul(nquad_e,normals[0],1,curv,1,curvfromnormaltimesnormalx,1);
							Vmath::Vmul(nquad_e,normals[1],1,curv,1,curvfromnormaltimesnormaly,1);


							// MOVE edge, i.e. calculate X^{n+1} at free surface boundary
							// TODO: Implement time integration instead}
							// in: physical values of txx and txy, out:weak values of forcing function
							m_fields[0]->GetExp(ElmtID[cnt])->AddEdgeNormBoundaryInt(edgeID,EdgeExp,curvfromnormaltimesnormalx,
						                                                e_outarrayu);
							m_fields[1]->GetExp(ElmtID[cnt])->AddEdgeNormBoundaryInt(edgeID,EdgeExp,
																		curvfromnormaltimesnormaly,
						                                                e_outarrayv);

						/*	m_fields[0]->GetExp(ElmtID[cnt])->AddEdgeNormBoundaryInt(edgeID,EdgeExp,dUdxedge,
																									dVdxedge,
													                                                e_outarrayu);
							m_fields[1]->GetExp(ElmtID[cnt])->AddEdgeNormBoundaryInt(edgeID,EdgeExp,dUdyedge,
																		dVdyedge,
																		e_outarrayv); */

							// calcuate (phi, dp/dn = [N-kinvis curl x curl v].n)
						/*	Uvals = BndExp[n]->UpdateCoeffs()+BndExp[n]->GetCoeff_Offset(el);
							Vvals = BndExpv[n]->UpdateCoeffs()+BndExpv[n]->GetCoeff_Offset(el);
							// Decide if normals facing outwards
							NegateNormals = (ElExp->GetEorient(edgeID) == StdRegions::eForwards)? false:true;

							EdgeExp->NormVectorIProductWRTBase(dUdxedge,dVdxedge,Uvals,NegateNormals);
							EdgeExpv->NormVectorIProductWRTBase(dUdyedge,dVdyedge,Vvals,NegateNormals); */
						}
				 }
				 else if(type == "" || type == "TimeDependent")  // for all other boundary conditi
				 {
					 cnt += BndExp[n]->GetExpSize();
				 }
			 }

      }

    	void MovingMesh::SetFreeSurfaceVelocityBC()
         {
   			Array<OneD, const SpatialDomains::BoundaryConditionShPtr > BndConds;
   			Array<OneD, MultiRegions::ExpListSharedPtr>  BndExp, BndExpv;
   			MultiRegions::ExpList1DSharedPtr  FreeSurfaceFct;
   			Array<OneD, Array<OneD, NekDouble> > normals, locnormals;
   			Array<OneD, NekDouble> e_outarrayu, e_outarrayv;
   			Array<OneD, NekDouble> curv;
   			bool NegateNormals;

   			BndConds   = m_fields[0]->GetBndConditions();
   			BndExp     = m_fields[0]->GetBndCondExpansions();
   			BndExpv    = m_fields[1]->GetBndCondExpansions();

   			StdRegions::StdExpansionSharedPtr ElExp;
   			StdRegions::StdExpansion1DSharedPtr EdgeExp,EdgeExpv;
   			Array<OneD, int> ElmtID,EdgeID;
   			int cnt,n,el,edgeID,elmtid;
   			Array<OneD, const NekDouble> U,V;
   			Array<OneD, NekDouble> Uvals,Vvals;


   			//Fills ElmtID and EdgeID with global ids of elements id number and edges id number of the boundary expansion
   			m_fields[0]->GetBoundaryToElmtMap(ElmtID,EdgeID);

   			for(cnt = n = 0; n < m_fields[0]->GetBndConditions().num_elements(); ++n)
   			{
   				 string type = m_fields[0]->GetBndConditions()[n]->GetUserDefined().GetEquation();
   		 //for the ones marked as Drag calculate drag along boundary edge
   				 if(type == "FreeSurface")
   				 {
   			  	   if(!(FreeSurfaceFct =  boost::dynamic_pointer_cast<MultiRegions::ExpList1D>(BndExp[n])))
					{
						ASSERTL0(false,"Dynamics cast failed");
					}

   			 	   int nq= BndExp[n]->GetTotPoints();
				   Array<OneD,NekDouble> x(nq,0.0);
					Array<OneD,NekDouble> y(nq,0.0);
					Array<OneD,NekDouble> z(nq,0.0);
					Array<OneD,NekDouble> dfdx(nq,0.0);
					BndExp[n]->GetCoords(x,y,z);

					// f(x,t) = y
					FreeSurfaceFct->SetPhys(y);
					FreeSurfaceFct->PhysDeriv(0,FreeSurfaceFct->GetPhys(),dfdx);


   						for(el = 0; el < BndExp[n]->GetExpSize(); ++el,cnt++)
   						{
   							EdgeExp =  boost::dynamic_pointer_cast<StdRegions::StdExpansion1D> (BndExp[n]->GetExp(el));
   							EdgeExpv =  boost::dynamic_pointer_cast<StdRegions::StdExpansion1D> (BndExpv[n]->GetExp(el));
   							int nquad_e = EdgeExp->GetNumPoints(0);
   							int id1  = BndExp[n]->GetPhys_Offset(el);
   							edgeID = EdgeID[cnt];
   							ElExp = m_fields[0]->GetExp(ElmtID[cnt]);
   							int nqel = ElExp->GetTotPoints();

   							Array<OneD,NekDouble> dUdx(nqel,0.0);
   							Array<OneD,NekDouble> dUdy(nqel,0.0);
   							Array<OneD,NekDouble> dVdx(nqel,0.0);
   							Array<OneD,NekDouble> dVdy(nqel,0.0);


   						    U = (m_fields[0]->GetPhys())+ m_fields[0]->GetPhys_Offset(ElmtID[cnt]);
   						    V = (m_fields[1]->GetPhys())+ m_fields[1]->GetPhys_Offset(ElmtID[cnt]);


   						    // Calculating deformation tensor
   						    ElExp->PhysDeriv(0,U,dUdx);
   						    ElExp->PhysDeriv(1,U,dUdy);
   						    ElExp->PhysDeriv(0,V,dVdx);
   						    ElExp->PhysDeriv(1,V,dVdy);

   							Array<OneD,NekDouble> dUdxedge(nquad_e,0.0);
   							Array<OneD,NekDouble> dUdyedge(nquad_e,0.0);
   							Array<OneD,NekDouble> dVdxedge(nquad_e,0.0);
   							Array<OneD,NekDouble> dVdyedge(nquad_e,0.0);
   							Array<OneD,NekDouble> dUdnedge(nquad_e,0.0);
   							Array<OneD,NekDouble> dVdnedge(nquad_e,0.0);

   							//normals	=  EdgeExp->GetMetricInfo()->GetNormal();

   							// Multiply with outward normals
   							//if(ElExp->GetEorient(edgeID) == StdRegions::eBackwards)
   							//{
   							//	 Vmath::Neg(nquad_e,normals[0],1);
   							//	 Vmath::Neg(nquad_e,normals[1],1);
   							//}

   							ElExp->GetEdgePhysVals(edgeID,EdgeExp,dUdx,dUdxedge);
   							ElExp->GetEdgePhysVals(edgeID,EdgeExp,dUdy,dUdyedge);
   							ElExp->GetEdgePhysVals(edgeID,EdgeExp,dVdx,dVdxedge);
   							ElExp->GetEdgePhysVals(edgeID,EdgeExp,dVdy,dVdyedge);

   							Array<OneD,NekDouble> xedge(nquad_e,0.0);
							Array<OneD,NekDouble> yedge(nquad_e,0.0);
							Array<OneD,NekDouble> zedge(nquad_e,0.0);

							EdgeExp->GetCoords(xedge,yedge,zedge);
							NekDouble dfdx,nx,ny;

							// compute wx/wy on the edge
							for(int i=0;i<nquad_e;i++)
							{
								m_freesurfacesplines[0]->GetSplineNormal(xedge[i],nx,ny);
								dfdx = -nx/ny;
								dUdnedge[i] = (2.0*dVdyedge[i]-dUdxedge[i]-dUdyedge[i]*dfdx)*nx-dVdxedge[i]*(dfdx*nx+ny);
								dVdnedge[i] = (2.0*dUdxedge[i]-dVdyedge[i])*ny+(dVdxedge[i]*nx)/(dfdx*dfdx)+ dUdyedge[i]*nx*(1.0/(dfdx*dfdx)-1.0);
							}

   							// calcuate (phi, dp/dn = [N-kinvis curl x curl v].n)
   							Uvals = BndExp[n]->UpdateCoeffs()+BndExp[n]->GetCoeff_Offset(el);
   							Vvals = BndExpv[n]->UpdateCoeffs()+BndExpv[n]->GetCoeff_Offset(el);
   							// Decide if normals facing outwards
   							// REMARK: Produced DieSwell with dUdxedge, dVdyedge
   							EdgeExp->IProductWRTBase(dUdnedge,Uvals);
   							EdgeExpv->IProductWRTBase(dUdnedge,Vvals);
   						}
   				 }
   				 else if(type == "" || type == "TimeDependent")  // for all other boundary conditi
   				 {
   					 cnt += BndExp[n]->GetExpSize();
   				 }
   			 }

         }

    void MovingMesh::SetFreeSurfaceVelocityBC2(
    		MultiRegions::ExpListSharedPtr pPressure, NekDouble tension)
	 {
		Array<OneD, const SpatialDomains::BoundaryConditionShPtr > BndConds;
		Array<OneD, MultiRegions::ExpListSharedPtr>  BndExp, BndExpv;
		MultiRegions::ExpList1DSharedPtr  FreeSurfaceFct;
		Array<OneD, Array<OneD, NekDouble> > normals, locnormals;
		Array<OneD, NekDouble> e_outarrayu, e_outarrayv;
		Array<OneD, NekDouble> curv;
		bool NegateNormals;

		BndConds   = m_fields[0]->GetBndConditions();
		BndExp     = m_fields[0]->GetBndCondExpansions();
		BndExpv    = m_fields[1]->GetBndCondExpansions();

		StdRegions::StdExpansionSharedPtr ElExp;
		StdRegions::StdExpansion1DSharedPtr EdgeExp,EdgeExpv;
		Array<OneD, int> ElmtID,EdgeID;
		int cnt,n,el,edgeID,elmtid;
		Array<OneD, const NekDouble> U,V,P;
		Array<OneD, NekDouble> Uvals,Vvals;

		Array<OneD, NekDouble> fieldcoeffs = Array<OneD, NekDouble>(m_fields[0]->GetNcoeffs());
		Array<OneD, NekDouble> Pphys = Array<OneD, NekDouble>(m_fields[0]->GetTotPoints());

		// Project Pressure onto velocity space
		pPressure->BwdTrans(pPressure->GetCoeffs(), pPressure->UpdatePhys());

		m_fields[0]->FwdTrans_IterPerExp(pPressure->GetPhys(),fieldcoeffs);
		m_fields[0]->BwdTrans_IterPerExp(fieldcoeffs,Pphys);

		//Fills ElmtID and EdgeID with global ids of elements id number and edges id number of the boundary expansion
		m_fields[0]->GetBoundaryToElmtMap(ElmtID,EdgeID);

		for(cnt = n = 0; n < m_fields[0]->GetBndConditions().num_elements(); ++n)
		{
			 string type = m_fields[0]->GetBndConditions()[n]->GetUserDefined().GetEquation();
	 //for the ones marked as Drag calculate drag along boundary edge
			 if(type == "FreeSurface")
			 {
			   if(!(FreeSurfaceFct =  boost::dynamic_pointer_cast<MultiRegions::ExpList1D>(BndExp[n])))
				{
					ASSERTL0(false,"Dynamics cast failed");
				}

			   int nq= BndExp[n]->GetTotPoints();
			   Array<OneD,NekDouble> x(nq,0.0);
				Array<OneD,NekDouble> y(nq,0.0);
				Array<OneD,NekDouble> z(nq,0.0);
				Array<OneD,NekDouble> dfdx(nq,0.0);
				BndExp[n]->GetCoords(x,y,z);


					for(el = 0; el < BndExp[n]->GetExpSize(); ++el,cnt++)
					{
						EdgeExp =  boost::dynamic_pointer_cast<StdRegions::StdExpansion1D> (BndExp[n]->GetExp(el));
						EdgeExpv =  boost::dynamic_pointer_cast<StdRegions::StdExpansion1D> (BndExpv[n]->GetExp(el));
						int nquad_e = EdgeExp->GetNumPoints(0);
						int id1  = BndExp[n]->GetPhys_Offset(el);
						edgeID = EdgeID[cnt];
						ElExp = m_fields[0]->GetExp(ElmtID[cnt]);
						int nqel = ElExp->GetTotPoints();

						Array<OneD,NekDouble> dUdx(nqel,0.0);
						Array<OneD,NekDouble> dUdy(nqel,0.0);
						Array<OneD,NekDouble> dVdx(nqel,0.0);
						Array<OneD,NekDouble> dVdy(nqel,0.0);


						U = (m_fields[0]->GetPhys())+ m_fields[0]->GetPhys_Offset(ElmtID[cnt]);
						V = (m_fields[1]->GetPhys())+ m_fields[1]->GetPhys_Offset(ElmtID[cnt]);
						P = Pphys + m_fields[0]->GetPhys_Offset(ElmtID[cnt]);


						// Calculating deformation tensor
						ElExp->PhysDeriv(0,U,dUdx);
						ElExp->PhysDeriv(1,U,dUdy);
						ElExp->PhysDeriv(0,V,dVdx);
						ElExp->PhysDeriv(1,V,dVdy);

						Array<OneD,NekDouble> dUdxedge(nquad_e,0.0);
						Array<OneD,NekDouble> dUdyedge(nquad_e,0.0);
						Array<OneD,NekDouble> dVdxedge(nquad_e,0.0);
						Array<OneD,NekDouble> dVdyedge(nquad_e,0.0);
						Array<OneD,NekDouble> Pedge(nquad_e,0.0);
						Array<OneD,NekDouble> dUdnedge(nquad_e,0.0);
						Array<OneD,NekDouble> dVdnedge(nquad_e,0.0);

						//normals	=  EdgeExp->GetMetricInfo()->GetNormal();

						// Multiply with outward normals
						//if(ElExp->GetEorient(edgeID) == StdRegions::eBackwards)
						//{
						//	 Vmath::Neg(nquad_e,normals[0],1);
						//	 Vmath::Neg(nquad_e,normals[1],1);
						//}

						ElExp->GetEdgePhysVals(edgeID,EdgeExp,dUdx,dUdxedge);
						ElExp->GetEdgePhysVals(edgeID,EdgeExp,dUdy,dUdyedge);
						ElExp->GetEdgePhysVals(edgeID,EdgeExp,dVdx,dVdxedge);
						ElExp->GetEdgePhysVals(edgeID,EdgeExp,dVdy,dVdyedge);
						ElExp->GetEdgePhysVals(edgeID,EdgeExp,P,Pedge);

						Array<OneD,NekDouble> xedge(nquad_e,0.0);
						Array<OneD,NekDouble> yedge(nquad_e,0.0);
						Array<OneD,NekDouble> zedge(nquad_e,0.0);

						EdgeExp->GetCoords(xedge,yedge,zedge);
						NekDouble dfdx,nx,ny,curvature;

						// compute wx/wy on the edge
						for(int i=0;i<nquad_e;i++)
						{
							m_freesurfacesplines[0]->GetSplineNormal(xedge[i],nx,ny);
							m_freesurfacesplines[0]->GetSplineCurvature(xedge[i],curvature);
							curvature = abs(curvature);
							//dUdnedge[i] = Pedge[i]*nx-dUdxedge[i]*nx-dVdxedge[i]*ny+tension*curvature*nx;
							//dVdnedge[i] = Pedge[i]*ny-dUdyedge[i]*nx-dVdyedge[i]*ny+tension*curvature*ny;
							dUdnedge[i] = Pedge[i]*nx-dUdxedge[i]*nx-dVdxedge[i]*ny-tension*curvature*nx;
							dVdnedge[i] = Pedge[i]*ny-dUdyedge[i]*nx-dVdyedge[i]*ny-tension*curvature*ny;
							//cout << "(x,y)=(" << xedge[i] << "," << yedge[i] << ")";
							//cout << ", P= " << Pedge[i] << ", dUdx=" << dUdxedge[i];
							//cout << ", dVdx= " << dVdxedge[i] << ", dVdy=" << dVdyedge[i] << endl;
						}

						// calcuate (phi, dp/dn = [N-kinvis curl x curl v].n)
						Uvals = BndExp[n]->UpdateCoeffs()+BndExp[n]->GetCoeff_Offset(el);
						Vvals = BndExpv[n]->UpdateCoeffs()+BndExpv[n]->GetCoeff_Offset(el);
						// Decide if normals facing outwards
						EdgeExp->IProductWRTBase(dUdnedge,Uvals);
						EdgeExpv->IProductWRTBase(dUdnedge,Vvals);

						// Copy result in boundary expansion
						//Vmath::Vcopy(nquad_e,&dUdnedge[0], 1,&(BndExp[n]->UpdatePhys())[id1],1);
						//Vmath::Vcopy(nquad_e,&dVdnedge[0], 1,&(BndExpv[n]->UpdatePhys())[id1],1);

					}
					//BndExp[n]->FwdTrans_BndConstrained(BndExp[n]->GetPhys(),BndExp[n]->UpdateCoeffs());
					//BndExpv[n]->FwdTrans_BndConstrained(BndExpv[n]->GetPhys(),BndExpv[n]->UpdateCoeffs());
			 }
			 else if(type == "" || type == "TimeDependent")  // for all other boundary conditi
			 {
				 cnt += BndExp[n]->GetExpSize();
			 }
		 }

	 }

    void MovingMesh::AddFreeSurfaceTension2(
    		const NekDouble tension,
    		MultiRegions::ExpListSharedPtr pPressure,
    		Array<OneD, Array<OneD, NekDouble> > &Weakoutarray)
    {

    	// Reshuffle wx time vector
    	static int ncallsFBC = 1;
		Array<OneD, NekDouble> tmpx,tmpy;
		int  nlevels = m_FBCwx.num_elements();
		int  nint    = min(ncallsFBC++,nlevels);
		int cnt,n,el,edgeID,elmtid;
		// Reshuffle Bc Storage vector
		tmpx = m_FBCwx[nlevels-1];
		tmpy = m_FBCwy[nlevels-1];
		for(n = nlevels-1; n > 0; --n)
		{
			m_FBCwx[n] = m_FBCwx[n-1];
			m_FBCwy[n] = m_FBCwy[n-1];
		}
		m_FBCwx[0] = tmpx;
		m_FBCwy[0] = tmpy;

		Array<OneD, const SpatialDomains::BoundaryConditionShPtr > BndConds;
		Array<OneD, MultiRegions::ExpListSharedPtr>  BndExp, BndExpv;
		Array<OneD, NekDouble> e_outarrayu, e_outarrayv;
		bool NegateNormals;
		Array<OneD, Array<OneD, NekDouble> > normals;

   			BndConds   = m_fields[0]->GetBndConditions();
   			BndExp     = m_fields[0]->GetBndCondExpansions();
   			BndExpv    = m_fields[1]->GetBndCondExpansions();

		StdRegions::StdExpansionSharedPtr ElExp;
		StdRegions::StdExpansion1DSharedPtr EdgeExp,EdgeExpv;
		Array<OneD, int> ElmtID,EdgeID;

		Array<OneD, const NekDouble> U,V,P;
		Array<OneD, NekDouble> Uvals,Vvals;

		Array<OneD, NekDouble> Pphys = Array<OneD, NekDouble>(m_fields[0]->GetTotPoints());
		Array<OneD, NekDouble> fieldcoeffs = Array<OneD, NekDouble>(m_fields[0]->GetNcoeffs());

		// Project Pressure onto velocity space
		pPressure->BwdTrans(pPressure->GetCoeffs(), pPressure->UpdatePhys());

		m_fields[0]->FwdTrans_IterPerExp(pPressure->GetPhys(),fieldcoeffs);
		m_fields[0]->BwdTrans_IterPerExp(fieldcoeffs,Pphys);

		//Fills ElmtID and EdgeID with global ids of elements id number and edges id number of the boundary expansion
		m_fields[0]->GetBoundaryToElmtMap(ElmtID,EdgeID);

		int nqbc;
		for(nqbc = n = 0; n < m_fields[0]->GetBndConditions().num_elements(); ++n)
		{
			  string type = m_fields[0]->GetBndConditions()[n]->GetUserDefined().GetEquation();
			  if(type == "FreeSurface")
			  {
				 nqbc += BndExp[n]->GetTotPoints();
			  }
		}


		for(cnt = n = 0; n < m_fields[0]->GetBndConditions().num_elements(); ++n)
		{
			 string type = m_fields[0]->GetBndConditions()[n]->GetUserDefined().GetEquation();
			 if(type == "FreeSurface")
			 {
				for(el = 0; el < BndExp[n]->GetExpSize(); ++el,cnt++)
				{
					EdgeExp =  boost::dynamic_pointer_cast<StdRegions::StdExpansion1D> (BndExp[n]->GetExp(el));
					EdgeExpv =  boost::dynamic_pointer_cast<StdRegions::StdExpansion1D> (BndExpv[n]->GetExp(el));
					int nquad_e = EdgeExp->GetNumPoints(0);
					edgeID = EdgeID[cnt];
					ElExp = m_fields[0]->GetExp(ElmtID[cnt]);
					int nqel = ElExp->GetTotPoints();

					Array<OneD,NekDouble> dUdx(nqel,0.0);
					Array<OneD,NekDouble> dUdy(nqel,0.0);
					Array<OneD,NekDouble> dVdx(nqel,0.0);
					Array<OneD,NekDouble> dVdy(nqel,0.0);

					U = (m_fields[0]->GetPhys())+ m_fields[0]->GetPhys_Offset(ElmtID[cnt]);
					V = (m_fields[1]->GetPhys())+ m_fields[1]->GetPhys_Offset(ElmtID[cnt]);
					P = Pphys + m_fields[0]->GetPhys_Offset(ElmtID[cnt]);

					// Calculating deformation tensor
					ElExp->PhysDeriv(0,U,dUdx);
					ElExp->PhysDeriv(1,U,dUdy);
					ElExp->PhysDeriv(0,V,dVdx);
					ElExp->PhysDeriv(1,V,dVdy);


					// TODO:test that right curvature values are taken for comparison

					Array<OneD,NekDouble> dUdxedge(nquad_e,0.0);
					Array<OneD,NekDouble> dUdyedge(nquad_e,0.0);
					Array<OneD,NekDouble> dVdxedge(nquad_e,0.0);
					Array<OneD,NekDouble> dVdyedge(nquad_e,0.0);
					Array<OneD,NekDouble> Pedge(nquad_e,0.0);

					ElExp->GetEdgePhysVals(edgeID,EdgeExp,dUdx,dUdxedge);
					ElExp->GetEdgePhysVals(edgeID,EdgeExp,dUdy,dUdyedge);
					ElExp->GetEdgePhysVals(edgeID,EdgeExp,dVdx,dVdxedge);
					ElExp->GetEdgePhysVals(edgeID,EdgeExp,dVdy,dVdyedge);
					ElExp->GetEdgePhysVals(edgeID,EdgeExp,P,Pedge);

					Array<OneD,NekDouble> xedge(nquad_e,0.0);
					Array<OneD,NekDouble> yedge(nquad_e,0.0);
					Array<OneD,NekDouble> zedge(nquad_e,0.0);

					EdgeExp->GetCoords(xedge,yedge,zedge);
					NekDouble curvature,nx,ny;
					int id = BndExp[n]->GetPhys_Offset(el);
				/*	normals	=  EdgeExp->GetMetricInfo()->GetNormal();

					if(ElExp->GetEorient(edgeID) == StdRegions::eBackwards)
					{
						 Vmath::Neg(nquad_e,normals[0],1);
						 Vmath::Neg(nquad_e,normals[1],1);
					} */

					for(int i=0;i<nquad_e;i++)
					{
						m_freesurfacesplines[0]->GetSplineNormal(xedge[i],nx,ny);
						m_freesurfacesplines[0]->GetSplineCurvature(xedge[i],curvature);
						//curvature = abs(curvature);
						//nx = normals[0][i];
						//ny = normals[1][i];
						m_FBCwx[0][id+i] = -tension*curvature*nx;//-dUdxedge[i]*nx-dVdxedge[i]*ny;//+ Pedge[i]*nx;
						m_FBCwy[0][id+i] = -tension*curvature*ny;//-dUdyedge[i]*nx-dVdyedge[i]*ny;//+ Pedge[i]*ny;
						//curvfromnormaltimesnormalx[i] = tension*curvature*nx-dUdxedge[i]*nx-dVdxedge[i]*ny;//+ Pedge[i]*nx;
						//curvfromnormaltimesnormaly[i] = tension*curvature*ny-dUdyedge[i]*nx-dVdyedge[i]*ny;//+ Pedge[i]*ny;

					}
				}
			 }
			 else // for all other boundary conditi
			 {
				 cnt += BndExp[n]->GetExpSize();
			 }
		 }

		// Extrapolate to n+1
		Vmath::Smul(nqbc,kFreeSurfaceBCsExtrapolation[nint-1][nint-1],
				m_FBCwx[nint-1],1,m_FBCwx[nlevels-1],1);
		Vmath::Smul(nqbc,kFreeSurfaceBCsExtrapolation[nint-1][nint-1],
				m_FBCwy[nint-1],1,m_FBCwy[nlevels-1],1);
		for(int n = 0; n < nint-1; ++n)
		{
			Vmath::Svtvp(nqbc,kFreeSurfaceBCsExtrapolation[nint-1][n],
					m_FBCwx[n],1,m_FBCwx[nlevels-1],1,
					m_FBCwx[nlevels-1],1);
			Vmath::Svtvp(nqbc,kFreeSurfaceBCsExtrapolation[nint-1][n],
							m_FBCwy[n],1,m_FBCwy[nlevels-1],1,
							m_FBCwy[nlevels-1],1);
		}

		for(cnt = n = 0; n < m_fields[0]->GetBndConditions().num_elements(); ++n)
		{
			 string type = m_fields[0]->GetBndConditions()[n]->GetUserDefined().GetEquation();
	 //for the ones marked as Drag calculate drag along boundary edge
			 if(type == "FreeSurface")
			 {
					for(el = 0; el < BndExp[n]->GetExpSize(); ++el,cnt++)
					{
						EdgeExp =  boost::dynamic_pointer_cast<StdRegions::StdExpansion1D> (BndExp[n]->GetExp(el));
						EdgeExpv =  boost::dynamic_pointer_cast<StdRegions::StdExpansion1D> (BndExpv[n]->GetExp(el));
						int nquad_e = EdgeExp->GetNumPoints(0);
						int id  = BndExp[n]->GetPhys_Offset(el);
						edgeID = EdgeID[cnt];
						Array<OneD,NekDouble> curvfromnormaltimesnormalx(nquad_e,0.0);
						Array<OneD,NekDouble> curvfromnormaltimesnormaly(nquad_e,0.0);

						Vmath::Vcopy(nquad_e,&(m_FBCwx[nlevels-1])[id],1,&(curvfromnormaltimesnormalx)[0],1);
						Vmath::Vcopy(nquad_e,&(m_FBCwy[nlevels-1])[id],1,&(curvfromnormaltimesnormaly)[0],1);

						e_outarrayu = Weakoutarray[0] + m_fields[0]->GetCoeff_Offset(ElmtID[cnt]);
						e_outarrayv = Weakoutarray[1] + m_fields[1]->GetCoeff_Offset(ElmtID[cnt]);

					/*	 for(int i=0;i<nquad_e;i++)
														{
							 cout << "BCx:" << curvfromnormaltimesnormalx[i] << endl;
							cout << "BCy:" << curvfromnormaltimesnormaly[i] << endl;

														} */

						m_fields[0]->GetExp(ElmtID[cnt])->AddEdgeNormBoundaryInt(edgeID,EdgeExp,curvfromnormaltimesnormalx,
																	e_outarrayu);
						m_fields[1]->GetExp(ElmtID[cnt])->AddEdgeNormBoundaryInt(edgeID,EdgeExpv,
																	curvfromnormaltimesnormaly,
																	e_outarrayv);
					}
			 }
			 else // for all other boundary conditi
			 {
				 cnt += BndExp[n]->GetExpSize();
			 }
		}
    }

    void MovingMesh::WriteValuesAlongFreeSurfaceToFile( string fname, int outputcount)
    {
    	Array<OneD, const SpatialDomains::BoundaryConditionShPtr > BndConds;
		Array<OneD, MultiRegions::ExpListSharedPtr>  BndExp;
		//static int outputcnt2 = 1;
		BndConds   = m_fields[m_meshvelocity[0]]->GetBndConditions();
		BndExp     = m_fields[m_meshvelocity[0]]->GetBndCondExpansions();

		StdRegions::StdExpansionSharedPtr ElExp;
		StdRegions::StdExpansion1DSharedPtr EdgeExp;
		Array<OneD, int> ElmtID,EdgeID;
		int cnt, n,el,edgeID,elmID;
		Array<OneD, const NekDouble> Uphyselement, Vphyselement,wxel,wyel;

		//Fills ElmtID and EdgeID with global ids of elements id number and edges id number of the boundary expansion
		m_fields[m_meshvelocity[0]]->GetBoundaryToElmtMap(ElmtID,EdgeID);

		stringstream convert;
		convert << "ValuesFS" << fname << "_" << outputcount << ".txt";
		string filename = convert.str();
		ofstream outfile(filename.c_str());

		int ispline=0;
		for(cnt = n = 0; n < m_fields[m_meshvelocity[0]]->GetBndConditions().num_elements(); ++n)
		{
			 string type = m_fields[m_meshvelocity[0]]->GetBndConditions()[n]->GetUserDefined().GetEquation();
	 //for the ones marked as Drag calculate drag along boundary edge
			 if(type == "FreeSurface")
			 {
				for(el = 0; el < BndExp[n]->GetExpSize(); ++el,cnt++)
				{
					Uphyselement = (m_fields[0]->GetPhys())+ m_fields[0]->GetPhys_Offset(ElmtID[cnt]);
					Vphyselement = (m_fields[1]->GetPhys())+ m_fields[1]->GetPhys_Offset(ElmtID[cnt]);
					wxel = (m_fields[m_meshvelocity[0]]->GetPhys())+ m_fields[m_meshvelocity[0]]->GetPhys_Offset(ElmtID[cnt]);
					wyel = (m_fields[m_meshvelocity[1]]->GetPhys())+ m_fields[m_meshvelocity[1]]->GetPhys_Offset(ElmtID[cnt]);

					ElExp = m_fields[m_meshvelocity[0]]->GetExp(ElmtID[cnt]);
					EdgeExp =  boost::dynamic_pointer_cast<StdRegions::StdExpansion1D> (BndExp[n]->GetExp(el));
					int nquad_e = EdgeExp->GetNumPoints(0);
					Array<OneD,NekDouble> xedge(nquad_e,0.0);
					Array<OneD,NekDouble> yedge(nquad_e,0.0);
					Array<OneD,NekDouble> zedge(nquad_e,0.0);
					Array<OneD,NekDouble> Uphysedge(nquad_e,0.0);
					Array<OneD,NekDouble> Vphysedge(nquad_e,0.0);
					Array<OneD,NekDouble> wxedge(nquad_e,0.0);
					Array<OneD,NekDouble> wyedge(nquad_e,0.0);

					EdgeExp->GetCoords(xedge,yedge,zedge);
					edgeID = EdgeID[cnt];
					ElExp->GetEdgePhysVals(edgeID,EdgeExp,Uphyselement,Uphysedge);
					ElExp->GetEdgePhysVals(edgeID,EdgeExp,Vphyselement,Vphysedge);
					ElExp->GetEdgePhysVals(edgeID,EdgeExp,wxel,wxedge);
					ElExp->GetEdgePhysVals(edgeID,EdgeExp,wyel,wyedge);

					EdgeExp->GetCoords(xedge,yedge,zedge);


					for(int i=0;i<nquad_e;i++)
					{
						outfile << xedge[i] << "\t"<< yedge[i] <<  "\t";
						outfile <<  Uphysedge[i] <<"\t"<< Vphysedge[i] <<  "\t";
						outfile << wxedge[i] << "\t" << wyedge[i] << endl;
					}
				}
				ispline++;
			 }
			 else if(type == "" || type == "TimeDependent")  // for all other boundary conditi
			 {
				 cnt += BndExp[n]->GetExpSize();
			 }
		 }

    }

    void MovingMesh::WriteStressAlongFreeSurfaceToFile(string fname, int outputcount, MultiRegions::ExpListSharedPtr pPressure)
    {
    	Array<OneD, const SpatialDomains::BoundaryConditionShPtr > BndConds;
		Array<OneD, MultiRegions::ExpListSharedPtr>  BndExp;
		//static int outputcnt2 = 1;
		BndConds   = m_fields[m_meshvelocity[0]]->GetBndConditions();
		BndExp     = m_fields[m_meshvelocity[0]]->GetBndCondExpansions();

		StdRegions::StdExpansionSharedPtr ElExp;
		StdRegions::StdExpansion1DSharedPtr EdgeExp;
		Array<OneD, int> ElmtID,EdgeID;
		int i,cnt, n,el,edgeID,elmID;
		Array<OneD, const NekDouble> Uphyselement, Vphyselement,Pphyselement,wxel,wyel;
		int maxpts =0;
		NekDouble nx,ny;
        Array<OneD, NekDouble> fieldcoeffs = Array<OneD, NekDouble>(m_fields[0]->GetNcoeffs());
        Array<OneD, NekDouble> Pphys = Array<OneD, NekDouble>(m_fields[0]->GetTotPoints());

		//Fills ElmtID and EdgeID with global ids of elements id number and edges id number of the boundary expansion
		m_fields[m_meshvelocity[0]]->GetBoundaryToElmtMap(ElmtID,EdgeID);

	    for(cnt = n = 0; n < BndConds.num_elements(); ++n)
	    {
	//Loop over all elements in boundary region
	    	//cnt starts with element 0 then going through elements
	        for(i = 0; i < BndExp[n]->GetExpSize(); ++i)
	        {
	            maxpts = max(maxpts, m_fields[m_meshvelocity[0]]->GetExp(ElmtID[cnt++])->GetTotPoints());
	        }
	    }

	    ASSERTL0(m_spacedim == 2,"Not set up for 3D expansions");

	    // Project Pressure onto velocity space
        pPressure->BwdTrans(pPressure->GetCoeffs(), pPressure->UpdatePhys());

        m_fields[0]->FwdTrans_IterPerExp(pPressure->GetPhys(),fieldcoeffs);
        m_fields[0]->BwdTrans_IterPerExp(fieldcoeffs,Pphys);

	    Array<OneD, NekDouble> dUdx(maxpts);
	    Array<OneD, NekDouble> dUdy(maxpts);
	    Array<OneD, NekDouble> dVdx(maxpts);
	    Array<OneD, NekDouble> dVdy(maxpts);
	    Array<OneD, NekDouble> cauchyxx(maxpts);
	    Array<OneD, NekDouble> cauchyxy(maxpts);
	    Array<OneD, NekDouble> cauchyyy(maxpts);


		int ispline=0;
		for(cnt = n = 0; n < m_fields[m_meshvelocity[0]]->GetBndConditions().num_elements(); ++n)
		{
			 string type = m_fields[m_meshvelocity[0]]->GetBndConditions()[n]->GetUserDefined().GetEquation();
	 //for the ones marked as Drag calculate drag along boundary edge
			 if(type == "FreeSurface")
			 {
				for(el = 0; el < BndExp[n]->GetExpSize(); ++el,cnt++)
				{
					int offset = m_fields[0]->GetPhys_Offset(ElmtID[cnt]);
					Uphyselement = (m_fields[0]->GetPhys())+ offset;
					Vphyselement = (m_fields[1]->GetPhys())+ m_fields[1]->GetPhys_Offset(ElmtID[cnt]);
					Pphyselement = Pphys + offset;
					wxel = (m_fields[m_meshvelocity[0]]->GetPhys())+ m_fields[m_meshvelocity[0]]->GetPhys_Offset(ElmtID[cnt]);
					wyel = (m_fields[m_meshvelocity[1]]->GetPhys())+ m_fields[m_meshvelocity[1]]->GetPhys_Offset(ElmtID[cnt]);

					ElExp = m_fields[m_meshvelocity[0]]->GetExp(ElmtID[cnt]);
					EdgeExp =  boost::dynamic_pointer_cast<StdRegions::StdExpansion1D> (BndExp[n]->GetExp(el));

					// Calculate Stress Tensor
		            // Calculating deformation tensor
					ElExp->PhysDeriv(0,Uphyselement,dUdx);
					ElExp->PhysDeriv(1,Uphyselement,dUdy);
					ElExp->PhysDeriv(0,Vphyselement,dVdx);
					ElExp->PhysDeriv(1,Vphyselement,dVdy);

					int nq     = ElExp->GetTotPoints();

					Vmath::Smul(nq,2.0,dUdx,1,dUdx,1);
				    Vmath::Vsub(nq,dUdx,1,Pphyselement,1,cauchyxx,1);

				    Vmath::Smul(nq,2.0,dVdy,1,dVdy,1);
				    Vmath::Vsub(nq,dVdy,1,Pphyselement,1,cauchyyy,1);

				    Vmath::Vadd(nq,dUdy,1,dVdx,1,cauchyxy,1);

					int nquad_e = EdgeExp->GetNumPoints(0);
					Array<OneD,NekDouble> xedge(nquad_e,0.0);
					Array<OneD,NekDouble> yedge(nquad_e,0.0);
					Array<OneD,NekDouble> zedge(nquad_e,0.0);
					Array<OneD,NekDouble> Uphysedge(nquad_e,0.0);
					Array<OneD,NekDouble> Vphysedge(nquad_e,0.0);
					Array<OneD,NekDouble> wxedge(nquad_e,0.0);
					Array<OneD,NekDouble> wyedge(nquad_e,0.0);

					EdgeExp->GetCoords(xedge,yedge,zedge);
					edgeID = EdgeID[cnt];
					ElExp->GetEdgePhysVals(edgeID,EdgeExp,cauchyxx,dUdx);
					ElExp->GetEdgePhysVals(edgeID,EdgeExp,cauchyxy,dUdy);
					ElExp->GetEdgePhysVals(edgeID,EdgeExp,cauchyyy,dVdy);

					EdgeExp->GetCoords(xedge,yedge,zedge);

					NekDouble tSn,nSn;
					// compute wx/wy on the edge
					stringstream convert;
					convert << "StressValues" << fname << "_" << outputcount << ".txt";
					string filename = convert.str();

					ofstream outfile(filename.c_str(),ios::app);

					for(int i=0;i<nquad_e;i++)
					{
						// wx = 0
						m_freesurfacesplines[0]->GetSplineNormal(xedge[i],nx,ny);
						tSn = ny*dUdx[i]*nx+dUdy[i]*(ny*ny-nx*nx)-nx*dVdy[i]*ny;
						nSn = nx*dUdx[i]*nx+ 2*nx*ny*dUdy[i] + ny*dVdy[i]*ny;
						outfile << xedge[i] << "\t" << tSn << "\t" << nSn << "\n";

					}
				}
				ispline++;
			 }
			 else if(type == "" || type == "TimeDependent")  // for all other boundary conditi
			 {
				 cnt += BndExp[n]->GetExpSize();
			 }
		 }

    }

    void MovingMesh::WriteArrayAlongFreeSurfaceToFile( string fname, int outputcount,Array<OneD, Array<OneD, NekDouble> > &inarray)
       {
       	Array<OneD, const SpatialDomains::BoundaryConditionShPtr > BndConds;
   		Array<OneD, MultiRegions::ExpListSharedPtr>  BndExp;
   		//static int outputcnt2 = 1;
   		BndConds   = m_fields[m_meshvelocity[0]]->GetBndConditions();
   		BndExp     = m_fields[m_meshvelocity[0]]->GetBndCondExpansions();

   		StdRegions::StdExpansionSharedPtr ElExp;
   		StdRegions::StdExpansion1DSharedPtr EdgeExp;
   		Array<OneD, int> ElmtID,EdgeID;
   		int i,cnt, n,el,edgeID,elmID;
   		int ninarray = inarray.num_elements();
   		Array<OneD, Array<OneD, const NekDouble> > inarrayelement(ninarray);
   		int maxpts =0;
   		NekDouble nx,ny;

   		for(int i=0;i<ninarray;i++)
   		{
   			inarrayelement[i] = Array<OneD, const NekDouble>(m_phystot);
   		}

   		//Fills ElmtID and EdgeID with global ids of elements id number and edges id number of the boundary expansion
   		m_fields[m_meshvelocity[0]]->GetBoundaryToElmtMap(ElmtID,EdgeID);
   		int ispline=0;
   		for(cnt = n = 0; n < m_fields[m_meshvelocity[0]]->GetBndConditions().num_elements(); ++n)
   		{
   			 string type = m_fields[m_meshvelocity[0]]->GetBndConditions()[n]->GetUserDefined().GetEquation();
   	 //for the ones marked as Drag calculate drag along boundary edge
   			 if(type == "FreeSurface")
   			 {
   				for(el = 0; el < BndExp[n]->GetExpSize(); ++el,cnt++)
   				{
   					int offset = m_fields[0]->GetPhys_Offset(ElmtID[cnt]);
   					for(int i=0;i<ninarray;i++)
					{
						inarrayelement[i] = inarray[i]+offset;
					}

   					ElExp = m_fields[m_meshvelocity[0]]->GetExp(ElmtID[cnt]);
   					EdgeExp =  boost::dynamic_pointer_cast<StdRegions::StdExpansion1D> (BndExp[n]->GetExp(el));

   					int nquad_e = EdgeExp->GetNumPoints(0);
   					Array<OneD,NekDouble> xedge(nquad_e,0.0);
   					Array<OneD,NekDouble> yedge(nquad_e,0.0);
   					Array<OneD,NekDouble> zedge(nquad_e,0.0);
   					Array<OneD, Array<OneD, NekDouble> > inarrayedge(ninarray);

   					for(int i=0;i<ninarray;i++)
					{
   						inarrayedge[i] = Array<OneD,NekDouble> (nquad_e,0.0);
					}

   					EdgeExp->GetCoords(xedge,yedge,zedge);
   					edgeID = EdgeID[cnt];
   					for(int i=0;i<ninarray;i++)
					{
   						ElExp->GetEdgePhysVals(edgeID,EdgeExp,inarrayelement[i],inarrayedge[i]);
					}

   					EdgeExp->GetCoords(xedge,yedge,zedge);

   					// compute wx/wy on the edge
   					stringstream convert;
   					convert << "ArrayValues" << fname << "_" << outputcount << ".txt";
   					string filename = convert.str();

   					ofstream outfile(filename.c_str(),ios::app);

   					for(int i=0;i<nquad_e;i++)
   					{
   						// wx = 0
   						for(int j=0;j<ninarray;j++)
   						{
   							outfile << inarrayedge[j][i] << "\t";
   						}
   						outfile << "\n";
   					}
   				}
   				ispline++;
   			 }
   			 else if(type == "" || type == "TimeDependent")  // for all other boundary conditi
   			 {
   				 cnt += BndExp[n]->GetExpSize();
   			 }
   		 }

       }





    void MovingMesh::SetInitialMeshSloshing()
    {
    	Array<OneD, const SpatialDomains::BoundaryConditionShPtr > BndConds;
		Array<OneD, MultiRegions::ExpListSharedPtr>  BndExp;

		// Grab boundary regions and expansions of velocity u
		BndConds   = m_fields[m_meshvelocity[0]]->GetBndConditions();
		BndExp     = m_fields[m_meshvelocity[0]]->GetBndCondExpansions();

		StdRegions::StdExpansion1DSharedPtr EdgeExp;
		StdRegions::StdExpansionSharedPtr ElExp;
		SpatialDomains::SegGeomSharedPtr edge;
		int vertID1,vertID2;
		SpatialDomains::VertexComponentSharedPtr vertex1,vertex2;
		SpatialDomains::MeshGraph2DSharedPtr mesh2D;

		Array<OneD, int> ElmtID,EdgeID;
		int cnt,n,el,edgeID,elmID;
		Array<OneD, const NekDouble> Uphyselement, Vphyselement, wxel,wyel;

		//Fills ElmtID and EdgeID with global ids of elements id number and edges id number of the boundary expansion
		m_fields[m_meshvelocity[0]]->GetBoundaryToElmtMap(ElmtID,EdgeID);

		for(cnt = n = 0; n < m_fields[0]->GetBndConditions().num_elements(); ++n)
		{
			  string type = m_fields[0]->GetBndConditions()[n]->GetUserDefined().GetEquation();
	  //for the ones marked as Drag calculate drag along boundary edge
			  if(type == "FreeSurface")
			  {
				 //loop over all elements along the boundary region
				for(el = 0; el < BndExp[n]->GetExpSize(); ++el,cnt++)
				{
					ElExp = m_fields[m_meshvelocity[0]]->GetExp(ElmtID[cnt]);
					EdgeExp =  boost::dynamic_pointer_cast<StdRegions::StdExpansion1D> (BndExp[n]->GetExp(el));

					int nquad_e = EdgeExp->GetNumPoints(0);
					Array<OneD,NekDouble> xedge(nquad_e,0.0);
					Array<OneD,NekDouble> yedge(nquad_e,0.0);
					Array<OneD,NekDouble> zedge(nquad_e,0.0);

					// EdgeID in within element
					edgeID = EdgeID[cnt];

					EdgeExp->GetCoords(xedge,yedge,zedge);
					NekDouble function,a;
					NekDouble PI = 3.14159265;

					for(int i=0;i<nquad_e;i++)
					{
						a = 0.05;
						//-0.05*sin(2*PI/2*x);
						function = -a*sin(PI/2*xedge[i]);
						yedge[i]= yedge[i] + function;
					}

					edgeID= m_mesh2D->GetEidFromElmt(ElExp->DetExpansionType(),
									edgeID, ElmtID[cnt]);

					edge = m_mesh2D->GetEdge(edgeID);


					vertID1=edge->GetVid(0);
					vertID2=edge->GetVid(1);

					// Make changes to m_graph

					 vertex1 = m_graph->GetVertex(vertID1);
					 vertex2 = m_graph->GetVertex(vertID2);

					vertex1->UpdatePosition(xedge[0],yedge[0],zedge[0]);
					vertex2->UpdatePosition(xedge[nquad_e-1],yedge[nquad_e-1],zedge[nquad_e-1]);

					 m_graph->CreateCurvedEdge(edgeID,xedge,yedge,zedge);
				}

			  }
			  else if(type == "" || type == "TimeDependent")  // for all other boundary conditi
			 {
				 cnt += BndExp[n]->GetExpSize();
			 }
		  }
    }

    void MovingMesh::SetInitialMeshDieSwell()
       {
       	Array<OneD, const SpatialDomains::BoundaryConditionShPtr > BndConds;
   		Array<OneD, MultiRegions::ExpListSharedPtr>  BndExp;

   		// Grab boundary regions and expansions of velocity u
   		BndConds   = m_fields[m_meshvelocity[0]]->GetBndConditions();
   		BndExp     = m_fields[m_meshvelocity[0]]->GetBndCondExpansions();

   		StdRegions::StdExpansion1DSharedPtr EdgeExp;
   		StdRegions::StdExpansionSharedPtr ElExp;
   		SpatialDomains::SegGeomSharedPtr edge;
   		int vertID1,vertID2;
   		SpatialDomains::VertexComponentSharedPtr vertex1,vertex2;
   		SpatialDomains::MeshGraph2DSharedPtr mesh2D;
   		NekDouble x;

   		Array<OneD, int> ElmtID,EdgeID;
   		int cnt,n,el,edgeID,elmID;
   		Array<OneD, const NekDouble> Uphyselement, Vphyselement, wxel,wyel;

   		//Fills ElmtID and EdgeID with global ids of elements id number and edges id number of the boundary expansion
   		m_fields[m_meshvelocity[0]]->GetBoundaryToElmtMap(ElmtID,EdgeID);

   		for(cnt = n = 0; n < m_fields[m_meshvelocity[0]]->GetBndConditions().num_elements(); ++n)
   		{
   			  string type = m_fields[m_meshvelocity[0]]->GetBndConditions()[n]->GetUserDefined().GetEquation();
   	  //for the ones marked as Drag calculate drag along boundary edge
   			  if(type == "FreeSurface")
   			  {
   				 //loop over all elements along the boundary region
   				for(el = 0; el < BndExp[n]->GetExpSize(); ++el,cnt++)
   				{
   					ElExp = m_fields[m_meshvelocity[0]]->GetExp(ElmtID[cnt]);
   					EdgeExp =  boost::dynamic_pointer_cast<StdRegions::StdExpansion1D> (BndExp[n]->GetExp(el));

   					int nquad_e = EdgeExp->GetNumPoints(0);
   					Array<OneD,NekDouble> xedge(nquad_e,0.0);
   					Array<OneD,NekDouble> yedge(nquad_e,0.0);
   					Array<OneD,NekDouble> zedge(nquad_e,0.0);

   					// EdgeID in within element
   					edgeID = EdgeID[cnt];

   					EdgeExp->GetCoords(xedge,yedge,zedge);
   					NekDouble function,a;
   					NekDouble PI = 3.14159265;

   					for(int i=0;i<nquad_e;i++)
   					{
   						x = xedge[i]-10.0; //if die exit is at x=10
   						//-0.05*sin(2*PI/2*x);
   						function = 0.4*(x)/((x)+exp(-x));
   						yedge[i]= yedge[i] + function;
   						cout << "(x,y)=(" << xedge[i] << "," << yedge[i] << ")" << endl;
   					}

   					edgeID= m_mesh2D->GetEidFromElmt(ElExp->DetExpansionType(),
   									edgeID, ElmtID[cnt]);

   					edge = m_mesh2D->GetEdge(edgeID);


   					vertID1=edge->GetVid(0);
   					vertID2=edge->GetVid(1);

   					// Make changes to m_graph

   					 vertex1 = m_graph->GetVertex(vertID1);
   					 vertex2 = m_graph->GetVertex(vertID2);

   					vertex1->UpdatePosition(xedge[0],yedge[0],zedge[0]);
   					vertex2->UpdatePosition(xedge[nquad_e-1],yedge[nquad_e-1],zedge[nquad_e-1]);

   					 m_graph->CreateCurvedEdge(edgeID,xedge,yedge,zedge);
   				}

   			  }
   			  else // for all other boundary conditi
   			 {
   				 cnt += BndExp[n]->GetExpSize();
   			 }
   		  }
       }

    void MovingMesh::Subdivwtu(  const Array<OneD, const Array<OneD, NekDouble> > &inarray,
    	                                                 Array<OneD, Array<OneD, NekDouble> > &outarray)
    {
    	int VelDim     = m_velocity.num_elements();
		Array<OneD, Array<OneD, NekDouble> > velocity(VelDim);
		Array<OneD, Array<OneD, NekDouble> > meshvelocity(VelDim);
		Array<OneD, Array<OneD, NekDouble> > divwtu(VelDim);
    	Array<OneD, NekDouble > dwxdx(m_phystot,0.0);
    	Array<OneD, NekDouble > dwydy(m_phystot,0.0);
    	Array<OneD, NekDouble > divw(m_phystot,0.0);

    	for(int i = 0; i < VelDim; ++i)
		{
			velocity[i] = inarray[m_velocity[i]];
			meshvelocity[i] = m_fields[m_meshvelocity[i]]->GetPhys();
			divwtu[i] =  Array<OneD, NekDouble > (m_phystot,0.0);
		}

    	m_fields[m_velocity[0]]->PhysDeriv(0,meshvelocity[0],dwxdx);
    	m_fields[m_velocity[0]]->PhysDeriv(1,meshvelocity[1],dwydy);

    	Vmath::Vadd(m_phystot,dwxdx,1,dwydy,1,divw,1);

    	for(int i=0;i<VelDim;i++)
    	{
    		Vmath::Vmul(m_phystot,divw,1,velocity[i],1,divwtu[i],1);
    		Vmath::Vsub(m_phystot,outarray[m_velocity[i]],1,divwtu[i],1,outarray[m_velocity[i]],1);
    	}

    }

	void MovingMesh::SetDirichletUnzeroBCALE()
	{
		int i,el,n,cnt, offset, phys_offset,nq,cnt2;
		int edgeID, elmtid;
		Array<OneD, NekDouble> e_outarray;
		Array<OneD, const NekDouble> U,V;
		Array<OneD, Array<OneD, NekDouble> > normals(2);
		StdRegions::StdExpansion1DSharedPtr EdgeExp;
		StdRegions::StdExpansionSharedPtr ElExp;
		MultiRegions::ExpList1DSharedPtr exp1d;
		Array<OneD, const NekDouble> Uphyselement, Vphyselement, WXel, WYel;
		int ispline = 0;
		NekDouble nx,ny;

		Array<OneD, int> ElmtID,EdgeID;

		//Fills ElmtID and EdgeID with global ids of elements id number and edges id number of the boundary expansion
		m_fields[m_velocity[1]]->GetBoundaryToElmtMap(ElmtID,EdgeID);

	  	Array<OneD, const SpatialDomains::BoundaryConditionShPtr > BndConds;
		Array<OneD, MultiRegions::ExpListSharedPtr>  BndExp;

		// Grab boundary regions and expansions of velocity u
		BndConds   = m_fields[m_velocity[1]]->GetBndConditions();
		BndExp     = m_fields[m_velocity[1]]->GetBndCondExpansions();

    	// Reshuffle ubc time vector
    	static int ncallsFBCu = 1;
		Array<OneD, NekDouble> tmpx;
		int  nint    = min(ncallsFBCu++,1);
		int  nlevels = m_FBCu.num_elements();
		// Reshuffle Bc Storage vector
		tmpx = m_FBCu[nlevels-1];
		for(int n = nlevels-1; n > 0; --n)
		{
			m_FBCu[n] = m_FBCu[n-1];
		}
		m_FBCu[0] = tmpx;

		cnt2=0;

		for(cnt = n = 0; n < BndConds.num_elements(); ++n)
		{
			if((BndConds[n]->GetUserDefined().GetEquation()== "FreeSurface"))
					//&& (BndConds[n]->GetBoundaryConditionType() == SpatialDomains::eDirichlet))
			{
				nq=BndExp[n]->GetTotPoints();
				Array<OneD,NekDouble> UBC(nq,0.0);

				for(el = 0; el < BndExp[n]->GetExpSize(); ++el,cnt++)
				{
					Uphyselement = (m_fields[0]->GetPhys())+ m_fields[0]->GetPhys_Offset(ElmtID[cnt]);
					WXel = (m_fields[m_meshvelocity[0]]->GetPhys())+ m_fields[m_meshvelocity[0]]->GetPhys_Offset(ElmtID[cnt]);
					WYel = (m_fields[m_meshvelocity[1]]->GetPhys())+ m_fields[m_meshvelocity[1]]->GetPhys_Offset(ElmtID[cnt]);
					ElExp = m_fields[1]->GetExp(ElmtID[cnt]);
					EdgeExp =  boost::dynamic_pointer_cast<StdRegions::StdExpansion1D> (BndExp[n]->GetExp(el));
					int nquad_e = EdgeExp->GetNumPoints(0);
					Array<OneD,NekDouble> Uphysedge(nquad_e,0.0);
					Array<OneD,NekDouble> Vphysedge(nquad_e,0.0);
					Array<OneD,NekDouble> Wxphysedge(nquad_e,0.0);
					Array<OneD,NekDouble> Wyphysedge(nquad_e,0.0);
					Array<OneD,NekDouble> nxny(nquad_e,0.0);
					Array<OneD,NekDouble> xedge(nquad_e,0.0);
					Array<OneD,NekDouble> yedge(nquad_e,0.0);
					Array<OneD,NekDouble> zedge(nquad_e,0.0);
					EdgeExp->GetCoords(xedge,yedge,zedge);

					// EdgeID in within element
					edgeID = EdgeID[cnt];
					ElExp->GetEdgePhysVals(edgeID,EdgeExp,Uphyselement,Uphysedge);
					ElExp->GetEdgePhysVals(edgeID,EdgeExp,WXel,Wxphysedge);
					ElExp->GetEdgePhysVals(edgeID,EdgeExp,WYel,Wyphysedge);

					int id1  = BndExp[n]->GetPhys_Offset(el);
//					EdgeExp->SetUpPhysNormals(ElExp, edgeID);
//					normals	=  EdgeExp->GetMetricInfo()->GetNormal();
					//Vmath::Vdiv(nquad_e,normals[0],1,normals[1],1,nxny,1);

					for(i=0;i<nquad_e;i++)
					{
						// wx = 0
						m_freesurfacesplines[ispline]->GetSplineNormal(xedge[i],nx,ny);
						UBC[id1+i] = Wyphysedge[i] - Uphysedge[i]*nx/ny;
					}

				/*	for(i=0;i<nquad_e;i++)
					{
						cout << "vedge=" << UBC[id1+i] << endl;
					} */
				}

				Vmath::Vcopy(nq,&UBC[0],1,&(m_FBCu[0])[cnt2],1);

		        // Extrapolate to n+1
		        Vmath::Smul(nq,kFreeSurfaceBCsExtrapolation[nint-1][nint-1],
		        		&(m_FBCu[nint-1])[cnt2],1,&(m_FBCu[nlevels-1])[cnt2],1);
		        for(int k = 0; k < nint-1; ++k)
		        {
		            Vmath::Svtvp(nq,kFreeSurfaceBCsExtrapolation[nint-1][k],
		            		&(m_FBCu[k])[cnt2],1,&(m_FBCu[nlevels-1])[cnt2],1,
		            		&(m_FBCu[nlevels-1])[cnt2],1);
		        }

		        Vmath::Vcopy(nq,&(m_FBCu[nlevels-1])[cnt2], 1,&(BndExp[n]->UpdatePhys())[0],1);
				BndExp[n]->FwdTrans_BndConstrained(BndExp[n]->GetPhys(),BndExp[n]->UpdateCoeffs());
				ispline++;
				cnt2 +=nq;
			}
            else  // for all other boundary conditi
			{
				 cnt += BndExp[n]->GetExpSize();
			}
		}

	}

	void MovingMesh::AddwnutoFreeSurface(  Array<OneD, Array<OneD, NekDouble> > &Weakoutarray)
	{

	    int i,el,n,cnt, offset, phys_offset;
	    int edgeID, elmtid;
	    Array<OneD, NekDouble> e_outarray;
	    Array<OneD, const NekDouble> U,V,WX,WY;

	   // Array<OneD, Array<OneD, StdRegions::StdExpansion1DSharedPtr> >
	    //    elmtToTrace = m_traceMap->GetElmtToTrace();
	    StdRegions::StdExpansion1DSharedPtr EdgeExp;
	    StdRegions::StdExpansionSharedPtr ElExp;

	    Array<OneD, int> ElmtID,EdgeID;

	    //Fills ElmtID and EdgeID with global ids of elements id number and edges id number of the boundary expansion
	    m_fields[m_velocity[0]]->GetBoundaryToElmtMap(ElmtID,EdgeID);

	    int maxpts = 0;

	     // find the maximum values of points
	     for(cnt = n = 0; n < m_fields[m_velocity[0]]->GetBndConditions().num_elements(); ++n)
	     {
	 //Loop over all elements in boundary region
	     	//cnt starts with element 0 then going through elements
	         for(i = 0; i < m_fields[m_velocity[0]]->GetBndCondExpansions()[n]->GetExpSize(); ++i)
	         {
	             maxpts = max(maxpts, m_fields[m_velocity[0]]->GetExp(ElmtID[cnt++])->GetTotPoints());
	         }
	     }

		Array<OneD,NekDouble> dUdx(4*maxpts,0.0);
		Array<OneD,NekDouble> dUdy = dUdx + maxpts;
		Array<OneD,NekDouble> dVdx = dUdy + maxpts;
		Array<OneD,NekDouble> dVdy = dVdx + maxpts;


	    // Neumann BC for velocity component u
		for(cnt = n = 0; n < m_fields[m_velocity[0]]->GetBndConditions().num_elements(); ++n)
		{
			// Waters solution
			if(((m_fields[m_velocity[0]]->GetBndConditions()[n])->GetBoundaryConditionType() == SpatialDomains::eNeumann)
					&& (m_fields[m_velocity[0]]->GetBndConditions()[n]->GetUserDefined().GetEquation()!= "FreeSurface"))
			//(m_fields[m_velocity[0]]->GetBndConditions()[n]->GetUserDefined().GetEquation()== "Outflow")
			{
				// loop over elements along boundary
				for(el = 0; el < m_fields[m_velocity[0]]->GetBndCondExpansions()[n]->GetExpSize(); ++el,cnt++)
				{
					EdgeExp =  boost::dynamic_pointer_cast<StdRegions::StdExpansion1D> (m_fields[m_velocity[0]]->GetBndCondExpansions()[n]->GetExp(el));

					// grab edge on the boundary
					edgeID = EdgeID[cnt];
					elmtid = ElmtID[cnt];
					ElExp   = m_fields[m_velocity[0]]->GetExp(elmtid);

					int nqel = ElExp->GetTotPoints();
					int nquad_e = EdgeExp->GetNumPoints(0);

					U = (m_fields[m_velocity[0]]->GetPhys())+ m_fields[m_velocity[0]]->GetPhys_Offset(ElmtID[cnt]);
					V = (m_fields[m_velocity[1]]->GetPhys())+ m_fields[m_velocity[1]]->GetPhys_Offset(ElmtID[cnt]);
					WX = (m_fields[m_meshvelocity[0]]->GetPhys())+ m_fields[m_meshvelocity[0]]->GetPhys_Offset(ElmtID[cnt]);
					WY = (m_fields[m_meshvelocity[1]]->GetPhys())+ m_fields[m_meshvelocity[1]]->GetPhys_Offset(ElmtID[cnt]);

					Array<OneD,NekDouble> Uedge(nquad_e,0.0);
					Array<OneD,NekDouble> Vedge(nquad_e,0.0);
					Array<OneD,NekDouble> Wxedge(nquad_e,0.0);
					Array<OneD,NekDouble> Wyedge(nquad_e,0.0);
					Array<OneD,NekDouble> xedge(nquad_e,0.0);
					Array<OneD,NekDouble> yedge(nquad_e,0.0);
					Array<OneD,NekDouble> zedge(nquad_e,0.0);
					NekDouble nx,ny;

					EdgeExp->GetCoords(xedge,yedge,zedge);

					ElExp->GetEdgePhysVals(edgeID,EdgeExp,U,Uedge);
					ElExp->GetEdgePhysVals(edgeID,EdgeExp,V,Vedge);
					ElExp->GetEdgePhysVals(edgeID,EdgeExp,WX,Wxedge);
					ElExp->GetEdgePhysVals(edgeID,EdgeExp,WY,Wyedge);

					for(int i=0;i<nquad_e;i++)
					{
						// wx = 0
						m_freesurfacesplines[0]->GetSplineNormal(xedge[i],nx,ny);

					}

					/*normals	=  EdgeExp->GetMetricInfo()->GetNormal();
					NegateNormals = (ElExp->GetEorient(edgeID) == StdRegions::eBackwards)? false:true;

					if(NegateNormals == true)
					{
						Vmath::Neg(nquad_e,normals[0],1);
						Vmath::Neg(nquad_e,normals[1],1);
					}

					Vmath::Vmul(nquad_e,normals[0],1,dUdxedge,1, dUdxedge,1);
					Vmath::Vvtvp(nquad_e,normals[1],1,dVdxedge,1,dUdxedge,1,dUdxedge,1);*/

					e_outarray = Weakoutarray[0] + m_fields[m_velocity[0]]->GetCoeff_Offset(elmtid);

					// in: physical values of txx and txy, out:weak values of forcing function
				/*	m_fields[m_velocity[0]]->GetExp(elmtid)->AddEdgeNormBoundaryInt(edgeID,EdgeExp,dUdxedge,
				                                                dVdxedge,
				                                                e_outarray); */
				}
			}
			else
			{
				cnt +=m_fields[m_velocity[0]]->GetBndCondExpansions()[n]->GetExpSize();
			}


		}

		   // Neumann BC for velocity component v
			for(cnt = n = 0; n < m_fields[m_velocity[1]]->GetBndConditions().num_elements(); ++n)
			{
				if(((m_fields[m_velocity[1]]->GetBndConditions()[n])->GetBoundaryConditionType() == SpatialDomains::eNeumann)
						&& (m_fields[m_velocity[1]]->GetBndConditions()[n]->GetUserDefined().GetEquation()!= "FreeSurface"))
				//(m_fields[m_velocity[1]]->GetBndConditions()[n]->GetUserDefined().GetEquation()== "Outflow")
				{
					// loop over elements along boundary
					for(el = 0; el < m_fields[m_velocity[1]]->GetBndCondExpansions()[n]->GetExpSize(); ++el,cnt++)
					{
						EdgeExp =  boost::dynamic_pointer_cast<StdRegions::StdExpansion1D> (m_fields[m_velocity[1]]->GetBndCondExpansions()[n]->GetExp(el));

						// grab edge on the boundary
						edgeID = EdgeID[cnt];
						elmtid = ElmtID[cnt];
						ElExp   = m_fields[m_velocity[1]]->GetExp(elmtid);

						int nqel = ElExp->GetTotPoints();
						int nquad_e = EdgeExp->GetNumPoints(0);

						U = (m_fields[m_velocity[0]]->GetPhys())+ m_fields[m_velocity[0]]->GetPhys_Offset(ElmtID[cnt]);
						V = (m_fields[m_velocity[1]]->GetPhys())+ m_fields[m_velocity[1]]->GetPhys_Offset(ElmtID[cnt]);

						// Calculating deformation tensor
						ElExp->PhysDeriv(1,U,dUdy);
						ElExp->PhysDeriv(1,V,dVdy);

						Array<OneD,NekDouble> dUdyedge(nquad_e,0.0);
						Array<OneD,NekDouble> dVdyedge(nquad_e,0.0);

						ElExp->GetEdgePhysVals(edgeID,EdgeExp,dUdy,dUdyedge);
						ElExp->GetEdgePhysVals(edgeID,EdgeExp,dVdy,dVdyedge);

						e_outarray = Weakoutarray[1] + m_fields[m_velocity[1]]->GetCoeff_Offset(elmtid);

						// in: physical values of txx and txy, out:weak values of forcing function
						m_fields[m_velocity[1]]->GetExp(elmtid)->AddEdgeNormBoundaryInt(edgeID,EdgeExp,dUdyedge,
					                                                dVdyedge,
					                                                e_outarray);
					}
				}
				else
				{
					cnt +=m_fields[m_velocity[1]]->GetBndCondExpansions()[n]->GetExpSize();
				}

			}
	}

	void MovingMesh::SetEdgePhysValsToZero(int edge, Array<OneD, NekDouble> &inarray, StdRegions::StdExpansionSharedPtr &ElExp)
	  {

		  int nquad0 = ElExp->GetBase()[0]->GetNumPoints();
	      int nquad1 = ElExp->GetBase()[1]->GetNumPoints();
		Array<OneD,NekDouble> e_tmp;

		StdRegions::EdgeOrientation edgedir = ElExp->GetEorient(edge);
		switch(edge)
		{
		case 0:
			if(edgedir == StdRegions::eForwards)
			{
				Vmath::Zero(nquad0,inarray,1);
			}
			else
			{
				e_tmp = inarray+(nquad0-1);
				Vmath::Zero(nquad0,e_tmp,-1);
			}

			break;
		case 1:
			if(edgedir == StdRegions::eForwards)
			{
				e_tmp = inarray+(nquad0-1);
				Vmath::Zero(nquad1,e_tmp,nquad0);
			}
			else
			{
				e_tmp = inarray+(nquad0*nquad1-1);
				Vmath::Zero(nquad1,e_tmp,-nquad0);
			}
			break;
		case 2:
			if(edgedir == StdRegions::eForwards)
			{
				e_tmp = inarray+(nquad0*nquad1-1);
				Vmath::Zero(nquad0,e_tmp,-1);
			}
			else
			{
				e_tmp = inarray+nquad0*(nquad1-1);
				Vmath::Zero(nquad0,e_tmp,1);
			}
			break;
		case 3:
			if(edgedir == StdRegions::eForwards)
			{
				e_tmp = inarray + nquad0*(nquad1-1);
				Vmath::Zero(nquad1,e_tmp,-nquad0);
			}
			else
			{
				Vmath::Zero(nquad1,inarray,nquad0);
			}
			break;
		default:
			ASSERTL0(false,"edge value (< 3) is out of range");
			break;
		}
	}


void MovingMesh::SetValueAlongFreeSurfaceToZero(  Array<OneD, Array<OneD, NekDouble> > &inarray)
    {
    	Array<OneD, const SpatialDomains::BoundaryConditionShPtr > BndConds;
		Array<OneD, MultiRegions::ExpListSharedPtr>  BndExp;
		//static int outputcnt2 = 1;
		BndConds   = m_fields[m_meshvelocity[0]]->GetBndConditions();
		BndExp     = m_fields[m_meshvelocity[0]]->GetBndCondExpansions();

		StdRegions::StdExpansionSharedPtr ElExp;
		StdRegions::StdExpansion1DSharedPtr EdgeExp;
		Array<OneD, int> ElmtID,EdgeID;
		int i,cnt, n,el,edgeID,elmID;
		int ninarray = inarray.num_elements();
		Array<OneD, Array<OneD, NekDouble> > inarrayelement(ninarray);
		int maxpts =0;
		NekDouble nx,ny;

		for(int i=0;i<ninarray;i++)
		{
			inarrayelement[i] = Array<OneD, NekDouble>(m_phystot);
		}

		//Fills ElmtID and EdgeID with global ids of elements id number and edges id number of the boundary expansion
		m_fields[m_meshvelocity[0]]->GetBoundaryToElmtMap(ElmtID,EdgeID);
		int ispline=0;
		for(cnt = n = 0; n < m_fields[m_meshvelocity[0]]->GetBndConditions().num_elements(); ++n)
		{
			 string type = m_fields[m_meshvelocity[0]]->GetBndConditions()[n]->GetUserDefined().GetEquation();
	 //for the ones marked as Drag calculate drag along boundary edge
			 if(type == "FreeSurface")
			 {
				for(el = 0; el < BndExp[n]->GetExpSize(); ++el,cnt++)
				{
					int offset = m_fields[0]->GetPhys_Offset(ElmtID[cnt]);
					for(int i=0;i<ninarray;i++)
					{
						inarrayelement[i] = inarray[i]+offset;
					}

					ElExp = m_fields[m_meshvelocity[0]]->GetExp(ElmtID[cnt]);
					EdgeExp =  boost::dynamic_pointer_cast<StdRegions::StdExpansion1D> (BndExp[n]->GetExp(el));

					int nquad_e = EdgeExp->GetNumPoints(0);
					Array<OneD,NekDouble> xedge(nquad_e,0.0);
					Array<OneD,NekDouble> yedge(nquad_e,0.0);
					Array<OneD,NekDouble> zedge(nquad_e,0.0);
					Array<OneD, Array<OneD, NekDouble> > inarrayedge(ninarray);

					for(int i=0;i<ninarray;i++)
					{
						inarrayedge[i] = Array<OneD,NekDouble> (nquad_e,0.0);
					}

					EdgeExp->GetCoords(xedge,yedge,zedge);
					edgeID = EdgeID[cnt];
					for(int i=0;i<ninarray;i++)
					{
						//ElExp->GetEdgePhysVals(edgeID,EdgeExp,inarrayelement[i],inarrayedge[i]);
						SetEdgePhysValsToZero(edgeID, inarrayelement[i], ElExp);
					}

					//EdgeExp->GetCoords(xedge,yedge,zedge);


				}
				ispline++;
			 }
			 else  // for all other boundary conditi
			 {
				 cnt += BndExp[n]->GetExpSize();
			 }
		 }

    }

void MovingMesh::SetOutflowBCALE()
	{

	    int i,el,n,cnt, offset, phys_offset;
	    int edgeID, elmtid;
	    Array<OneD, NekDouble> e_outarray;
	    Array<OneD, const NekDouble> U;
	    Array<OneD, Array<OneD, NekDouble> > normals;

	    Array<OneD, MultiRegions::ExpListSharedPtr>  BndExp;
	    Array<OneD, const SpatialDomains::BoundaryConditionShPtr > BndConds;
		BndConds   = m_fields[m_meshvelocity[0]]->GetBndConditions();
		BndExp     = m_fields[m_meshvelocity[0]]->GetBndCondExpansions();

	   // Array<OneD, Array<OneD, StdRegions::StdExpansion1DSharedPtr> >
	    //    elmtToTrace = m_traceMap->GetElmtToTrace();
	    StdRegions::StdExpansion1DSharedPtr EdgeExp;
	    StdRegions::StdExpansionSharedPtr ElExp;

	    Array<OneD, int> ElmtID,EdgeID;

	    //Fills ElmtID and EdgeID with global ids of elements id number and edges id number of the boundary expansion
	    m_fields[m_velocity[0]]->GetBoundaryToElmtMap(ElmtID,EdgeID);

	    // Neumann BC for velocity component u
		for(cnt = n = 0; n < m_fields[m_velocity[0]]->GetBndConditions().num_elements(); ++n)
		{
			// Waters solution
			if(((m_fields[m_velocity[0]]->GetBndConditions()[n])->GetBoundaryConditionType() == SpatialDomains::eNeumann)
					&& (m_fields[m_velocity[0]]->GetBndConditions()[n]->GetUserDefined().GetEquation()== "Outflow"))
			{
				// loop over elements along boundary
				for(el = 0; el < m_fields[m_velocity[0]]->GetBndCondExpansions()[n]->GetExpSize(); ++el,cnt++)
				{
					EdgeExp =  boost::dynamic_pointer_cast<StdRegions::StdExpansion1D> (m_fields[m_velocity[0]]->GetBndCondExpansions()[n]->GetExp(el));

					// grab edge on the boundary
					edgeID = EdgeID[cnt];
					elmtid = ElmtID[cnt];
					ElExp   = m_fields[m_velocity[0]]->GetExp(elmtid);

					int nqel = ElExp->GetTotPoints();
					int nquad_e = EdgeExp->GetNumPoints(0);

					U = (m_fields[m_velocity[0]]->GetPhys())+ m_fields[m_velocity[0]]->GetPhys_Offset(ElmtID[cnt]);

					Array<OneD,NekDouble> Uedge(nquad_e,0.0);
					Array<OneD,NekDouble> wxedge(nquad_e,0.0);

					ElExp->GetEdgePhysVals(edgeID,EdgeExp,U,Uedge);

					for(int i=0;i<nquad_e;i++)
					{
						wxedge[i]= Uedge[i];
					}

					// Copy result in boundary expansion
					int id1  = BndExp[n]->GetPhys_Offset(el);
					Vmath::Vcopy(nquad_e,&wxedge[0], 1,&(BndExp[n]->UpdatePhys())[id1],1);
				}
				BndExp[n]->FwdTrans_BndConstrained(BndExp[n]->GetPhys(),BndExp[n]->UpdateCoeffs());
			}
			else
			{
				cnt +=m_fields[m_velocity[0]]->GetBndCondExpansions()[n]->GetExpSize();
			}

		}
	}

void MovingMesh::SetPlateBCALE()
	{

	    int i,el,n,cnt, offset, phys_offset;
	    int edgeID, elmtid;
	    Array<OneD, NekDouble> e_outarray;
	    Array<OneD, const NekDouble> U;
	    Array<OneD, Array<OneD, NekDouble> > normals;

	    Array<OneD, MultiRegions::ExpListSharedPtr>  BndExp;
	    Array<OneD, const SpatialDomains::BoundaryConditionShPtr > BndConds;
		BndConds   = m_fields[m_meshvelocity[0]]->GetBndConditions();
		BndExp     = m_fields[m_meshvelocity[0]]->GetBndCondExpansions();
		NekDouble PlateVelocity = 1;

	   // Array<OneD, Array<OneD, StdRegions::StdExpansion1DSharedPtr> >
	    //    elmtToTrace = m_traceMap->GetElmtToTrace();
	    StdRegions::StdExpansion1DSharedPtr EdgeExp;
	    StdRegions::StdExpansionSharedPtr ElExp;

	    Array<OneD, int> ElmtID,EdgeID;

	    //Fills ElmtID and EdgeID with global ids of elements id number and edges id number of the boundary expansion
	    m_fields[m_meshvelocity[0]]->GetBoundaryToElmtMap(ElmtID,EdgeID);

	    // Neumann BC for velocity component u
		for(cnt = n = 0; n < m_fields[m_meshvelocity[0]]->GetBndConditions().num_elements(); ++n)
		{
			// Waters solution
			if((m_fields[m_meshvelocity[0]]->GetBndConditions()[n]->GetUserDefined().GetEquation()== "Plate"))
			{
				// loop over elements along boundary
				for(el = 0; el < m_fields[m_meshvelocity[0]]->GetBndCondExpansions()[n]->GetExpSize(); ++el,cnt++)
				{
					EdgeExp =  boost::dynamic_pointer_cast<StdRegions::StdExpansion1D> (m_fields[m_meshvelocity[0]]->GetBndCondExpansions()[n]->GetExp(el));

					// grab edge on the boundary
					edgeID = EdgeID[cnt];
					elmtid = ElmtID[cnt];
					ElExp   = m_fields[m_meshvelocity[0]]->GetExp(elmtid);

					int nqel = ElExp->GetTotPoints();
					int nquad_e = EdgeExp->GetNumPoints(0);

					U = (m_fields[m_velocity[0]]->GetPhys())+ m_fields[m_velocity[0]]->GetPhys_Offset(ElmtID[cnt]);

					Array<OneD,NekDouble> Uedge(nquad_e,0.0);
					Array<OneD,NekDouble> wxedge(nquad_e,0.0);

					ElExp->GetEdgePhysVals(edgeID,EdgeExp,U,Uedge);


					for(int i=0;i<nquad_e;i++)
					{
						wxedge[i]= Uedge[i];
					}

					// Copy result in boundary expansion
					int id1  = BndExp[n]->GetPhys_Offset(el);
					Vmath::Vcopy(nquad_e,&wxedge[0], 1,&(BndExp[n]->UpdatePhys())[id1],1);
					//cout << "Set Plate velocity" << endl;
				}
				BndExp[n]->FwdTrans_BndConstrained(BndExp[n]->GetPhys(),BndExp[n]->UpdateCoeffs());
			}
			else
			{
				cnt +=m_fields[m_meshvelocity[0]]->GetBndCondExpansions()[n]->GetExpSize();
			}

		}
	}
/*void ALE::WriteSwellRatioToFileFromSpline()
    {
    	int phystot = m_freesurfacesplines[0]->GetTotPoints();
    	Array<OneD, Array<OneD, NekDouble > > coords(2);
    	stringstream filenamestream;

    	filenamestream <<  "SwellRatio.txt";

    	string filename = filenamestream.str();

    	for(int i=0;i<2;i++)
    	{
    		coords[i] = Array<OneD, NekDouble >(phystot,0.0);
    	}

    	coords = m_freesurfacesplines[0]->GetCoords();

    	// Calculate Swell Ratio
    	//Determine quadrature point with Minimum x and maximum y
    	NekDouble xmin=0,ymaxtoxmin=0;
    	NekDouble xmaxtoymax=0,ymax=0;
    	NekDouble xmaxexit=0,ymaxtoxmaxexit=0;
    	int idxminymax, idxmaxymax;
    	xmin = coords[0][0];
    	xmaxexit = coords[0][phystot-1];
    	ymaxtoxmin = coords[1][0];
    	ymaxtoxmaxexit = coords[0][phystot-1];

    	for(int i=0;i<phystot;i++)
    	{
    		if(coords[1][i]>=ymax)
    		{
    			ymax = coords[1][i];
    		}
    	}

    	cout << "Height at Entry: at x=" << xmin <<"y= "<< ymaxtoxmin << endl;
    	cout << "Height at Maximum: at x="  << xmax <<"y= "<< ymaxtoxmax << endl;
    	cout << "Height at Exit: at x= "  << xmaxexit <<", y= "<< ymaxtoxmaxexit << endl;

    	//Compute Swell Ratio
    	NekDouble ratio = ymaxtoxmax/ymaxtoxmin;
    	NekDouble ratioexit = ymaxtoxmaxexit/ymaxtoxmin;

    	// Check if file is empty
    	int length;
    	ifstream filestr;

    	filestr.open(filename.c_str(), ios::binary); // open your file
    	filestr.seekg(0, ios::end); // put the "cursor" at the end of the file
    	length = filestr.tellg(); // find the position of the cursor
    	filestr.close(); // close your file

    	if ( length == 0 )
    	{
    		// Create file and write header
    		ofstream outfile(filename.c_str());
        	outfile.width(10);
        	outfile.precision(8);
        	outfile << "Re" << "\t";
        	outfile.width(10);
        	outfile.precision(8);
           outfile << "Swell" << "\t";
           cout << "Creating " << filename << endl;

         	cout << "Writing Swell ratio into " << filename << endl;

         	outfile.width(10);
         	outfile.precision(8);
         	outfile << m_reynolds << "\t";
         	outfile.width(10);
         	outfile.precision(8);
           outfile << ratio << "\t";
        	outfile.width(10);
        	outfile.precision(8);
           outfile << ratioexit << "\t";
           outfile << "\n";

           outfile.close();
    	}
    	else
    	{
    		//Open File to append data
    		ofstream outfile(filename.c_str(),ios::app);

        	cout << "Writing Swell ratio into " << filename << endl;

        	outfile.width(10);
        	outfile.precision(8);
        	outfile << m_reynolds << "\t";
        	outfile.width(10);
        	outfile.precision(8);
           outfile << ratio << "\t";
        	outfile.width(10);
        	outfile.precision(8);
           outfile << ratioexit << "\t";
           outfile << "\n";

           outfile.close();
    	}

    } */



	  } //end of namespace






