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

#include <LocalRegions/Expansion2D.h>
#include <MultiRegions/ContField1D.h>

#include <IncNavierStokesSolver/EquationSystems/ALE.h>
#include <cstdio>
#include <cstdlib>

#include <cmath>

#include <string>
namespace Nektar
{
    /**
     * Constructor. Creates ...
     *
     * \param
     * \param
     */

    ALE::ALE(LibUtilities::SessionReaderSharedPtr &pSession):
            EquationSystem(pSession)
    {
          int i,j;
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
                  std::string var = m_boundaryConditions->GetVariable(j);
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
        m_ExpField = MemoryManager<MultiRegions::ExpList2D>::AllocateSharedPtr(m_comm,*m_mesh2D);

        int intSteps = 1;
        m_meshvelwx =  Array<OneD, Array<OneD, NekDouble> >  (intSteps);
        m_meshvelwy =  Array<OneD, Array<OneD, NekDouble> >  (intSteps);
        for(int i=0;i<intSteps;i++)
        {
             m_meshvelwx[i] = Array<OneD, NekDouble>(m_phystot);
             m_meshvelwy[i] = Array<OneD, NekDouble>(m_phystot);
        }

        //Setting up free surface spline
        SetUpFreeSurfaceSplines();


    }


    ALE::~ALE()
    {
    }

    void ALE::SetUpFreeSurfaceSplines()
    {
        int cnt,n;
        Array<OneD, MultiRegions::ExpListSharedPtr>  BndExp     = m_fields[m_meshvelocity[0]]->GetBndCondExpansions();
        Array<OneD, const SpatialDomains::BoundaryConditionShPtr > BndConds = m_fields[m_meshvelocity[0]]->GetBndConditions();

        // Count number of free surface regions
        for(cnt = n = 0; n < BndConds.num_elements(); ++n)
         {
             string type = BndConds[n]->GetUserDefined();
     //for the ones marked as Drag calculate drag along boundary edge
             if(type == "FreeSurface")
             {
                 cnt++;
             }
         }

        m_freesurfacesplines = Array<OneD, LibUtilities::CubicSplineSharedPtr>(cnt);
        LibUtilities::SplineBoundaryType bcleft,bcright;
        //bcleft = LibUtilities::eClampedQuadraticLagrange;
        bcleft = LibUtilities::eNotaKnot;
        bcright = LibUtilities::eNotaKnot;
        int bcl=0,bcr=0;

        for(cnt = n = 0; n < BndConds.num_elements(); ++n)
       {
           string type = BndConds[n]->GetUserDefined();
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
/*    void ALE::GenerateBoundaryConditionsMeshVelocity(Array<OneD, MultiRegions::ExpListSharedPtr > m_fields)
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

    void ALE::DoALE(const NekDouble aii_Dt, NekDouble time)
    {
        CreateFreeSurfaceSpline(aii_Dt, time);
        Array<OneD,NekDouble> forcing(m_phystot,0.0);
        StdRegions::ConstFactorMap factors;
        factors[StdRegions::eFactorLambda] = 0.0;

        // At this place we could solve a poisson equation with RHS=0
        for(int i=0;i<m_spacedim;i++)
        {
            m_fields[m_meshvelocity[i]]->HelmSolve(forcing,m_fields[m_meshvelocity[i]]->UpdateCoeffs(),NullFlagList,factors);
            m_fields[m_meshvelocity[i]]->BwdTrans(m_fields[m_meshvelocity[i]]->GetCoeffs(),m_fields[m_meshvelocity[i]]->UpdatePhys());
        }

        DeformMesh2(aii_Dt);
        DeformMeshBoundary(aii_Dt);

        // Updates Mesh, Fields and Coordinates of Mesh
        UpdateMesh();
        // Compute Mesh Velocity from Coordinates of Mesh
        ComputeMeshVelocity(aii_Dt);
    }



/*
 * Updates m_mesh2D with new mesh information
 */
    void ALE::SetMeshVelocityatBoundary(const NekDouble aii_Dt)
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
        SpatialDomains::PointGeomSharedPtr vertex1,vertex2;
        SpatialDomains::MeshGraph2DSharedPtr mesh2D;

        MultiRegions::ExpList1DSharedPtr  FreeSurfaceFct;

        Array<OneD, int> ElmtID,EdgeID;
        int cnt,n,el,edgeID;
        Array<OneD, const NekDouble> Uphyselement, Vphyselement, wxel,wyel;

        //Fills ElmtID and EdgeID with global ids of elements id number and edges id number of the boundary expansion
        m_fields[m_meshvelocity[0]]->GetBoundaryToElmtMap(ElmtID,EdgeID);

        int cnt2=0;
        for(cnt = n = 0; n < m_fields[0]->GetBndConditions().num_elements(); ++n)
         {
             string type = m_fields[0]->GetBndConditions()[n]->GetUserDefined();
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
                   EdgeExp->SetUpPhysNormals(EdgeID[cnt2]);
                   //locnormals    =  EdgeExp->GetMetricInfo()->GetNormal();
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

                if(!(FreeSurfaceFct =  boost::dynamic_pointer_cast<MultiRegions::ExpList1D>(BndExp[n])))
                    {
                        ASSERTL0(false,"Dynamics cast failed");
                    }
                     MultiRegions::ContField1DSharedPtr contfield =
                                                               MemoryManager<MultiRegions::ContField1D>
                                                               ::AllocateSharedPtr(m_session,*FreeSurfaceFct);
                     contfield->FwdTrans(normals[0],contfield->UpdateCoeffs());
                     contfield->BwdTrans(contfield->GetCoeffs(),normals[0]);
                     contfield->FwdTrans(normals[1],contfield->UpdateCoeffs());
                     contfield->BwdTrans(contfield->GetCoeffs(),normals[1]);

                 //loop over all elements along the boundary region
                for(el = 0; el < BndExp[n]->GetExpSize(); ++el,cnt++)
                {
                    Uphyselement = (m_fields[0]->GetPhys())+ m_fields[0]->GetPhys_Offset(ElmtID[cnt]);
                    Vphyselement = (m_fields[1]->GetPhys())+ m_fields[1]->GetPhys_Offset(ElmtID[cnt]);
                    wxel = (m_fields[m_meshvelocity[0]]->GetPhys())+ m_fields[m_meshvelocity[0]]->GetPhys_Offset(ElmtID[cnt]);
                    wyel = (m_fields[m_meshvelocity[1]]->GetPhys())+ m_fields[m_meshvelocity[1]]->GetPhys_Offset(ElmtID[cnt]);

                    ElExp = m_fields[m_meshvelocity[0]]->GetExp(ElmtID[cnt]);
                /*    int nq = locExp->GetTotPoints();
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

                    EdgeExp->SetUpPhysNormals(edgeID);
                    //EdgeExp->SetUpPhysNormals(ElExp,edgeID);

                    //locnormals    =  EdgeExp->GetMetricInfo()->GetNormal();
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

                /*    cout<< "After smoothing" << endl;
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

                /*    if(first)
                    { */

                    /*    for(int i=0;i<nquad_e;i++)
                        {
                            //if(xedge[i]==11.0 && yedge[i]<=0.0)
                                                            {
                            cout << i << ": (" << xedge[i] << "," << yedge[i] << ")"  << "\t";
                            cout << "(n1,n2) = (" << locnormals[0][i] << ", " << locnormals[1][i] << ")" << endl;
                            cout << "Edge: (u,v)=(" << Uphysedge[i] <<","<< Vphysedge[i] << ")"<< "\t";
                            cout << "(wx,wy)=(" << wxedge[i] <<","<< wyedge[i] << ")"<< endl;
                                                            }
                        } */
                //    }


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

    void ALE::DeformMesh2(const NekDouble aii_Dt)
   {
    // number of elements
    int  nlevels = m_meshvelwx.num_elements();
    NekDouble x,y,z;
    int eid, el,j;
    int edgeId;
    StdRegions::StdExpansionSharedPtr ElExp;
       SpatialDomains::SegGeomSharedPtr edge;
       StdRegions::StdExpansion1DSharedPtr EdgeExp;
    Array<OneD, Array<OneD, NekDouble> > normals, locnormals,weaknormals;
    Array<OneD, const NekDouble> wxel,wyel;
    int vertID1;
     SpatialDomains::PointGeomSharedPtr vertex1,vertex2;
     SpatialDomains::Geometry1DSharedPtr ElmtSegGeom;

     // local expansion vector
     const LocalRegions::ExpansionVector &locExpVector = *(m_ExpField->GetExp());

     int nel = locExpVector.size();

       // Loop over elemente and collect forward expansion
    int nquad_e;
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
               edgeId = (locExpVector[eid]->as<LocalRegions::Expansion2D>())->GetGeom2D()->GetEid(j);
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

                   x = x + (aii_Dt)*wxedge[0];
                   y = y + (aii_Dt)*wyedge[0];

                   vertex1->UpdatePosition(x,y,z);
               }

           }
       }
   }

    void ALE::DeformMeshBoundary(const NekDouble aii_Dt)
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
                 SpatialDomains::PointGeomSharedPtr vertex1,vertex2;
                 SpatialDomains::MeshGraph2DSharedPtr mesh2D;

                Array<OneD, int> ElmtID,EdgeID;
                int cnt,n,el,edgeID;
                Array<OneD, const NekDouble> Uphyselement, Vphyselement, wxel,wyel;

                //Fills ElmtID and EdgeID with global ids of elements id number and edges id number of the boundary expansion
                m_fields[m_meshvelocity[0]]->GetBoundaryToElmtMap(ElmtID,EdgeID);

                for(cnt = n = 0; n < m_fields[m_meshvelocity[0]]->GetBndConditions().num_elements(); ++n)
                      {
                          string type = m_fields[m_meshvelocity[0]]->GetBndConditions()[n]->GetUserDefined();
                  //for the ones marked as Drag calculate drag along boundary edge
                          if(type == "FreeSurface")
                          {
                              //loop over all elements along the boundary region
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

                                 for(int i=0;i<nquad_e;i++)
                                 {
                                     xedge[i]= xedge[i] + (aii_Dt)*wxedge[i];
                                     yedge[i]= yedge[i] + (aii_Dt)*wyedge[i];
                                 }


                                 edgeID= m_mesh2D->GetEidFromElmt(ElExp->DetShapeType(),
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

    void ALE::ComputeNewCoordsatFreeSurfaceBoundary(Array<OneD, Array<OneD, NekDouble> > &newcoords, const NekDouble aii_Dt)
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
         SpatialDomains::PointGeomSharedPtr vertex1,vertex2;
         SpatialDomains::MeshGraph2DSharedPtr mesh2D;

        Array<OneD, int> ElmtID,EdgeID;
        int cnt,n,el,edgeID;
        Array<OneD, const NekDouble> Uphyselement, Vphyselement, wxel,wyel;

        //Fills ElmtID and EdgeID with global ids of elements id number and edges id number of the boundary expansion
        m_fields[m_meshvelocity[0]]->GetBoundaryToElmtMap(ElmtID,EdgeID);

        for(cnt = n = 0; n < m_fields[m_meshvelocity[0]]->GetBndConditions().num_elements(); ++n)
          {
              string type = m_fields[m_meshvelocity[0]]->GetBndConditions()[n]->GetUserDefined();
      //for the ones marked as Drag calculate drag along boundary edge
              if(type == "FreeSurface")
              {
                 //loop over all elements along the boundary region
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
                        newcoords[0][id1+i]= xedge[i] + (aii_Dt)*wxedge[i];
                        newcoords[1][id1+i]= yedge[i] + (aii_Dt)*wyedge[i];
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

    void ALE::UpdateMeshFreeSurfaceBoundary(const NekDouble aii_Dt)
    {
        int cnt,n,i;
        int el,edgeID;
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
        SpatialDomains::PointGeomSharedPtr vertex1,vertex2;
        SpatialDomains::MeshGraph2DSharedPtr mesh2D;

        Array<OneD, int> ElmtID,EdgeID;

        m_fields[m_meshvelocity[0]]->GetBoundaryToElmtMap(ElmtID,EdgeID);

        // For every free surface region compute gll points along new free surface boundary
        for(cnt = n = 0; n < BndConds.num_elements(); ++n)
         {
             string type = BndConds[n]->GetUserDefined();
     //for the ones marked as Drag calculate drag along boundary edge
             if(type == "FreeSurface")
             {
                 int nq = BndExp[n]->GetTotPoints();
                 for(i=0;i<newcoords.num_elements();i++)
                 {
                     newcoords[i] = Array<OneD, NekDouble>(nq,0.0);
                 }
                 /* Compute new coordinates of free surface boundary */
                 ComputeNewCoordsatFreeSurfaceBoundary(newcoords,aii_Dt);
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

                    /*    cout << "Coords before= " << endl;
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
                        edgeID= m_mesh2D->GetEidFromElmt(ElExp->DetShapeType(),
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

    void ALE::UpdateFields()
    {
        // velocities plus pressure
        int nvariables=m_fields.num_elements();
        int i;

        Array<OneD, Array<OneD, NekDouble> >   fieldsphys(nvariables);
        Array<OneD, Array<OneD, NekDouble> >   boundaryphys(nvariables);

        for(i = 0 ; i < m_fields.num_elements(); i++)
        {
            MultiRegions::ContField2DSharedPtr field;

            if(!(field = boost::dynamic_pointer_cast<
                    MultiRegions::ContField2D>(m_fields[i])))
                {
                    ASSERTL0(false,"Dynamics cast failed");
                }

            field->UpdateContField2D(m_comm,*m_mesh2D,*m_boundaryConditions,
                 m_boundaryConditions->GetVariable(i),
                                 false);
        }
    }

    void ALE::UpdateMesh()
    {
        m_mesh2D = m_mesh2D->UpdateMesh(m_graph);
        SpatialDomains::MeshGraph *meshptr2 = m_mesh2D.get();

        m_boundaryConditions = m_boundaryConditions->UpdateBoundaryConditions(meshptr2,*m_boundaryConditions);

        if(!(m_graph = boost::dynamic_pointer_cast<
                               SpatialDomains::MeshGraph>(m_mesh2D)))
        {
          ASSERTL0(false,"Dynamics cast failed");
        }
         // Dummy to deal with mesh properties
        /*  m_DisField->UpdateDisContField2D(m_comm,*m_mesh2D,*m_boundaryConditions,
                         m_boundaryConditions->GetVariable(0),
                                         m_solnType,
                                         true); */
        m_ExpField->UpdateExpList2D(m_comm,*m_mesh2D,m_boundaryConditions->GetVariable(m_meshvelocity[0]));

        for(int i = 0; i < 3; ++i)
        {
            Vmath::Vcopy(m_phystot,m_meshcoords[i],1,m_meshcoordsold[i],1);
        }

        // new coordinates
        // TODO: we need a field that saves the new coords
        //m_DisField->GetCoords(m_meshcoords[0],m_meshcoords[1],m_meshcoords[2]);
        m_ExpField->GetCoords(m_meshcoords[0],m_meshcoords[1],m_meshcoords[2]);
    }

    void ALE::ComputeMeshVelocity(const NekDouble aii_Dt)
    {
        int VelDim = m_spacedim;
        for(int i=0;i<VelDim;i++)
        {
            for(int j=0;j<m_phystot;j++)
            {
                m_fields[m_meshvelocity[i]]->UpdatePhys()[j]=(m_meshcoords[i][j]-m_meshcoordsold[i][j])/aii_Dt;
            }
        }
    }

    void ALE::ComputeVelocityMinusMeshVelocity(Array<OneD, Array<OneD, NekDouble> > &velocity)
    {
        int VelDim=velocity.num_elements();

        for(int i=0;i<VelDim;i++)
        {
            Vmath::Vsub(m_phystot,velocity[i],1,m_fields[m_meshvelocity[i]]->GetPhys(),1,velocity[i],1);

        }
    }

    void ALE::CreateFreeSurfaceSpline(const NekDouble aii_Dt, NekDouble time)
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
         int cnt, n,el,edgeID;
         Array<OneD, const NekDouble> Uphyselement, Vphyselement;

         //Fills ElmtID and EdgeID with global ids of elements id number and edges id number of the boundary expansion
         m_fields[m_meshvelocity[0]]->GetBoundaryToElmtMap(ElmtID,EdgeID);

         NekDouble xmax = 0.0;
         int l;
         l =0;
         // Determine xmax because its outflow value that needs to be set to value from before
        for(cnt = n = 0; n < m_fields[m_meshvelocity[0]]->GetBndConditions().num_elements(); ++n)
        {
             string type = m_fields[m_meshvelocity[0]]->GetBndConditions()[n]->GetUserDefined();
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
             string type = m_fields[m_meshvelocity[0]]->GetBndConditions()[n]->GetUserDefined();
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
                        if(xedge[i]==0)
                        {
                            wyedge[i]=0.0;
                        }

                    }

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

    void ALE::SetFreeSurfaceVelocityBC()
    {
    Array<OneD, const SpatialDomains::BoundaryConditionShPtr > BndConds;
    Array<OneD, MultiRegions::ExpListSharedPtr>  BndExp, BndExpv;
    MultiRegions::ExpList1DSharedPtr  FreeSurfaceFct;
    Array<OneD, Array<OneD, NekDouble> > normals, locnormals;
    Array<OneD, NekDouble> e_outarrayu, e_outarrayv;
    Array<OneD, NekDouble> curv;

    BndConds   = m_fields[0]->GetBndConditions();
    BndExp     = m_fields[0]->GetBndCondExpansions();
    BndExpv    = m_fields[1]->GetBndCondExpansions();

    StdRegions::StdExpansionSharedPtr ElExp;
    StdRegions::StdExpansion1DSharedPtr EdgeExp,EdgeExpv;
    Array<OneD, int> ElmtID,EdgeID;
    int cnt,n,el,edgeID;
    Array<OneD, const NekDouble> U,V;
    Array<OneD, NekDouble> Uvals,Vvals;


    //Fills ElmtID and EdgeID with global ids of elements id number and edges id number of the boundary expansion
    m_fields[0]->GetBoundaryToElmtMap(ElmtID,EdgeID);

    for(cnt = n = 0; n < m_fields[0]->GetBndConditions().num_elements(); ++n)
    {
         string type = m_fields[0]->GetBndConditions()[n]->GetUserDefined();
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

                    //normals    =  EdgeExp->GetMetricInfo()->GetNormal();

                    // Multiply with outward normals
                    //if(ElExp->GetEorient(edgeID) == StdRegions::eBackwards)
                    //{
                    //     Vmath::Neg(nquad_e,normals[0],1);
                    //     Vmath::Neg(nquad_e,normals[1],1);
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

    void ALE::WriteArrayAlongFreeSurfaceToFile(string fname, int outputcount,Array<OneD, Array<OneD, NekDouble> > &inarray)
       {
           Array<OneD, const SpatialDomains::BoundaryConditionShPtr > BndConds;
           Array<OneD, MultiRegions::ExpListSharedPtr>  BndExp;
           //static int outputcnt2 = 1;
           BndConds   = m_fields[m_meshvelocity[0]]->GetBndConditions();
           BndExp     = m_fields[m_meshvelocity[0]]->GetBndCondExpansions();

           StdRegions::StdExpansionSharedPtr ElExp;
           StdRegions::StdExpansion1DSharedPtr EdgeExp;
           Array<OneD, int> ElmtID,EdgeID;
           int cnt, n,el,edgeID;
           int ninarray = inarray.num_elements();
           Array<OneD, Array<OneD, const NekDouble> > inarrayelement(ninarray);

           for(int i=0;i<ninarray;i++)
           {
               inarrayelement[i] = Array<OneD, const NekDouble>(m_phystot);
           }

           //Fills ElmtID and EdgeID with global ids of elements id number and edges id number of the boundary expansion
           m_fields[m_meshvelocity[0]]->GetBoundaryToElmtMap(ElmtID,EdgeID);
           int ispline=0;
           for(cnt = n = 0; n < m_fields[m_meshvelocity[0]]->GetBndConditions().num_elements(); ++n)
           {
                string type = m_fields[m_meshvelocity[0]]->GetBndConditions()[n]->GetUserDefined();
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

NekDouble ALE::LagrangeInterpolant(NekDouble x, int npts, const Array<OneD, const NekDouble>& xpts,
           const Array<OneD, const NekDouble>& funcvals)
{
   NekDouble sum = 0.0;

   NekDouble dif=fabs(x-xpts[0]);
   NekDouble dift;
   int index=0;

   // find index of closest x-value on centre to x
   for(int i=0;i<npts;++i)
   {
       if ( (dift=fabs(x-xpts[i])) < dif)
       {
           index=i;
           dif=dift;
       }
   }

//   cout << "x=" << x << "id=" << index << "x_c=" << xpts[index] << "u_c=" << funcvals[index] << endl;

  /* Array<OneD, NekDouble> uc2(npts2);
   Array<OneD, NekDouble> xc2(npts2);
   int id2=0;

   if(((index-npts2/2)<=0))
   {
       id2 = index;
   }
   else if((index+npts2/2)>=npts)
   {
          id2 = index-npts2;
   }
   else
   {
       id2 = index-npts2/2;
   }

   for(int i=0;i<npts2;++i)
   {
       xc2[i] = xpts[id2];
       uc2[i] = funcvals[id2];
       id2++;
   } */

   sum = funcvals[index];
   //for(int i=0;i<npts2;++i)
   //{
//       sum += uc2[i]*LagrangePoly(x,i,npts2,xc2);
 //  }

 //  cout << "u_i=" << sum << endl;

   return sum;
}

NekDouble ALE::LinearInterpolation(NekDouble x, NekDouble y, int npts, const Array<OneD, const NekDouble>& xfs,
                                const Array<OneD, const NekDouble>& yfs, const Array<OneD, const NekDouble>& vfs)
{
   NekDouble sum = 0.0;

   NekDouble dif=fabs(x-xfs[0]);
   NekDouble dift;

   int index=0;

   // find index of closest x-value on centre to x
   for(int i=0;i<npts;++i)
   {
       if ( (dift=fabs(x-xfs[i])) < dif)
       {
           index=i;
           dif=dift;
       }
   }

   sum = y*(vfs[index]/yfs[index]);

  // cout << "x=" << x << "id=" << index << "x_fs=" << xfs[index] << "y_fs=" << yfs[index] << "u_c=" << vfs[index] << "wy=" << sum << endl;

   return sum;
}


NekDouble ALE::LagrangePoly(NekDouble x, int pt, int npts, const Array<OneD, const NekDouble>& xpts)
{
   NekDouble h=1.0;

   for(int i=0;i<pt; ++i)
   {
       h = h * (x - xpts[i])/(xpts[pt]-xpts[i]);
   }

   for(int i=pt+1;i<npts;++i)
   {
       h = h * (x - xpts[i])/(xpts[pt]-xpts[i]);
   }

   return h;
}

} //end of namespace

