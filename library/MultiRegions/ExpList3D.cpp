///////////////////////////////////////////////////////////////////////////////
//
// File ExpList3D.cpp
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
// Description: Expansion list 3D definition
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/ExpList3D.h>

namespace Nektar
{
  namespace MultiRegions
  {

    ExpList3D::ExpList3D(): ExpList()
    {
    }
        
    ExpList3D::ExpList3D(const ExpList3D &In): ExpList(In)
    {
    }

    ExpList3D::~ExpList3D()
    {
    }

    namespace
        {
            // Adds up the number of cells in a truncated Nc by Nc by Nc pyramid, 
            // where the longest Na rows and longest Nb columns are kept.
            // Example: (Na, Nb, Nc) = (3, 4, 5); The number of coefficients is the 
            // sum of the elements of the following matrix:
            //     |5  4  3  2  0|
            //     |4  3  2  0   |
            //     |3  2  0      |
            //     |0  0         |
            //     |0            |
            // Sum = 28 = number of tet coefficients
            inline int getNumberOfCoefficients( int Na, int Nb, int Nc ) 
            {
                int nCoef = 0;
                for( int a = 0; a < Na; ++a )
                {
                    for( int b = 0; b < Nb - a; ++b )
                    {
                        for( int c = 0; c < Nc - a - b; ++c )
                        {
                            ++nCoef;
                        }
                    }
                }
                return nCoef;
            }
        }

       
        ExpList3D::ExpList3D(const LibUtilities::BasisKey &Ba,
                             const LibUtilities::BasisKey &Bb,
                             const LibUtilities::BasisKey &Bc,
                             const SpatialDomains::MeshGraph3D &graph3D,
                             const LibUtilities::PointsType TetNb):
            ExpList()
        {

            LocalRegions::TetExpSharedPtr tet;
            LocalRegions::HexExpSharedPtr hex;
            LocalRegions::PrismExpSharedPtr prism;
            LocalRegions::PyrExpSharedPtr pyramid;

            const SpatialDomains::ExpansionVector &expansions = graph3D.GetExpansions();
                        
            for(int i = 0; i < expansions.size(); ++i)
            {
                SpatialDomains::TetGeomSharedPtr TetGeom;
                SpatialDomains::HexGeomSharedPtr HexGeom;
                SpatialDomains::PrismGeomSharedPtr PrismGeom;
                SpatialDomains::PyrGeomSharedPtr PyrGeom;
                
                if(TetGeom = boost::dynamic_pointer_cast<SpatialDomains::TetGeom>(expansions[i]->m_GeomShPtr)) // Tetrahedron
                {
                    if(TetNb < LibUtilities::SIZE_PointsType)
                    {
                           
//                         Ntet = MemoryManager<LocalRegions::NodalTetExp>::AllocateSharedPtr(TetBa,TetBb,TetBc,TetNb,TetGeom);
//                         (*m_exp).push_back(Ntet);
                    }
                    else
                    {
                        tet = MemoryManager<LocalRegions::TetExp>::AllocateSharedPtr(Ba,Bb,Bc,TetGeom);
                        (*m_exp).push_back(tet);
                    }

                       m_ncoeffs += getNumberOfCoefficients(Ba.GetNumModes(), Bb.GetNumModes(), Bc.GetNumModes());
                        
                       m_npoints += Ba.GetNumPoints()*Bb.GetNumPoints()*Bc.GetNumPoints();
                }
                else if(PrismGeom = boost::dynamic_pointer_cast<SpatialDomains::PrismGeom>(expansions[i]->m_GeomShPtr)) // Prism
                {
                      prism = MemoryManager<LocalRegions::PrismExp>::AllocateSharedPtr(Ba,Bb,Bc,PrismGeom);
                      (*m_exp).push_back(prism);

                      m_ncoeffs += getNumberOfCoefficients(Ba.GetNumModes(), Bb.GetNumModes(), Bc.GetNumModes());
                      m_npoints +=  Ba.GetNumPoints()*Bb.GetNumPoints()*Bc.GetNumPoints();
                
                }
                else if(PyrGeom = boost::dynamic_pointer_cast<SpatialDomains::PyrGeom>(expansions[i]->m_GeomShPtr)) // Pyramid
                {
                     pyramid = MemoryManager<LocalRegions::PyrExp>::AllocateSharedPtr(Ba,Bb,Bc,PyrGeom);
                     (*m_exp).push_back(pyramid);

                      m_ncoeffs += getNumberOfCoefficients(Ba.GetNumModes(), Bb.GetNumModes(), Bc.GetNumModes());
                      m_npoints +=  Ba.GetNumPoints()*Bb.GetNumPoints()*Bc.GetNumPoints();

                }
                else if(HexGeom = boost::dynamic_pointer_cast<SpatialDomains::HexGeom>(expansions[i]->m_GeomShPtr)) // Hexahedron
                {
                    hex = MemoryManager<LocalRegions::HexExp>::AllocateSharedPtr(Ba,Bb,Bc, HexGeom);
                    (*m_exp).push_back(hex);
                    
                    m_ncoeffs += Ba.GetNumModes()*Bb.GetNumModes()*Bc.GetNumModes();
                    m_npoints += Ba.GetNumPoints()*Bb.GetNumPoints()*Bc.GetNumPoints();
                }
                else
                {
                    ASSERTL0(false,"dynamic cast to a proper Geometry3D failed");
                }  
                
            }
            
            m_coeffs = Array<OneD, NekDouble>(m_ncoeffs);
            m_phys   = Array<OneD, NekDouble>(m_npoints);
        }


        ExpList3D::ExpList3D(SpatialDomains::MeshGraph3D &graph3D):ExpList()
        {

            LocalRegions::TetExpSharedPtr tet;
            LocalRegions::HexExpSharedPtr hex;
            LocalRegions::PrismExpSharedPtr prism;
            LocalRegions::PyrExpSharedPtr pyramid;
            LibUtilities::PointsType TetNb;

            const SpatialDomains::ExpansionVector &expansions = graph3D.GetExpansions();
            
            for(int i = 0; i < expansions.size(); ++i)
            {
                SpatialDomains::TetGeomSharedPtr TetGeom;
                SpatialDomains::HexGeomSharedPtr HexGeom;
                SpatialDomains::PrismGeomSharedPtr PrismGeom;
                SpatialDomains::PyrGeomSharedPtr PyrGeom;
                
                if(TetGeom = boost::dynamic_pointer_cast<SpatialDomains::TetGeom>(expansions[i]->m_GeomShPtr))
                {
                    LibUtilities::BasisKey TetBa = graph3D.GetBasisKey(expansions[i],0);
                    LibUtilities::BasisKey TetBb = graph3D.GetBasisKey(expansions[i],1);
                    LibUtilities::BasisKey TetBc = graph3D.GetBasisKey(expansions[i],2);
                    
                    if(expansions[i]->m_ExpansionType == SpatialDomains::eNodal)
                    {
                      ASSERTL0(false,"LocalRegions::NodalTetExp is not implemented yet");
                    }
                    else
                    {
                        tet = MemoryManager<LocalRegions::TetExp>::AllocateSharedPtr(TetBa,TetBb,TetBc,TetGeom);
                        (*m_exp).push_back(tet);
                    }
                        
                    m_ncoeffs += getNumberOfCoefficients(TetBa.GetNumModes(), TetBb.GetNumModes(), TetBc.GetNumModes());
                    m_npoints += TetBa.GetNumPoints()*TetBb.GetNumPoints()*TetBc.GetNumPoints();
                }
                else if(PrismGeom = boost::dynamic_pointer_cast<SpatialDomains::PrismGeom>(expansions[i]->m_GeomShPtr))
                {
                    LibUtilities::BasisKey PrismBa = graph3D.GetBasisKey(expansions[i],0);
                    LibUtilities::BasisKey PrismBb = graph3D.GetBasisKey(expansions[i],0);
                    LibUtilities::BasisKey PrismBc = graph3D.GetBasisKey(expansions[i],1);

                    prism = MemoryManager<LocalRegions::PrismExp>::AllocateSharedPtr(PrismBa,PrismBb,PrismBc,PrismGeom);
                    (*m_exp).push_back(prism);

                    m_ncoeffs += getNumberOfCoefficients(PrismBa.GetNumModes(), PrismBb.GetNumModes(), PrismBc.GetNumModes());
                    m_npoints += PrismBa.GetNumPoints()*PrismBb.GetNumPoints()*PrismBc.GetNumPoints();
                }
                else if(PyrGeom = boost::dynamic_pointer_cast<SpatialDomains::PyrGeom>(expansions[i]->m_GeomShPtr))
                {
                    LibUtilities::BasisKey PyrBa = graph3D.GetBasisKey(expansions[i],0);
                    LibUtilities::BasisKey PyrBb = graph3D.GetBasisKey(expansions[i],0);
                    LibUtilities::BasisKey PyrBc = graph3D.GetBasisKey(expansions[i],2);

                    pyramid = MemoryManager<LocalRegions::PyrExp>::AllocateSharedPtr(PyrBa,PyrBb,PyrBc,PyrGeom);
                    (*m_exp).push_back(pyramid);

                    m_ncoeffs += getNumberOfCoefficients(PyrBa.GetNumModes(), PyrBb.GetNumModes(), PyrBc.GetNumModes());
                    m_npoints += PyrBa.GetNumPoints()*PyrBb.GetNumPoints()*PyrBc.GetNumPoints();
                }
                else if(HexGeom = boost::dynamic_pointer_cast<SpatialDomains::HexGeom>(expansions[i]->m_GeomShPtr))
                {
                    LibUtilities::BasisKey HexBa = graph3D.GetBasisKey(expansions[i],0);
                    LibUtilities::BasisKey HexBb = graph3D.GetBasisKey(expansions[i],0);
                    LibUtilities::BasisKey HexBc = graph3D.GetBasisKey(expansions[i],0);
                    
                    hex = MemoryManager<LocalRegions::HexExp>::AllocateSharedPtr(HexBa,HexBb,HexBc,HexGeom);
                    (*m_exp).push_back(hex);
                    
                    m_ncoeffs += HexBa.GetNumModes()*HexBb.GetNumModes()*HexBc.GetNumModes();
                    m_npoints += HexBa.GetNumPoints()*HexBb.GetNumPoints()*HexBc.GetNumPoints();
                }
                else
                {
                    ASSERTL0(false,"dynamic cast to a proper Geometry3D failed");
                }  
                
            }            
            m_coeffs = Array<OneD, NekDouble>(m_ncoeffs);
            m_phys   = Array<OneD, NekDouble>(m_npoints);
        }

        void ExpList3D::SetBoundaryConditionExpansion(SpatialDomains::MeshGraph3D &graph3D,
                                                      SpatialDomains::BoundaryConditions &bcs, 
                                                      const std::string variable,
                                                      Array<OneD, ExpList2DSharedPtr> &bndCondExpansions,
                                                      Array<OneD, SpatialDomains::BoundaryConditionShPtr> &bndConditions)
        {
            int i;
            int cnt  = 0;
            
            SpatialDomains::BoundaryRegionCollection    &bregions = bcs.GetBoundaryRegions();
            SpatialDomains::BoundaryConditionCollection &bconditions = bcs.GetBoundaryConditions();   
            
            MultiRegions::ExpList2DSharedPtr       locExpList;
            SpatialDomains::BoundaryConditionShPtr locBCond; 

            int nbnd = bregions.size();
          
            cnt=0;
            // list Dirichlet boundaries first
            for(i = 0; i < nbnd; ++i)
            {  
                locBCond = (*(bconditions[i]))[variable];  
                if(locBCond->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                {
                   //TODO : implement -- just commented out to avoid the error                
                   // locExpList = MemoryManager<MultiRegions::ExpList2D>::AllocateSharedPtr(*(bregions[i]),graph3D);
                    bndCondExpansions[cnt]  = locExpList;
                    bndConditions[cnt++]    = locBCond;
                } // end if Dirichlet
            }
            // then, list the other (non-periodic) boundaries
            for(i = 0; i < nbnd; ++i)
            {        
                locBCond = (*(bconditions[i]))[variable];  
                if(locBCond->GetBoundaryConditionType() == SpatialDomains::eNeumann)
                {
                   //TODO: implement below                
                   // locExpList = MemoryManager<MultiRegions::ExpList2D>::AllocateSharedPtr(*(bregions[i]),graph3D);
                    bndCondExpansions[cnt]  = locExpList;
                    bndConditions[cnt++]    = locBCond;
                }  
                else if((locBCond->GetBoundaryConditionType() != SpatialDomains::eDirichlet) && 
                        (locBCond->GetBoundaryConditionType() != SpatialDomains::ePeriodic))
                {
                    ASSERTL0(false,"This type of BC not implemented yet");
                }                  
            }
        }

        void ExpList3D::EvaluateBoundaryConditions(const NekDouble time,
                                                   Array<OneD, ExpList2DSharedPtr> &bndCondExpansions,
                                                   Array<OneD, SpatialDomains::BoundaryConditionShPtr> &bndConditions)
        {         
            int i,j;
            int npoints;
            int nbnd = bndCondExpansions.num_elements();
            MultiRegions::ExpList2DSharedPtr locExpList;
            
            for(i = 0; i < nbnd; ++i)
            {                 
                locExpList = bndCondExpansions[i];                  
                npoints = locExpList->GetPointsTot();
                
                Array<OneD,NekDouble> x0(npoints,0.0);
                Array<OneD,NekDouble> x1(npoints,0.0);
                Array<OneD,NekDouble> x2(npoints,0.0);  
                
                locExpList->GetCoords(x0,x1,x2);

                if(bndConditions[i]->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                {            
                    for(j = 0; j < npoints; j++)
                    {
                        (locExpList->UpdatePhys())[j] = (boost::static_pointer_cast<SpatialDomains::DirichletBoundaryCondition>(bndConditions[i])->m_DirichletCondition).Evaluate(x0[j],x1[j],x2[j],time);
                    }
                    
                    locExpList->FwdTrans_BndConstrained(*locExpList);
                }
                else if(bndConditions[i]->GetBoundaryConditionType() == SpatialDomains::eNeumann)
                {          
                    for(j = 0; j < npoints; j++)
                    {
                        (locExpList->UpdatePhys())[j] = (boost::static_pointer_cast<SpatialDomains::NeumannBoundaryCondition>(bndConditions[i])->m_NeumannCondition).Evaluate(x0[j],x1[j],x2[j],time);
                    }

                    locExpList->IProductWRTBase(*locExpList); 
                }
                else
                {
                    ASSERTL0(false,"This type of BC not implemented yet");
                }
            }           
        }


      //TODO: implement
      void ExpList3D::GetPeriodicFaces(SpatialDomains::MeshGraph3D &graph3D,
                                       SpatialDomains::BoundaryConditions &bcs,
                                       const std::string variable,
                                       map<int,int>& periodicVertices,
                                       map<int,int>& periodicEdges,
                                       map<int,int>& periodicFaces)
        {
            int i,j,k;
            
            SpatialDomains::BoundaryRegionCollection    &bregions = bcs.GetBoundaryRegions();
            SpatialDomains::BoundaryConditionCollection &bconditions = bcs.GetBoundaryConditions();

            int region1ID;
            int region2ID;

            SpatialDomains::Composite comp1;
            SpatialDomains::Composite comp2;

            SpatialDomains::SegGeomSharedPtr segmentGeom1;
            SpatialDomains::SegGeomSharedPtr segmentGeom2;

            SpatialDomains::ElementEdgeVectorSharedPtr element1;
            SpatialDomains::ElementEdgeVectorSharedPtr element2;

            StdRegions::EdgeOrientation orient1;
            StdRegions::EdgeOrientation orient2;
            
            SpatialDomains::BoundaryConditionShPtr locBCond; 

            // This std::map is a check so that the periodic pairs
            // are not treated twice
            map<int, int> doneBndRegions;

            int nbnd = bregions.size();
          
            for(i = 0; i < nbnd; ++i)
            {        
                locBCond = (*(bconditions[i]))[variable];  
                if(locBCond->GetBoundaryConditionType() == SpatialDomains::ePeriodic)
                {    
                   region1ID = i;
                   region2ID = (boost::static_pointer_cast<SpatialDomains::PeriodicBoundaryCondition>(locBCond))->m_ConnectedBoundaryRegion;

                    if(doneBndRegions.count(region1ID)==0)
                    {                    
                        ASSERTL0(bregions[region1ID]->size() == bregions[region2ID]->size(),
                                 "Size of the 2 periodic boundary regions should be equal");
                    
                        for(j = 0; j < bregions[region1ID]->size(); j++)
                        {
                            comp1 = (*(bregions[region1ID]))[j];
                            comp2 = (*(bregions[region2ID]))[j];
                            
                            ASSERTL0(comp1->size() == comp2->size(),
                                     "Size of the 2 periodic composites should be equal");
                            
                            for(k = 0; k < comp1->size(); k++)
                            {                                      
                                if(!(segmentGeom1 = boost::dynamic_pointer_cast<SpatialDomains::SegGeom>((*comp1)[k]))||
                                   !(segmentGeom2 = boost::dynamic_pointer_cast<SpatialDomains::SegGeom>((*comp2)[k])))
                                {
                                    ASSERTL0(false,"dynamic cast to a SegGeom failed");
                                } 

                                // Extract the periodic edges
                                periodicEdges[segmentGeom1->GetEid()] = segmentGeom2->GetEid();
                                periodicEdges[segmentGeom2->GetEid()] = segmentGeom1->GetEid();

                                // Extract the periodic vertices
//                                 element1 = graph3D.GetElementsFromEdge(segmentGeom1);
//                                 element2 = graph3D.GetElementsFromEdge(segmentGeom2);

                                ASSERTL0(element1->size()==1,"The periodic boundaries belong to more than one element of the mesh");
                                ASSERTL0(element2->size()==1,"The periodic boundaries belong to more than one element of the mesh");

                                orient1 = (boost::dynamic_pointer_cast<SpatialDomains::Geometry2D>((*element1)[0]->m_Element))->
                                    GetEorient((*element1)[0]->m_EdgeIndx);
                                orient2 = (boost::dynamic_pointer_cast<SpatialDomains::Geometry2D>((*element2)[0]->m_Element))->
                                    GetEorient((*element2)[0]->m_EdgeIndx);

                                if(orient1!=orient2)
                                {
                                    periodicVertices[segmentGeom1->GetVid(0)] = segmentGeom2->GetVid(0);
                                    periodicVertices[segmentGeom1->GetVid(1)] = segmentGeom2->GetVid(1);
                                }
                                else
                                {
                                    periodicVertices[segmentGeom1->GetVid(0)] = segmentGeom2->GetVid(1);
                                    periodicVertices[segmentGeom1->GetVid(1)] = segmentGeom2->GetVid(0);
                                }
                            }
                        }
                    }
                    else
                    {
                        ASSERTL0(doneBndRegions[region1ID]==region2ID,
                                 "Boundary regions are not mutually periodic");
                    }
                    doneBndRegions[region2ID] = region1ID;
                }                  
            }
        } // end of GetPeriodicEdge()


  } //end of namespace
} //end of namespace

