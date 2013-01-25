#include <cstdio>
#include <cstdlib>

#include <MultiRegions/ExpList.h>
#include <MultiRegions/ExpList2D.h>
#include <MultiRegions/ExpList3D.h>

using namespace Nektar;

int main(int argc, char *argv[])
{
    int i,j;

    if(argc != 3)
    {
        fprintf(stderr,"Usage: CalcEffPerm  meshfile first_fld second_fld\n");
        exit(1);
    }

    LibUtilities::SessionReaderSharedPtr vSession
            = LibUtilities::SessionReader::CreateInstance(argc, argv);

    //----------------------------------------------
    // Read in mesh from input file
    string meshfile(argv[argc-2]);
    SpatialDomains::MeshGraphSharedPtr graphShPt = SpatialDomains::MeshGraph::Read(vSession);//meshfile);
    //----------------------------------------------

    //----------------------------------------------
    // Import field file.
    string fieldfile(argv[argc-1]);
    vector<SpatialDomains::FieldDefinitionsSharedPtr> fielddef;
    vector<vector<NekDouble> > fielddata;
    graphShPt->Import(fieldfile,fielddef,fielddata);
	bool useFFT = false;
	bool dealiasing = false;
    //----------------------------------------------

    //----------------------------------------------
    // Define Expansion
    int expdim  = graphShPt->GetMeshDimension();
    int nfields = fielddef[0]->m_fields.size();
    int addfields = (nfields == 3)? 3:1;           //<-- nfields in 3D is 4, since pressure is field as well
    Array<OneD, MultiRegions::ExpListSharedPtr> Exp(nfields + addfields);
	
    switch(expdim)
    {
    case 1:
        {
            ASSERTL0(false,"Calculation of the effective permeability is not set up for 1D problems");
        }
        break;
    case 2:
        {
            ASSERTL0(fielddef[0]->m_numHomogeneousDir <= 1,"NumHomogeneousDir is only set up for 1");

            if(fielddef[0]->m_numHomogeneousDir == 1)
            {
                ASSERTL0(fielddef[0]->m_numHomogeneousDir <= 1,"Only full 2D");
            }
            else
            {
                MultiRegions::ExpList2DSharedPtr Exp2D;
                Exp2D = MemoryManager<MultiRegions::ExpList2D>
                                                        ::AllocateSharedPtr(vSession,graphShPt);
                Exp[0] =  Exp2D;
                
                for(i = 1; i < nfields + addfields; ++i)
                {
                    Exp[i] = MemoryManager<MultiRegions::ExpList2D>
                        ::AllocateSharedPtr(*Exp2D);
                }
            }
        }
        break;
    case 3:
        {
            MultiRegions::ExpList3DSharedPtr Exp3D;
            Exp3D = MemoryManager<MultiRegions::ExpList3D>
                                                    ::AllocateSharedPtr(vSession,graphShPt);
            Exp[0] =  Exp3D;

            for(i = 1; i < nfields + addfields; ++i)
            {
                Exp[i] = MemoryManager<MultiRegions::ExpList3D>
                    ::AllocateSharedPtr(*Exp3D);
            }
        }
        break;
    default:
        ASSERTL0(false,"Expansion dimension not recognised");
        break;
    }
    //----------------------------------------------

    //----------------------------------------------
    // Copy data from field file
    for(j = 0; j < nfields; ++j)
    {
        for(int i = 0; i < fielddata.size(); ++i)
        {
            Exp[j]->ExtractDataToCoeffs(fielddef [i],
                                        fielddata[i],
                                        fielddef [i]->m_fields[j]);
        }
        Exp[j]->BwdTrans(Exp[j]->GetCoeffs(),Exp[j]->UpdatePhys());
    }
    //----------------------------------------------

    //----------------------------------------------
    // Initialise fields 
    ASSERTL0(nfields >= 2, "Need two fields (u,v) to add reentricity");
    int nq = Exp[0]->GetNpoints();
    Array<OneD, Array<OneD, NekDouble> > grad(nfields*expdim);
    Array<OneD, Array<OneD, NekDouble> > forc(expdim);
    Array<OneD, NekDouble> avvel(expdim);
    Array<OneD, NekDouble> avforc(expdim);
    
    for(i = 0; i < expdim; ++i)
    {
        forc[i] = Array<OneD, NekDouble>(nq);
    }
    
    for(i = 0; i < nfields*expdim; ++i)
    {
        grad[i] = Array<OneD, NekDouble>(nq);
    }
    //----------------------------------------------

    //----------------------------------------------
    // Extracting extend for reduced "cube"
    
    int err;
    SpatialDomains::Composite comp_it;
    Array <OneD, vector<NekDouble> > extend(expdim);
    Array <OneD, Array <OneD, NekDouble> > tmp_Vcoords(expdim);
    for (i=0; i<expdim; ++i)
    {
        tmp_Vcoords[i] = Array <OneD, NekDouble>(3);
    }
    NekDouble tol = 1e-06;
    NekDouble vol = 0;
    
    // looping over all 6 walls
    // assumming:
    // C[0]:   tets
    // C[1]:   prisms
    // C[2-7]: outer walls
    // C[8]:   inner walls (SMC walls)
    // C[9]:   all faces, that have to be assigned
    for (j=1; j<7; ++j) // <-- needs to be changed back to "for (j=2; j<8; ++j)"
    {
        comp_it = graphShPt->GetComposite(j);
        
        for (i=0; i<3; ++i)
        {
            graphShPt
                ->GetVertex((*comp_it)[0]->GetVid(i))
                ->GetCoords(tmp_Vcoords[0][i],
                            tmp_Vcoords[1][i],
                            tmp_Vcoords[2][i]);
        }

        err = 0;
        for (i=0; i<3; ++i)
        {
            if ( abs( tmp_Vcoords[i][0] - tmp_Vcoords[i][1] ) < tol )
            {
                if ( abs( tmp_Vcoords[i][0] - tmp_Vcoords[i][2] ) < tol )
                {
                    extend[i].push_back(tmp_Vcoords[i][0]);
                }
                else
                {
                    err++;
                }
            }
            else
            {
                err++;
            }
        }
        ASSERTL0(err < 3, "Extend of the geometry could not be calculated. Outer walls appear not to be in spatial planes.");
    }

    vol = abs( ( extend[0][0] - extend[0][1] )
              *( extend[1][0] - extend[1][1] )
              *( extend[2][0] - extend[2][1] ) );
    //----------------------------------------------

    //----------------------------------------------
    // Calculate average Velocity

    for (i=0; i<expdim; ++i)
    {
        avvel[i] = (Exp[i]->PhysIntegral())/vol;
    }
    //----------------------------------------------

    //----------------------------------------------
    // Calculate average forcing term

    NekDouble m_kinvis;
    vSession->LoadParameter("Kinvis", m_kinvis);

    // calculating first deriv of velocity and pressure
    for(i = 0; i < nfields; ++i)
    {
        Exp[i]->PhysDeriv(Exp[i]->GetPhys(), grad[i*expdim], grad[i*expdim+1], grad[i*expdim+2]);
    }
    
    // calculating second deriv of verlocity
    for(i = 0; i < expdim; ++i)
    {
        for(j = 0; j < expdim; ++j)
        {
            Exp[i]->PhysDeriv(MultiRegions::DirCartesianMap[j], grad[expdim*i+j], grad[expdim*i+j]);
        }
    }
    
    // evaluate force function
    std::string forc_var[3] = {
        "u",
        "v",
        "w"
    };
// Read in of 'Body Force' needs to be done correctly... or at all
//    for(i = 0; i < expdim; ++i)
//    {
//        EvaluateFunction(forc_var[i], forc[i], "BodyForce");
//    }
    // f
    Vmath::Fill(nq, 2.0, forc[0], 1);
    Vmath::Fill(nq, 1.0, forc[1], 1);
    Vmath::Fill(nq, 0.0, forc[2], 1);

    for (i=0; i<expdim; ++i)
    {
        // f - grad*p
        Vmath::Vsub(nq, forc[i], 1, grad[expdim*expdim+i], 1, forc[i], 1);
        for (j=0; j<expdim; ++j)
        {
            // f - grad*p + kinvis*grad^2*u
            Vmath::Svtvp(nq, m_kinvis, grad[expdim*i+j], 1, forc[i], 1, forc[i], 1);
        }
    }
 
    // writing forcing term to m_phys
    for (i=0; i<expdim; ++i)
    {
        Vmath::Vcopy(nq, forc[i], 1, Exp[i]->UpdatePhys(), 1);
    }
    
    // calculating average over domain
    for (i=0; i<expdim; ++i)
    {
        avforc[i] = (Exp[i]->PhysIntegral())/(vol*m_kinvis);
    }
    //----------------------------------------------

    //----------------------------------------------
    // basic calculation for Permeability Matrix
    // without off-diagonal entries
    
    NekDouble K_xx, K_yy, K_zz;
    K_xx = avvel[0]/avforc[0];
    K_yy = avvel[1]/avforc[1];
    K_zz = avvel[2]/avforc[2];
    
    // output
    cout << "K_xx = " << K_xx << endl;
    cout << "K_yy = " << K_yy << endl;
    cout << "K_zz = " << K_zz << endl;
    //----------------------------------------------

    return 0;
}

