#include <cstdio>
#include <cstdlib>

#include <MultiRegions/ExpList.h>
#include <MultiRegions/ExpList1D.h>
#include <MultiRegions/ExpList2D.h>
#include <MultiRegions/ExpList3D.h>
#include <MultiRegions/ExpList2DHomogeneous1D.h>
#include <MultiRegions/ExpList3DHomogeneous1D.h>
#include <MultiRegions/ExpList3DHomogeneous2D.h>

using namespace Nektar;

int main(int argc, char *argv[])
{
    int i,j;

    if(argc != 4)
    {
        fprintf(stderr,"Usage: FldAddVort  meshfile infld outfld\n");
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
    int addfields = (nfields == 3)? 3:1;
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
    Array<OneD, Array<OneD, NekDouble> > integral(nfields);
    Array<OneD, Array<OneD, NekDouble> > grad(nfields*nfields);
    Array<OneD, Array<OneD, NekDouble> > forc(nfields);
    Array<OneD, NekDouble> avvel(nfields);
    Array<OneD, NekDouble> avforc(nfields);
    
    for(i = 0; i < nfields; ++i)
    {
        integral[i] = Array<OneD, NekDouble>(nq);
        grad[i] = Array<OneD, NekDouble>(nq);
        forc[i] = Array<OneD, NekDouble>(nq);
    }

    // Calculate average Velocity
    

/*    // Calculate Gradient & Vorticity
    if(nfields == 2)
    {
        for(i = 0; i < nfields; ++i)
        {
            Exp[i]->PhysDeriv(Exp[i]->GetPhys(), grad[i*nfields], grad[i*nfields+1]);
            Exp[i]->PhysDeriv(MultiRegions::DirCartesianMap[0], grad[i], grad[i]);
            Exp[i]->PhysDeriv(MultiRegions::DirCartesianMap[1], grad[i], grad[i]);
        }
        // Ux.Vy - Uy.Vx
        Vmath::Vmul (nq,grad[1],1,grad[nfields],1,outfield[0],1);
        Vmath::Vvtvm(nq,grad[0],1,grad[nfields+1],1,outfield[0],1,outfield[0],1);
    }
    else
    {
        for(i = 0; i < nfields; ++i)
        {
            Exp[i]->PhysDeriv(Exp[i]->GetPhys(), grad[i*nfields],grad[i*nfields+1],grad[i*nfields+2]);
        }

        // W_x = Vy.Wz - Vz.Wy
        Vmath::Vmul (nq,grad[nfields+2],1,grad[2*nfields+1],1,outfield[0],1);
        Vmath::Vvtvm(nq,grad[nfields+1],1,grad[2*nfields+2],1,outfield[0],1,outfield[0],1);
        // W_y = Wx.Uz - Ux.Wz
        Vmath::Vmul (nq,grad[0],1,grad[2*nfields+2],1,outfield[1],1);
        Vmath::Vvtvm(nq,grad[2*nfields],1,grad[2],1,outfield[1],1,outfield[1],1);
        // W_z = Ux.Vy - Uy.Vx
        Vmath::Vmul (nq,grad[1],1,grad[nfields],1,outfield[2],1);
        Vmath::Vvtvm(nq,grad[0],1,grad[nfields+1],1,outfield[2],1,outfield[2],1);
    }
    
    for (i = 0; i < addfields; ++i)
    {
        Exp[nfields + i]->FwdTrans(outfield[i], Exp[nfields+i]->UpdateCoeffs());
    }
*/    
    //-----------------------------------------------
    // Write solution to file with additional computed fields
/*    string   out(argv[argc-1]);
    std::vector<SpatialDomains::FieldDefinitionsSharedPtr> FieldDef
                                                = Exp[0]->GetFieldDefinitions();
    std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());
    
    vector<string > outname;
    
    if(addfields == 1)
    {
        outname.push_back("W_z");
    }
    else
    {
        outname.push_back("W_x");
        outname.push_back("W_y");
        outname.push_back("W_z");
    }

    for(j = 0; j < nfields + addfields; ++j)
    {
        for(i = 0; i < FieldDef.size(); ++i)
        {
            if (j >= nfields)
            {
                FieldDef[i]->m_fields.push_back(outname[j-nfields]);
            }
            else
            {
                FieldDef[i]->m_fields.push_back(fielddef[i]->m_fields[j]);
            }
            Exp[j]->AppendFieldData(FieldDef[i], FieldData[i]);
        }
    }
    graphShPt->Write(out, FieldDef, FieldData);
*/    //-----------------------------------------------

    return 0;
}

