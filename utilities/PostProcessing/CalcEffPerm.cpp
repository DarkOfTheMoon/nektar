#include <cstdio>
#include <cstdlib>

#include <MultiRegions/ExpList.h>
#include <MultiRegions/ExpList2D.h>
#include <MultiRegions/ExpList3D.h>

#include <SolverUtils/EquationSystem.h>
#include <SolverUtils/UnsteadySystem.h>

#include <LibUtilities/LinearAlgebra/Lapack.hpp>
#include <LibUtilities/LinearAlgebra/NekLinSys.hpp>

#include <tinyxml/tinyxml.h>

using namespace Nektar;

void showme(
    int rows,
    int cols,
    DNekMatSharedPtr Matrix);

void showme(
    int size,
    Array<OneD, NekDouble> Vektor);

Array<OneD, NekDouble> proj(
    int size,
    Array<OneD, NekDouble> a,
    Array<OneD, NekDouble> e);

void get_a(
    int size,
    Array<OneD, NekDouble> a,
    DNekMatSharedPtr A,
    int index);

int main(int argc, char *argv[])
{
    int i,j,n;

    if(argc < 4)
    {
        fprintf(stderr,"Usage for 2D: CalcEffPerm meshfile first_fld second_fld\nUsage for 3D: CalcEffPerm meshfile first_fld second_fld third_fld\n");
        exit(1);
    }

    LibUtilities::SessionReaderSharedPtr vSession
            = LibUtilities::SessionReader::CreateInstance(argc, argv);

    //----------------------------------------------
    // Read in mesh from input file
    string meshfile(argv[1]);
    SpatialDomains::MeshGraphSharedPtr graphShPt = SpatialDomains::MeshGraph::Read(vSession);//meshfile);
    //----------------------------------------------

    //----------------------------------------------
    // Import field file.
    int num_fld = argc-2;
    Array<OneD, string> fieldfile(num_fld);
    for(i = 0; i < num_fld; ++i)
    {
        fieldfile[i] = string(argv[2+i]);
    }
    Array<OneD, vector<LibUtilities::FieldDefinitionsSharedPtr> > fielddef(num_fld);
    Array<OneD, vector<vector<NekDouble> > > fielddata(num_fld);
    for(i = 0; i < num_fld; ++i)
    {
        LibUtilities::Import(fieldfile[i],fielddef[i],fielddata[i]);
    }
    //----------------------------------------------

    //----------------------------------------------
    // some more initialations

    // dimension
    int expdim  = graphShPt->GetMeshDimension();
    ASSERTL0(expdim != 1, "Calculation of the effective permeability is not set up for 1D problems");
    ASSERTL0(num_fld >= expdim, "Need two fld files for 2D problems and three for 3D problems to calculate effective permeability");
    int ncols   = 3*(expdim-1);
    int nrows   = num_fld*expdim;
    int nfields = fielddef[0][0]->m_fields.size(); // assuming both fields are the same size

    // Initialise fields 
    Array<OneD, Array<OneD, NekDouble> > grad(nfields*expdim);
    //Array<OneD, Array<OneD, NekDouble> > grad2(2*expdim);
    Array<OneD, Array<OneD, NekDouble> > forc(expdim);
    //Array<OneD, Array<OneD, NekDouble> > visc(expdim);
    Array<OneD, NekDouble> avvel(nrows);
    Array<OneD, NekDouble> avforc(nrows);

    // Extracting extend for reduced "cube"
    int err;
    SpatialDomains::Composite comp_it;
    Array <OneD, vector<NekDouble> > extend(expdim);
    Array <OneD, Array <OneD, NekDouble> > tmp_Vcoords(3);
    for (i=0; i<3; ++i)
    {
        tmp_Vcoords[i] = Array <OneD, NekDouble>(3);
    }
    NekDouble tol = 1e-10; 
    NekDouble vol;

    // Calculate average forcing term
    NekDouble m_kinvis;
    vSession->LoadParameter("Kinvis", m_kinvis);
    // evaluate force function
    std::string forc_var[3] = {
        "u",
        "v",
        "w"
    };
    std::string body_forces[3] = {
        "BodyForce1",
        "BodyForce2",
        "BodyForce3"
    };

    Array<OneD, NekDouble > x;
    Array<OneD, NekDouble > y;
    Array<OneD, NekDouble > z;

    //----------------------------------------------

    //----------------------------------------------
    // Define Expansion
    Array<OneD, Array<OneD, MultiRegions::ExpListSharedPtr> > Exp(num_fld);
    for(i = 0; i < num_fld; ++i)
    {
        Exp[i] = Array<OneD, MultiRegions::ExpListSharedPtr>(nfields);
    }	

    // looping over all fld-files
    for(n = 0; n < num_fld; ++n)
    {
        switch(expdim)
        {
            case 1:
            {
                ASSERTL0(false,"Calculation of the effective permeability is not set up for 1D problems");
            }
            break;
            case 2:
            {
                MultiRegions::ExpList2DSharedPtr Exp2D;
                Exp2D = MemoryManager<MultiRegions::ExpList2D>
                    ::AllocateSharedPtr(vSession,graphShPt);
                Exp[n][0] =  Exp2D;
                
                for(i = 1; i < nfields; ++i)
                {
                    Exp[n][i] = MemoryManager<MultiRegions::ExpList2D>
                        ::AllocateSharedPtr(*Exp2D);
                }
            }
            break;
            case 3:
            {
                MultiRegions::ExpList3DSharedPtr Exp3D;
                Exp3D = MemoryManager<MultiRegions::ExpList3D>
                    ::AllocateSharedPtr(vSession,graphShPt);
                Exp[n][0] =  Exp3D;

                for(i = 1; i < nfields; ++i)
                {
                    Exp[n][i] = MemoryManager<MultiRegions::ExpList3D>
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
            for(i = 0; i < fielddata[n].size(); ++i)
            {
                Exp[n][j]->ExtractDataToCoeffs(fielddef [n][i],
                                               fielddata[n][i],
                                               fielddef [n][i]->m_fields[j],
                                               Exp[n][j]->UpdateCoeffs());
            }
            Exp[n][j]->BwdTrans(Exp[n][j]->GetCoeffs(),Exp[n][j]->UpdatePhys());
        }
        //----------------------------------------------

        //----------------------------------------------
        // Initialisation
        int nq = Exp[0][0]->GetNpoints();

        for(i = 0; i < expdim; ++i)
        {
            forc[i] = Array<OneD, NekDouble>(nq);
        }
    
        for(i = 0; i < nfields*expdim; ++i)
        {
            grad[i] = Array<OneD, NekDouble>(nq);
            //grad2[i] = Array<OneD, NekDouble>(nq);
        }
        //----------------------------------------------
        
        //----------------------------------------------
        // Calculate average Velocity
    
        //NekDouble x_tmp, y_tmp, z_tmp;
        for (i=0; i<expdim; ++i)
        {
            // for (j=0; j<nq; ++j)
            // {
            //     Exp[i]->GetCoords(x_tmp, y_tmp, z_tmp);
            //     //if ( (*Exp[i])[j]
            // }
            avvel[n*expdim+i] = (Exp[n][i]->PhysIntegral());
        }
        //----------------------------------------------

        //----------------------------------------------
        // Calculate average forcing term

        // calculating first deriv of velocity and pressure
        if(expdim == 2)
        {
            for(i = 0; i < nfields; ++i)
            {
                Exp[n][i]->PhysDeriv(Exp[n][i]->GetPhys(), grad[i*expdim], grad[i*expdim+1]);
            }
        }
        else
        {
            for(i = 0; i < nfields; ++i)
            {
                Exp[n][i]->PhysDeriv(Exp[n][i]->GetPhys(), grad[i*expdim], grad[i*expdim+1], grad[i*expdim+2]);
            }
        }
    
        // calculating second deriv of velocity
        for(i = 0; i < expdim; ++i)
        {
            for(j = 0; j < expdim; ++j)
            {
                Exp[n][i]->PhysDeriv(MultiRegions::DirCartesianMap[j], grad[expdim*i+j], grad[expdim*i+j]);
            }
        }

        // added
        // ux,0,vy,0,wz,0
        // for(i = 0; i < nfields; ++i)
        // {
        //     Exp[n][i]->PhysDeriv(MultiRegions::DirCartesianMap[i], Exp[n][i]->GetPhys(), grad2[2*i]);
        // }
        
        // Exp[n][0]->PhysDeriv(MultiRegions::DirCartesianMap[2], grad2[0], grad2[1]);
        // Exp[n][0]->PhysDeriv(MultiRegions::DirCartesianMap[1], grad2[0], grad2[0]);
        // Exp[n][1]->PhysDeriv(MultiRegions::DirCartesianMap[2], grad2[2], grad2[3]);
        // Exp[n][1]->PhysDeriv(MultiRegions::DirCartesianMap[0], grad2[2], grad2[2]);
        // Exp[n][2]->PhysDeriv(MultiRegions::DirCartesianMap[1], grad2[4], grad2[5]);
        // Exp[n][2]->PhysDeriv(MultiRegions::DirCartesianMap[0], grad2[4], grad2[4]);
        // uxy,uxz,vyx,vyz,wzx,wzy
        
        //Vmath::Vcopy(nq, grad[], 1, viscExp[n][i]->UpdatePhys(), 1);

        // added
    
        // evaluate force function
        x = Array<OneD, NekDouble >(nq);
        y = Array<OneD, NekDouble >(nq);
        z = Array<OneD, NekDouble >(nq);
        Exp[n][0]->GetCoords(x,y,z);

        for(i = 0; i < expdim; ++i)
        {
            // f
            vSession->GetFunction(body_forces[n], forc_var[i])->Evaluate(x,y,z,forc[i]);
        }

        for (i=0; i<expdim; ++i)
        {
            // f - grad*p
            //Vmath::Vsub(nq, forc[i], 1, grad[expdim*expdim+i], 1, forc[i], 1);
            for (j=0; j<expdim; ++j)
            {
                // f - grad*p + kinvis*grad^2*u
                Vmath::Svtvp(nq, m_kinvis, grad[expdim*i+j], 1, forc[i], 1, forc[i], 1);
            }
        }
 
        // writing forcing term to m_phys
        for (i=0; i<expdim; ++i)
        {
            Vmath::Vcopy(nq, forc[i], 1, Exp[n][i]->UpdatePhys(), 1);
        }
    
        // calculating average over domain
        for (i=0; i<expdim; ++i)
        {
            avforc[n*expdim+i] = (Exp[n][i]->PhysIntegral())/m_kinvis;
        }
        //----------------------------------------------
    }
    //----------------------------------------------

    //----------------------------------------------
    // Extracting extend for reduced "cube"/"rectangle"
    TiXmlDocument doc(argv[1]);
    bool loadOkay = doc.LoadFile();
    ASSERTL0(loadOkay, "Unable to load input meshfile.");

    TiXmlElement* comps = doc.FirstChildElement("NEKTAR")
        ->FirstChildElement("GEOMETRY")
        ->FirstChildElement("COMPOSITE");
    ASSERTL0(comps, "Tag COMPOSITE could not be found in meshfile.");

    TiXmlElement* it_comp = comps->FirstChild("C")->ToElement();
    string check_body;
    int start;

    int n_comps = atoi(comps->LastChild("C")->ToElement()->Attribute("ID"));

    for(i=0; i<n_comps; ++i)
    {
        check_body = it_comp->GetText();
        if((check_body[0] == 'E') || (check_body[0] == 'F'))
        {
            start = atoi(it_comp->Attribute("ID"));
            break;
        }
        else
        {
            it_comp = it_comp->NextSiblingElement();
        }
    }

    // looping over all 6 walls or 4 edges
    // assumming:
    // C[0]:   tets
    // C[1]:   prisms (not necessary)
    // C[2-7]: outer walls
    // C[8]:   inner walls (SMC walls)
    // C[9]:   all faces, that have to be assigned
    // or:
    // C[0]:   tri
    // C[1]:   rectangle (not necesaary)
    // C[2-5]: outer edges
    // C[6]:   inner edges (SMC walls)

    for (j=start; j<(start+2*expdim); ++j)
    {
        comp_it = graphShPt->GetComposite(j);
        
        for (i=0; i<expdim; ++i)
        {
            graphShPt
                ->GetVertex((*comp_it)[0]->GetVid(i))
                ->GetCoords(tmp_Vcoords[0][i],
                            tmp_Vcoords[1][i],
                            tmp_Vcoords[2][i]);
        }

        err = 0;
        for (i=0; i<expdim; ++i)
        {
            if ( abs( tmp_Vcoords[i][0] - tmp_Vcoords[i][1] ) < tol )
            {
                if(expdim == 2)
                {
                    extend[i].push_back(tmp_Vcoords[i][0]);
                }
                else
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
            }
            else
            {
                err++;
            }
        }
        ASSERTL0(err < 3, "Extend of the geometry could not be calculated. Outer walls appear not to be in spatial planes.");
    }
    
    // rearranging extend entries and calculating area/volume
    vol = 1;
    for (i=0; i<expdim; ++i)
    {
        if (extend[i][0] > extend[i][1])
        {
            NekDouble tmp = extend[i][0];
            extend[i][0]  = extend[i][1];
            extend[i][1]  = tmp;
        }
        vol *= (extend[i][1] - extend[i][0]);
    }

    // vol = ( extend[0][1] - extend[0][0] )
    //      *( extend[1][1] - extend[1][0] )
    //      *( extend[2][1] - extend[2][0] );

    Vmath::Smul(nrows,1/vol,avvel,1,avvel,1);
    Vmath::Smul(nrows,1/vol,avforc,1,avforc,1);
    //----------------------------------------------

    //----------------------------------------------
    // calculation of the effective permeability matrix
    // by solving A*K=b
    
    //int v,m;
    NekDouble zero = 0.0;
    MatrixStorage storage = eFULL;

    DNekMatSharedPtr tmp_mat = 
        MemoryManager<DNekMat>::AllocateSharedPtr
        (nrows,ncols,zero,storage);

    NekDouble MatValue;

    // Assembling Matrix A
    for(i=0; i<expdim; ++i)
    {
        for(n=0; n<num_fld; ++n)
        {
            MatValue = avforc[n*expdim+i];
            tmp_mat->SetValue(n*expdim+i,i,MatValue);
        }
    }
    if(expdim == 2)
    {
        for(n=0; n<num_fld; ++n)
        {
            MatValue = avforc[n*expdim+1];
            tmp_mat->SetValue(n*expdim+0,expdim,MatValue);
            MatValue = avforc[n*expdim+0];
            tmp_mat->SetValue(n*expdim+1,expdim,MatValue);
        }
    }
    else
    {
        for(i=0; i<expdim; ++i)
        {
            for(j=0; j<expdim; ++j)
            {
                if( i == j && i != 1 )
                {
                    for(n=0; n<num_fld; ++n)
                    {
                        MatValue = avforc[n*expdim+1];
                        tmp_mat->SetValue(n*expdim+i,expdim+i,MatValue);
                    }
                }
                if( i == j+1 )
                {
                    for(n=0; n<num_fld; ++n)
                    {
                        MatValue = avforc[n*expdim];
                        tmp_mat->SetValue(n*expdim+i,expdim+j,MatValue);
                    }
                }
                if( i == j-1 )
                {
                    for(n=0; n<num_fld; ++n)
                    {
                        MatValue = avforc[n*expdim+2];
                        tmp_mat->SetValue(n*expdim+i,expdim+j,MatValue);
                    }
                }  
            }
        }
    }

    cout << "A" << endl;
    showme(nrows,ncols,tmp_mat);
    cout << "avvel" << endl;
    showme(nrows,avvel);
    cout << "avforc" << endl;
    showme(nrows,avforc);

    // QR decomposition
    Array<OneD, Array<OneD, NekDouble> > Q_tmp(ncols);
    for(i=0; i<ncols; ++i)
    {
        Q_tmp[i] = Array<OneD, NekDouble>(nrows);
        Vmath::Zero(nrows,Q_tmp[i],1);
    }

    // constructing orthonormal basis using Gram-Schmidt process
    int cnt = 0;
    NekDouble alpha;
    Array<OneD, NekDouble> a(nrows),u(nrows),proj_tmp(nrows);
//    if(expdim == 3)
//    {
        for(i=0; i<ncols; ++i)
        {
            //Vmath::Vcopy(nrows,test_case[i],1,a,1);
            get_a(nrows,a,tmp_mat,i);
            Vmath::Vcopy(nrows,a,1,u,1);
            
            for(j=0; j<cnt; ++j)
            {
                proj_tmp = proj(nrows,a,Q_tmp[j]);
                Vmath::Vsub(nrows,u,1,proj_tmp,1,u,1);
            }
            alpha = sqrt(Vmath::Dot(nrows,u,u));
            Vmath::Smul(nrows,1/alpha,u,1,u,1);
            Vmath::Vcopy(nrows,u,1,Q_tmp[i],1);
            cnt++;
        }
    // }
    // else
    // {
    //     // not implemented yet
    // }

    // creating matrix Q
    DNekMatSharedPtr Q = 
        MemoryManager<DNekMat>::AllocateSharedPtr
        (nrows,ncols,zero,storage);

    for(i=0; i<nrows; ++i)
    {
        for(j=0; j<ncols; ++j)
        {
            Q->SetValue(i,j,Q_tmp[j][i]);
        }
    }    
    
    cout << "Q" << endl;
    showme(nrows,ncols,Q);

    // calculating matrix R
    Q->Transpose();

    DNekMatSharedPtr R = 
        MemoryManager<DNekMat>::AllocateSharedPtr
        (ncols,ncols,zero,storage);

    DNekMat &Q_wrap = (*Q);
    DNekMat &A_wrap = (*tmp_mat);
    DNekMat &R_wrap = (*R);

    R_wrap = Q_wrap * A_wrap;

    cout << "R" << endl;
    showme(ncols,ncols,R);

    // setting values of R smaller than the tolerance to zero
    // in order to create upper triangle matrix
    for(i=0; i<ncols; ++i)
    {
        for(j=0; j<ncols; ++j)
        {
            if( abs(R->GetValue(i,j)) < tol )
            {
                R->SetValue(i,j,0);
            }
        }
    }  

    showme(ncols,ncols,R);

    // calculating effective permeability
    Array<OneD, NekDouble> K(ncols);
    NekVector<NekDouble> K_wrap(ncols,K,eWrapper);
    NekVector<NekDouble> b_wrap(nrows,avvel,eWrapper);

    R->Invert();

    cout << "R^-1" << endl;
    showme(ncols,ncols,R);

    K_wrap = R_wrap * Q_wrap * b_wrap;

    // setting values of K smaller than the tolerance to zero    
    for(i=0; i<ncols; ++i)
    {
        if( abs(K[i]) < tol )
        {
            K[i] = 0;
        }
    }
    //----------------------------------------------
    
    //----------------------------------------------
    // Check if permeability matrix is positive definite
    if (expdim == 2)
    {
        ASSERTL0(K[0] > 0,"Permeability Matrix is not positive definite");
        NekDouble detTemp = K[0]*K[1] - K[2]*K[2];
        ASSERTL0(detTemp > 0,"Permeability Matrix is not positive definite");
    }
    else
    {
        ASSERTL0(K[0] > 0,"Permeability Matrix is not positive definite");
        NekDouble detTemp = K[0]*K[1] - K[3]*K[3];
        ASSERTL0(detTemp > 0,"Permeability Matrix is not positive definite");
        detTemp = K[0]*(K[1]*K[2]-K[5]*K[5])
                 -K[3]*(K[2]*K[3]-K[4]*K[5])
                 +K[4]*(K[3]*K[5]-K[1]*K[4]);
        ASSERTL0(detTemp > 0,"Permeability Matrix is not positive definite");
    }
    //----------------------------------------------

    //----------------------------------------------
    // basic output of K
    for(i=0; i<ncols; ++i)
    {
        cout << "K[" << i << "] = " << K[i] << endl;
    }
    //----------------------------------------------

    return 0;
}

void showme(
    int rows,
    int cols,
    DNekMatSharedPtr Matrix)
{
    for(int i=0; i<rows; ++i)
    {
        for(int j=0; j<cols; ++j)
        {
            cout << Matrix->GetValue(i,j) << "  ";
        }
        cout << endl;
    }
    cout << endl;
    cout << endl;
}

void showme(
    int size,
    Array<OneD, NekDouble> Vektor)
{
    for(int i=0; i<size; ++i)
    {
        cout << Vektor[i] << endl;
    }
    cout << endl;
}


Array<OneD, NekDouble> proj(
    int size,
    Array<OneD, NekDouble> a,
    Array<OneD, NekDouble> e)
{
    Array<OneD, NekDouble> return_vec(size);
    NekDouble dotea, dotee;
    NekDouble tol = 1e-10;
    
    dotea = Vmath::Dot(size,e,a);
    dotee = Vmath::Dot(size,e,e);

    ASSERTL0( dotee-1 < tol, "Base vector in Gram-Schmidt process not of length 1");

    Vmath::Smul(size,dotea/dotee,e,1,return_vec,1);

    return return_vec;
}
    
void get_a(
    int size,
    Array<OneD, NekDouble> a,
    DNekMatSharedPtr A,
    int index)
{
    for(int i=0; i<size; ++i)
    {
        a[i] = A->GetValue(i,index);
    }
}
