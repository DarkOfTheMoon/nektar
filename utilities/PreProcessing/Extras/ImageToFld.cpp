#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/Communication/Comm.h>
#include <SpatialDomains/MeshGraph2D.h>
#include <MultiRegions/ExpList2D.h>
#include <SpatialDomains/MeshGraph3D.h>
#include <MultiRegions/ExpList3D.h>

#include <vtkTIFFReader.h>
#include <vtkImageData.h>
#include <vtkSmartPointer.h>
// Usage: VtkToFld session.xml singleimage000 (not tif extension) out.fld

int main(int argc, char* argv[])
{
    MultiRegions::ExpListSharedPtr Exp;

    std::vector<std::string> vFilenames;
    vFilenames.push_back(std::string(argv[1]));

    LibUtilities::SessionReaderSharedPtr vSession
            = LibUtilities::SessionReader::CreateInstance(2, argv, vFilenames);
    try
    {
        //----------------------------------------------
        // Read in mesh from input file
        SpatialDomains::MeshGraphSharedPtr graphShPt = SpatialDomains::MeshGraph::Read(vSession);
        int coordim  = graphShPt->GetMeshDimension();

        //----------------------------------------------
        switch(coordim)
        {
        case 2:
            //MultiRegions::ExpList2DSharedPtr Exp2D;
            Exp = MemoryManager<MultiRegions::ExpList2D>
                ::AllocateSharedPtr(vSession,graphShPt);
            break;
        case 3:
            //MultiRegions::ExpList2DSharedPtr Exp2D;
            Exp = MemoryManager<MultiRegions::ExpList3D>
                ::AllocateSharedPtr(vSession,graphShPt);
            break;
        default:
            ASSERTL0(false,"Coordim not valid");
            break;
        }
        //----------------------------------------------

        //Get number of quadrature points
        int nq      = Exp->GetNpoints();
        
        //Coordinate arrays
        Array<OneD, NekDouble> xc0(nq,0.0);
        Array<OneD, NekDouble> xc1(nq,0.0);
        Array<OneD, NekDouble> xc2(nq,0.0);

        //coordinates of quadrature points
        Exp->GetCoords(xc0,xc1,xc2);

        //----------------------------------------------
        vtkSmartPointer<vtkImageData> imageData = 
            vtkSmartPointer<vtkImageData>::New(); 

        vtkTIFFReader *vtkImageReader = vtkTIFFReader::New();
        //vtkImageReader->SetFileName ( argv[2] );

        NekDouble spacing[3];
        vtkImageReader->SetFilePrefix(argv[2]);
        vtkImageReader->SetDataExtent(0,100,0,100,0,249);
        vtkImageReader->GetDataSpacing(spacing);
        vtkImageReader->SetFilePattern("%s%03d.tif");
        vtkImageReader->UpdateWholeExtent();
        imageData= vtkImageReader->GetOutput();

        int* dim = imageData->GetDimensions(); 

        NekDouble origin[3];
        imageData->GetOrigin(origin);

        std::cout << "Dims: " << " x: " << dim[0] << " y: " << dim[1] << " z: " << dim[2] << std::endl; 
        std::cout << "Number of points: " << imageData->GetNumberOfPoints() << std::endl; 
        std::cout << "Number of cells: " << imageData->GetNumberOfCells() << std::endl; 
        std::cout << "Origin: " << " x: " << origin[0]<< " y: " <<origin[1]<< " z: " <<origin[2]  << std::endl; 
  
        NekDouble intensity;
        NekDouble max_val=imageData->GetScalarTypeMax();
        NekDouble min_val=imageData->GetScalarTypeMin();

        Array<OneD, Array<OneD, NekDouble> >  m_spatialperm;
        m_spatialperm = Array<OneD, Array<OneD, NekDouble> > (coordim);

        for (int i = 0; i < coordim; ++i)
        {
            //Array of permeability values
            m_spatialperm[i] = Array<OneD, NekDouble>(nq);
        }

        for(int i = 0; i < nq; ++i)
        {
            NekDouble p[3];
            NekDouble whl[3];
            NekDouble pcoords[3];
            int ijk[3];
                
            //width, length and height of pixel/voxel
            imageData->GetSpacing(whl);
            
            //coordinates of the quadrature points (scaled)
            /*p[0]=xc0[i]*(dim[0]-1)*whl[0];
            p[1]=xc1[i]*(dim[1]-1)*whl[1];
            p[2]=xc2[i]*(dim[2]-1);*/

            p[0]=xc0[i]*whl[0];
            p[1]=xc1[i]*whl[1];
            p[2]=xc2[i]*whl[2];
            
            cout<<"coodinates of quadrature point: "<<p[0]<<" "<<p[1]<<" "<<p[2]<<endl;
            //this is not need (just for testing)
            //imageData->GetPoint(i,p);
            //std::cout << "Point " << i << " : (" << p[0] << " " << p[1] << " " << p[2] << ")" << std::endl;
            
            //find the closest point in the image to the quadrature point
            int iD = imageData->FindPoint(p[0], p[1], p[2]);
            
            //Get the structure coordinate location
            NekDouble closestPoint[3];
            imageData->GetPoint(iD, closestPoint);
            std::cout << "Coordinates found: " << closestPoint[0] << " " << closestPoint[1] << " " << closestPoint[2] << std::endl;
            imageData->ComputeStructuredCoordinates( closestPoint, ijk, pcoords );
            cout << "ijk: " << ijk[0] << " " << ijk[1] << " " << ijk[2] << endl;
                
            //Extract the intensity value
            intensity=imageData->GetScalarComponentAsDouble(ijk[0],ijk[1],ijk[2],0);
            cout<< "intensity value: " <<intensity<<endl;
            //intensity=intensity/max_val;

            if(intensity < 100)
            {
                for (int j = 0; j < coordim; ++j)
                {
                    m_spatialperm[j][i]=0.005;
                }
            }
            else
            {
                for (int j = 0; j < coordim; ++j)
                {
                    m_spatialperm[j][i]=1/intensity;
                }
            }
        }

        Exp->FwdTrans(m_spatialperm[0], Exp->UpdateCoeffs());
        //-----------------------------------------------
        // Write solution to file
        string   out(argv[3]);
        if (vSession->GetComm()->GetSize() > 1)
        {
            out += "." + boost::lexical_cast<string>(vSession->GetComm()->GetRank());
        }
        std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef
                                                    = Exp->GetFieldDefinitions();
        std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());

        for(int i = 0; i < FieldDef.size(); ++i)
        {
            FieldDef[i]->m_fields.push_back("kxx");
            Exp->AppendFieldData(FieldDef[i], FieldData[i]);
        }
        LibUtilities::Write(out, FieldDef, FieldData);
        //-----------------------------------------------

    }
    catch (...) {
        cout << "An error occurred." << endl;
    }
}
