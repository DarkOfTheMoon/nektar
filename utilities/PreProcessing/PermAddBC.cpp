#include <cstdio>
#include <cstdlib>
#include <string>

#include <MultiRegions/ContField3D.h>

using namespace std;

int main(int argc, char *argv[])
{
    if(argc != 2)
    {
        fprintf(stderr,"Usage: PermAddBC in_meshfile\n");
        exit(1);
    }

    //----------------------------------------------
    // Read in mesh from input file
    string meshfile(argv[argc-1]);
    SpatialDomains::MeshGraphSharedPtr graphShPt = SpatialDomains::MeshGraph::Read(meshfile);//meshfile);
    //----------------------------------------------

    //----------------------------------------------
    // Initialisation and necessary values
    int i, j, k;
    int cnt_edge = 0;
    int dim  = graphShPt->GetMeshDimension();
    ASSERTL0(dim == 3, "The geometry should be 3D.");
    //----------------------------------------------

    //----------------------------------------------
    // determining faces that need to be added to
    // existing composite

    // getting egdes of composite to compare
    SpatialDomains::Composite comp_it;
    comp_it = graphShPt->GetComposite(1);
    int nedge = 3*comp_it->size();
    cnt_edge = 0;
    Array <OneD, int> comp_edge(nedge);
    for (i=0; i<comp_it->size(); ++i)
    {
        for (j=0; j<3; ++j)
        {
            comp_edge[cnt_edge++] = (*comp_it)[i]->GetEid(j);
        }
    }

    // getting edges of composite to split
    SpatialDomains::Composite comp_prism;
    comp_prism = graphShPt->GetComposite(2);
    int nedge_prism = 3*comp_prism->size();
    cnt_edge = 0;
    Array <OneD, int> comp_prism_edge(nedge_prism);
    for (i=0; i<comp_prism->size(); ++i)
    {
        for (j=0; j<3; ++j)
        {
            comp_prism_edge[cnt_edge++] = (*comp_prism)[i]->GetEid(j);
        }
    }    

    // determining matching edges
    vector<int> matching_edges;
    for (i=0; i<nedge_prism; ++i)
    {
        for (j=0; j<nedge; ++j)
        {
            if ( comp_prism_edge[i] == comp_edge[j] )
            {
                matching_edges.push_back(comp_prism_edge[i]);
            }
        }
    }

    // determining faces of composite to splt
    // from matching edges
    vector<int> faces;
    for (i=0; i<comp_prism->size(); ++i)
    {
        for (j=0; j<matching_edges.size(); ++j)
        {
            for (k=0; k<3; ++k)
            {
                if ( (*comp_prism)[i]->GetEid(k) == matching_edges[j] )
                {
                    faces.push_back((*comp_prism)[i]->GetGlobalID());
                }
            }
        }
    }

    // output
    cout << "number of boundary faces: " << faces.size() << endl;
    for (i=0; i<faces.size(); ++i)
    {
        cout << "face " << i+1 << ": " << faces[i] << endl;
    }

    return 0;
}
