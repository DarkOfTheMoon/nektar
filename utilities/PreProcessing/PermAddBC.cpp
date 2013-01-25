#include <cstdio>
#include <cstdlib>
#include <string>
#include <sstream>

#include <MultiRegions/ContField3D.h>

#include <tinyxml/tinyxml.h>

int main(int argc, char *argv[])
{
    if(argc != 3)
    {
        fprintf(stderr,"Usage: PermAddBC in_meshfile out_meshfile\n");
        exit(1);
    }

    //----------------------------------------------
    // Read in mesh from input file
    string meshfile(argv[argc-2]);
    SpatialDomains::MeshGraphSharedPtr graphShPt = SpatialDomains::MeshGraph::Read(meshfile);

    TiXmlDocument doc(argv[argc-2]);
    bool loadOkay = doc.LoadFile();
    ASSERTL0(loadOkay, "Unable to load input meshfile.");
    //----------------------------------------------

    //----------------------------------------------
    // Initialisation and necessary values
    int n, i, j, k, nedge, nedge_prism, lc_id, err;
    vector<int> comp_edge, comp_prism_edge, matching_edges, faces;
    stringstream body;

    int dim  = graphShPt->GetMeshDimension();
    ASSERTL0(dim == 3, "The geometry should be 3D.");

    Array <OneD, vector<NekDouble> > extend(dim);
    Array <OneD, Array <OneD, NekDouble> > tmp_Vcoords(dim);
    for (i=0; i<dim; ++i)
    {
        tmp_Vcoords[i] = Array <OneD, NekDouble>(3);
    }
    NekDouble tol = 1e-06;

    TiXmlElement* comps = doc.FirstChildElement("NEKTAR")
        ->FirstChildElement("GEOMETRY")
        ->FirstChildElement("COMPOSITE");
    ASSERTL0(comps, "Tag COMPOSITE could not be found in meshfile.");
    
    TiXmlElement* last_comp = comps->LastChild("C")->ToElement();
    TiXmlComment* comment;
    TiXmlElement* newC;
    TiXmlText* body_text;

    SpatialDomains::Composite comp_it, comp_prism;
    //----------------------------------------------

    //----------------------------------------------
    // determining faces that need to be added to
    // existing composite

    // looping over all 6 walls
    // assumming:
    // C[0]:   tets
    // C[1]:   prisms
    // C[2-7]: outer walls
    // C[8]:   inner walls (SMC walls)
    // C[9]:   all faces, that have to be assigned
    for (n=2; n<8; ++n)
    {
        // getting egdes of composite to compare
        comp_edge.clear();
        comp_it = graphShPt->GetComposite(n);
        nedge = 3*comp_it->size();
        for (i=0; i<comp_it->size(); ++i)
        {
            for (j=0; j<3; ++j)
            {
                comp_edge.push_back((*comp_it)[i]->GetEid(j));
            }
        }

        // getting edges of composite to split
        comp_prism = graphShPt->GetComposite(9);
        nedge_prism = 4*comp_prism->size();
        comp_prism_edge.clear();
        for (i=0; i<comp_prism->size(); ++i)
        {
            for (j=0; j<4; ++j)
            {
                comp_prism_edge.push_back((*comp_prism)[i]->GetEid(j));
            }
        }    

        // determining matching edges
        matching_edges.clear();
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
        faces.clear();
        for (i=0; i<comp_prism->size(); ++i)
        {
            for (j=0; j<matching_edges.size(); ++j)
            {
                for (k=0; k<4; ++k)
                {
                    if ( (*comp_prism)[i]->GetEid(k) == matching_edges[j] )
                    {
                        faces.push_back((*comp_prism)[i]->GetGlobalID());
                    }
                }
            }
        }
        //----------------------------------------------
    

        // output
        cout << "number of boundary faces: " << faces.size() << endl;
        for (i=0; i<faces.size(); ++i)
        {
            cout << "face " << i+1 << ": " << faces[i] << endl;
        }

        //----------------------------------------------
        // writing out identified faces to a new compound

        // Only create new composite if matching faces
        // have been found
        if (faces.size() > 0)
        {
            // Getting ID of the last composite
            lc_id = atoi(last_comp->Attribute("ID"));

            // Creating comment
            comment = new TiXmlComment();
            body.str("");
            body << "Following Composite is associated to Composite with ID:" << n;
            comment->SetValue(body.str());  
            comps->LinkEndChild(comment);

            // Creating the body of the new compound
            body.str("");
            body << "F[";
            for (i=0; i<faces.size(); ++i)
            {
                body << faces[i];
                if (i < faces.size()-1)
                {
                    body << ",";
                }
            }
            body << "]";

            // Creating new composite entry
            newC = new TiXmlElement("C");
            body_text = new TiXmlText(body.str());
            newC->LinkEndChild(body_text);
            comps->LinkEndChild(newC);

            // Setting correct ID of new composite
            last_comp = last_comp->NextSiblingElement();
            last_comp->SetAttribute("ID", lc_id+1);
        }
        //----------------------------------------------

/*      //----------------------------------------------
        // calculating extend of "cube"
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
        //----------------------------------------------
*/
    }

/*  //----------------------------------------------
    // Writing extend of "cube" to function "extend"
    TiXmlElement* cond = doc.FirstChildElement("NEKTAR")
        ->FirstChildElement("CONDITIONS");
    
    TiXmlElement* tmp = new TiXmlElement("FUNCTION");
    cond->LinkEndChild(tmp);
    TiXmlElement* extend_func = cond->FirstChildElement("FUNCTION");
    extend_func->SetAttribute("NAME", "Extend");

    string extend_name[6] = {
        "x_min",
        "x_max",
        "y_min",
        "y_max",
        "z_min",
        "z_max"
    };

    //TiXmlElement* tmp;
    for (i=0; i<3; ++i)
    {
        tmp = new TiXmlElement("E");
        extend_func->LinkEndChild(tmp);
        tmp = new TiXmlElement("E");
        extend_func->LinkEndChild(tmp);
        
        tmp = extend_func->LastChild()->ToElement();
        tmp->SetAttribute("VAR", extend_name[2*i+1]);
        tmp = tmp->PreviousSibling()->ToElement();
        tmp->SetAttribute("VAR", extend_name[2*i]);

        if (extend[i][0] > extend[i][1])
        {
            tmp->SetDoubleAttribute("VALUE", extend[i][1]);
            tmp = tmp->NextSiblingElement();
            tmp->SetDoubleAttribute("VALUE", extend[i][0]);
        }
        else
        {
            tmp->SetDoubleAttribute("VALUE", extend[i][0]);
            tmp = tmp->NextSiblingElement();
            tmp->SetDoubleAttribute("VALUE", extend[i][1]);
        }
    }
    //----------------------------------------------
*/    
    //----------------------------------------------
    // Writing out original meshfile and changes to
    // new meshfile
    doc.SaveFile(argv[argc-1]);
    //----------------------------------------------

    // todo:
    // (- adding the additional composites to BCs)
    // - adding calculation of the simplyfied domain size
    //   with more digits...

    return 0;
}

