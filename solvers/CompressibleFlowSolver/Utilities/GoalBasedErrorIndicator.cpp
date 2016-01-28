#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>

#include <LibUtilities/Foundations/Interp.h>
#include <iomanip>
#include <MultiRegions/ExpList.h>
#include <MultiRegions/ContField1D.h>
#include <MultiRegions/ContField2D.h>
#include <MultiRegions/ContField3D.h>
#include <MultiRegions/ContField3DHomogeneous1D.h>
#define PI 3.14159265
//#include <tinyxml/tinyxml.h>

#include <boost/math/special_functions/fpclassify.hpp>

using namespace Nektar;

int main(int argc, char *argv[])
{
    void SetFields(
                   SpatialDomains::MeshGraphSharedPtr              &mesh,
                   vector<LibUtilities::FieldDefinitionsSharedPtr> fielddef,
                   LibUtilities::SessionReaderSharedPtr            &session,
                   Array<OneD,MultiRegions::ExpListSharedPtr>      &Exp,
                   int                                             nvariables,
                   const vector<std::string>                       &variables,
                   bool                                            homogeneous);
    
    void Writefield(
                    LibUtilities::SessionReaderSharedPtr        vSession,
                    const vector<std::string>                   &variables,
                    string                                      fieldfile,
                    SpatialDomains::MeshGraphSharedPtr          &graph,
                    Array<OneD, MultiRegions::ExpListSharedPtr> &outfield);
    
    void GetSensor(
                   Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                   Array<OneD,                   NekDouble>   &Sensor,
                   Array<OneD,                   NekDouble>   &SensorKappa,
                   NekDouble Skappa,
                   NekDouble Kappa,
                   NekDouble mu0);
    
    if (argc != 9)
    {
        fprintf(stderr, "Usage: ./FieldToField  meshfileResiduals fieldfileResiduals  "
                "meshfile1  fieldfile1\n");
        exit(1);
    }
    
    
    // exact.xml exact.fld ap.fld adjoint.xml adjoint.fld output.fld
    //----------------------------------------------
    string meshfilePlowOnPhigh(argv[argc-8]);
    string fieldfilePlowOnPhigh(argv[argc-7]);
    
    
    
    std::vector<std::string> filenamesPlowOnPhigh;
    filenamesPlowOnPhigh.push_back(meshfilePlowOnPhigh);
    LibUtilities::SessionReaderSharedPtr vSession
    = LibUtilities::SessionReader::CreateInstance(argc, argv, filenamesPlowOnPhigh);
    
    // Read in mesh from input file0
    SpatialDomains::MeshGraphSharedPtr graphShPtPlowOnPhigh
    = SpatialDomains::MeshGraph::Read(meshfilePlowOnPhigh);
    int expdim = graphShPtPlowOnPhigh->GetMeshDimension();
    
    //----------------------------------------------
    // Import fieldfilePlowOnPhigh.
    vector<LibUtilities::FieldDefinitionsSharedPtr> fielddefPlowOnPhigh;
    vector<vector<NekDouble> > fielddataPlowOnPhigh;
    LibUtilities::Import(fieldfilePlowOnPhigh,
                         fielddefPlowOnPhigh,
                         fielddataPlowOnPhigh);
    
    //----------------------------------------------
    
    
    std::string                         m_ViscosityType;
    
    NekDouble                           m_gamma;
    NekDouble                           m_pInf;
    NekDouble                           m_rhoInf;
    NekDouble                           m_uInf;
    NekDouble                           m_vInf;
    NekDouble                           m_wInf;
    NekDouble                           m_gasConstant;
    NekDouble                           m_Skappa;
    NekDouble                           m_Kappa;
    NekDouble                           m_mu0;
    NekDouble                           m_Twall;
    NekDouble                           m_mu;
    NekDouble                           m_tol1;
    NekDouble                           m_tol2;
    NekDouble                           m_tol3;
    NekDouble                           m_tol4;
    NekDouble                           m_tol5;
    NekDouble                           m_tol6;
    NekDouble                           m_tol7;
    NekDouble                           m_tol8;
    NekDouble                           m_P1;
    NekDouble                           m_P2;
    NekDouble                           m_P3;
    NekDouble                           m_P4;
    NekDouble                           m_P5;
    NekDouble                           m_P6;
    NekDouble                           m_P7;
    NekDouble                           m_P8;
    NekDouble                           m_thermalConductivity;
    
    ASSERTL0(vSession->DefinesParameter("Gamma"),
             "Compressible flow sessions must define a Gamma parameter.");
    vSession->LoadParameter("Gamma", m_gamma, 1.4);
    
    // Get E0 parameter from session file.
    ASSERTL0(vSession->DefinesParameter("pInf"),
             "Compressible flow sessions must define a pInf parameter.");
    vSession->LoadParameter("pInf", m_pInf, 101325);
    
    // Get rhoInf parameter from session file.
    ASSERTL0(vSession->DefinesParameter("rhoInf"),
             "Compressible flow sessions must define a rhoInf parameter.");
    vSession->LoadParameter("rhoInf", m_rhoInf, 1.225);
    
    // Get uInf parameter from session file.
    ASSERTL0(vSession->DefinesParameter("uInf"),
             "Compressible flow sessions must define a uInf parameter.");
    vSession->LoadParameter("uInf", m_uInf, 0.1);
    
    // Get vInf parameter from session file.
    if (expdim == 2 || expdim == 3)
    {
        ASSERTL0(vSession->DefinesParameter("vInf"),
                 "Compressible flow sessions must define a vInf parameter"
                 "for 2D/3D problems.");
        vSession->LoadParameter("vInf", m_vInf, 0.0);
    }
    
    // Get wInf parameter from session file.
    if (expdim == 3)
    {
        ASSERTL0(vSession->DefinesParameter("wInf"),
                 "Compressible flow sessions must define a wInf parameter"
                 "for 3D problems.");
        vSession->LoadParameter("wInf", m_wInf, 0.0);
    }
    vSession->LoadParameter ("Skappa",        m_Skappa,        -2.048);
    vSession->LoadParameter ("Kappa",         m_Kappa,         0.0);
    vSession->LoadParameter ("mu0",           m_mu0,           1.0);
    vSession->LoadParameter ("GasConstant",   m_gasConstant,   287.058);
    vSession->LoadParameter ("Twall",         m_Twall,         300.15);
    vSession->LoadSolverInfo("ViscosityType", m_ViscosityType, "Constant");
    vSession->LoadParameter ("mu",            m_mu,            1.78e-05);
    vSession->LoadParameter ("TOL1",          m_tol1,            0.0);
    vSession->LoadParameter ("TOL2",          m_tol2,            0.0);
    vSession->LoadParameter ("TOL3",          m_tol3,            0.0);
    vSession->LoadParameter ("TOL4",          m_tol4,            0.0);
    vSession->LoadParameter ("TOL5",          m_tol5,            0.0);
    vSession->LoadParameter ("TOL6",          m_tol6,            0.0);
    vSession->LoadParameter ("TOL7",          m_tol7,            0.0);
    vSession->LoadParameter ("TOL8",          m_tol8,            0.0);
    vSession->LoadParameter ("P1",          m_P1,            0.0);
    vSession->LoadParameter ("P2",          m_P2,            0.0);
    vSession->LoadParameter ("P3",          m_P3,            0.0);
    vSession->LoadParameter ("P4",          m_P4,            0.0);
    vSession->LoadParameter ("P5",          m_P5,            0.0);
    vSession->LoadParameter ("P6",          m_P6,            0.0);
    vSession->LoadParameter ("P7",          m_P7,            0.0);
    vSession->LoadParameter ("P8",          m_P8,            0.0);


    vSession->LoadParameter ("thermalConductivity",
                             m_thermalConductivity, 0.0257);
    
    
    
    //read info from fldfile
    // const std::vector<std::string> variables = fielddef[0]->m_fields;
    const std::vector<std::string> variables = vSession->GetVariables();
    bool homo = (fielddefPlowOnPhigh[0]->m_numHomogeneousDir > 0)? true: false;
    
    // Define Expansion
    int nfields;
    nfields = variables.size();
    Array<OneD, MultiRegions::ExpListSharedPtr> fieldsPlowOnPhigh;
    fieldsPlowOnPhigh = Array<OneD, MultiRegions::ExpListSharedPtr>(nfields);
    SetFields(graphShPtPlowOnPhigh,fielddefPlowOnPhigh, vSession,fieldsPlowOnPhigh,nfields,variables,homo);
    int nq; //pointsper plane
    
    //-----------------------------------------------
    // Copy data from file:fill fields with the fielddata
    if (fielddefPlowOnPhigh[0]->m_numHomogeneousDir == 1)
    {
        nq = fieldsPlowOnPhigh[0]->GetPlane(1)->GetTotPoints();
        
        //THE IM PHYS VALUES ARE WRONG USING bwdTrans !!!
        for (int j = 0; j < nfields; ++j)
        {
            for (int i = 0; i < fielddataPlowOnPhigh.size(); ++i)
            {
                fieldsPlowOnPhigh[j]->ExtractDataToCoeffs(fielddefPlowOnPhigh[i],
                                                          fielddataPlowOnPhigh[i],
                                                          fielddefPlowOnPhigh[i]->m_fields[j],
                                                          fieldsPlowOnPhigh[j]->UpdateCoeffs());
            }
            
            //bwd plane 0
            fieldsPlowOnPhigh[j]->GetPlane(0)->BwdTrans_IterPerExp(
                                                                   fieldsPlowOnPhigh[j]->GetPlane(0)->GetCoeffs(),
                                                                   fieldsPlowOnPhigh[j]->GetPlane(0)->UpdatePhys());
            
            //bwd plane 1
            fieldsPlowOnPhigh[j]->GetPlane(1)->BwdTrans_IterPerExp(
                                                                   fieldsPlowOnPhigh[j]->GetPlane(1)->GetCoeffs(),
                                                                   fieldsPlowOnPhigh[j]->GetPlane(1)->UpdatePhys());
        }
    }
    else
    {
        nq = fieldsPlowOnPhigh[0]->GetTotPoints();
        for (int j = 0; j < nfields; ++j)
        {
            for (int i = 0; i < fielddataPlowOnPhigh.size(); ++i)
            {
                fieldsPlowOnPhigh[j]->ExtractDataToCoeffs(fielddefPlowOnPhigh[i],
                                                          fielddataPlowOnPhigh[i],
                                                          fielddefPlowOnPhigh[i]->m_fields[j],
                                                          fieldsPlowOnPhigh[j]->UpdateCoeffs());
            }
            fieldsPlowOnPhigh[j]->BwdTrans_IterPerExp(
                                                      fieldsPlowOnPhigh[j]->GetCoeffs(),
                                                      fieldsPlowOnPhigh[j]->UpdatePhys());
        }
    }
    
    // store mesh0 quadrature points
    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> y0(nq);
    Array<OneD, NekDouble> z0(nq);
    
    if (fielddefPlowOnPhigh[0]->m_numHomogeneousDir == 1)
    {
        fieldsPlowOnPhigh[0]->GetPlane(1)->GetCoords(x0, y0, z0);
    }
    else
    {
        if (expdim == 2)
        {
            fieldsPlowOnPhigh[0]->GetCoords(x0, y0);
        }
        else if (expdim == 3)
        {
            fieldsPlowOnPhigh[0]->GetCoords(x0, y0, z0);
        }
    }
    
    string meshfilePhigh(argv[argc-6]);
    string fieldfilePhigh(argv[argc-5]);
    
    std::vector<std::string> filenamesPhigh;
    filenamesPhigh.push_back(meshfilePhigh);
    //LibUtilities::SessionReaderSharedPtr vSession
    //= LibUtilities::SessionReader::CreateInstance(argc, argv, filenamesPhigh);
    
    // Read in mesh from input file0
    SpatialDomains::MeshGraphSharedPtr graphShPtPhigh
    = SpatialDomains::MeshGraph::Read(meshfilePhigh);
    //int expdim = graphShPtPhigh->GetMeshDimension();
    
    //----------------------------------------------
    // Import fieldfilePhigh.
    vector<LibUtilities::FieldDefinitionsSharedPtr> fielddefPhigh;
    vector<vector<NekDouble> > fielddataPhigh;
    LibUtilities::Import(fieldfilePhigh,
                         fielddefPhigh,
                         fielddataPhigh);
    
    //----------------------------------------------
    
    //read info from fldfile
    // const std::vector<std::string> variables = fielddef[0]->m_fields;
    //const std::vector<std::string> variables = vSession->GetVariables();
    //bool homo = (fielddefPhigh[0]->m_numHomogeneousDir > 0)? true: false;
    
    // Define Expansion
    //int nfields;
    nfields = variables.size();
    Array<OneD, MultiRegions::ExpListSharedPtr> fieldsPhigh;
    fieldsPhigh = Array<OneD, MultiRegions::ExpListSharedPtr>(nfields);
    SetFields(graphShPtPhigh,fielddefPhigh, vSession,fieldsPhigh,nfields,variables,homo);
    //int nq; //pointsper plane
    
    //-----------------------------------------------
    // Copy data from file:fill fields with the fielddata
    if (fielddefPhigh[0]->m_numHomogeneousDir == 1)
    {
        nq = fieldsPhigh[0]->GetPlane(1)->GetTotPoints();
        
        //THE IM PHYS VALUES ARE WRONG USING bwdTrans !!!
        for (int j = 0; j < nfields; ++j)
        {
            for (int i = 0; i < fielddataPhigh.size(); ++i)
            {
                fieldsPhigh[j]->ExtractDataToCoeffs(fielddefPhigh[i],
                                                    fielddataPhigh[i],
                                                    fielddefPhigh[i]->m_fields[j],
                                                    fieldsPhigh[j]->UpdateCoeffs());
            }
            
            //bwd plane 0
            fieldsPhigh[j]->GetPlane(0)->BwdTrans_IterPerExp(
                                                             fieldsPhigh[j]->GetPlane(0)->GetCoeffs(),
                                                             fieldsPhigh[j]->GetPlane(0)->UpdatePhys());
            
            //bwd plane 1
            fieldsPhigh[j]->GetPlane(1)->BwdTrans_IterPerExp(
                                                             fieldsPhigh[j]->GetPlane(1)->GetCoeffs(),
                                                             fieldsPhigh[j]->GetPlane(1)->UpdatePhys());
        }
    }
    else
    {
        nq = fieldsPhigh[0]->GetTotPoints();
        for (int j = 0; j < nfields; ++j)
        {
            for (int i = 0; i < fielddataPhigh.size(); ++i)
            {
                fieldsPhigh[j]->ExtractDataToCoeffs(fielddefPhigh[i],
                                                    fielddataPhigh[i],
                                                    fielddefPhigh[i]->m_fields[j],
                                                    fieldsPhigh[j]->UpdateCoeffs());
            }
            fieldsPhigh[j]->BwdTrans_IterPerExp(fieldsPhigh[j]->GetCoeffs(),
                                                fieldsPhigh[j]->UpdatePhys());
        }
    }
    
    // store mesh0 quadrature points
    Array<OneD, NekDouble> x0Phigh(nq);
    Array<OneD, NekDouble> y0Phigh(nq);
    Array<OneD, NekDouble> z0Phigh(nq);
    
    if (fielddefPhigh[0]->m_numHomogeneousDir == 1)
    {
        fieldsPhigh[0]->GetPlane(1)->GetCoords(x0Phigh, y0Phigh, z0Phigh);
    }
    else
    {
        if (expdim == 2)
        {
            fieldsPhigh[0]->GetCoords(x0Phigh, y0Phigh);
        }
        else if (expdim == 3)
        {
            fieldsPhigh[0]->GetCoords(x0Phigh, y0Phigh, z0Phigh);
        }
    }
    
    Array<OneD,                   NekDouble>   Sensor(nq, 0.0);
    Array<OneD,                   NekDouble>   SensorKappa(nq, 0.0);
    
    GetSensor(fieldsPhigh, Sensor, SensorKappa, m_Skappa, m_Kappa, m_mu0);
    
    //----------------------------------------------
    //----------------------------------------------
    string meshfileAdjoint(argv[argc-4]);
    string fieldfileAdjoint(argv[argc-3]);
    
    std::vector<std::string> filenamesAdjoint;
    filenamesAdjoint.push_back(meshfileAdjoint);
    
    // Read in mesh from input file0
    SpatialDomains::MeshGraphSharedPtr graphShPtAdjoint
    = SpatialDomains::MeshGraph::Read(meshfileAdjoint);
    
    //----------------------------------------------
    // Import fieldfileAdjoint.
    vector<LibUtilities::FieldDefinitionsSharedPtr> fielddefAdjoint;
    vector<vector<NekDouble> > fielddataAdjoint;
    LibUtilities::Import(fieldfileAdjoint,
                         fielddefAdjoint,
                         fielddataAdjoint);
    
    //----------------------------------------------
    
    //read info from fldfile
    // const std::vector<std::string> variables = fielddef[0]->m_fields;
    homo = (fielddefAdjoint[0]->m_numHomogeneousDir > 0)? true: false;
    
    // Define Expansion
    Array<OneD, MultiRegions::ExpListSharedPtr> fieldsAdjoint;
    fieldsAdjoint = Array<OneD, MultiRegions::ExpListSharedPtr>(nfields);
    SetFields(graphShPtAdjoint,fielddefAdjoint, vSession,fieldsAdjoint,nfields,variables,homo);
    
    //-----------------------------------------------
    // Copy data from file:fill fields with the fielddata
    if (fielddefAdjoint[0]->m_numHomogeneousDir == 1)
    {
        nq = fieldsAdjoint[0]->GetPlane(1)->GetTotPoints();
        
        for (int j = 0; j < nfields; ++j)
        {
            for (int i = 0; i < fielddataAdjoint.size(); ++i)
            {
                fieldsAdjoint[j]->ExtractDataToCoeffs(fielddefAdjoint[i],
                                                      fielddataAdjoint[i],
                                                      fielddefAdjoint[i]->m_fields[j],
                                                      fieldsAdjoint[j]->UpdateCoeffs());
            }
            
            //bwd plane 0
            fieldsAdjoint[j]->GetPlane(0)->BwdTrans_IterPerExp(
                                                               fieldsAdjoint[j]->GetPlane(0)->GetCoeffs(),
                                                               fieldsAdjoint[j]->GetPlane(0)->UpdatePhys());
            
            //bwd plane 1
            fieldsAdjoint[j]->GetPlane(1)->BwdTrans_IterPerExp(
                                    fieldsAdjoint[j]->GetPlane(1)->GetCoeffs(),
                                    fieldsAdjoint[j]->GetPlane(1)->UpdatePhys());
        }
    }
    else
    {
        nq = fieldsAdjoint[0]->GetTotPoints();
        for (int j = 0; j < nfields; ++j)
        {
            for (int i = 0; i < fielddataAdjoint.size(); ++i)
            {
                fieldsAdjoint[j]->ExtractDataToCoeffs(fielddefAdjoint[i],
                                                      fielddataAdjoint[i],
                                                      fielddefAdjoint[i]->m_fields[j],
                                                       fieldsAdjoint[j]->UpdateCoeffs());
            }
            fieldsAdjoint[j]->BwdTrans_IterPerExp(fieldsAdjoint[j]->GetCoeffs(),
                                                  fieldsAdjoint[j]->UpdatePhys());
        }
    }
    
    string fieldfileOutPut(argv[argc-2]);
    
    Array<OneD, MultiRegions::ExpListSharedPtr> fieldsOutPut;
    fieldsOutPut = Array<OneD, MultiRegions::ExpListSharedPtr>(nfields);
    
    SetFields(graphShPtPhigh, fielddefPhigh, vSession, fieldsOutPut, nfields,
              variables, homo);
    
    string fieldfileResiduals(argv[argc-1]);
    
    Array<OneD, MultiRegions::ExpListSharedPtr> fieldsResiduals;
    fieldsResiduals = Array<OneD, MultiRegions::ExpListSharedPtr>(nfields);
    
    SetFields(graphShPtPhigh, fielddefPhigh, vSession, fieldsResiduals, nfields,
              variables, homo);
    
    
    for (int j = 0; j < nfields; ++j)
    {
        
        Vmath::Vsub(nq,
                    fieldsPhigh[j]->GetPhys(), 1,
                    fieldsPlowOnPhigh[j]->GetPhys(), 1,
                    fieldsResiduals[j]->UpdatePhys(), 1);
        
        Vmath::Vabs(nq,
                    fieldsResiduals[j]->GetPhys(), 1,
                    fieldsResiduals[j]->UpdatePhys(), 1);
    }
    
    NekDouble WeightedRes = 0.0;
    
    for (int j = 0; j < nfields; ++j)
    {
        
        Vmath::Vabs(nq,
                    fieldsAdjoint[j]->GetPhys(), 1,
                    fieldsAdjoint[j]->UpdatePhys(), 1);
        
        Vmath::Vmul(nq,
                    fieldsResiduals[j]->GetPhys(), 1,
                    fieldsAdjoint[j]->GetPhys(), 1,
                    fieldsOutPut[j]->UpdatePhys(), 1);
        
        WeightedRes += Vmath::Vsum(nq, fieldsOutPut[j]->GetPhys(), 1);
    }
    
    
    int nElements       = fieldsOutPut[0]->GetExpSize();
    int offset          = 0;
    
    int print = 0;
    std::ofstream m_file( "VariablePComposites.txt", std::ios_base::app);
    for (int e = 0; e < nElements; e++)
    {
        m_file << "<C ID=\"" << e+1 << "\"> T[" << e << "] </C>"<< endl;
    }
    m_file.close();
    
    std::ofstream m_file2( "VariablePExpansions.txt", std::ios_base::app);
    
    for (int e = 0; e < nElements; e++)
    {
        
        int nqel      = fieldsOutPut[0]->GetExp(e)->GetTotPoints();
        int ncoeffsel = fieldsOutPut[0]->GetExp(e)->GetNcoeffs();
        
        NekDouble sum0 = Vmath::Vsum(nqel, &fieldsOutPut[0]->GetPhys()[offset], 1);
        NekDouble sum1 = Vmath::Vsum(nqel, &fieldsOutPut[1]->GetPhys()[offset], 1);
        NekDouble sum2 = Vmath::Vsum(nqel, &fieldsOutPut[2]->GetPhys()[offset], 1);
        NekDouble sum3 = Vmath::Vsum(nqel, &fieldsOutPut[3]->GetPhys()[offset], 1);
        
        NekDouble sum = abs(sum0)+abs(sum1)+abs(sum2)+abs(sum3);
        //sum = abs(sum0)+abs(sum1)+abs(sum2)+abs(sum3);
        Array<OneD, NekDouble> error_sum(nqel, sum);
        
        //cout << Vmath::Vmax(nqel, &fieldsOutPut[0]->UpdatePhys()[offset], 1) << endl;
        
        Vmath::Vcopy(nqel,
                     &error_sum[0], 1,
                     &fieldsOutPut[0]->UpdatePhys()[offset], 1);
        
        //cout << abs(sum) << endl;
        
        print = 3;
        
        /*
        if (abs(sum) >= tol7 && abs(sum) < tol6)
        {
         print = 4;
         }
         if (abs(sum) >= tol6 && abs(sum) < tol5)
         {
         print = 5;
         }
         if (abs(sum) >= tol5 && abs(sum) < tol4)
         {
         print = 6;
         }
         if (abs(sum) >= tol4 && abs(sum) < tol3)
         {
         print = 7;
        }*/
        //cout << m_tol1 << " " << m_tol2 << " "  << m_tol3 << " "  << m_tol4 << endl;
        /*if (abs(sum) >= m_tol8 && abs(sum) < m_tol7)
        {
            print = m_P8;
        }
        if (abs(sum) >= m_tol7 && abs(sum) < m_tol6)
        {
            print = m_P7;
        }
        if (abs(sum) >= m_tol6 && abs(sum) < m_tol5)
        {
            print = m_P6;
        }
        if (abs(sum) >= m_tol5 && abs(sum) < m_tol4)
        {
            print = m_P5;
        }
        if (abs(sum) >= m_tol4 && abs(sum) < m_tol3)
        {
            print = m_P4;
        }*/
        print = 3;
        /*if (abs(sum) >= m_tol3 && abs(sum) < m_tol2)
        {
            print = m_P3;
        }*/
        if (abs(sum) >= m_tol2 && abs(sum) <= m_tol1)
        {
            print = m_P2;
        }
        if (abs(sum) > m_tol1)
        {
            print = m_P1;
        }
        //if (Sensor[offset] > m_Skappa)
        //{
        //    print = 4;
        //}
        
        Array<OneD, NekDouble> p_order(nqel, print-1);
        
        Vmath::Vcopy(nqel,
                     &p_order[0], 1,
                     &fieldsOutPut[3]->UpdatePhys()[offset], 1);
        
        m_file2 << "<E COMPOSITE= \"C[" << e+1
        << "]\" NUMMODES=\"" << print
        << "\" TYPE=\"MODIFIED\" FIELDS=\"rho,rhou,rhov,E\" />"
        << endl;
        
        offset += nqel;
    }
    
    m_file2.close();
    
    Writefield(vSession, variables, fieldfileOutPut, graphShPtAdjoint, fieldsOutPut);
    
    Writefield(vSession, variables, fieldfileResiduals, graphShPtAdjoint, fieldsResiduals);
    
    NekDouble dynpres = 0.5*m_rhoInf*(m_uInf*m_uInf+m_vInf*m_vInf);
    cout <<
    cout << " ========================================= " << endl;
    cout << " ======== dF_{l/d} = " << WeightedRes << " =========" <<  endl;
    cout << " ======== dC_{l/d} = " << WeightedRes/dynpres << " =========" <<  endl;
    cout << " ========================================= " << endl;
}



// Define Expansion
void SetFields(
               SpatialDomains::MeshGraphSharedPtr              &mesh,
               vector<LibUtilities::FieldDefinitionsSharedPtr> fielddef,
               LibUtilities::SessionReaderSharedPtr            &session,
               Array<OneD,MultiRegions::ExpListSharedPtr>      &Exp,
               int                                             nvariables,
               const vector<std::string>                       &variables,
               bool                                            homogeneous)
{
    // Session reader stuff has to be evaluated only from the
    // first session which refers to mesh0
    static int cnt = 0;
    
    // Setting parameteres for homogenous problems
    NekDouble static LhomY;///< physical length in Y direction (if homogeneous)
    NekDouble static LhomZ;///< physical length in Z direction (if homogeneous)
    
    bool static DeclareCoeffPhysArrays = true;
    
    int static npointsY;///< number of points in Y direction (if homogeneous)
    int static npointsZ;///< number of points in Z direction (if homogeneous)
    
    bool static useFFT = false;
    bool static deal   = false;
    
    cnt++;
    int i;
    int expdim = mesh->GetMeshDimension();
    
    switch (expdim)
    {
        case 1:
        {
            if (fielddef[0]->m_numHomogeneousDir == 1)
            {
                const LibUtilities::PointsKey PkeyY(
                                                    npointsY,
                                                    LibUtilities::eFourierEvenlySpaced);
                const LibUtilities::BasisKey  BkeyY(
                                                    LibUtilities::eFourier, npointsY, PkeyY);
                
                for (i = 0 ; i < nvariables; i++)
                {
                    Exp[i] = MemoryManager<MultiRegions::
                    ContField3DHomogeneous1D>::AllocateSharedPtr(
                                                                 session, BkeyY, LhomY,
                                                                 useFFT, deal, mesh,
                                                                 variables[i]);
                }
            }
            else
            {
                for (i = 0 ; i < nvariables; i++)
                {
                    Exp[i] = MemoryManager<MultiRegions::ContField1D>
                    ::AllocateSharedPtr(session, mesh, variables[i]);
                }
            }
            break;
        }
        case 2:
        {
            if (fielddef[0]->m_numHomogeneousDir == 1)
            {
                const LibUtilities::PointsKey PkeyZ(
                                                    npointsZ,
                                                    LibUtilities::eFourierEvenlySpaced);
                const LibUtilities::BasisKey  BkeyZ(
                                                    LibUtilities::eFourier, npointsZ, PkeyZ);
                
                MultiRegions::ContField3DHomogeneous1DSharedPtr firstfield =
                MemoryManager<MultiRegions::ContField3DHomogeneous1D>::
                AllocateSharedPtr(session, BkeyZ, LhomZ, useFFT,
                                  deal, mesh, variables[0]);
                
                Exp[0] = firstfield;
                for (i = 1 ; i < nvariables; i++)
                {
                    Exp[i] = MemoryManager<MultiRegions::
                    ContField3DHomogeneous1D>::
                    AllocateSharedPtr(*firstfield);
                }
            }
            else
            {
                i = 0;
                MultiRegions::DisContField2DSharedPtr firstfield;
                firstfield = MemoryManager<MultiRegions::DisContField2D>::
                AllocateSharedPtr(session, mesh, variables[i],
                                  DeclareCoeffPhysArrays);
                
                Exp[0] = firstfield;
                for (i = 1 ; i < nvariables; i++)
                {
                    Exp[i] = MemoryManager<MultiRegions::DisContField2D>
                    ::AllocateSharedPtr(*firstfield, mesh, variables[i],
                                        DeclareCoeffPhysArrays);
                }
            }
            break;
        }
        case 3:
        {
            if (fielddef[0]->m_numHomogeneousDir == 1)
            {
                ASSERTL0(false,
                         "3D fully periodic problems not implemented yet");
            }
            else
            {
                i = 0;
                MultiRegions::ExpList3DSharedPtr firstfield =
                MemoryManager<MultiRegions::ExpList3D>::
                AllocateSharedPtr(session, mesh);
                
                Exp[0] = firstfield;
                for (i = 1 ; i < nvariables; i++)
                {
                    Exp[i] = MemoryManager<MultiRegions::ExpList3D>
                    ::AllocateSharedPtr(*firstfield);
                }
            }
            break;
        }
        default:
            ASSERTL0(false, "Expansion dimension not recognised");
            break;
    }
    
}

void GetSensor(
               Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
               Array<OneD,                   NekDouble>   &Sensor,
               Array<OneD,                   NekDouble>   &SensorKappa,
               NekDouble Skappa,
               NekDouble Kappa,
               NekDouble mu0)
{
    
    int e, NumModesElement, nQuadPointsElement;
    int nTotQuadPoints  = fields[0]->GetTotPoints();
    int nElements       = fields[0]->GetExpSize();
    int nfields = fields.num_elements();
    Array<OneD, Array<OneD, NekDouble> > tmp(nfields);
    
    for (int i = 0; i < nfields; ++i)
    {
        tmp[i] = fields[i]->GetPhys();
    }
    
    // Find solution (SolP) at p = P;
    // The input array (physarray) is the solution at p = P;
    
    Array<OneD, NekDouble> SolP(nTotQuadPoints,0.0);
    Array<OneD, NekDouble> SolPmOne(nTotQuadPoints,0.0);
    Array<OneD, NekDouble> SolNorm(nTotQuadPoints,0.0);
    Array<OneD, LibUtilities::BasisSharedPtr> base;
    Vmath::Vcopy(nTotQuadPoints,tmp[0],1,SolP,1);
    
    int CoeffsCount = 0;
    
    for (e = 0; e < nElements; e++)
    {
        base         = fields[0]->GetExp(e)->GetBase();
        
        NumModesElement = base[0]->GetNumModes();
        
        int nQuadPointsElement = fields[0]->GetExp(e)->GetTotPoints();
        int nCoeffsElement = fields[0]->GetExp(e)->GetNcoeffs();
        int numCutOff = NumModesElement - 1;
        
        // Set-up of the Orthogonal basis for a Quadrilateral element which is
        // needed to obtain thesolution at P =  p - 1;
        
        Array<OneD, NekDouble> SolPElementPhys(nQuadPointsElement,0.0);
        Array<OneD, NekDouble> SolPElementCoeffs(nCoeffsElement,0.0);
        
        Array<OneD, NekDouble> SolPmOneElementPhys(nQuadPointsElement,0.0);
        Array<OneD, NekDouble> SolPmOneElementCoeffs(nCoeffsElement,0.0);
        
        // create vector the save the solution points per element at P = p;
        
        for (int i = 0; i < nQuadPointsElement; i++)
        {
            SolPElementPhys[i] = SolP[CoeffsCount+i];
        }
        
        fields[0]->GetExp(e)->FwdTrans(SolPElementPhys,
                                       SolPElementCoeffs);
        
        // ReduceOrderCoeffs reduces the polynomial order of the solution that
        // is represented by the coeffs given as an inarray. This is done by
        // projecting the higher order solution onto the orthogonal basis and
        // padding the higher order coefficients with zeros.
        
        fields[0]->GetExp(e)->ReduceOrderCoeffs(numCutOff,
                                                SolPElementCoeffs,
                                                SolPmOneElementCoeffs);
        
        fields[0]->GetExp(e)->BwdTrans(SolPmOneElementCoeffs,
                                       SolPmOneElementPhys);
        
        for (int i = 0; i < nQuadPointsElement; i++)
        {
            SolPmOne[CoeffsCount+i] = SolPmOneElementPhys[i];
        }
        
        NekDouble SolPmeanNumerator = 0.0;
        NekDouble SolPmeanDenumerator = 0.0;
        
        // Determining the norm of the numerator of the Sensor
        
        Vmath::Vsub(nQuadPointsElement,
                    SolPElementPhys, 1,
                    SolPmOneElementPhys, 1,
                    SolNorm, 1);
        
        Vmath::Vmul(nQuadPointsElement,
                    SolNorm, 1,
                    SolNorm, 1,
                    SolNorm, 1);
        
        for (int i = 0; i < nQuadPointsElement; i++)
        {
            SolPmeanNumerator   += SolNorm[i];
            SolPmeanDenumerator += SolPElementPhys[i];
        }
        
        for (int i = 0; i < nQuadPointsElement; ++i)
        {
            Sensor[CoeffsCount+i] = sqrt(SolPmeanNumerator/nQuadPointsElement)
            /sqrt(SolPmeanDenumerator/nQuadPointsElement);
            
            Sensor[CoeffsCount+i] = log10(Sensor[CoeffsCount+i]);
        }
        CoeffsCount += nQuadPointsElement;
    }
    
    CoeffsCount = 0.0;
    
    for (e = 0; e < nElements; e++)
    {
        base         = fields[0]->GetExp(e)->GetBase();
        
        NumModesElement = base[0]->GetNumModes();
        
        NekDouble ThetaS        = mu0;
        NekDouble Phi0          = Skappa;
        NekDouble DeltaPhi      = Kappa;
        nQuadPointsElement      = fields[0]->GetExp(e)->GetTotPoints();
        
        for (int i = 0; i < nQuadPointsElement; i++)
        {
            if (Sensor[CoeffsCount+i] <= (Phi0 - DeltaPhi))
            {
                SensorKappa[CoeffsCount+i] = 0;
            }
            else if(Sensor[CoeffsCount+i] >= (Phi0 + DeltaPhi))
            {
                SensorKappa[CoeffsCount+i] = ThetaS;
            }
            else if(abs(Sensor[CoeffsCount+i]-Phi0) < DeltaPhi)
            {
                SensorKappa[CoeffsCount+i] = ThetaS/2*(1+sin(M_PI*
                                        (Sensor[CoeffsCount+i]-Phi0)
                                        /(2*DeltaPhi)));
            }
        }
        
        CoeffsCount += nQuadPointsElement;
    }
    
}

void Writefield(
                LibUtilities::SessionReaderSharedPtr        vSession,
                const std::vector<std::string>              &variables,
                string                                      fieldfile,
                SpatialDomains::MeshGraphSharedPtr          &graph,
                Array<OneD, MultiRegions::ExpListSharedPtr> &outfield)
{
    string var;
    std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef
    = outfield[0]->GetFieldDefinitions();
    std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());
    Array<OneD, Array<OneD, NekDouble> > fieldcoeffs(outfield.num_elements());
    
    for (int j = 0; j < fieldcoeffs.num_elements(); ++j)
    {
        if (FieldDef[0]->m_numHomogeneousDir == 1)
        {
            // plane 0
            outfield[j]->GetPlane(0)->FwdTrans_IterPerExp(
                                                          outfield[j]->GetPlane(0)->GetPhys(),
                                                          outfield[j]->GetPlane(0)->UpdateCoeffs());
            
            // plane 1
            outfield[j]->GetPlane(1)->FwdTrans_IterPerExp(
                                                          outfield[j]->GetPlane(1)->GetPhys(),
                                                          outfield[j]->GetPlane(1)->UpdateCoeffs());
        }
        else
        {  
            outfield[j]->FwdTrans_IterPerExp(
                                             outfield[j]->GetPhys(), 
                                             outfield[j]->UpdateCoeffs());
        }
        
        fieldcoeffs[j] = outfield[j]->UpdateCoeffs();	
        
        for (int i = 0; i < FieldDef.size(); i++)
        {		     	     
            FieldDef[i]->m_fields.push_back(variables[j]);   
            outfield[0]->AppendFieldData(FieldDef[i], FieldData[i], 
                                         fieldcoeffs[j]);  
        }
    }
    LibUtilities::Write(fieldfile, FieldDef, FieldData);
}    		

