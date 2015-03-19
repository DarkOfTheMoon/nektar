///////////////////////////////////////////////////////////////////////////////
//
// File SessionReader.cpp
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
// Description: Session reader
//
///////////////////////////////////////////////////////////////////////////////

#ifndef TIXML_USE_STL
#define TIXML_USE_STL
#endif

#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicConst/GitRevision.h>

#include <iostream>
#include <fstream>
#include <string>
using namespace std;

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/algorithm/string.hpp>
#include <tinyxml.h>
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/BasicUtils/Equation.h>
#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <LibUtilities/BasicUtils/ParseUtils.hpp>
#include <LibUtilities/BasicUtils/Thread.h>
#include <LibUtilities/Communication/ThreadedComm.h>
#include <LibUtilities/BasicUtils/FileSystem.h>

#include <boost/program_options.hpp>
#include <boost/format.hpp>

#ifndef NEKTAR_VERSION
#define NEKTAR_VERSION "Unknown"
#endif

namespace po = boost::program_options;
namespace io = boost::iostreams;

namespace Nektar
{
    namespace LibUtilities
    {
        /**
         * @class SessionReader
         *
         * This class provides an interface to Nektar++-specific content in a
         * supplied XML document. It also initialises a Nektar++ session
         * including setting up communication for parallel execution and where
         * necessary partitioning the supplied mesh for running across multiple
         * processes.
         *
         * A session should be initialised at the beginning of a user's
         * application by passing the command-line arguments. This not only
         * allows the SessionReader to extract the name of the XML document to
         * load containing Nektar++ session information, but also supplies the
         * MPI arguments necessary for setting up parallel communication. The
         * SessionReader should be initialised using the #CreateInstance
         * function:
         * @code
         * LibUtilities::SessionReaderSharedPtr vSession
         *          = LibUtilities::SessionReader::CreateInstance(argc, argv);
         * @endcode
         * The instance \c vSession can now be passed to other key Nektar++
         * components during their construction.
         * @note At the end of the user application, it is important to call the
         * #Finalise routine in order to finalise any MPI communication and
         * correctly free resources.
         *
         * The SessionReader class provides streamlined, validated access to
         * session parameters, solver information and functions defined within a
         * Nektar++ XML document. The available routines and their usage is
         * documented below.
         *
         * In the case of solver information properties, the classes to which
         * these parameters are pertinent may register with the SessionReader
         * class the set of valid values for a given property. Such values may
         * also be associated with an enumeration value for more transparent use
         * of the property values in code.
         */

        /**
         * This map of maps stores the list of valid string values for a number
         * of solver information parameters. The top level map connects
         * different parameter names to their list of possible values. The list
         * of possible values is also a map, mapping a valid string to a
         * corresponding enum value.
         *
         * This list is populated through the #RegisterEnumValue static member
         * function which is called statically from various classes to register
         * the valid values for solver info parameters associated with them. The
         * map is therefore fully populated before the SessionReader class is
         * instantiated and a file is read in and parsed.
         */
        EnumMapList& SessionReader::GetSolverInfoEnums()
        {
            static EnumMapList solverInfoEnums;
            return solverInfoEnums;
        }


        /**
         * List of default values for solver information parameters to be used
         * in the case of them not being provided.
         *
         * This list is populated through the #RegisterDefaultSolverInfo static
         * member variable which is called statically from various classes to
         * register the default value for a given parameter.
         */
        SolverInfoMap& SessionReader::GetSolverInfoDefaults()
        {
            static SolverInfoMap solverInfoMap;
            return solverInfoMap;
        }


        /**
         * List of values for GlobalSysSoln parameters to be used to override
         * details given in SolverInfo
         *
         * This list is populated by ReadGlobalSysSoln if the
         * GLOBALSYSSOLNINFO section is defined in the input file.
         * This List allows for details to define for the Global Sys
         * solver for each variable. 
         */
        GloSysSolnInfoList& SessionReader::GetGloSysSolnList()
        {
            static GloSysSolnInfoList gloSysSolnInfoList;
            return gloSysSolnInfoList;
        }

        /**
         * Lists the possible command-line argument which can be specified for
         * this executable.
         *
         * This list is populated through the #RegisterCmdLineArgument static
         * member function which is called statically from various classes to
         * register command-line arguments they need.
         */
        CmdLineArgMap& SessionReader::GetCmdLineArgMap()
        {
            static CmdLineArgMap cmdLineArguments;
            return cmdLineArguments;
        }


        /**
         * This constructor parses the command-line arguments given to the user
         * application to set up any MPI communication, read supplied XML
         * session files, and partition meshes where necessary.
         *
         * @param   argc        Number of command-line arguments
         * @param   argv        Array of command-line arguments
         * @param   pMainFunc   The "main" part of the application using the library.
         */
        SessionReader::SessionReader(int argc, char *argv[],
            void (*pMainFunc)(SessionReaderSharedPtr)) :
        		m_xmlDoc(1), m_parameters(1), m_solverInfo(1), m_geometricInfo(1),
        		m_expressions(1), m_exprEvaluator(1), m_functions(1), m_variables(1), m_tags(1),
        		m_filters(1), m_mainFunc(pMainFunc), m_bndRegOrder(1)
        {
            m_xmlDoc[0] = 0;
            m_filenames = ParseCommandLineArguments(argc, argv);

            ASSERTL0(m_filenames.size() > 0, "No session file(s) given.");

            m_sessionName = ParseSessionName(m_filenames);
            m_exprEvaluator[0] = new AnalyticExpressionEvaluator();

            Nektar::InitMemoryPools(1, false);

            // Create communicator
            CreateComm(argc, argv);

            // If running in parallel change the default global sys solution
            // type.
            if (m_comm->GetSize() > 1)
            {
                GetSolverInfoDefaults()["GLOBALSYSSOLN"] = 
                    "IterativeStaticCond";
            }
        }


        /**
         *
         */
        SessionReader::SessionReader(
            int                             argc, 
            char                           *argv[], 
            const std::vector<std::string> &pFilenames, 
            const CommSharedPtr            &pComm,
            void (*pMainFunc)(SessionReaderSharedPtr)) :
    			m_xmlDoc(1), m_parameters(1), m_solverInfo(1), m_geometricInfo(1),
    			m_expressions(1), m_exprEvaluator(1), m_functions(1), m_variables(1), m_tags(1),
        		m_filters(1), m_mainFunc(pMainFunc), m_bndRegOrder(1)
        {
            ASSERTL0(pFilenames.size() > 0, "No filenames specified.");

            ParseCommandLineArguments(argc, argv);
            m_xmlDoc[0]   = 0;
            m_filenames   = pFilenames;

            m_sessionName = ParseSessionName(m_filenames);
            m_exprEvaluator[0] = new AnalyticExpressionEvaluator();

            Nektar::InitMemoryPools(1, false);

            // Create communicator
            if (!pComm.get())
            {
                CreateComm(argc, argv);
            }
            else
            {
                m_comm = pComm;

            }

            // If running in parallel change the default global sys solution
            // type.
            if (m_comm->GetSize() > 1)
            {
                GetSolverInfoDefaults()["GLOBALSYSSOLN"] = 
                    "IterativeStaticCond";
            }
        }


        /**
         *
         */
        SessionReader::~SessionReader()
        {
        	unsigned int vNumW = m_threadManager->GetMaxNumWorkers();
        	for (unsigned int i=0; i < vNumW; i++)
        	{
        		delete m_exprEvaluator[i];
                delete m_xmlDoc[i];
        	}
        }

        SessionReader::SessionJob::SessionJob(SessionReaderSharedPtr psession) :
					m_session(psession)
        {
        	// empty
        }
        SessionReader::SessionJob::~SessionJob()
        {
        	// empty
        }
        void SessionReader::SessionJob::Run()
        {
            unsigned int vThr = GetWorkerNum();

            // Spin until all threads are running and main thread has
            // entered Wait().  If we don't do this then calls to Hold()
            // will interfere with Wait()'s resizing the worker pool.
            while (m_session->m_threadManager->GetNumWorkers() !=
                    m_session->m_threadManager->GetMaxNumWorkers() )
            {
            }

        	// Partition mesh
        	m_session->PartitionMesh();

        	// Parse the XML data in #m_xmlDoc
        	m_session->ParseDocument();
            m_session->m_exprEvaluator[vThr]->
                        SetRandomSeed((m_session->m_comm->GetRank() + 1) * time(NULL));

        	// Override SOLVERINFO and parameters with any specified on the
        	// command line.
        	m_session->CmdLineOverride();

            // In verbose mode, print out parameters and solver info sections
            if (m_session->m_verbose && m_session->m_comm)
            {
                if (m_session->m_comm->GetRank() == 0 &&
                        m_session->m_parameters[vThr].size() > 0)
                {
                    cout << "Parameters:" << endl;
                    ParameterMap::iterator x;
                    for (x = m_session->m_parameters[vThr].begin();
                            x != m_session->m_parameters[vThr].end(); ++x)
                    {
                        cout << "\t" << x->first << " = " << x->second << endl;
                    }
                    cout << endl;
                }

                if (m_session->m_comm->GetRank() == 0 &&
                        m_session->m_solverInfo[vThr].size() > 0)
                {
                    cout << "Solver Info:" << endl;
                    SolverInfoMap::iterator x;
                    for (x = m_session->m_solverInfo[vThr].begin();
                            x != m_session->m_solverInfo[vThr].end(); ++x)
                    {
                        cout << "\t" << x->first << " = " << x->second << endl;
                    }
                    cout << endl;
                }
            }

            // Run the main function passed in.  If it's not been set
            // then we just return.  Control will be returned to whatever
            // just constructed the SessionReader which will run _serially_.
            // A check has already been made that if nthreads > 1 then
            // there is a mainFunc set.
            if (m_session->m_mainFunc) m_session->m_mainFunc(m_session);

        }

        /**
         * Performs the main initialisation of the object. The XML file provided
         * on the command-line is loaded and any mesh partitioning is done. The
         * resulting process-specific XML file (containing the process's
         * geometry partition) is then reloaded and parsed.
         */
        void SessionReader::InitSession()
        {
            // Session not yet loaded, so load parameters section
            m_xmlDoc[0] = MergeDoc(m_filenames);
            // Check we actually have a document loaded.
            ASSERTL0(m_xmlDoc[0], "No XML document loaded.");
            TiXmlHandle docHandle(m_xmlDoc[0]);
            TiXmlElement* e;
            e = docHandle.FirstChildElement("NEKTAR").
                FirstChildElement("CONDITIONS").Element();

            ReadParameters(e);

                cout << "InitSession printing parameters" << endl;
                if (m_comm->GetRank() == 0)
                {
                    cout << "Parameters:" << endl;
                    ParameterMap::iterator x;
                    for (x = m_parameters[0].begin();
                            x != m_parameters[0].end(); ++x)
                    {
                        cout << "\t" << x->first << " = " << x->second << endl;
                    }
                    cout << endl;
                }

        	// Start threads
        	StartThreads();

        	// Split up the communicator
        	PartitionComm(e);

            // If running in parallel change the default global sys solution
            // type.
            if (m_comm->GetSize() > 1 )
            {
            	GetSolverInfoDefaults()["GLOBALSYSSOLN"] =
            			"IterativeStaticCond";
            }

            unsigned int vNumWorkers = m_threadManager->GetMaxNumWorkers();
            for (unsigned int i=0; i < vNumWorkers; i++)
            {
            	m_threadManager->QueueJob(new SessionReader::SessionJob(GetSharedThisPtr()));
            }

            m_threadManager->Wait();
			m_comm->Block();

        }

        /**
         * @brief Parses the command-line arguments for known options and
         * filenames.
         */
        std::vector<std::string> SessionReader::ParseCommandLineArguments(
            int argc, char *argv[])
        {
            // List the publically visible options (listed using --help)
            po::options_description desc("Allowed options");
            desc.add_options()
                ("verbose,v",    "be verbose")
                ("version,V",    "print version information")
                ("help,h",       "print this help message")
                ("solverinfo,I", po::value<vector<std::string> >(), 
                                 "override a SOLVERINFO property")
                ("parameter,P",  po::value<vector<std::string> >(),
                                 "override a parameter")
                ("shared-filesystem,s", "Using shared filesystem.")
                ("npx",          po::value<int>(),
                                 "number of procs in X-dir")
                ("npy",          po::value<int>(),
                                 "number of procs in Y-dir")
                ("npz",          po::value<int>(),
                                 "number of procs in Z-dir")
                ("nsz",          po::value<int>(),
                                 "number of slices in Z-dir")
                ("part-only",    po::value<int>(),
                                 "only partition mesh into N partitions.")
                ("part-info",    "Output partition information")
            ;
            
            CmdLineArgMap::const_iterator cmdIt;
            for (cmdIt  = GetCmdLineArgMap().begin();
                 cmdIt != GetCmdLineArgMap().end(); ++cmdIt)
            {
                std::string names = cmdIt->first;
                if (cmdIt->second.shortName != "")
                {
                    names += "," + cmdIt->second.shortName;
                }
                if (cmdIt->second.isFlag)
                {
                    desc.add_options()
                        (names.c_str(), cmdIt->second.description.c_str())
                    ;
                }
                else
                {
                    desc.add_options()
                        (names.c_str(), po::value<std::string>(),
                         cmdIt->second.description.c_str())
                    ;
                }
            }

            // List hidden options (e.g. session file arguments are not actually
            // specified using the input-file option by the user).
            po::options_description hidden("Hidden options");
            hidden.add_options()
                    ("input-file", po::value< vector<string> >(), 
                                   "input filename")
            ;

            // Combine all options for the parser
            po::options_description all("All options");
            all.add(desc).add(hidden);

            // Session file is a positional option
            po::positional_options_description p;
            p.add("input-file", -1);

            // Parse the command-line options
            po::parsed_options parsed = po::command_line_parser(argc, argv).
                                                options(all).
                                                positional(p).
                                                allow_unregistered().
                                                run();

            // Extract known options to map and update
            po::store(parsed, m_cmdLineOptions);
            po::notify(m_cmdLineOptions);

            // Help message
            if (m_cmdLineOptions.count("help"))
            {
                cout << desc;
                exit(0);
            }

            // Version information
            if (m_cmdLineOptions.count("version"))
            {
                cout << "Nektar++ version " << NEKTAR_VERSION;

                if (NekConstants::kGitSha1 != "GITDIR-NOTFOUND")
                {
                    string sha1(NekConstants::kGitSha1);
                    string branch(NekConstants::kGitBranch);
                    boost::replace_all(branch, "refs/heads/", "");

                    cout << " (git changeset " << sha1.substr(0, 8) << ", ";

                    if (branch == "")
                    {
                        cout << "detached head";
                    }
                    else
                    {
                        cout << "head " << branch;
                    }

                    cout << ")";
                }

                cout << endl;
                exit(0);
            }

            // Enable verbose mode
            if (m_cmdLineOptions.count("verbose"))
            {
                m_verbose = true;
            }
            else
            {
                m_verbose = false;
            }
            
            // Print a warning for unknown options
            std::vector< po::basic_option<char> >::iterator x;
            for (x = parsed.options.begin(); x != parsed.options.end(); ++x)
            {
                if (x->unregistered)
                {
                    cout << "Warning: Unknown option: " << x->string_key 
                         << endl;
                }
            }

            // Return the vector of filename(s) given as positional options
            if (m_cmdLineOptions.count("input-file"))
            {
                return m_cmdLineOptions["input-file"].as<
                    std::vector<std::string> >();
            }
            else
            {
                return std::vector<std::string>();
            }
        }


        /**
         *
         */
        std::string SessionReader::ParseSessionName(
                std::vector<std::string> &filenames)
        {
            ASSERTL0(!filenames.empty(),
                     "At least one filename expected.");

            std::string retval = "";

            // First input file defines the session name
            std::string fname = filenames[0];

            // If loading a pre-partitioned mesh, remove _xml extension
            if (fname.size() > 4 &&
                fname.substr(fname.size() - 4, 4) == "_xml")
            {
                retval = fname.substr(0, fname.find_last_of("_"));
            }
            // otherwise remove the .xml extension
            else if (fname.size() > 4 &&
                fname.substr(fname.size() - 4, 4) == ".xml")
            {
                retval = fname.substr(0, fname.find_last_of("."));
            }
            // If compressed .xml.gz, remove both extensions
            else if (fname.size() > 7 &&
                fname.substr(fname.size() - 7, 7) == ".xml.gz")
            {
                retval = fname.substr(0, fname.find_last_of("."));
                retval = retval.substr(0, retval.find_last_of("."));
            }

            return retval;
        }


        /**
         *
         */
        TiXmlDocument& SessionReader::GetDocument()
        {
            ASSERTL1(m_xmlDoc[m_threadManager->GetWorkerNum()], "XML Document not defined.");
            return *m_xmlDoc[m_threadManager->GetWorkerNum()];
        }


        /**
         * The single parameter specifies a path to the requested element in a
         * similar format to the filesystem path. Given the following XML:
         * @code
         * <NEKTAR>
         *   <CONDITIONS>
         *     <PARAMETERS>
         *     ...
         *     </PARAMETERS>
         *   </CONDITIONS>
         * </NEKTAR>
         * @endcode
         * the PARAMETERS element would be retrieved by requesting the path:
         * @code
         * Nektar/Conditions/Parameters
         * @endcode
         * @note Paths are case-insensitive.
         *
         * @param   pPath       Path to requested element.
         *
         * @return Direct pointer to requested XML Element.
         */
        TiXmlElement* SessionReader::GetElement(const string& pPath)
        {
            std::string vPath = boost::to_upper_copy(pPath);
            std::vector<std::string> st;
            boost::split(st, vPath, boost::is_any_of("\\/ "));
            ASSERTL0(st.size() > 0, "No path given in XML element request.");

            TiXmlElement* vReturn = m_xmlDoc[m_threadManager->GetWorkerNum()]->FirstChildElement(st[0].c_str());
            ASSERTL0(vReturn, std::string("Cannot find element '")
                              + st[0] + std::string("'."));
            for (int i = 1; i < st.size(); ++i)
            {
                vReturn = vReturn->FirstChildElement(st[i].c_str());
                ASSERTL0(vReturn, std::string("Cannot find element '")
                                  + st[i] + std::string("'."));
            }
            return vReturn;
        }


        /**
         *
         */
        bool SessionReader::DefinesElement(const std::string &pPath) const
        {
            std::string vPath = boost::to_upper_copy(pPath);
            std::vector<std::string> st;
            boost::split(st, vPath, boost::is_any_of("\\/ "));
            ASSERTL0(st.size() > 0, "No path given in XML element request.");

            TiXmlElement* vReturn = m_xmlDoc[m_threadManager->GetWorkerNum()]->FirstChildElement(st[0].c_str());
            ASSERTL0(vReturn, std::string("Cannot find element '")
                              + st[0] + std::string("'."));
            for (int i = 1; i < st.size(); ++i)
            {
                vReturn = vReturn->FirstChildElement(st[i].c_str());
                if (!vReturn) return false;
            }
            return true;
        }


        /**
         *
         */
        const std::vector<std::string>& SessionReader::GetFilenames() const
        {
            return m_filenames;
        }


        /**
         *
         */
        const std::string& SessionReader::GetSessionName() const
        {
            return m_sessionName;
        }


        /**
         * Output is of the form [sessionName]_P[idx] where idx is the rank
         * of the process.
         */
        const std::string SessionReader::GetSessionNameRank() const
        {
            std::string  dirname = m_sessionName + "_xml"; 
            fs::path     pdirname(dirname);
            
            std::string vFilename = "P" + boost::lexical_cast<std::string>(m_comm->GetRowComm()->GetRank());
            fs::path    pFilename(vFilename);            

            fs::path fullpath = pdirname / pFilename;

            return PortablePath(fullpath);
        }

        /**
         *
         */
        CommSharedPtr& SessionReader::GetComm()
        {
            return m_comm;
        }


        /**
         * This routine finalises any parallel communication.
         *
         * @note This routine should be called at the very end of a users
         * application.
         */
        void SessionReader::Finalise()
        {
            m_comm->Finalise();
        }


        /**
         *
         */
        bool SessionReader::DefinesParameter(const std::string& pName) const
        {
            std::string vName = boost::to_upper_copy(pName);
            unsigned int vThr = Nektar::Thread::GetThreadMaster().GetInstance(Nektar::Thread::ThreadMaster::SessionJob)->GetWorkerNum();
            ParameterMap::const_iterator paramIter = m_parameters[vThr].find(vName);
            return (paramIter != m_parameters[vThr].end());
        }


        /**
         * If the parameter is not defined, termination occurs. Therefore, the
         * parameters existence should be tested for using #DefinesParameter
         * before calling this function.
         *
         * @param   pName       The name of a floating-point parameter.
         * @returns The value of the floating-point parameter.
         */
        const NekDouble& SessionReader::GetParameter(
            const std::string& pName) const
        {
            std::string vName = boost::to_upper_copy(pName);
            unsigned int vThr = Nektar::Thread::GetThreadMaster().GetInstance(Nektar::Thread::ThreadMaster::SessionJob)->GetWorkerNum();
        	ParameterMap::const_iterator paramIter = m_parameters[vThr].find(vName);

            ASSERTL0(paramIter != m_parameters[vThr].end(),
                     "Unable to find requested parameter: " + pName);

            return paramIter->second;
        }


        /**
         *
         */
        void SessionReader::LoadParameter(
            const std::string &pName, int &pVar) const
        {
            std::string vName = boost::to_upper_copy(pName);
            unsigned int vThr = Nektar::Thread::GetThreadMaster().GetInstance(Nektar::Thread::ThreadMaster::SessionJob)->GetWorkerNum();
            ParameterMap::const_iterator paramIter = m_parameters[vThr].find(vName);
            ASSERTL0(paramIter != m_parameters[vThr].end(), "Required parameter '" +
                     pName + "' not specified in session.");
            pVar = (int)floor(paramIter->second);
        }


        /**
         *
         */
        void SessionReader::LoadParameter(
            const std::string &pName, int &pVar, const int &pDefault) const
        {
            std::string vName = boost::to_upper_copy(pName);
            unsigned int vThr = Nektar::Thread::GetThreadMaster().GetInstance(Nektar::Thread::ThreadMaster::SessionJob)->GetWorkerNum();
            ParameterMap::const_iterator paramIter = m_parameters[vThr].find(vName);
            if(paramIter != m_parameters[vThr].end())
            {
                pVar = (int)floor(paramIter->second);
            }
            else
            {
                pVar  = pDefault;
            }
        }


        /**
         *
         */
        void SessionReader::LoadParameter(
            const std::string &pName, NekDouble& pVar) const
        {
            std::string vName = boost::to_upper_copy(pName);
            unsigned int vThr = Nektar::Thread::GetThreadMaster().GetInstance(Nektar::Thread::ThreadMaster::SessionJob)->GetWorkerNum();
            ParameterMap::const_iterator paramIter = m_parameters[vThr].find(vName);
            ASSERTL0(paramIter != m_parameters[vThr].end(), "Required parameter '" +
                     pName + "' not specified in session.");
            pVar = paramIter->second;
        }


        /**
         *
         */
        void SessionReader::LoadParameter(
            const std::string &pName, 
                  NekDouble   &pVar, 
            const NekDouble   &pDefault) const
        {
            std::string vName = boost::to_upper_copy(pName);
            unsigned int vThr = Nektar::Thread::GetThreadMaster().GetInstance(Nektar::Thread::ThreadMaster::SessionJob)->GetWorkerNum();
            ParameterMap::const_iterator paramIter = m_parameters[vThr].find(vName);
            if(paramIter != m_parameters[vThr].end())
            {
                pVar = paramIter->second;
            }
            else
            {
                pVar = pDefault;
            }
        }



        /**
         *
         */
        void SessionReader::SetParameter(const std::string &pName, int &pVar) 
        {
            std::string vName = boost::to_upper_copy(pName);
            unsigned int vThr = Nektar::Thread::GetThreadMaster().GetInstance(Nektar::Thread::ThreadMaster::SessionJob)->GetWorkerNum();
            m_parameters[vThr][vName] = pVar;
        }


        /**
         *
         */
        void SessionReader::SetParameter(
            const std::string &pName, NekDouble& pVar) 
        {
            std::string vName = boost::to_upper_copy(pName);
            unsigned int vThr = Nektar::Thread::GetThreadMaster().GetInstance(Nektar::Thread::ThreadMaster::SessionJob)->GetWorkerNum();
            m_parameters[vThr][vName] = pVar;
        }



        /**
         *
         */
        bool SessionReader::DefinesSolverInfo(const std::string &pName) const
        {
            std::string vName = boost::to_upper_copy(pName);
            unsigned int vThr = Nektar::Thread::GetThreadMaster().GetInstance(Nektar::Thread::ThreadMaster::SessionJob)->GetWorkerNum();
            SolverInfoMap::const_iterator infoIter = m_solverInfo[vThr].find(vName);
            return (infoIter != m_solverInfo[vThr].end());
        }


        /**
         *
         */
        const std::string& SessionReader::GetSolverInfo(
            const std::string &pProperty) const
        {
            std::string vProperty = boost::to_upper_copy(pProperty);
            unsigned int vThr = Nektar::Thread::GetThreadMaster().GetInstance(Nektar::Thread::ThreadMaster::SessionJob)->GetWorkerNum();
            SolverInfoMap::const_iterator iter = m_solverInfo[vThr].find(vProperty);

            ASSERTL1(iter != m_solverInfo[vThr].end(),
                     "Unable to find requested property: " + pProperty);

            return iter->second;
        }

        /**
         *
         */
        void SessionReader::SetSolverInfo(
            const std::string &pProperty, const std::string &pValue) 
        {
            std::string vProperty = boost::to_upper_copy(pProperty);
            unsigned int vThr = Nektar::Thread::GetThreadMaster().GetInstance(Nektar::Thread::ThreadMaster::SessionJob)->GetWorkerNum();
            SolverInfoMap::iterator iter = m_solverInfo[vThr].find(vProperty);

            ASSERTL1(iter != m_solverInfo[vThr].end(),
                     "Unable to find requested property: " + pProperty);

            iter->second = pValue;
        }

        /**
         *
         */
        void SessionReader::LoadSolverInfo(
            const std::string &pName, 
                  std::string &pVar, 
            const std::string &pDefault) const
        {
            std::string vName = boost::to_upper_copy(pName);
            unsigned int vThr = Nektar::Thread::GetThreadMaster().GetInstance(Nektar::Thread::ThreadMaster::SessionJob)->GetWorkerNum();
            SolverInfoMap::const_iterator infoIter = m_solverInfo[vThr].find(vName);
            if(infoIter != m_solverInfo[vThr].end())
            {
                pVar = infoIter->second;
            }
            else
            {
                pVar = pDefault;
            }
        }


        /**
         *
         */
        void SessionReader::MatchSolverInfo(
            const std::string &pName,
            const std::string &pTrueVal,
                  bool        &pVar,
            const bool        &pDefault) const
        {
            std::string vName = boost::to_upper_copy(pName);
            unsigned int vThr = Nektar::Thread::GetThreadMaster().GetInstance(Nektar::Thread::ThreadMaster::SessionJob)->GetWorkerNum();
            SolverInfoMap::const_iterator infoIter = m_solverInfo[vThr].find(vName);
            if(infoIter != m_solverInfo[vThr].end())
            {
                pVar = boost::iequals(infoIter->second, pTrueVal);
            }
            else
            {
                pVar = pDefault;
            }
        }


        /**
         *
         */
        bool SessionReader::MatchSolverInfo(
            const std::string &pName,
            const std::string &pTrueVal) const
        {
            if (DefinesSolverInfo(pName))
            {
                std::string vName = boost::to_upper_copy(pName);
                unsigned int vThr = Nektar::Thread::GetThreadMaster().GetInstance(Nektar::Thread::ThreadMaster::SessionJob)->GetWorkerNum();
                SolverInfoMap::const_iterator iter = m_solverInfo[vThr].find(vName);
                if(iter != m_solverInfo[vThr].end())
                {
                    return true;
                }
            }
            return false;
        }


        /**
         *
         */
        bool SessionReader::DefinesGlobalSysSolnInfo(const std::string &pVariable, 
                                                     const std::string &pProperty) const
        {

            GloSysSolnInfoList::const_iterator iter =
                    GetGloSysSolnList().find(pVariable);
            if(iter == GetGloSysSolnList().end())
            {
                return false;
            }

            std::string vProperty = boost::to_upper_copy(pProperty);
            
            GloSysInfoMap::const_iterator iter1 = iter->second.find(vProperty);
            if(iter1 == iter->second.end())
            {
                return false;
            }
            
            return true;
        }

        
        /**
         *
         */
        const std::string &SessionReader::GetGlobalSysSolnInfo(const std::string &pVariable, const std::string &pProperty) const
        {
            GloSysSolnInfoList::const_iterator iter; 

            ASSERTL0( (iter = GetGloSysSolnList().find(pVariable)) !=
                              GetGloSysSolnList().end(),
                      "Failed to find variable in GlobalSysSolnInfoList");

            std::string vProperty = boost::to_upper_copy(pProperty);
            GloSysInfoMap::const_iterator iter1; 

            ASSERTL0( (iter1 = iter->second.find(vProperty)) != iter->second.end(),
                      "Failed to find property: " + vProperty + " in GlobalSysSolnInfoList");
            
            return iter1->second;
        }
        
        /**
         *
         */
        bool SessionReader::DefinesGeometricInfo(const std::string &pName) const
        {
            std::string vName = boost::to_upper_copy(pName);
            unsigned int vThr = m_threadManager->GetWorkerNum();
            GeometricInfoMap::const_iterator iter = m_geometricInfo[vThr].find(vName);
            return (iter != m_geometricInfo[vThr].end());
        }


        /**
         *
         */
        void SessionReader::LoadGeometricInfo(
            const std::string &pName,
                  std::string &pVar,
            const std::string &pDefault) const
        {
            std::string vName = boost::to_upper_copy(pName);
            unsigned int vThr = m_threadManager->GetWorkerNum();
            GeometricInfoMap::const_iterator iter = m_geometricInfo[vThr].find(vName);
            if(iter != m_geometricInfo[vThr].end())
            {
                pVar = iter->second;
            }
            else
            {
                pVar = pDefault;
            }
        }


        /**
         *
         */
        void SessionReader::LoadGeometricInfo(
            const std::string &pName,
                  bool        &pVar,
            const bool        &pDefault) const
        {
            std::string vName = boost::to_upper_copy(pName);
            unsigned int vThr = m_threadManager->GetWorkerNum();
            GeometricInfoMap::const_iterator iter = m_geometricInfo[vThr].find(vName);
            if(iter != m_geometricInfo[vThr].end())
            {
                if (iter->second == "TRUE")
                {
                    pVar = true;
                }
                else
                {
                    pVar = false;
                }
            }
            else
            {
                pVar = pDefault;
            }
        }


        /**
         *
         */
        void SessionReader::LoadGeometricInfo(
            const std::string &pName,
                  NekDouble   &pVar,
            const NekDouble   &pDefault) const
        {
            std::string vName = boost::to_upper_copy(pName);
            unsigned int vThr = m_threadManager->GetWorkerNum();
            GeometricInfoMap::const_iterator iter = m_geometricInfo[vThr].find(vName);
            if(iter != m_geometricInfo[vThr].end())
            {
                pVar = std::atoi(iter->second.c_str());
            }
            else
            {
                pVar = pDefault;
            }
        }


        /**
         *
         */
        void SessionReader::MatchGeometricInfo(
            const std::string &pName,
            const std::string &pTrueVal,
                  bool        &pVar,
            const bool        &pDefault) const
        {
            std::string vName = boost::to_upper_copy(pName);
            unsigned int vThr = m_threadManager->GetWorkerNum();
            GeometricInfoMap::const_iterator iter = m_geometricInfo[vThr].find(vName);
            if(iter != m_geometricInfo[vThr].end())
            {
                pVar = boost::iequals(iter->second, pTrueVal);
            }
            else
            {
                pVar  = pDefault;
            }
        }


        /**
         *
         */
        const std::string& SessionReader::GetVariable(
            const unsigned int &idx) const
        {
            unsigned int vThr = m_threadManager->GetWorkerNum();
            ASSERTL0(idx < m_variables[vThr].size(), "Variable index out of range.");
            return m_variables[vThr][idx];
        }



        /**
         *
         */
        void SessionReader::SetVariable(const unsigned int &idx, 
                                        std::string newname) 
        {
            unsigned int vThr = m_threadManager->GetWorkerNum();
            ASSERTL0(idx < m_variables[vThr].size(), "Variable index out of range.");
            m_variables[vThr][idx] = newname;
        }


        /**
         *
         */
        std::vector<std::string> SessionReader::GetVariables() const
        {
            unsigned int vThr = m_threadManager->GetWorkerNum();
            return m_variables[vThr];
        }


        /**
         *
         */
        bool SessionReader::DefinesFunction(const std::string &pName) const
        {
            FunctionMap::const_iterator it1;
            std::string vName = boost::to_upper_copy(pName);
            unsigned int vThr = m_threadManager->GetWorkerNum();

            if ((it1 = m_functions[vThr].find(vName)) != m_functions[vThr].end())
            {
                return true;
            }
            return false;
        }


        /**
         *
         */
        bool SessionReader::DefinesFunction(
            const std::string &pName,
            const std::string &pVariable,
            const int pDomain) const
        {
            FunctionMap::const_iterator it1;
            FunctionVariableMap::const_iterator it2;
            std::string vName = boost::to_upper_copy(pName);
            unsigned int vThr = m_threadManager->GetWorkerNum();

            // Check function exists
            if ((it1 = m_functions[vThr].find(vName))     != m_functions[vThr].end())
            {
                pair<std::string, int> key(pVariable,pDomain);
                pair<std::string, int> defkey("*",pDomain);
                bool varExists =
                    (it2 = it1->second.find(key)) != it1->second.end() ||
                    (it2 = it1->second.find(defkey)) != it1->second.end();
                return varExists;
            }
            return false;
        }


        /**
         *
         */
        EquationSharedPtr SessionReader::GetFunction(
            const std::string &pName,
            const std::string &pVariable,
            const int pDomain) const
        {
            FunctionMap::const_iterator it1;
            FunctionVariableMap::const_iterator it2, it3;
            std::string vName = boost::to_upper_copy(pName);
            unsigned int vThr = m_threadManager->GetWorkerNum();

            ASSERTL0((it1 = m_functions[vThr].find(vName)) != m_functions[vThr].end(),
                     std::string("No such function '") + pName
                     + std::string("' has been defined in the session file."));

            // Check for specific and wildcard definitions
            pair<std::string,int> key(pVariable,pDomain);
            pair<std::string,int> defkey("*",pDomain);
            bool specific = (it2 = it1->second.find(key)) !=
                it1->second.end();
            bool wildcard = (it3 = it1->second.find(defkey)) !=
                it1->second.end();

            // Check function is defined somewhere
            ASSERTL0(specific || wildcard,
                     "No such variable " + pVariable
                     + " in domain " + boost::lexical_cast<string>(pDomain) 
                     + " defined for function " + pName
                     + " in session file.");

            // If not specific, must be wildcard
            if (!specific)
            {
                it2 = it3;
            }

            ASSERTL0((it2->second.m_type == eFunctionTypeExpression),
                    std::string("Function is defined by a file."));
            return it2->second.m_expression;
        }


        /**
         *
         */
        EquationSharedPtr SessionReader::GetFunction(
            const std::string  &pName,
            const unsigned int &pVar,
            const int pDomain) const
        {
            unsigned int vThr = m_threadManager->GetWorkerNum();
            ASSERTL0(pVar < m_variables[vThr].size(), "Variable index out of range.");
            return GetFunction(pName, m_variables[vThr][pVar],pDomain);
        }


        /**
         *
         */
        enum FunctionType SessionReader::GetFunctionType(
            const std::string &pName,
            const std::string &pVariable,
            const int pDomain) const
        {
            FunctionMap::const_iterator it1;
            FunctionVariableMap::const_iterator it2, it3;
            std::string vName = boost::to_upper_copy(pName);
            unsigned int vThr = m_threadManager->GetWorkerNum();

            it1 = m_functions[vThr].find(vName);
            ASSERTL0 (it1 != m_functions[vThr].end(),
                      std::string("Function '") + pName
                      + std::string("' not found."));

            // Check for specific and wildcard definitions
            pair<std::string,int> key(pVariable,pDomain);
            pair<std::string,int> defkey("*",pDomain);
            bool specific = (it2 = it1->second.find(key)) !=
                            it1->second.end();
            bool wildcard = (it3 = it1->second.find(defkey)) !=
                            it1->second.end();

            // Check function is defined somewhere
            ASSERTL0(specific || wildcard,
                     "No such variable " + pVariable
                     + " in domain " + boost::lexical_cast<string>(pDomain) 
                     + " defined for function " + pName
                     + " in session file.");

            // If not specific, must be wildcard
            if (!specific)
            {
                it2 = it3;
            }

            return it2->second.m_type;
        }


        /**
         *
         */
        enum FunctionType SessionReader::GetFunctionType(
            const std::string  &pName,
            const unsigned int &pVar,
            const int pDomain) const
        {
            unsigned int vThr = m_threadManager->GetWorkerNum();
            ASSERTL0(pVar < m_variables[vThr].size(), "Variable index out of range.");
            return GetFunctionType(pName, m_variables[vThr][pVar],pDomain);
        }


        /**
         *
         */
        std::string SessionReader::GetFunctionFilename(
            const std::string &pName, 
            const std::string &pVariable,
            const int pDomain) const
        {
            FunctionMap::const_iterator it1;
            FunctionVariableMap::const_iterator it2, it3;
            std::string vName = boost::to_upper_copy(pName);
            unsigned int vThr = m_threadManager->GetWorkerNum();

            it1 = m_functions[vThr].find(vName);
            ASSERTL0 (it1 != m_functions[vThr].end(),
                      std::string("Function '") + pName
                      + std::string("' not found."));

            // Check for specific and wildcard definitions
            pair<std::string,int> key(pVariable,pDomain);
            pair<std::string,int> defkey("*",pDomain);
            bool specific = (it2 = it1->second.find(key)) !=
                            it1->second.end();
            bool wildcard = (it3 = it1->second.find(defkey)) !=
                            it1->second.end();

            // Check function is defined somewhere
            ASSERTL0(specific || wildcard,
                     "No such variable " + pVariable
                     + " in domain " + boost::lexical_cast<string>(pDomain) 
                     + " defined for function " + pName
                     + " in session file.");
            
            // If not specific, must be wildcard
            if (!specific)
            {
                it2 = it3;
            }

            return it2->second.m_filename;
        }


        /**
         *
         */
        std::string SessionReader::GetFunctionFilename(
            const std::string  &pName, 
            const unsigned int &pVar,
            const int pDomain) const
        {
            unsigned int vThr = m_threadManager->GetWorkerNum();
            ASSERTL0(pVar < m_variables[vThr].size(), "Variable index out of range.");
            return GetFunctionFilename(pName, m_variables[vThr][pVar],pDomain);
        }


        /**
         *
         */
        std::string SessionReader::GetFunctionFilenameVariable(
            const std::string &pName,
            const std::string &pVariable,
            const int pDomain) const
        {
            unsigned int vThr = m_threadManager->GetWorkerNum();
            FunctionMap::const_iterator it1;
            FunctionVariableMap::const_iterator it2, it3;
            std::string vName = boost::to_upper_copy(pName);

            it1 = m_functions[vThr].find(vName);
            ASSERTL0 (it1 != m_functions[vThr].end(),
                      std::string("Function '") + pName
                      + std::string("' not found."));

            // Check for specific and wildcard definitions
            pair<std::string,int> key(pVariable,pDomain);
            pair<std::string,int> defkey("*",pDomain);
            bool specific = (it2 = it1->second.find(key)) !=
                            it1->second.end();
            bool wildcard = (it3 = it1->second.find(defkey)) !=
                            it1->second.end();

            // Check function is defined somewhere
            ASSERTL0(specific || wildcard,
                     "No such variable " + pVariable
                     + " in domain " + boost::lexical_cast<string>(pDomain)
                     + " defined for function " + pName
                     + " in session file.");

            // If not specific, must be wildcard
            if (!specific)
            {
                it2 = it3;
            }

            return it2->second.m_fileVariable;
        }


        /**
         *
         */
        AnalyticExpressionEvaluator& SessionReader::GetExpressionEvaluator()
        {
            unsigned int vThr = Nektar::Thread::GetThreadMaster().GetInstance(Nektar::Thread::ThreadMaster::SessionJob)->GetWorkerNum();
            return *m_exprEvaluator[vThr];
        }


        /**
         *
         */
        bool SessionReader::DefinesTag(const std::string &pName) const
        {
            unsigned int vThr = m_threadManager->GetWorkerNum();
            std::string vName = boost::to_upper_copy(pName);
            TagMap::const_iterator vTagIterator = m_tags[vThr].find(vName);
            return (vTagIterator != m_tags[vThr].end());
        }


        /**
         *
         */
        void SessionReader::SetTag(
            const std::string &pName, 
            const std::string &pValue)
        {
        	unsigned int vThr = m_threadManager->GetWorkerNum();
        	std::string vName = boost::to_upper_copy(pName);
            m_tags[vThr][vName] = pValue;
        }


        /**
         *
         */
        const std::string &SessionReader::GetTag(const std::string& pName) const
        {
        	unsigned int vThr = m_threadManager->GetWorkerNum();
        	std::string vName = boost::to_upper_copy(pName);
            TagMap::const_iterator vTagIterator = m_tags[vThr].find(vName);
            ASSERTL0(vTagIterator != m_tags[vThr].end(),
                     "Requested tag does not exist.");
            return vTagIterator->second;
        }


        /**
         *
         */
        const FilterMap &SessionReader::GetFilters() const
        {
            unsigned int vThr = m_threadManager->GetWorkerNum();
            return m_filters[vThr];
        }


        /**
         *
         */
        bool SessionReader::DefinesCmdLineArgument(
            const std::string& pName) const
        {
            return (m_cmdLineOptions.find(pName) != m_cmdLineOptions.end());
        }


        /**
         *
         */
        void SessionReader::SubstituteExpressions(std::string& pExpr)
        {
            ExpressionMap::iterator exprIter;
            unsigned int vThr = m_threadManager->GetWorkerNum();
            for (exprIter  = m_expressions[vThr].begin();
                 exprIter != m_expressions[vThr].end(); ++exprIter)
            {
                boost::replace_all(pExpr, exprIter->first, exprIter->second);
            }
        }

        CompositeOrdering SessionReader::GetCompositeOrdering() const
        {
            return m_compOrder;
        }

        BndRegionOrdering SessionReader::GetBndRegionOrdering() const
        {
            unsigned int vThr = m_threadManager->GetWorkerNum();
            return m_bndRegOrder[vThr];
        }

        /**
         *
         */
        void SessionReader::LoadDoc(
            const std::string &pFilename,
            TiXmlDocument* pDoc) const
        {
            if (pFilename.size() > 3 &&
                pFilename.substr(pFilename.size() - 3, 3) == ".gz")
            {
                ifstream file(pFilename.c_str(),
                              ios_base::in | ios_base::binary);
                ASSERTL0(file.good(), "Unable to open file: " + pFilename);
                stringstream ss;
                io::filtering_streambuf<io::input> in;
                in.push(io::gzip_decompressor());
                in.push(file);
                try
                {
                    io::copy(in, ss);
                    ss >> (*pDoc);
                }
                catch (io::gzip_error& e)
                {
                    ASSERTL0(false,
                             "Error: File '" + pFilename + "' is corrupt.");
                }
            }
            else if (pFilename.size() > 4 &&
                    pFilename.substr(pFilename.size() - 4, 4) == "_xml")
            {
                fs::path    pdirname(pFilename);
                boost::format pad("P%1$07d.xml");
                pad % m_comm->GetRank();
                fs::path    pRankFilename(pad.str());
                fs::path fullpath = pdirname / pRankFilename;

                ifstream file(PortablePath(fullpath).c_str());
                ASSERTL0(file.good(), "Unable to open file: " + fullpath.string());
                file >> (*pDoc);
            }
            else
            {
                ifstream file(pFilename.c_str());
                ASSERTL0(file.good(), "Unable to open file: " + pFilename);
                file >> (*pDoc);
            }
        }

        /**
         *
         */
        TiXmlDocument *SessionReader::MergeDoc(
            const std::vector<std::string> &pFilenames) const
        {
            ASSERTL0(pFilenames.size() > 0, "No filenames for merging.");

            // Read the first document
            TiXmlDocument *vMainDoc = new TiXmlDocument;
            LoadDoc(pFilenames[0], vMainDoc);

            TiXmlHandle vMainHandle(vMainDoc);
            TiXmlElement* vMainNektar = 
                vMainHandle.FirstChildElement("NEKTAR").Element();

            // Read all subsequent XML documents.
            // For each element within the NEKTAR tag, use it to replace the
            // version already present in the loaded XML data.
            for (int i = 1; i < pFilenames.size(); ++i)
            {
                if((pFilenames[i].compare(pFilenames[i].size()-3,3,"xml") == 0)
                   ||(pFilenames[i].compare(pFilenames[i].size()-6,6,"xml.gz") == 0))
                {
                    TiXmlDocument* vTempDoc = new TiXmlDocument;
                    LoadDoc(pFilenames[i], vTempDoc);
                    
                    TiXmlHandle docHandle(vTempDoc);
                    TiXmlElement* vTempNektar;
                    vTempNektar = docHandle.FirstChildElement("NEKTAR").Element();
                    ASSERTL0(vTempNektar, "Unable to find NEKTAR tag in file.");
                    TiXmlElement* p = vTempNektar->FirstChildElement();
                    
                    while (p)
                    {
                        TiXmlElement *vMainEntry = 
                            vMainNektar->FirstChildElement(p->Value());
                        TiXmlElement *q = new TiXmlElement(*p);
                        if (vMainEntry)
                        {
                            vMainNektar->RemoveChild(vMainEntry);
                        }
                        vMainNektar->LinkEndChild(q);
                        p = p->NextSiblingElement();
                    }
                    
                    delete vTempDoc;
                }
            }
            return vMainDoc;
        }


        /**
         *
         */
        void SessionReader::ParseDocument()
        {
            ASSERTL0(m_threadManager->IsInitialised(),
                "ThreadManager not initialised in ParseDocument.");
            unsigned int vThr = m_threadManager->GetWorkerNum();

            m_threadManager->Hold();

            // Check we actually have a document loaded.
            ASSERTL0(m_xmlDoc[vThr], "No XML document loaded.");

            // Look for all data in CONDITIONS block.
            TiXmlHandle docHandle(m_xmlDoc[vThr]);
            TiXmlElement* e;
            e = docHandle.FirstChildElement("NEKTAR").
                FirstChildElement("CONDITIONS").Element();

            // Read the various sections of the CONDITIONS block
            ReadParameters        (e);
            ReadSolverInfo        (e);
            if (vThr == 0)
            {
            	ReadGlobalSysSolnInfo (e);
            }
            ReadExpressions       (e);
            ReadVariables         (e);
            ReadFunctions         (e);

            e = docHandle.FirstChildElement("NEKTAR").
                FirstChildElement("FILTERS").Element();

            ReadFilters(e);

        }


        /**
         *
         */
        void SessionReader::CreateComm(
            int               &argc, 
            char*              argv[])
        {
            TiXmlHandle docHandle(m_xmlDoc[0]); // threads not initialised yet
            TiXmlElement* e;
            e = docHandle.FirstChildElement("NEKTAR").
                FirstChildElement("CONDITIONS").Element();

            ReadSolverInfo(e);

            int nthreads;
            LoadParameter("NThreads", nthreads, 1);

            if (argc == 0 && nthreads == 1)
            {
                m_comm = GetCommFactory().CreateInstance("Serial", 0, 0);
            }
            else
            {

                string vCommModule("Serial");
                if (nthreads > 1)
                {
                    vCommModule = "ParallelMPI";
                }
                if (GetCommFactory().ModuleExists("ParallelMPI"))
                {
                    vCommModule = "ParallelMPI";
                }

                m_comm = GetCommFactory().CreateInstance(vCommModule,argc,argv);

                // Code altering default solver if parallel job moved to InitSession
            }
        }


        /**
         *
         */
        void SessionReader::PartitionMesh()
        {
            ASSERTL0(m_comm.get(), "Communication not initialised.");

            // Get row of comm, or the whole comm if not split
            CommSharedPtr vCommMesh = m_comm->GetRowComm();
            const bool isRoot = (m_comm->GetRank() == 0);
            unsigned int vThr = m_threadManager->GetWorkerNum();

            // Delete any existing loaded mesh
            if (m_xmlDoc[vThr])
            {
                delete m_xmlDoc[vThr];
            }

            // Load file for root process only (since this is always needed)
            // and determine if the provided geometry has already been
            // partitioned. This will be the case if the user provides the
            // directory of mesh partitions as an input. Partitioned geometries
            // have the attribute
            //    PARTITION=X
            // where X is the number of the partition (and should match the
            // process rank). The result is shared with all other processes.
            int isPartitioned = 0;
            if (isRoot)
            {
                m_xmlDoc[0] = MergeDoc(m_filenames);
                if (DefinesElement("Nektar/Geometry"))
                {
                    if (GetElement("Nektar/Geometry")->Attribute("PARTITION"))
                    {
                        cout << "Using pre-partitioned mesh." << endl;
                        isPartitioned = 1;
                    }
                }
            }
            GetComm()->AllReduce(isPartitioned, LibUtilities::ReduceMax);

            // If the mesh is already partitioned, we are done. Remaining
            // processes must load their partitions.
            if (isPartitioned) {
                if (!isRoot)
                {
                    m_xmlDoc[vThr] = MergeDoc(m_filenames);
                }
                return;
            }

            // Default partitioner to use is Metis. Use Scotch as default
            // if it is installed. Override default with command-line flags
            // if they are set.
            string vPartitionerName = "Metis";
            if (GetMeshPartitionFactory().ModuleExists("Scotch"))
            {
                vPartitionerName = "Scotch";
            }
            if (DefinesCmdLineArgument("use-metis"))
            {
                vPartitionerName = "Metis";
            }
            if (DefinesCmdLineArgument("use-scotch"))
            {
                vPartitionerName = "Scotch";
            }

            // Mesh has not been partitioned so do partitioning if required.
            // Note in the serial case nothing is done as we have already loaded
            // the mesh.
            if (DefinesCmdLineArgument("part-only"))
            {
                // Perform partitioning of the mesh only. For this we insist
                // the code is run in serial (parallel execution is pointless).
                ASSERTL0(GetComm()->GetSize() == 1,
                        "The 'part-only' option should be used in serial.");

                // Number of partitions is specified by the parameter.
                int nParts = GetCmdLineArgument<int>("part-only");
                SessionReaderSharedPtr vSession     = GetSharedThisPtr();
                MeshPartitionSharedPtr vPartitioner = 
                        GetMeshPartitionFactory().CreateInstance(
                                            vPartitionerName, vSession);
                vPartitioner->PartitionMesh(nParts, true);
                vPartitioner->WriteAllPartitions(vSession);
                vPartitioner->GetCompositeOrdering(m_compOrder);
                vPartitioner->GetBndRegionOrdering(m_bndRegOrder[0]);

                if (isRoot && DefinesCmdLineArgument("part-info"))
                {
                    vPartitioner->PrintPartInfo(std::cout);
                }

                Finalise();
                exit(0);
            }
            else if (vCommMesh->GetSize() > 1)
            {
                SessionReaderSharedPtr vSession     = GetSharedThisPtr();
                int nParts = vCommMesh->GetSize();
                if (DefinesCmdLineArgument("shared-filesystem"))
                {
                    ASSERTL0(1,
                        "I ain't done this bit cos it looks weird.");
                    /* SJCFIXME
                    CommSharedPtr vComm = GetComm();
                    vector<unsigned int> keys, vals;
                    int i;

                    if (vComm->GetRank() == 0)
                    {
                        m_xmlDoc[0] = MergeDoc(m_filenames);

                        MeshPartitionSharedPtr vPartitioner =
                                GetMeshPartitionFactory().CreateInstance(
                                                    vPartitionerName, vSession);
                        vPartitioner->PartitionMesh(nParts, true);
                        vPartitioner->WriteAllPartitions(vSession);
                        vPartitioner->GetCompositeOrdering(m_compOrder);
                        for (unsigned int i = 0; i < vNumThr; ++i) {
                            vPartitioner->GetBndRegionOrdering(m_bndRegOrder[i], i);
                        }

                        // Communicate orderings to the other processors/threads.

                        // First send sizes of the orderings and boundary
                        // regions to allocate storage on the remote end.
                        keys.resize(2);
                        keys[0] = m_compOrder.size();
                        keys[1] = m_bndRegOrder.size();

                        for (i = 1; i < vComm->GetSize(); ++i)
                        {
                            vComm->Send(i, keys);
                        }

                        // Construct the keys and sizes of values for composite
                        // ordering
                        CompositeOrdering::iterator cIt;
                        keys.resize(m_compOrder.size());
                        vals.resize(m_compOrder.size());

                        for (cIt  = m_compOrder.begin(), i = 0;
                             cIt != m_compOrder.end(); ++cIt, ++i)
                        {
                            keys[i] = cIt->first;
                            vals[i] = cIt->second.size();
                        }

                        // Send across data.
                        for (i = 1; i < vComm->GetSize(); ++i)
                        {
                            vComm->Send(i, keys);
                            vComm->Send(i, vals);

                            for (cIt  = m_compOrder.begin();
                                 cIt != m_compOrder.end(); ++cIt)
                            {
                                vComm->Send(i, cIt->second);
                            }
                        }

                        // Construct the keys and sizes of values for composite
                        // ordering
                        BndRegionOrdering::iterator bIt;
                        keys.resize(m_bndRegOrder.size());
                        vals.resize(m_bndRegOrder.size());

                        for (bIt  = m_bndRegOrder.begin(), i = 0;
                             bIt != m_bndRegOrder.end(); ++bIt, ++i)
                        {
                            keys[i] = bIt->first;
                            vals[i] = bIt->second.size();
                        }

                        // Send across data.
                        for (i = 1; i < vComm->GetSize(); ++i)
                        {
                            vComm->Send(i, keys);
                            vComm->Send(i, vals);

                            for (bIt  = m_bndRegOrder.begin();
                                 bIt != m_bndRegOrder.end(); ++bIt)
                            {
                                vComm->Send(i, bIt->second);
                            }
                        }

                        if (DefinesCmdLineArgument("part-info"))
                        {
                            vPartitioner->PrintPartInfo(std::cout);
                        }
                    }
                    else
                    {
                        keys.resize(2);
                        vComm->Recv(0, keys);

                        int cmpSize = keys[0];
                        int bndSize = keys[1];

                        keys.resize(cmpSize);
                        vals.resize(cmpSize);
                        vComm->Recv(0, keys);
                        vComm->Recv(0, vals);

                        for (int i = 0; i < keys.size(); ++i)
                        {
                            vector<unsigned int> tmp(vals[i]);
                            vComm->Recv(0, tmp);
                            m_compOrder[keys[i]] = tmp;
                        }

                        keys.resize(bndSize);
                        vals.resize(bndSize);
                        vComm->Recv(0, keys);
                        vComm->Recv(0, vals);

                        for (int i = 0; i < keys.size(); ++i)
                        {
                            vector<unsigned int> tmp(vals[i]);
                            vComm->Recv(0, tmp);
                            m_bndRegOrder[keys[i]] = tmp;
                        }
                    }
                    */
                }
                else
                {
                    // Need to load mesh on non-root processes.
                    if (!isRoot)
                    {
                        m_xmlDoc[vThr] = MergeDoc(m_filenames);
                    }

                    // Partitioner now operates in parallel
                    // Each process receives partitioning over interconnect
                    // and writes its own session file to the working directory.
                    MeshPartitionSharedPtr vPartitioner =
                                GetMeshPartitionFactory().CreateInstance(
                                                    vPartitionerName, vSession);
                    vPartitioner->PartitionMesh(nParts, false);
                    vPartitioner->WriteLocalPartition(vSession);
                    vPartitioner->GetCompositeOrdering(m_compOrder);
                    vPartitioner->GetBndRegionOrdering(m_bndRegOrder[vThr]);

                    if (DefinesCmdLineArgument("part-info") && isRoot)
                    {
                        vPartitioner->PrintPartInfo(std::cout);
                    }
                }
                m_comm->Block();

                std::string  dirname = GetSessionName() + "_xml";
                fs::path    pdirname(dirname);
                boost::format pad("P%1$07d.xml");
                pad % m_comm->GetRowComm()->GetRank();
                fs::path    pFilename(pad.str());
                fs::path fullpath = pdirname / pFilename;

                std::string vFilename = PortablePath(fullpath);

                if (m_xmlDoc[vThr])
                {
                    delete m_xmlDoc[vThr];
                }
                m_xmlDoc[vThr] = new TiXmlDocument(vFilename);

                ASSERTL0(m_xmlDoc[vThr], "Failed to create XML document object.");

            	bool loadOkay = m_xmlDoc[vThr]->LoadFile();
            	ASSERTL0(loadOkay, "Unable to load file: " + vFilename      +
            			". Check XML standards compliance. Error on line: " +
            			boost::lexical_cast<std::string>(m_xmlDoc[vThr]->Row()));

            }
            else
            {
                m_xmlDoc[vThr] = MergeDoc(m_filenames);
            }
        }


        /**
         * Splits the processes into a cartesian grid and creates communicators
         * for each row and column of the grid. The grid is defined by the
         * PROC_X parameter which, if specified, gives the number of processes
         * spanned by the Fourier direction. PROC_X must exactly divide the
         * total number of processes or an error is thrown.
         *
         * Also attaches the threading decorator if threading is being used.
         */
        void SessionReader::PartitionComm(TiXmlElement* e)
        {
        	bool threadedCommDone = false;
            unsigned int vNumThr = m_threadManager->GetMaxNumWorkers();
            if (m_comm->GetSize() > 1 || vNumThr > 1)
            {
                int nProcZ  = 1;
                int nProcY  = 1;
                int nProcX  = 1;
                int nStripZ = 1;
                if (DefinesCmdLineArgument("npx")) {
                    nProcX = GetCmdLineArgument<int>("npx");
                }
                if (DefinesCmdLineArgument("npy")) {
                    nProcY = GetCmdLineArgument<int>("npy");
                }
                if (DefinesCmdLineArgument("npz")) {
                    nProcZ = GetCmdLineArgument<int>("npz");
                }
                if (DefinesCmdLineArgument("nsz")) {
                    nStripZ = GetCmdLineArgument<int>("nsz");
                }
                ASSERTL0(m_comm->GetSize() % (nProcZ*nProcY*nProcX) == 0,
                         "Cannot exactly partition using PROC_Z value.");
                ASSERTL0(nProcZ % nProcY == 0,
                         "Cannot exactly partition using PROC_Y value.");
                ASSERTL0(nProcY % nProcX == 0,
                         "Cannot exactly partition using PROC_X value.");

                // Number of processes associated with the spectral method
                int nProcSm  = nProcZ * nProcY * nProcX;

                // Number of processes associated with the spectral element
                // method.
                int nProcSem = m_comm->GetSize() / nProcSm;
                
                if (vNumThr > 1)
                {
                	CommSharedPtr vTmpComm(new ThreadedComm(m_comm, m_threadManager));
                	m_comm = vTmpComm;
                	threadedCommDone = true;
                }

                m_comm->SplitComm(nProcSm,nProcSem);
                m_comm->GetColumnComm()->SplitComm(nProcZ/nStripZ,nStripZ);
                m_comm->GetColumnComm()->GetColumnComm()->SplitComm(
                                            (nProcY*nProcX),nProcZ/nStripZ);
                m_comm->GetColumnComm()->GetColumnComm()->GetColumnComm()
                                            ->SplitComm(nProcX,nProcY);
            }

            /*
            if (!threadedCommDone && vNumThr > 1)
            {
            	CommSharedPtr vTmpComm(new ThreadedComm(m_comm, m_threadManager));
            	m_comm = vTmpComm;
            }
            */
        }


        /**
         *
         */
        void SessionReader::ReadParameters(TiXmlElement *conditions)
        {
            unsigned int vThr = Nektar::Thread::GetThreadMaster().GetInstance(Nektar::Thread::ThreadMaster::SessionJob)->GetWorkerNum();


            if (!conditions)
            {
                return;
            }

            TiXmlElement *parametersElement = conditions->FirstChildElement(
                "PARAMETERS");

            // See if we have parameters defined.  They are optional so we go on
            // if not.
            if (parametersElement)
            {
                TiXmlElement *parameter = 
                    parametersElement->FirstChildElement("P");

                ParameterMap caseSensitiveParameters;

                // Multiple nodes will only occur if there is a comment in
                // between definitions.
                while (parameter)
                {
                    stringstream tagcontent;
                    tagcontent << *parameter;
                    TiXmlNode *node = parameter->FirstChild();

                    while (node && node->Type() != TiXmlNode::TINYXML_TEXT)
                    {
                        node = node->NextSibling();
                    }

                    if (node)
                    {
                        // Format is "paramName = value"
                        std::string line = node->ToText()->Value(), lhs, rhs;

                        try {
                            ParseEquals(line, lhs, rhs);
                        }
                        catch (...)
                        {
                            ASSERTL0(false, "Syntax error in parameter "
                                     "expression '" + line 
                                     + "' in XML element: \n\t'"
                                     + tagcontent.str() + "'");
                        }

                        // We want the list of parameters to have their RHS
                        // evaluated, so we use the expression evaluator to do
                        // the dirty work.
                        if (!lhs.empty() && !rhs.empty())
                        {
                            NekDouble value=0.0;
                            try
                            {
                                LibUtilities::Equation expession(
                                    GetSharedThisPtr(), rhs);
                                value = expession.Evaluate();
                            }
                            catch (const std::runtime_error &)
                            {
                                ASSERTL0(false, 
                                         "Error evaluating parameter expression"
                                         " '" + rhs + "' in XML element: \n\t'"
                                         + tagcontent.str() + "'");
                            }
                            m_exprEvaluator[vThr]->SetParameter(lhs, value);
                            caseSensitiveParameters[lhs] = value;
                            boost::to_upper(lhs);
                            m_parameters[vThr][lhs] = value;
                        }
                    }

                    parameter = parameter->NextSiblingElement();
                }
            }
        }


        /**
         *
         */
        void SessionReader::ReadSolverInfo(TiXmlElement *conditions)
        {
            unsigned int vThr;
            vThr     = Nektar::Thread::GetThreadMaster().GetInstance(Nektar::Thread::ThreadMaster::SessionJob)->GetWorkerNum();

            m_solverInfo.clear();
            m_solverInfo[vThr] = GetSolverInfoDefaults();

            if (!conditions)
            {
                return;
            }

            TiXmlElement *solverInfoElement = 
                conditions->FirstChildElement("SOLVERINFO");

            if (solverInfoElement)
            {
                TiXmlElement *solverInfo = 
                    solverInfoElement->FirstChildElement("I");

                while (solverInfo)
                {
                    std::stringstream tagcontent;
                    tagcontent << *solverInfo;
                    // read the property name
                    ASSERTL0(solverInfo->Attribute("PROPERTY"),
                             "Missing PROPERTY attribute in solver info "
                             "XML element: \n\t'" + tagcontent.str() + "'");
                    std::string solverProperty = 
                        solverInfo->Attribute("PROPERTY");
                    ASSERTL0(!solverProperty.empty(),
                             "PROPERTY attribute must be non-empty in XML "
                             "element: \n\t'" + tagcontent.str() + "'");

                    // make sure that solver property is capitalised
                    std::string solverPropertyUpper =
                        boost::to_upper_copy(solverProperty);

                    // read the value
                    ASSERTL0(solverInfo->Attribute("VALUE"),
                            "Missing VALUE attribute in solver info "
                            "XML element: \n\t'" + tagcontent.str() + "'");
                    std::string solverValue    = solverInfo->Attribute("VALUE");
                    ASSERTL0(!solverValue.empty(),
                             "VALUE attribute must be non-empty in XML "
                             "element: \n\t'" + tagcontent.str() + "'");

                    EnumMapList::const_iterator propIt = 
                        GetSolverInfoEnums().find(solverPropertyUpper);
                    if (propIt != GetSolverInfoEnums().end())
                    {
                        EnumMap::const_iterator valIt = 
                            propIt->second.find(solverValue);
                        ASSERTL0(valIt != propIt->second.end(),
                                 "Value '" + solverValue + "' is not valid for "
                                 "property '" + solverProperty + "'");
                    }

                    // Set Variable
                    m_solverInfo[vThr][solverPropertyUpper] = solverValue;
                    solverInfo = solverInfo->NextSiblingElement("I");
                }
            }
            
            if (m_comm && m_comm->GetRowComm()->GetSize() > 1)
            {
                ASSERTL0(
                    m_solverInfo[vThr]["GLOBALSYSSOLN"] == "IterativeFull"       ||
                    m_solverInfo[vThr]["GLOBALSYSSOLN"] == "IterativeStaticCond" ||
                    m_solverInfo[vThr]["GLOBALSYSSOLN"] ==
                        "IterativeMultiLevelStaticCond"                          ||
                    m_solverInfo[vThr]["GLOBALSYSSOLN"] == "XxtFull"             ||
                    m_solverInfo[vThr]["GLOBALSYSSOLN"] == "XxtStaticCond"       ||
                    m_solverInfo[vThr]["GLOBALSYSSOLN"] ==
                        "XxtMultiLevelStaticCond"                                ||
                    m_solverInfo[vThr]["GLOBALSYSSOLN"] == "PETScFull"           ||
                    m_solverInfo[vThr]["GLOBALSYSSOLN"] == "PETScStaticCond"     ||
                    m_solverInfo[vThr]["GLOBALSYSSOLN"] ==
                        "PETScMultiLevelStaticCond",
                    "A parallel solver must be used when run in parallel.");
            }
        }



        /**
         *
         */
        void SessionReader::ReadGlobalSysSolnInfo(TiXmlElement *conditions)
        {
        	// NB, this is only supposed to run on the master thread as it
        	// is setting a static.
            GetGloSysSolnList().clear();

            if (!conditions)
            {
                return;
            }

            TiXmlElement *GlobalSys =
                            conditions->FirstChildElement("GLOBALSYSSOLNINFO");

            if(!GlobalSys)
            {
                return;
            }

            TiXmlElement *VarInfo   = GlobalSys->FirstChildElement("V");

            while (VarInfo)
            {
                std::stringstream tagcontent;
                tagcontent << *VarInfo;
                ASSERTL0(VarInfo->Attribute("VAR"),
                         "Missing VAR attribute in GobalSysSolnInfo XML "
                         "element: \n\t'" + tagcontent.str() + "'");

                std::string VarList = VarInfo->Attribute("VAR");
                ASSERTL0(!VarList.empty(),
                         "VAR attribute must be non-empty in XML element:\n\t'"
                         + tagcontent.str() + "'");

                // generate a list of variables.
                std::vector<std::string> varStrings;
                bool valid = ParseUtils::GenerateOrderedStringVector(
                                                VarList.c_str(),varStrings);

                ASSERTL0(valid,"Unable to process list of variable in XML "
                               "element \n\t'" + tagcontent.str() + "'");

                if(varStrings.size())
                {
                    TiXmlElement *SysSolnInfo = VarInfo->FirstChildElement("I");

                    while (SysSolnInfo)
                    {
                        tagcontent.clear();
                        tagcontent << *SysSolnInfo;
                        // read the property name
                        ASSERTL0(SysSolnInfo->Attribute("PROPERTY"),
                                 "Missing PROPERTY attribute in "
                                 "GlobalSysSolnInfo for variable(s) '"
                                 + VarList + "' in XML element: \n\t'" 
                                 + tagcontent.str() + "'");

                        std::string SysSolnProperty =
                            SysSolnInfo->Attribute("PROPERTY");

                        ASSERTL0(!SysSolnProperty.empty(),
                                 "GlobalSysSolnIno properties must have a "
                                 "non-empty name for variable(s) : '"
                                 + VarList + "' in XML element: \n\t'"
                                 + tagcontent.str() + "'");

                        // make sure that solver property is capitalised
                        std::string SysSolnPropertyUpper =
                                        boost::to_upper_copy(SysSolnProperty);

                        // read the value
                        ASSERTL0(SysSolnInfo->Attribute("VALUE"),
                                 "Missing VALUE attribute in GlobalSysSolnInfo "
                                 "for variable(s) '" + VarList
                                 + "' in XML element: \n\t"
                                 + tagcontent.str() + "'");

                        std::string SysSolnValue =
                                            SysSolnInfo->Attribute("VALUE");
                        ASSERTL0(!SysSolnValue.empty(),
                                 "GlobalSysSolnInfo properties must have a "
                                 "non-empty value for variable(s) '"
                                 + VarList + "' in XML element: \n\t'"
                                 + tagcontent.str() + "'");

                        // Store values under variable map.
                        for(int i = 0; i < varStrings.size(); ++i)
                        {
                            GloSysSolnInfoList::iterator x;
                            if ((x = GetGloSysSolnList().find(varStrings[i])) ==
                                    GetGloSysSolnList().end())
                            {
                                (GetGloSysSolnList()[varStrings[i]])[
                                        SysSolnPropertyUpper] = SysSolnValue;
                            }
                            else
                            {
                                x->second[SysSolnPropertyUpper] = SysSolnValue;
                            }
                        }

                        SysSolnInfo = SysSolnInfo->NextSiblingElement("I");
                    }
                    VarInfo = VarInfo->NextSiblingElement("V");
                }
            }

            if (m_verbose && GetGloSysSolnList().size() > 0 && m_comm)
            {
                if(m_comm->GetRank() == 0)
                {
                    cout << "GlobalSysSoln Info:" << endl;

                    GloSysSolnInfoList::iterator x;
                    for (x = GetGloSysSolnList().begin();
                         x != GetGloSysSolnList().end();
                         ++x)
                    {
                        cout << "\t Variable: " << x->first <<  endl;

                        GloSysInfoMap::iterator y;
                        for (y = x->second.begin(); y != x->second.end(); ++y)
                        {
                            cout << "\t\t " << y->first  << " = " << y->second
                                 << endl;
                        }
                    }
                    cout << endl;
                }
            }
        }


        /**
         *
         */
        void SessionReader::ReadExpressions(TiXmlElement *conditions)
        {
            unsigned int vThr = m_threadManager->GetWorkerNum();

            if (!conditions)
            {
                return;
            }

            TiXmlElement *expressionsElement = 
                conditions->FirstChildElement("EXPRESSIONS");

            if (expressionsElement)
            {
                TiXmlElement *expr = expressionsElement->FirstChildElement("E");

                while (expr)
                {
                    stringstream tagcontent;
                    tagcontent << *expr;
                    ASSERTL0(expr->Attribute("NAME"),
                             "Missing NAME attribute in expression "
                             "definition: \n\t'" + tagcontent.str() + "'");
                    std::string nameString = expr->Attribute("NAME");
                    ASSERTL0(!nameString.empty(),
                             "Expressions must have a non-empty name: \n\t'"
                             + tagcontent.str() + "'");

                    ASSERTL0(expr->Attribute("VALUE"),
                             "Missing VALUE attribute in expression "
                             "definition: \n\t'" + tagcontent.str() + "'");
                    std::string valString = expr->Attribute("VALUE");
                    ASSERTL0(!valString.empty(),
                             "Expressions must have a non-empty value: \n\t'"
                            + tagcontent.str() + "'");

                    ExpressionMap::iterator exprIter
                                            = m_expressions[vThr].find(nameString);
                    ASSERTL0(exprIter == m_expressions[vThr].end(),
                             std::string("Expression '") + nameString
                             + std::string("' already specified."));

                    m_expressions[vThr][nameString] = valString;
                    expr = expr->NextSiblingElement("E");
                }
            }
        }


        /**
         *
         */
        void SessionReader::ReadVariables(TiXmlElement *conditions)
        {
            unsigned int vThr = m_threadManager->GetWorkerNum();

            if (!conditions)
            {
                return;
            }

            TiXmlElement *variablesElement = 
                conditions->FirstChildElement("VARIABLES");

            // See if we have parameters defined. They are optional so we go on
            // if not.
            if (variablesElement)
            {
                TiXmlElement *varElement = 
                    variablesElement->FirstChildElement("V");

                // Sequential counter for the composite numbers.
                int nextVariableNumber = -1;

                while (varElement)
                {
                    stringstream tagcontent;
                    tagcontent << *varElement;

                    /// All elements are of the form: "<V ID="#"> name = value
                    /// </V>", with ? being the element type.
                    nextVariableNumber++;

                    int i;
                    int err = varElement->QueryIntAttribute("ID", &i);
                    ASSERTL0(err == TIXML_SUCCESS,
                             "Variables must have a unique ID number attribute "
                             "in XML element: \n\t'" + tagcontent.str() + "'");
                    ASSERTL0(i == nextVariableNumber,
                             "ID numbers for variables must begin with zero and"
                             " be sequential in XML element: \n\t'"
                             + tagcontent.str() + "'");

                    TiXmlNode* varChild = varElement->FirstChild();
                    // This is primarily to skip comments that may be present.
                    // Comments appear as nodes just like elements.  We are
                    // specifically looking for text in the body of the
                    // definition.
                    while(varChild && varChild->Type() != TiXmlNode::TINYXML_TEXT)
                    {
                        varChild = varChild->NextSibling();
                    }

                    ASSERTL0(varChild,
                             "Unable to read variable definition body for "
                             "variable with ID "
                             + boost::lexical_cast<string>(i) 
                             + " in XML element: \n\t'"
                             + tagcontent.str() + "'");
                    std::string variableName = varChild->ToText()->ValueStr();

                    std::istringstream variableStrm(variableName);
                    variableStrm >> variableName;

                    ASSERTL0(std::find(m_variables[vThr].begin(), m_variables[vThr].end(), 
                                       variableName) == m_variables[vThr].end(),
                             "Variable with ID "
                             + boost::lexical_cast<string>(i) 
                             + " in XML element \n\t'" + tagcontent.str()
                             + "'\nhas already been defined.");

                    m_variables[vThr].push_back(variableName);

                    varElement = varElement->NextSiblingElement("V");
                }

                ASSERTL0(nextVariableNumber > -1,
                         "Number of variables must be greater than zero.");
            }
        }


        /**
         *
         */
        void SessionReader::ReadFunctions(TiXmlElement *conditions)
        {
            unsigned int vThr = m_threadManager->GetWorkerNum();

            if (!conditions)
            {
                return;
            }

            // Scan through conditions section looking for functions.
            TiXmlElement *function = conditions->FirstChildElement("FUNCTION");
            while (function)
            {
                stringstream tagcontent;
                tagcontent << *function;

                // Every function must have a NAME attribute
                ASSERTL0(function->Attribute("NAME"),
                         "Functions must have a NAME attribute defined in XML "
                         "element: \n\t'" + tagcontent.str() + "'");
                std::string functionStr = function->Attribute("NAME");
                ASSERTL0(!functionStr.empty(),
                         "Functions must have a non-empty name in XML "
                         "element: \n\t'" + tagcontent.str() + "'");

                // Store function names in uppercase to remain case-insensitive.
                boost::to_upper(functionStr);

                // Retrieve first entry (variable, or file)
                TiXmlElement *variable  = function->FirstChildElement();

                // Create new function structure with default type of none.
                FunctionVariableMap functionVarMap;

                // Process all entries in the function block
                while (variable)
                {
                    FunctionVariableDefinition funcDef;
                    std::string conditionType = variable->Value();

                    // If no var is specified, assume wildcard
                    std::string variableStr;
                    if (!variable->Attribute("VAR"))
                    {
                        variableStr = "*";
                    }
                    else
                    {
                        variableStr = variable->Attribute("VAR");
                    }

                    // Parse list of variables
                    std::vector<std::string> variableList;
                    ParseUtils::GenerateOrderedStringVector(variableStr.c_str(),
                                                            variableList);

                    // If no domain string put to 0
                    std::string domainStr;
                    if (!variable->Attribute("DOMAIN"))
                    {
                        domainStr = "0";
                    }
                    else
                    {
                        domainStr = variable->Attribute("DOMAIN");
                    }

                    // Parse list of variables
                    std::vector<std::string> varSplit;
                    std::vector<unsigned int> domainList;
                    ParseUtils::GenerateSeqVector(domainStr.c_str(), domainList);

                    // Expressions are denoted by E
                    if (conditionType == "E")
                    {
                        funcDef.m_type = eFunctionTypeExpression;

                        // Expression must have a VALUE.
                        ASSERTL0(variable->Attribute("VALUE"),
                                 "Attribute VALUE expected for function '"
                                 + functionStr + "'.");
                        std::string fcnStr = variable->Attribute("VALUE");

                        ASSERTL0(!fcnStr.empty(),
                                 (std::string("Expression for var: ")
                                 + variableStr
                                 + std::string(" must be specified.")).c_str());

                        SubstituteExpressions(fcnStr);

                        // set expression
                        funcDef.m_expression = MemoryManager<Equation>
                            ::AllocateSharedPtr(GetSharedThisPtr(),fcnStr);
                    }

                    // Files are denoted by F
                    else if (conditionType == "F")
                    {
                        if (variable->Attribute("TIMEDEPENDENT") &&
                            boost::lexical_cast<bool>(variable->Attribute("TIMEDEPENDENT")))
                        {
                            funcDef.m_type = eFunctionTypeTransientFile;
                        }
                        else
                        {
                            funcDef.m_type = eFunctionTypeFile;
                        }

                        // File must have a FILE.
                        ASSERTL0(variable->Attribute("FILE"),
                                 "Attribute FILE expected for function '"
                                 + functionStr + "'.");
                        std::string filenameStr = variable->Attribute("FILE");

                        ASSERTL0(!filenameStr.empty(),
                                 "A filename must be specified for the FILE "
                                 "attribute of function '" + functionStr
                                 + "'.");

                        std::vector<std::string> fSplit;
                        boost::split(fSplit, filenameStr, boost::is_any_of(":"));

                        ASSERTL0(fSplit.size() == 1 || fSplit.size() == 2,
                                 "Incorrect filename specification in function "
                                 + functionStr + "'. "
                                 "Specify variables inside file as: "
                                 "filename:var1,var2");

                        // set the filename
                        funcDef.m_filename = fSplit[0];

                        if (fSplit.size() == 2)
                        {
                            ASSERTL0(variableList[0] != "*",
                                     "Filename variable mapping not valid "
                                     "when using * as a variable inside "
                                     "function '" + functionStr + "'.");

                            boost::split(
                                varSplit, fSplit[1], boost::is_any_of(","));
                            ASSERTL0(varSplit.size() == variableList.size(),
                                     "Filename variables should contain the "
                                     "same number of variables defined in "
                                     "VAR in function " + functionStr + "'.");
                        }
                    }

                    // Nothing else supported so throw an error
                    else
                    {
                        stringstream tagcontent;
                        tagcontent << *variable;

                        ASSERTL0(false,
                                "Identifier " + conditionType + " in function "
                                + std::string(function->Attribute("NAME"))
                                + " is not recognised in XML element: \n\t'"
                                + tagcontent.str() + "'");
                    }


                    
                    // Add variables to function
                    for (unsigned int i = 0; i < variableList.size(); ++i)
                    {
                        for(unsigned int j = 0; j < domainList.size(); ++j)
                        {
                            // Check it has not already been defined
                            pair<std::string,int> key(variableList[i],domainList[j]);
                            FunctionVariableMap::iterator fcnsIter
                                = functionVarMap.find(key);
                            ASSERTL0(fcnsIter == functionVarMap.end(),
                                     "Error setting expression '" + variableList[i]
                                     + " in domain " 
                                     + boost::lexical_cast<std::string>(domainList[j]) 
                                     + "' in function '" + functionStr + "'. "
                                     "Expression has already been defined.");

                            if (varSplit.size() > 0)
                            {
                                FunctionVariableDefinition funcDef2 = funcDef;
                                funcDef2.m_fileVariable = varSplit[i];
                                functionVarMap[key] = funcDef2;
                            }
                            else
                            {
                                functionVarMap[key] = funcDef;
                            }
                        }
                    }
                    
                    variable = variable->NextSiblingElement();
                }
                // Add function definition to map
                m_functions[vThr][functionStr] = functionVarMap;
                function = function->NextSiblingElement("FUNCTION");
            }
        }


        /**
         *
         */
        void SessionReader::ReadFilters(TiXmlElement *filters)
        {
            if (!filters)
            {
                return;
            }

            unsigned int vThr = m_threadManager->GetWorkerNum();

            TiXmlElement *filter = filters->FirstChildElement("FILTER");
            while (filter)
            {
                ASSERTL0(filter->Attribute("TYPE"),
                        "Missing attribute 'TYPE' for filter.");
                std::string typeStr = filter->Attribute("TYPE");

                std::map<std::string, std::string> vParams;

                TiXmlElement *param = filter->FirstChildElement("PARAM");
                while (param)
                {
                    ASSERTL0(param->Attribute("NAME"),
                            "Missing attribute 'NAME' for parameter in filter "
                            + typeStr + "'.");
                    std::string nameStr = param->Attribute("NAME");

                    ASSERTL0(param->GetText(), "Empty value string for param.");
                    std::string valueStr = param->GetText();

                    vParams[nameStr] = valueStr;

                    param = param->NextSiblingElement("PARAM");
                }

                m_filters[vThr].push_back(
                    std::pair<std::string, FilterParams>(typeStr, vParams));

                filter = filter->NextSiblingElement("FILTER");
            }
        }
        
        void SessionReader::ParseEquals(
            const std::string &line,
                  std::string &lhs,
                  std::string &rhs)
        {
            /// Pull out lhs and rhs and eliminate any spaces.
            int beg = line.find_first_not_of(" ");
            int end = line.find_first_of("=");
            // Check for no parameter name
            if (beg == end) throw 1;
            // Check for no parameter value
            if (end != line.find_last_of("=")) throw 1;
            // Check for no equals sign
            if (end == std::string::npos) throw 1;
            
            lhs = line.substr(line.find_first_not_of(" "),
                              end-beg);
            lhs = lhs .substr(0, lhs.find_last_not_of(" ")+1);
            rhs = line.substr(line.find_last_of("=")+1);
            rhs = rhs .substr(rhs.find_first_not_of(" "));
            rhs = rhs .substr(0, rhs.find_last_not_of(" ")+1);
        }
        
        /**
         *
         */
        void SessionReader::CmdLineOverride()
        {
            unsigned int vThr = m_threadManager->GetWorkerNum();
            // Parse solver info overrides
            if (m_cmdLineOptions.count("solverinfo"))
            {
                std::vector<std::string> solverInfoList = 
                    m_cmdLineOptions["solverinfo"].as<
                        std::vector<std::string> >();
                
                for (int i = 0; i < solverInfoList.size(); ++i)
                {
                    std::string lhs, rhs;

                    try
                    {
                        ParseEquals(solverInfoList[i], lhs, rhs);
                    } 
                    catch (...)
                    {
                        ASSERTL0(false, "Parse error with command line "
                                 "option: "+solverInfoList[i]);
                    }

                    std::string lhsUpper = boost::to_upper_copy(lhs);
                    m_solverInfo[vThr][lhsUpper] = rhs;
                }
            }
            
            if (m_cmdLineOptions.count("parameter"))
            {
                std::vector<std::string> parametersList = 
                    m_cmdLineOptions["parameter"].as<
                        std::vector<std::string> >();
                
                for (int i = 0; i < parametersList.size(); ++i)
                {
                    std::string lhs, rhs;

                    try
                    {
                        ParseEquals(parametersList[i], lhs, rhs);
                    } 
                    catch (...)
                    {
                        ASSERTL0(false, "Parse error with command line "
                                 "option: "+parametersList[i]);
                    }

                    std::string lhsUpper = boost::to_upper_copy(lhs);
                    
                    try
                    {
                        m_parameters[vThr][lhsUpper] =
                            boost::lexical_cast<NekDouble>(rhs);
                    }
                    catch (...)
                    {
                        ASSERTL0(false, "Unable to convert string: "+rhs+
                                 "to double value.");
                    }
                }
            }
        }

        /**
         *
         */
        void SessionReader::StartThreads()
        {
            int nthreads;
            LoadParameter("NThreads", nthreads, 1);

            if (nthreads > 1 && m_mainFunc == 0)
            {
                cerr << "This solver has not been set up to do threaded"
                    " parallelism (no MainFunc() passed to SessionReader())."
                    << endl << "Number of threads set to 1." << endl;
                nthreads = 1;
            }

            // Decide what implementation of ThreadManager you want here.
            ThreadMaster &vTM = Thread::GetThreadMaster();
            vTM.SetThreadingType("ThreadManagerBoost"); // or whatever
            // This ThreadManager will run the SessionJob jobs.
            // It runs the ThreadedComm class.
            m_threadManager = vTM.CreateInstance(Nektar::Thread::ThreadMaster::SessionJob, nthreads);
            cerr << "Number of threads will be: " << nthreads << endl;

            Nektar::InitMemoryPools(nthreads, true);
            m_xmlDoc.resize(nthreads);
            m_parameters[0].clear();
            m_parameters.resize(nthreads);
            m_solverInfo[0].clear();
            m_solverInfo.resize(nthreads);
            m_expressions[0].clear();
            m_expressions.resize(nthreads);
            m_variables[0].clear();
            m_variables.resize(nthreads);
            m_functions[0].clear();
            m_functions.resize(nthreads);
            m_exprEvaluator.resize(nthreads);
            for (unsigned int i=1; i<nthreads; i++)
            {
                m_exprEvaluator[i] = new AnalyticExpressionEvaluator();
            }
            m_geometricInfo[0].clear();
            m_geometricInfo.resize(nthreads);
            m_filters[0].clear();
            m_filters.resize(nthreads);
            m_tags.resize(nthreads);
            m_bndRegOrder[0].clear();
            m_bndRegOrder.resize(nthreads);

        }

        Nektar::Thread::ThreadManagerSharedPtr SessionReader::GetThreadManager()
        {
        	return m_threadManager;
        }

        void SessionReader::SetUpXmlDoc(void)
        {
            m_xmlDoc[Nektar::Thread::GetThreadMaster().GetInstance(Nektar::Thread::ThreadMaster::SessionJob)->GetWorkerNum()]
                = MergeDoc(m_filenames);
        }
    }
}
