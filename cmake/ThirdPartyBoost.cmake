OPTION(THIRDPARTY_BUILD_BOOST "Build Boost libraries" OFF)
SET(Boost_USE_MULTITHREADED ON CACHE BOOL
    "Search for multithreaded boost libraries")
MARK_AS_ADVANCED(Boost_USE_MULTITHREADED)

IF (THIRDPARTY_BUILD_BOOST)
    INCLUDE(ExternalProject)
   
    IF (NOT WIN32)
        # Only build the libraries we need
        SET(BOOST_LIB_LIST --with-system --with-iostreams --with-filesystem 
                           --with-program_options --with-date_time --with-thread
                           --with-regex)
                           
        # We need -fPIC for 64-bit builds
        IF( CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )
            SET(BOOST_FLAGS cxxflags=-fPIC cflags=-fPIC linkflags=-fPIC)
        ENDIF ()
        
        # Build Boost
        IF (APPLE)
            SET(TOOLSET darwin)
        ELSE(APPLE)
            SET(TOOLSET gcc)
        ENDIF(APPLE)
        EXTERNALPROJECT_ADD(
            boost
            PREFIX ${TPSRC}
            URL ${TPURL}/boost_1_55_0.tar.bz2
            URL_MD5 "d6eef4b4cacb2183f2bf265a5a03a354"
            STAMP_DIR ${TPBUILD}/stamp
            DOWNLOAD_DIR ${TPSRC}
            SOURCE_DIR ${TPBUILD}/boost
            BINARY_DIR ${TPBUILD}/boost
            TMP_DIR ${TPBUILD}/boost-tmp
            INSTALL_DIR ${TPDIST}
            CONFIGURE_COMMAND CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} ./bootstrap.sh --prefix=${TPDIST}
            BUILD_COMMAND NO_BZIP2=1 ./b2
                            variant=release
                            link=shared 
                            include=${TPDIST}/include
                            linkflags="-L${TPDIST}/lib"
                            ${BOOST_FLAGS} ${BOOST_LIB_LIST} 
                            --layout=system toolset=${TOOLSET} install
            INSTALL_COMMAND ""
        )
    
        IF (APPLE)
            EXTERNALPROJECT_ADD_STEP(boost patch-install-path
                COMMAND sed -i ".bak" "s|-install_name \"|&${TPDIST}/lib/|" ${TPSRC}/boost/tools/build/v2/tools/darwin.jam
                DEPENDERS build
                DEPENDEES download)
        ENDIF (APPLE)

        # If building ThirdParty zlib, force zlib build before boost
        IF (THIRDPARTY_BUILD_ZLIB)
             ADD_DEPENDENCIES(boost zlib-1.2.7)
        ENDIF(THIRDPARTY_BUILD_ZLIB)

        # Set up CMake variables
        SET(Boost_DATE_TIME_LIBRARY boost_date_time)
        SET(Boost_DATE_TIME_LIBRARY_DEBUG boost_date_time)
        SET(Boost_DATE_TIME_LIBRARY_RELEASE boost_date_time)
        SET(Boost_FILESYSTEM_LIBRARY boost_filesystem)
        SET(Boost_FILESYSTEM_LIBRARY_DEBUG boost_filesystem)
        SET(Boost_FILESYSTEM_LIBRARY_RELEASE boost_filesystem)
        SET(Boost_IOSTREAMS_LIBRARY boost_iostreams)
        SET(Boost_IOSTREAMS_LIBRARY_DEBUG boost_iostreams)
        SET(Boost_IOSTREAMS_LIBRARY_RELEASE boost_iostreams)
        SET(Boost_PROGRAM_OPTIONS_LIBRARY boost_program_options)
        SET(Boost_PROGRAM_OPTIONS_LIBRARY_DEBUG boost_program_options)
        SET(Boost_PROGRAM_OPTIONS_LIBRARY_RELEASE boost_program_options)
        SET(Boost_REGEX_LIBRARY boost_regex)
        SET(Boost_REGEX_LIBRARY_DEBUG boost_regex)
        SET(Boost_REGEX_LIBRARY_RELEASE boost_regex)
        SET(Boost_SYSTEM_LIBRARY boost_system)
        SET(Boost_SYSTEM_LIBRARY_DEBUG boost_system)
        SET(Boost_SYSTEM_LIBRARY_RELEASE boost_system)
        SET(Boost_THREAD_LIBRARY boost_thread)
        SET(Boost_THREAD_LIBRARY_DEBUG boost_thread)
        SET(Boost_THREAD_LIBRARY_RELEASE boost_thread)
        SET(Boost_INCLUDE_DIRS ${TPSRC}/dist/include 
                               ${TPSRC}/dist/include/boost-1_49)
        SET(Boost_LIBRARY_DIRS ${TPSRC}/dist/lib)
        LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})
    ELSE ()
        EXTERNALPROJECT_ADD(
            boost
            PREFIX ${TPSRC}
            URL ${TPURL}/boost_1_55_0.tar.bz2
            URL_MD5 "d6eef4b4cacb2183f2bf265a5a03a354"
            STAMP_DIR ${TPSRC}/stamp
            DOWNLOAD_DIR ${TPSRC}
            SOURCE_DIR ${TPSRC}/boost
            BINARY_DIR ${TPBUILD}/boost
            TMP_DIR ${TPBUILD}/boost-tmp
            INSTALL_DIR ${TPDIST}
            CONFIGURE_COMMAND bootstrap.bat --prefix=${TPDIST}/boost
            BUILD_COMMAND b2 --layout=system install
            INSTALL_COMMAND ""
        )
    ENDIF ()
ELSE (THIRDPARTY_BUILD_BOOST)
    ADD_CUSTOM_TARGET(boost ALL)
    SET(Boost_DEBUG 0)
    SET(Boost_NO_BOOST_CMAKE ON)
    #If the user has not set BOOST_ROOT, look in a couple common places first.
    IF( NOT BOOST_ROOT )
        SET(TEST_ENV1 $ENV{BOOST_HOME})
        SET(TEST_ENV2 $ENV{BOOST_DIR})
        IF (DEFINED TEST_ENV1)
            SET(BOOST_ROOT $ENV{BOOST_HOME})
            FIND_PACKAGE( Boost QUIET COMPONENTS thread iostreams date_time
                filesystem system program_options regex )
        ELSEIF (DEFINED TEST_ENV2)
            SET(BOOST_ROOT $ENV{BOOST_DIR})
            FIND_PACKAGE( Boost QUIET COMPONENTS thread iostreams date_time
                filesystem system program_options regex )
        ELSE ()
            SET(BOOST_ROOT ${TPDIST})
            FIND_PACKAGE( Boost QUIET COMPONENTS thread iostreams date_time filesystem system program_options regex)
        ENDIF()
    ELSE()
        FIND_PACKAGE( Boost COMPONENTS thread iostreams zlib date_time filesystem system program_options regex)
    ENDIF()
    IF (Boost_IOSTREAMS_FOUND)
        MESSAGE(STATUS "Found Boost iostreams library: ${Boost_IOSTREAMS_LIBRARY}")
    ENDIF ()
    IF (Boost_THREAD_FOUND)
        MESSAGE(STATUS "Found Boost thread library: ${Boost_THREAD_LIBRARY}")
    ENDIF ()
    IF (Boost_DATE_TIME_FOUND)
        MESSAGE(STATUS "Found Boost date_time library: ${Boost_DATE_TIME_LIBRARY}")
    ENDIF ()
    IF (Boost_FILESYSTEM_FOUND)
        MESSAGE(STATUS "Found Boost filesystem library: ${Boost_FILESYSTEM_LIBRARY}")
    ENDIF ()
    IF (Boost_SYSTEM_FOUND)
        MESSAGE(STATUS "Found Boost system library: ${Boost_SYSTEM_LIBRARY}")
    ENDIF ()
    IF (Boost_PROGRAM_OPTIONS_FOUND)
        MESSAGE(STATUS "Found Boost program_options library: ${Boost_PROGRAM_OPTIONS_LIBRARY}")
    ENDIF ()
    IF (Boost_REGEX_FOUND)
        MESSAGE(STATUS "Found Boost regex library: ${Boost_REGEX_LIBRARY}")
    ENDIF ()
    IF (BOOST_THREAD_LIBRARY)
        MARK_AS_ADVANCED(BOOST_THREAD_LIBRARY)
    ENDIF (BOOST_THREAD_LIBRARY)
ENDIF (THIRDPARTY_BUILD_BOOST)

INCLUDE_DIRECTORIES(SYSTEM ${Boost_INCLUDE_DIRS})
