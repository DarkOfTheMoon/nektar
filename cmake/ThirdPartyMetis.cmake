SET(THIRDPARTY_BUILD_METIS ON CACHE BOOL
    "Build ModMetis library from ThirdParty")

IF (THIRDPARTY_BUILD_METIS)
    INCLUDE(ExternalProject)
    EXTERNALPROJECT_ADD(
        modmetis-5.1.0
        PREFIX ${TPSRC}
        URL ${TPURL}/modmetis-5.1.0_1.tar.bz2
        URL_MD5 "6c6816aea0f53db6c71b1d98ed4ad42b"
        STAMP_DIR ${TPBUILD}/stamp
        DOWNLOAD_DIR ${TPSRC}
        SOURCE_DIR ${TPSRC}/modmetis-5.1.0
        BINARY_DIR ${TPBUILD}/modmetis-5.1.0
        TMP_DIR ${TPBUILD}/modmetis-5.1.0-tmp
        INSTALL_DIR ${TPDIST}
        CONFIGURE_COMMAND ${CMAKE_COMMAND}
            -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
            -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
            -DCMAKE_INSTALL_PREFIX:PATH=${TPDIST}
            -DCMAKE_C_FLAGS:STRING=-fPIC
            -DGKLIB_PATH:PATH=${TPSRC}/modmetis-5.1.0/GKlib
            ${TPSRC}/modmetis-5.1.0
    )
    SET(METIS_LIB metis CACHE FILEPATH
        "METIS library" FORCE)
    MARK_AS_ADVANCED(METIS_LIB)
    LINK_DIRECTORIES(${TPDIST}/lib)
    INCLUDE_DIRECTORIES(${TPDIST}/include)
    MESSAGE(STATUS "Build Metis: ${TPDIST}/lib/lib${METIS_LIB}.a")
ELSE (THIRDPARTY_BUILD_METIS)
    INCLUDE (FindMetis)
ENDIF (THIRDPARTY_BUILD_METIS)

