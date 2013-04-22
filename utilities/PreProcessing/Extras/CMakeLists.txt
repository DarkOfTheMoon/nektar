SET(EnsiteToXmlSources EnsiteToXml.cpp)
ADD_UTILITIES_EXECUTABLE(EnsiteToXml util-extra ${EnsiteToXmlSources})

SET(CheckXmlFileSources CheckXmlFile.cpp)
ADD_UTILITIES_EXECUTABLE(CheckXmlFile util-extra ${CheckXmlFileSources})

IF (NEKTAR_USE_VTK)
    SET(VtkToFldSources VtkToFld.cpp)
    SET(ImageToFldSources ImageToFld.cpp)
    SET(VtkStripsToPolysSources VtkStripsToPolys.cpp)
    ADD_UTILITIES_EXECUTABLE(VtkToFld util-extra ${VtkToFldSources})
    ADD_UTILITIES_EXECUTABLE(ImageToFld util-extra ${ImageToFldSources})
    ADD_UTILITIES_EXECUTABLE(VtkStripsToPolys util-extra ${VtkStripsToPolysSources})
ENDIF (NEKTAR_USE_VTK)

