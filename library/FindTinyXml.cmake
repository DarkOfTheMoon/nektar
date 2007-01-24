
FIND_PATH(TINYXML_INCLUDE_DIR tinyxml/tinyxml.h /usr/include /usr/local/include 
${CMAKE_SOURCE_DIR}/../ThirdParty )


IF (TINYXML_INCLUDE_DIR)
  SET(TINYXML_FOUND TRUE)
ENDIF (TINYXML_INCLUDE_DIR)

SET (TINYXML_LIB_DIR ${TINYXML_INCLUDE_DIR}../lib )

IF (TINYXML_FOUND)
  IF (NOT TINYXML_FIND_QUIETLY)
     MESSAGE(STATUS "Found TinyXML: ${TINYXML_INCLUDE_DIR}")
  ENDIF (NOT TINYXML_FIND_QUIETLY)
ELSE(TINYXML_FOUND)
  IF (TINYXML_FIND_REQUIRED)
     MESSAGE(FATAL_ERROR "Could not find TinyXML")
  ENDIF (TINYXML_FIND_REQUIRED)
ENDIF (TINYXML_FOUND)
