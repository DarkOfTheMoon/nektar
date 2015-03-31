SET(LinkLibraries SpatialDomains)


SET(SpatialDomains3DSource  Graph2Dto3D.cpp)
ADD_NEKTAR_EXECUTABLE(Graph2Dto3D demos SpatialDomains3DSource)
TARGET_LINK_LIBRARIES(Graph2Dto3D ${LinkLibraries})

