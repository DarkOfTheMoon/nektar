SET(LinkLibraries MultiRegions)

SET(CollectionTimingSource  CollectionTiming.cpp)
ADD_NEKTAR_EXECUTABLE(CollectionTiming demos CollectionTimingSource)
TARGET_LINK_LIBRARIES(CollectionTiming ${LinkLibraries})

SET(CollectionTiming2Source  CollectionTiming2.cpp)
ADD_NEKTAR_EXECUTABLE(CollectionTiming2 demos CollectionTiming2Source)
TARGET_LINK_LIBRARIES(CollectionTiming2 ${LinkLibraries})

IF (NEKTAR_USE_MPI)
#    ADD_NEKTAR_TEST(Helmholtz2D_CG_P7_Modes_AllBCs_xxt_full)
ENDIF (NEKTAR_USE_MPI)
