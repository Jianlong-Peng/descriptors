ADD_EXECUTABLE(calcDescriptors main.cpp)
ADD_EXECUTABLE(runMopac runMopac.cpp)

INCLUDE_DIRECTORIES("../include")
TARGET_LINK_LIBRARIES(calcDescriptors descriptors ${OPENBABEL_LIB})
TARGET_LINK_LIBRARIES(runMopac ${OPENBABEL_LIB})

INSTALL(TARGETS calcDescriptors runMopac RUNTIME DESTINATION bin CONFIGURATIONS Release)
INSTALL(PROGRAMS runMopac.py fingerprints.py DESTINATION bin)
