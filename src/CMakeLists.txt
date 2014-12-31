SET(descriptors_src
    constitutional.cpp
    dmath.cpp
    electrostatic.cpp
    geometrical.cpp
    graph.cpp
    logP.cpp
    mymatrix.cpp
    names.cpp
    physicochemical.cpp
    qc.cpp
    tools.cpp
    topological.cpp)


ADD_LIBRARY(descriptors ${BUILD_TYPE} ${descriptors_src})

add_definitions(-DMAKE_DLL)
INCLUDE_DIRECTORIES("../include")
TARGET_LINK_LIBRARIES(descriptors ${OPENBABEL_LIB})

INSTALL(TARGETS descriptors 
    RUNTIME DESTINATION ${BIN_INSTALL_DIR} CONFIGURATIONS
    LIBRARY DESTINATION ${LIB_INSTALL_DIR} CONFIGURATIONS
    ARCHIVE DESTINATION ${LIB_INSTALL_DIR} CONFIGURATIONS
)
