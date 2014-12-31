#file(GLOB header_files *.h)

SET(header_files
    config.h
    constitutional.h
    electrostatic.h
    geometrical.h
    graph.h
    mymatrix.h
    names.h
    physicochemical.h
    qc.h
    tools.h
    topological.h)

INSTALL(FILES ${header_files} DESTINATION include)
