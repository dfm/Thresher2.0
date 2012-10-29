# Metadata for the build
CMAKE_MINIMUM_REQUIRED (VERSION 2.8.0)
PROJECT (Thresher C CXX)
SET (THRESHER_VERSION 0.0.1)

# Initialize the dependencies to an empty list.
SET (THRESHER_DEPS)

# Default search directories.
LIST (APPEND SEARCH_LIBS /usr/lib)
LIST (APPEND SEARCH_LIBS /usr/local/lib)
LIST (APPEND SEARCH_LIBS /usr/local/homebrew/lib)
LIST (APPEND SEARCH_LIBS /opt/local/lib)

LIST (APPEND SEARCH_HEADERS /usr/include)
LIST (APPEND SEARCH_HEADERS /usr/local/include)
LIST (APPEND SEARCH_HEADERS /usr/local/homebrew/include)
LIST (APPEND SEARCH_HEADERS /opt/local/include)


# ==========================================================================
#                                                                      CERES
# ==========================================================================

MESSAGE ("-- Check for Ceres")
SET (CERES_FOUND TRUE)

SET (CERES_SEARCH_LIBS SEARCH_LIBS)
LIST (APPEND CERES_SEARCH_LIBS /usr/lib/ceres)
LIST (APPEND CERES_SEARCH_LIBS /usr/local/lib/ceres)
LIST (APPEND CERES_SEARCH_LIBS /opt/local/lib)
FIND_LIBRARY (CERES_LIB NAMES ceres PATHS ${CERES_SEARCH_LIBS})
IF (NOT EXISTS ${CERES_LIB})
    MESSAGE(FATAL_ERROR "Can't find Ceres. Try passing -DCERES_LIB=...")
ENDIF (NOT EXISTS ${CERES_LIB})
MESSAGE("-- Found Ceres library: ${CERES_LIB}")
LIST (APPEND THRESHER_DEPS ${CERES_LIB})

SET (CERES_SEARCH_HEADERS SEARCH_HEADERS)
LIST (APPEND CERES_SEARCH_HEADERS /usr/include/ceres)
LIST (APPEND CERES_SEARCH_HEADERS /usr/local/include/ceres)
LIST (APPEND CERES_SEARCH_HEADERS /opt/local/include)
FIND_PATH (CERES_INCLUDE NAMES ceres/ceres.h PATHS ${CERES_SEARCH_HEADERS})
IF (NOT EXISTS ${CERES_INCLUDE})
    MESSAGE(FATAL_ERROR "Can't find Ceres headers. Try passing -DCERES_INCLUDE=...")
ENDIF (NOT EXISTS ${CERES_INCLUDE})
MESSAGE("-- Found Ceres headers: ${CERES_INCLUDE}")


# ==========================================================================
#                                                                    EIGEN 3
# ==========================================================================

SET (EIGEN_SEARCH_HEADERS ${SEARCH_HEADERS})
LIST (APPEND EIGEN_SEARCH_HEADERS /usr/include/eigen3) # Ubuntu 10.04
LIST (APPEND EIGEN_SEARCH_HEADERS /usr/local/include/eigen3)
LIST (APPEND EIGEN_SEARCH_HEADERS /usr/local/homebrew/include/eigen3)  # Mac
LIST (APPEND EIGEN_SEARCH_HEADERS /opt/local/var/macports/software/eigen3/opt/local/include/eigen3) # Mac

MESSAGE ("-- Check for Eigen 3.0")
FIND_PATH (EIGEN_INCLUDE NAMES Eigen/Core PATHS ${EIGEN_SEARCH_HEADERS})
IF (NOT EXISTS ${EIGEN_INCLUDE})
    MESSAGE(FATAL_ERROR "Can't find Eigen. Try passing -DEIGEN_INCLUDE=...")
ENDIF (NOT EXISTS ${EIGEN_INCLUDE})
MESSAGE("-- Found Eigen 3.0: ${EIGEN_INCLUDE}")


# ==========================================================================
#                                                             THE EXECUTABLE
# ==========================================================================

INCLUDE_DIRECTORIES (include
                     src
                     ${EIGEN_INCLUDE}
                     ${CERES_INCLUDE})

SET (THRESHER_SRC
     src/test.cc
     src/thresher.cc)

SET (THRESHER_EXE
     ${THRESHER_SRC}
     ${THRESHER_HDR})


ADD_EXECUTABLE (thresher ${THRESHER_EXE})
TARGET_LINK_LIBRARIES (thresher ${THRESHER_DEPS})
