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
LIST (APPEND EIGEN_SEARCH_HEADERS /usr/include/eigen3)
LIST (APPEND EIGEN_SEARCH_HEADERS /usr/local/include/eigen3)
LIST (APPEND EIGEN_SEARCH_HEADERS /usr/local/homebrew/include/eigen3)
LIST (APPEND EIGEN_SEARCH_HEADERS /opt/local/var/macports/software/eigen3/opt/local/include/eigen3)

MESSAGE ("-- Check for Eigen 3.0")
FIND_PATH (EIGEN_INCLUDE NAMES Eigen/Core PATHS ${EIGEN_SEARCH_HEADERS})
IF (NOT EXISTS ${EIGEN_INCLUDE})
    MESSAGE(FATAL_ERROR "Can't find Eigen. Try passing -DEIGEN_INCLUDE=...")
ENDIF (NOT EXISTS ${EIGEN_INCLUDE})
MESSAGE("-- Found Eigen 3.0: ${EIGEN_INCLUDE}")


# ==========================================================================
#                                                                    CFITSIO
# ==========================================================================

MESSAGE ("-- Check for Fitsio")

SET (FITSIO_SEARCH_LIBS ${SEARCH_LIBS})
LIST (APPEND FITSIO_SEARCH_LIBS /usr/lib/cfitsio)
LIST (APPEND FITSIO_SEARCH_LIBS /usr/local/lib/cfitsio)
LIST (APPEND FITSIO_SEARCH_LIBS /usr/local/homebrew/lib/cfitsio)
FIND_LIBRARY (FITSIO_LIB NAMES cfitsio PATHS ${FITSIO_SEARCH_LIBS})
IF (NOT EXISTS ${FITSIO_LIB})
    MESSAGE(FATAL_ERROR "Can't find Fitsio library. Try passing -DFITSIO_LIB=...")
ENDIF (NOT EXISTS ${FITSIO_LIB})
LIST (APPEND THRESHER_DEPS ${FITSIO_LIB})
MESSAGE("-- Found Fitsio library: ${FITSIO_LIB}")

SET (FITSIO_SEARCH_HEADERS ${SEARCH_HEADERS})
LIST (APPEND FITSIO_SEARCH_HEADERS /usr/include/cfitsio)
LIST (APPEND FITSIO_SEARCH_HEADERS /usr/local/include/cfitsio)
LIST (APPEND FITSIO_SEARCH_HEADERS /usr/local/homebrew/include/cfitsio)
FIND_PATH (FITSIO_INCLUDE NAMES fitsio.h PATHS ${FITSIO_SEARCH_HEADERS})
IF (NOT EXISTS ${FITSIO_INCLUDE})
    MESSAGE(FATAL_ERROR "Can't find Fitsio headers. Try passing -DFITSIO_INCLUDE=...")
ENDIF (NOT EXISTS ${FITSIO_INCLUDE})
MESSAGE("-- Found Fitsio headers: ${FITSIO_INCLUDE}")


# ==========================================================================
#                                                                     GFLAGS
# ==========================================================================

MESSAGE("-- Check for Gflags")

SET (GFLAGS_SEARCH_LIBS ${SEARCH_LIBS})
LIST (APPEND GFLAGS_SEARCH_LIBS /usr/lib/gflags)
LIST (APPEND GFLAGS_SEARCH_LIBS /usr/local/lib/gflags)
LIST (APPEND GFLAGS_SEARCH_LIBS /usr/local/homebrew/lib/gflags)
FIND_LIBRARY(GFLAGS_LIB NAMES gflags PATHS ${GFLAGS_SEARCH_LIBS})
IF (NOT EXISTS ${GFLAGS_LIB})
    MESSAGE(FATAL_ERROR "Can't find Gflags library. Try passing -DGFLAGS_LIB=...")
ENDIF (NOT EXISTS ${GFLAGS_LIB})
LIST (APPEND THRESHER_DEPS ${GFLAGS_LIB})
MESSAGE("-- Found Gflags library: ${GFLAGS_LIB}")

SET (GFLAGS_SEARCH_HEADERS ${SEARCH_HEADERS})
LIST (APPEND GFLAGS_SEARCH_HEADERS /usr/include/gflags)
LIST (APPEND GFLAGS_SEARCH_HEADERS /usr/local/include/gflags)
LIST (APPEND GFLAGS_SEARCH_HEADERS /usr/local/homebrew/include/gflags)
FIND_PATH(GFLAGS_INCLUDE NAMES gflags/gflags.h PATHS ${GFLAGS_SEARCH_HEADERS})
IF (NOT EXISTS ${GFLAGS_INCLUDE})
    MESSAGE(FATAL_ERROR "Can't find Gflags header. Please specify: -DGFLAGS_INCLUDE=...")
ENDIF (NOT EXISTS ${GFLAGS_INCLUDE})
MESSAGE("-- Found Gflags header in: ${GFLAGS_INCLUDE}")


# ==========================================================================
#                                                                     FFTW 3
# ==========================================================================

MESSAGE("-- Check for FFTW 3")

SET (FFTW_SEARCH_LIBS ${SEARCH_LIBS})
LIST (APPEND FFTW_SEARCH_LIBS /usr/lib/fftw3)
LIST (APPEND FFTW_SEARCH_LIBS /usr/local/lib/fftw3)
LIST (APPEND FFTW_SEARCH_LIBS /usr/local/homebrew/lib/fftw3)
FIND_LIBRARY(FFTW_LIB NAMES fftw3 PATHS ${FFTW_SEARCH_LIBS})
IF (NOT EXISTS ${FFTW_LIB})
    MESSAGE(FATAL_ERROR "Can't find FFTW 3 library. Try passing -DFFTW_LIB=...")
ENDIF (NOT EXISTS ${FFTW_LIB})
LIST (APPEND THRESHER_DEPS ${FFTW_LIB})
MESSAGE("-- Found FFTW 3 library: ${FFTW_LIB}")

SET (FFTW_SEARCH_HEADERS ${SEARCH_HEADERS})
LIST (APPEND FFTW_SEARCH_HEADERS /usr/include/fftw3)
LIST (APPEND FFTW_SEARCH_HEADERS /usr/local/include/fftw3)
LIST (APPEND FFTW_SEARCH_HEADERS /usr/local/homebrew/include/fftw3)
FIND_PATH(FFTW_INCLUDE NAMES fftw3.h PATHS ${FFTW_SEARCH_HEADERS})
IF (NOT EXISTS ${FFTW_INCLUDE})
    MESSAGE(FATAL_ERROR "Can't find FFTW 3 header. Please specify: -DFFTW_INCLUDE=...")
ENDIF (NOT EXISTS ${FFTW_INCLUDE})
MESSAGE("-- Found FFTW 3 header in: ${FFTW_INCLUDE}")


# ==========================================================================
#                                                            THE EXECUTABLES
# ==========================================================================

INCLUDE_DIRECTORIES (include
                     ${CERES_INCLUDE}
                     ${EIGEN_INCLUDE}
                     ${FITSIO_INCLUDE}
                     ${GFLAGS_INCLUDE}
                     ${FFTW_INCLUDE})


# ./thresh
SET (THRESHER_SRC
     src/fits.cc
     src/thresher.cc
     src/thresh.cc)

ADD_EXECUTABLE (thresh ${THRESHER_SRC})
TARGET_LINK_LIBRARIES (thresh ${THRESHER_DEPS})

# ./lucky
SET (LUCKY_SRC
     src/fits.cc
     src/thresher.cc
     src/lucky.cc)

ADD_EXECUTABLE (lucky ${LUCKY_SRC})
TARGET_LINK_LIBRARIES (lucky ${THRESHER_DEPS})

# Installation
INSTALL (TARGETS thresh lucky DESTINATION bin)
