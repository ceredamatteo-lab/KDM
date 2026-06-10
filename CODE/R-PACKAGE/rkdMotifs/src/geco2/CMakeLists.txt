cmake_minimum_required(VERSION 3.18)
project(geco2 VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# ── Shared library only ───────────────────────────────────────────────────────
set(BUILD_SHARED_LIBS ON)

# ── Options ───────────────────────────────────────────────────────────────────
option(GECO_ENABLE_CUDA "Use CUDA GPU acceleration if available" ON)

# ── Hardened-linker detection (Fedora / RHEL / CentOS) ───────────────────────
# These distros inject -specs=redhat-hardened-ld into every link command via
# the rpm build system, forcing PIE mode on all executables.  CMake's CUDA
# compiler test builds a plain executable; the CUDA object files are not
# compiled with -fPIE by default, so the linker rejects them with:
#   "relocation R_X86_64_32 can not be used when making a PIE object"
#
# We detect the spec file and apply two fixes:
#   1. CMAKE_CUDA_FLAGS  -Xcompiler -fPIE  → fixes the CMake compiler-test exe
#   2. Per-target        -Xcompiler -fPIC  → fixes the actual shared-lib objects
#
# Both flags must be set BEFORE enable_language(CUDA) so they are visible
# during the compiler-test link step.
if(EXISTS "/usr/lib/rpm/redhat/redhat-hardened-ld")
    message(STATUS "geco2: hardened linker detected → enabling PIE/PIC flags for CUDA")
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler -fPIE")
endif()

# ── CUDA detection ────────────────────────────────────────────────────────────
set(GECO_HAS_CUDA OFF)

if(GECO_ENABLE_CUDA)

    # ── Locate nvcc ───────────────────────────────────────────────────────────
    # Honour explicit env vars before probing PATH
    if(DEFINED ENV{CUDA_HOME} AND NOT CMAKE_CUDA_COMPILER)
        set(CMAKE_CUDA_COMPILER "$ENV{CUDA_HOME}/bin/nvcc"
            CACHE PATH "CUDA compiler (from CUDA_HOME)")
        message(STATUS "geco2: nvcc from CUDA_HOME → ${CMAKE_CUDA_COMPILER}")
    endif()

    if(DEFINED ENV{CUDA_PATH} AND NOT CMAKE_CUDA_COMPILER)
        set(CMAKE_CUDA_COMPILER "$ENV{CUDA_PATH}/bin/nvcc"
            CACHE PATH "CUDA compiler (from CUDA_PATH)")
        message(STATUS "geco2: nvcc from CUDA_PATH → ${CMAKE_CUDA_COMPILER}")
    endif()

    # Last resort: search the most common install locations on Linux
    if(NOT CMAKE_CUDA_COMPILER)
        find_program(CMAKE_CUDA_COMPILER
            NAMES nvcc
            PATHS
                /usr/local/cuda/bin
                /usr/local/cuda-12/bin
                /usr/local/cuda-11/bin
                /opt/cuda/bin
                /usr/bin
            NO_DEFAULT_PATH
        )
        if(CMAKE_CUDA_COMPILER)
            message(STATUS "geco2: nvcc found by path search → ${CMAKE_CUDA_COMPILER}")
        endif()
    endif()

    # check_language probes without hard-failing
    include(CheckLanguage)
    check_language(CUDA)

    if(CMAKE_CUDA_COMPILER)
        enable_language(CUDA)
        find_package(CUDAToolkit REQUIRED)

        set(GECO_HAS_CUDA ON)
        set(CMAKE_CUDA_STANDARD 17)
        set(CMAKE_CUDA_STANDARD_REQUIRED ON)

        message(STATUS "geco2: CUDA found")
        message(STATUS "       compiler  : ${CMAKE_CUDA_COMPILER}")
        message(STATUS "       toolkit   : ${CUDAToolkit_VERSION}")
    else()
        message(STATUS "geco2: CUDA not found — building CPU-only")
        message(STATUS "       Searched: PATH, CUDA_HOME, CUDA_PATH,")
        message(STATUS "                 /usr/local/cuda*/bin, /opt/cuda/bin")
        message(STATUS "       Fixes:")
        message(STATUS "         export PATH=/usr/local/cuda/bin:\$PATH && R CMD INSTALL .")
        message(STATUS "         CMAKE_CUDA_COMPILER=/usr/local/cuda/bin/nvcc R CMD INSTALL .")
        message(STATUS "         cmake .. -DCMAKE_CUDA_COMPILER=/usr/local/cuda/bin/nvcc")
        message(STATUS "         cmake .. -DGECO_ENABLE_CUDA=OFF  (to build CPU-only intentionally)")
    endif()

else()
    message(STATUS "geco2: CUDA disabled (GECO_ENABLE_CUDA=OFF)")
endif()

# ── OpenMP (REQUIRED) ─────────────────────────────────────────────────────────
# macOS: brew install libomp && cmake .. -DOpenMP_ROOT=$(brew --prefix libomp)
find_package(OpenMP REQUIRED)
message(STATUS "geco2: OpenMP ${OpenMP_CXX_VERSION}")

# ── OpenBLAS (REQUIRED) ───────────────────────────────────────────────────────
# cmake .. -DOpenBLAS_ROOT=/opt/OpenBLAS   or   -DOpenBLAS_ROOT=$CONDA_PREFIX
find_package(PkgConfig REQUIRED)
pkg_check_modules(OpenBLAS QUIET openblas)

if(NOT OpenBLAS_FOUND)
    find_path(OpenBLAS_INCLUDE_DIR
        NAMES cblas.h openblas/cblas.h
        HINTS ${OpenBLAS_ROOT} $ENV{OpenBLAS_ROOT}
              /usr/include/openblas /usr/local/include/openblas)
    find_library(OpenBLAS_LIB
        NAMES openblas
        HINTS ${OpenBLAS_ROOT} $ENV{OpenBLAS_ROOT}
              /usr/lib /usr/local/lib)
    if(OpenBLAS_INCLUDE_DIR AND OpenBLAS_LIB)
        add_library(OpenBLAS::OpenBLAS IMPORTED INTERFACE)
        set_target_properties(OpenBLAS::OpenBLAS PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${OpenBLAS_INCLUDE_DIR}"
            INTERFACE_LINK_LIBRARIES      "${OpenBLAS_LIB}")
        set(OpenBLAS_FOUND TRUE)
        message(STATUS "geco2: OpenBLAS found via find_library")
    endif()
endif()

if(NOT OpenBLAS_FOUND)
    message(FATAL_ERROR
        "\n  REQUIRED: OpenBLAS not found.\n"
        "  Ubuntu/Debian : sudo apt install libopenblas-dev\n"
        "  Fedora/RHEL   : sudo dnf install openblas-devel\n"
        "  macOS         : brew install openblas\n"
        "  Conda         : conda install -c conda-forge openblas\n"
        "  Custom prefix : cmake .. -DOpenBLAS_ROOT=/path/to/openblas\n")
endif()
message(STATUS "geco2: OpenBLAS ${OpenBLAS_VERSION}")

# ── NLopt (REQUIRED) ──────────────────────────────────────────────────────────
# cmake .. -DNLopt_ROOT=/opt/nlopt
find_package(NLopt QUIET CONFIG)

if(NOT NLopt_FOUND)
    pkg_check_modules(NLopt QUIET nlopt)
    if(NLopt_FOUND AND NOT TARGET NLopt::nlopt)
        add_library(NLopt::nlopt IMPORTED INTERFACE)
        set_target_properties(NLopt::nlopt PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${NLopt_INCLUDE_DIRS}"
            INTERFACE_LINK_LIBRARIES      "${NLopt_LINK_LIBRARIES}")
    endif()
endif()

if(NOT NLopt_FOUND)
    message(FATAL_ERROR
        "\n  REQUIRED: NLopt not found.\n"
        "  Ubuntu/Debian : sudo apt install libnlopt-dev\n"
        "  Fedora/RHEL   : sudo dnf install nlopt-devel\n"
        "  macOS         : brew install nlopt\n"
        "  Conda         : conda install -c conda-forge nlopt\n"
        "  Custom prefix : cmake .. -DNLopt_ROOT=/path/to/nlopt\n")
endif()
message(STATUS "geco2: NLopt ${NLopt_VERSION}")

# ── Boost (REQUIRED) — iostreams, filesystem ──────────────────────────────────
# CMake 3.30 removed the legacy FindBoost module (policy CMP0167).
# We explicitly request config-file mode (BoostConfig.cmake shipped with
# Boost >= 1.70) and fall back to the legacy module only on older CMake.
#
# cmake .. -DBOOST_ROOT=/opt/boost
# cmake .. -DBOOST_ROOT=$CONDA_PREFIX
cmake_policy(SET CMP0167 NEW)          # silence the deprecation warning
                                       # and use BoostConfig.cmake directly

set(Boost_USE_STATIC_LIBS   OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)

find_package(Boost REQUIRED
    CONFIG                             # config-file mode: uses BoostConfig.cmake
    COMPONENTS iostreams filesystem system
)

if(NOT Boost_FOUND)
    message(FATAL_ERROR
        "\n  REQUIRED: Boost (iostreams, filesystem) not found.\n"
        "  Ubuntu/Debian : sudo apt install libboost-iostreams-dev libboost-filesystem-dev\n"
        "  Fedora/RHEL   : sudo dnf install boost-iostreams boost-filesystem\n"
        "  macOS         : brew install boost\n"
        "  Conda         : conda install -c conda-forge boost\n"
        "  Custom prefix : cmake .. -DBOOST_ROOT=/path/to/boost\n")
endif()
message(STATUS "geco2: Boost ${Boost_VERSION}")

# ── Source list ───────────────────────────────────────────────────────────────
set(GECO2_SOURCES
    src/base.cpp
    src/array.cpp
    src/sequence.cpp
    src/element.cpp
    src/transcript.cpp
    src/genotype.cpp
    src/ioutils.cpp
    src/2bit_retriever.cpp
    src/gMethods.cpp
    src/fastcluster/fastcluster_cd.cpp
    src/gCluster.cpp
    src/gKmers.cpp
    src/gNMF.cpp
    src/TFM-Pvalue/Matrix.cpp
    src/gPWM.cpp
    src/csv_helpers.cpp
    src/gKDMotifs.cpp
)

if(GECO_HAS_CUDA)
    list(APPEND GECO2_SOURCES
        src/gKDMotifs.cu
        src/gKmers_cuda.cu
        src/gMethods_cuda.cu
        src/gNMF_cuda.cu
    )
endif()

# ── Library target ────────────────────────────────────────────────────────────
add_library(geco2 SHARED ${GECO2_SOURCES})
add_library(geco2::geco2 ALIAS geco2)

target_include_directories(geco2
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

set_target_properties(geco2 PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    # Do NOT set VERSION / SOVERSION here.
    # When VERSION is set CMake embeds a soname (libgeco2.so.1) in the shared
    # library.  R's linker then records that soname and the dynamic loader looks
    # for libgeco2.so.1 at runtime — but inst/lib/ only contains libgeco2.so,
    # causing "cannot open shared object file: No such file or directory".
    # Without VERSION the soname stays as libgeco2.so, which matches the file
    # we copy into inst/lib/ and the rpath finds it correctly.
)

# ── CUDA target settings ──────────────────────────────────────────────────────
if(GECO_HAS_CUDA)
    target_compile_definitions(geco2 PUBLIC GECO_HAS_CUDA)

    set_target_properties(geco2 PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
        CUDA_ARCHITECTURES         "${CMAKE_CUDA_ARCHITECTURES}"
    )

    if(NOT CMAKE_CUDA_ARCHITECTURES)
        set_target_properties(geco2 PROPERTIES CUDA_ARCHITECTURES "75;80;86")
        message(STATUS "geco2: CUDA_ARCHITECTURES not set — using default 75;80;86")
    endif()

    # Per-target PIC fix for hardened systems (companion to the global fPIE above)
    # -fPIC is correct for shared libraries; -fPIE is only for executables.
    if(EXISTS "/usr/lib/rpm/redhat/redhat-hardened-ld")
        target_compile_options(geco2 PRIVATE
            $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler -fPIC>
        )
    endif()
endif()

# ── Link dependencies (all PRIVATE — none exposed in geco2.h) ─────────────────
if(GECO_HAS_CUDA)
    # cudart: static preferred — baked into libgeco2.so, no extra runtime dep
    if(TARGET CUDA::cudart_static)
        target_link_libraries(geco2 PRIVATE CUDA::cudart_static)
    else()
        target_link_libraries(geco2 PRIVATE CUDA::cudart)
    endif()
    # cuBLAS / cuSOLVER / cuSPARSE: shared only (no static variant from NVIDIA)
    target_link_libraries(geco2 PRIVATE
        CUDA::cublas
        CUDA::cusolver
        CUDA::cusparse
    )
endif()

target_link_libraries(geco2 PRIVATE OpenMP::OpenMP_CXX)

if(TARGET OpenBLAS::OpenBLAS)
    target_link_libraries(geco2 PRIVATE OpenBLAS::OpenBLAS)
else()
    target_include_directories(geco2 PRIVATE ${OpenBLAS_INCLUDE_DIRS})
    target_link_libraries(geco2 PRIVATE ${OpenBLAS_LIBRARIES})
endif()

target_link_libraries(geco2 PRIVATE NLopt::nlopt)

target_link_libraries(geco2 PRIVATE
    Boost::iostreams
    Boost::filesystem
)

# ── Installation ──────────────────────────────────────────────────────────────
install(TARGETS geco2
    EXPORT  geco2Targets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/geco2
    FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp"
)

install(EXPORT geco2Targets
    FILE        geco2Targets.cmake
    NAMESPACE   geco2::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/geco2
)

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/geco2Config.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/geco2Config.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/geco2
)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/geco2ConfigVersion.cmake"
    VERSION       ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/geco2Config.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/geco2ConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/geco2
)

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/geco2.pc.in"
    "${CMAKE_CURRENT_BINARY_DIR}/geco2.pc"
    @ONLY
)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/geco2.pc"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig
)

# ── Summary ───────────────────────────────────────────────────────────────────
message(STATUS "")
message(STATUS "┌─ geco2 ${PROJECT_VERSION} configuration ───────────────────────────")
message(STATUS "│  Install prefix  : ${CMAKE_INSTALL_PREFIX}")
message(STATUS "│  Build type      : ${CMAKE_BUILD_TYPE}")
message(STATUS "│  Library type    : SHARED")
message(STATUS "│  CUDA (optional) : ${GECO_HAS_CUDA}")
if(GECO_HAS_CUDA)
message(STATUS "│    toolkit       : ${CUDAToolkit_VERSION}")
message(STATUS "│    architectures : ${CMAKE_CUDA_ARCHITECTURES}")
message(STATUS "│    libraries     : cudart cublas cusolver cusparse")
endif()
message(STATUS "│  OpenMP          : ${OpenMP_CXX_VERSION}  [required]")
message(STATUS "│  OpenBLAS        : ${OpenBLAS_VERSION}  [required]")
message(STATUS "│  NLopt           : ${NLopt_VERSION}  [required]")
message(STATUS "│  Boost           : ${Boost_VERSION}  [required — iostreams, filesystem]")
message(STATUS "└────────────────────────────────────────────────────────────────")
