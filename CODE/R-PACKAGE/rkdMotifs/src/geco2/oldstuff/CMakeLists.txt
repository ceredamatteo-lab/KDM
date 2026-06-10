cmake_minimum_required(VERSION 3.18)   # 3.18+ for native CUDA target support

# ── Project declaration ───────────────────────────────────────────────────────
# CUDA is listed as an *optional* language here; we enable it conditionally
# below after probing for the compiler.
project(geco2 VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)


# ── Build options ─────────────────────────────────────────────────────────────
option(BUILD_SHARED_LIBS
    "Build a shared library (.so/.dylib/.dll) instead of a static one" ON)

option(GECO_ENABLE_CUDA
    "Use CUDA GPU acceleration if a CUDA compiler is available" ON)

# ── CUDA detection ────────────────────────────────────────────────────────────
set(GECO_HAS_CUDA OFF)

if(GECO_ENABLE_CUDA)
    # Must be set BEFORE check_language(CUDA) / enable_language(CUDA)
    # so it applies to the CMake compiler-test executable too
    if(EXISTS "/usr/lib/rpm/redhat/redhat-hardened-ld")
        set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler -fPIE")
    endif()
    # include(CheckLanguage) probes for the compiler without hard-failing
    include(CheckLanguage)
    check_language(CUDA)

    if(CMAKE_CUDA_COMPILER)
        enable_language(CUDA)               # activate CUDA toolchain
        find_package(CUDAToolkit REQUIRED)  # brings in CUDA::cudart etc.

        set(GECO_HAS_CUDA ON)
        message(STATUS "geco2: CUDA found → building GPU backend")
        message(STATUS "       CUDA compiler : ${CMAKE_CUDA_COMPILER}")
        message(STATUS "       CUDA version  : ${CUDAToolkit_VERSION}")

        # Default to C++11 for device code too
        set(CMAKE_CUDA_STANDARD 17)
        set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    else()
        message(STATUS "geco2: CUDA compiler not found → building CPU-only")
    endif()
else()
    message(STATUS "geco2: CUDA disabled by GECO_ENABLE_CUDA=OFF")
endif()

# ── OpenMP (REQUIRED) ─────────────────────────────────────────────────────────
# On macOS, Apple Clang does not bundle OpenMP — install via Homebrew first:
#   brew install libomp
#   cmake .. -DOpenMP_ROOT=$(brew --prefix libomp)
find_package(OpenMP REQUIRED)
message(STATUS "geco2: OpenMP found → ${OpenMP_CXX_VERSION}")

# ── OpenBLAS (REQUIRED) ───────────────────────────────────────────────────────
# CMake has no built-in FindOpenBLAS module, so we use pkg-config.
# If not found, provide a clear fatal message with install instructions.
#
# Override search path:
#   cmake .. -DOpenBLAS_ROOT=/opt/OpenBLAS
#   cmake .. -DOpenBLAS_ROOT=$CONDA_PREFIX
find_package(PkgConfig REQUIRED)
pkg_check_modules(OpenBLAS QUIET openblas)

if(NOT OpenBLAS_FOUND)
    # pkg-config missed it — try a direct header+library search before failing
    find_path(OpenBLAS_INCLUDE_DIR
        NAMES cblas.h openblas/cblas.h
        HINTS ${OpenBLAS_ROOT} $ENV{OpenBLAS_ROOT}
              /usr/include/openblas /usr/local/include/openblas
    )
    find_library(OpenBLAS_LIB
        NAMES openblas
        HINTS ${OpenBLAS_ROOT} $ENV{OpenBLAS_ROOT}
              /usr/lib /usr/local/lib
    )
    if(OpenBLAS_INCLUDE_DIR AND OpenBLAS_LIB)
        # Normalise into an IMPORTED target so the linking block below is uniform
        add_library(OpenBLAS::OpenBLAS IMPORTED INTERFACE)
        set_target_properties(OpenBLAS::OpenBLAS PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${OpenBLAS_INCLUDE_DIR}"
            INTERFACE_LINK_LIBRARIES      "${OpenBLAS_LIB}"
        )
        set(OpenBLAS_FOUND TRUE)
        message(STATUS "geco2: OpenBLAS found via find_library → ${OpenBLAS_LIB}")
    endif()
endif()

if(NOT OpenBLAS_FOUND)
    message(FATAL_ERROR
        "\n"
        "  REQUIRED library OpenBLAS was not found.\n"
        "  Install it, then re-run cmake.\n"
        "\n"
        "  Ubuntu/Debian : sudo apt install libopenblas-dev\n"
        "  Fedora/RHEL   : sudo dnf install openblas-devel\n"
        "  macOS         : brew install openblas\n"
        "  Conda         : conda install -c conda-forge openblas\n"
        "\n"
        "  Non-standard prefix: cmake .. -DOpenBLAS_ROOT=/path/to/openblas\n"
    )
endif()
message(STATUS "geco2: OpenBLAS found → ${OpenBLAS_VERSION}")

# ── NLopt (REQUIRED) ──────────────────────────────────────────────────────────
# NLopt ships an official CMake config (NLoptConfig.cmake) since v2.6.
# Older installations fall back to pkg-config.
#
# Override search path:
#   cmake .. -DNLopt_ROOT=/opt/nlopt
#   cmake .. -DNLopt_ROOT=$CONDA_PREFIX
find_package(NLopt QUIET CONFIG)

if(NOT NLopt_FOUND)
    pkg_check_modules(NLopt QUIET nlopt)
    if(NLopt_FOUND AND NOT TARGET NLopt::nlopt)
        add_library(NLopt::nlopt IMPORTED INTERFACE)
        set_target_properties(NLopt::nlopt PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${NLopt_INCLUDE_DIRS}"
            INTERFACE_LINK_LIBRARIES      "${NLopt_LINK_LIBRARIES}"
        )
    endif()
endif()

if(NOT NLopt_FOUND)
    message(FATAL_ERROR
        "\n"
        "  REQUIRED library NLopt was not found.\n"
        "  Install it, then re-run cmake.\n"
        "\n"
        "  Ubuntu/Debian : sudo apt install libnlopt-dev\n"
        "  Fedora/RHEL   : sudo dnf install nlopt-devel\n"
        "  macOS         : brew install nlopt\n"
        "  Conda         : conda install -c conda-forge nlopt\n"
        "\n"
        "  Non-standard prefix: cmake .. -DNLopt_ROOT=/path/to/nlopt\n"
    )
endif()
message(STATUS "geco: NLopt found → ${NLopt_VERSION}")

# ── Boost (REQUIRED) — iostreams, filesystem ──────────────────────────────────
# CMake's FindBoost module locates headers and per-component libraries.
# Each requested component maps to an imported target: Boost::<Component>
#
# Boost::filesystem links Boost::system automatically on older Boost (<1.69);
# on Boost >=1.69 filesystem is self-contained, but listing system is harmless.
#
# Override search path:
#   cmake .. -DBOOST_ROOT=/opt/boost
#   cmake .. -DBOOST_ROOT=$CONDA_PREFIX
#
# Force header-only mode (no pre-built libs, not recommended here):
#   cmake .. -DBoost_USE_STATIC_LIBS=OFF -DBOOST_ROOT=...
set(Boost_USE_STATIC_LIBS    OFF)   # prefer shared Boost libs (matches system installs)
set(Boost_USE_MULTITHREADED  ON)    # thread-safe iostreams/filesystem
set(Boost_USE_STATIC_RUNTIME OFF)

find_package(Boost REQUIRED COMPONENTS
    iostreams
    filesystem
    system          # implicit dependency of filesystem on Boost < 1.69
)

if(NOT Boost_FOUND)
    message(FATAL_ERROR
        "\n"
        "  REQUIRED library Boost (iostreams, filesystem) was not found.\n"
        "  Install it, then re-run cmake.\n"
        "\n"
        "  Ubuntu/Debian : sudo apt install libboost-iostreams-dev libboost-filesystem-dev\n"
        "  Fedora/RHEL   : sudo dnf install boost-iostreams boost-filesystem\n"
        "  macOS         : brew install boost\n"
        "  Conda         : conda install -c conda-forge boost\n"
        "\n"
        "  Non-standard prefix: cmake .. -DBOOST_ROOT=/path/to/boost\n"
    )
endif()
message(STATUS "geco2: Boost found → ${Boost_VERSION}")
message(STATUS "       iostreams : ${Boost_IOSTREAMS_LIBRARY}")
message(STATUS "       filesystem: ${Boost_FILESYSTEM_LIBRARY}")





# ── Source lists ──────────────────────────────────────────────────────────────
set(GECO_SOURCES
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
    list(APPEND GECO_SOURCES src/gKDMotifs.cu src/gKmers_cuda.cu src/gMethods_cuda.cu src/gNMF_cuda.cu )
endif()

# ── Library target ────────────────────────────────────────────────────────────
add_library(geco2 ${GECO_SOURCES})
add_library(geco2::geco2 ALIAS geco2)

# Inject GECO_HAS_CUDA into every translation unit that includes geco.h,
# both inside this build and in installed downstream consumers.
if(GECO_HAS_CUDA)
    target_compile_definitions(geco2
        PUBLIC GECO_HAS_CUDA          # visible to consumers via the header
    )

    # Applied only to .cu compilation units of the mylib target
    target_compile_options(geco2 PRIVATE
        $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler -fPIC>
    )
endif()

target_include_directories(geco2
    PUBLIC
#        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/src>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/geco2>
)

set_target_properties(geco2 PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    VERSION   ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
)

# ── CUDA-specific target settings ─────────────────────────────────────────────
if(GECO_HAS_CUDA)
    set_target_properties(geco2 PROPERTIES
        # Allow .cu files to call __host__ functions across compilation units
        CUDA_SEPARABLE_COMPILATION ON
        # Target the GPU architecture of the current machine when possible;
        # fall back to a safe baseline (Kepler, sm_35) otherwise.
        # Override with: cmake -DCMAKE_CUDA_ARCHITECTURES="70;80" ..
        CUDA_ARCHITECTURES "${CMAKE_CUDA_ARCHITECTURES}"
    )

    if(NOT CMAKE_CUDA_ARCHITECTURES)
        set_target_properties(geco2 PROPERTIES CUDA_ARCHITECTURES "52;60;70;75;80;86")
        message(STATUS "geco2: CUDA_ARCHITECTURES not set — using safe multi-arch default")
    endif()

    # Link against the CUDA runtime (static preferred to avoid rpath issues
    # when the .a is later embedded in the R package .so).
    # PRIVATE: cuda_runtime.h is only included internally , never in
    # the public headers, so consumers don't need the CUDA runtime.
    if(TARGET CUDA::cudart_static)
        target_link_libraries(geco2 PRIVATE CUDA::cudart_static)
        message(STATUS "geco2: linking CUDA::cudart_static (PRIVATE)")
    else()
        target_link_libraries(geco2 PRIVATE CUDA::cudart)
        message(STATUS "geco2: linking CUDA::cudart (PRIVATE, shared runtime)")
    endif()
endif()

# ── Link all dependencies ─────────────────────────────────────────────────────
# All PRIVATE: none of these types appear in the public headers.

# CUDA — optional
if(GECO_HAS_CUDA)
    if(TARGET CUDA::cudart_static)
        target_link_libraries(geco2 PRIVATE CUDA::cudart_static)
        message(STATUS "geco2: linking CUDA::cudart_static (PRIVATE)")
    else()
        target_link_libraries(geco2 PRIVATE CUDA::cudart)
        message(STATUS "geco2: linking CUDA::cudart (PRIVATE, shared runtime)")
    endif()

    target_link_libraries(geco2 PRIVATE
        CUDA::cublas       # dense linear algebra  (libcublas.so)
        CUDA::cusolver     # eigenvalues, SVD, LU  (libcusolver.so)
        CUDA::cusparse     # sparse matrix ops     (libcusparse.so)
    )
    message(STATUS "geco2: linking CUDA::cublas, CUDA::cusolver, CUDA::cusparse")

endif()

# OpenMP, OpenBLAS, NLopt — all required, so unconditional
target_link_libraries(geco2 PRIVATE OpenMP::OpenMP_CXX)

if(TARGET OpenBLAS::OpenBLAS)
    target_link_libraries(geco2 PRIVATE OpenBLAS::OpenBLAS)
else()
    # pkg-config path: variables populated by pkg_check_modules
    target_include_directories(geco2 PRIVATE ${OpenBLAS_INCLUDE_DIRS})
    target_link_libraries(geco2 PRIVATE ${OpenBLAS_LIBRARIES})
endif()

target_link_libraries(geco2 PRIVATE NLopt::nlopt)

target_link_libraries(geco2 PRIVATE
    Boost::iostreams
    Boost::filesystem
)

# ── Installation: library binary ──────────────────────────────────────────────
install(TARGETS geco2
    EXPORT  geco2Targets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

# ── Installation: public headers ──────────────────────────────────────────────
install(DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/geco2
    FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp"
)

# ── CMake package config files ────────────────────────────────────────────────
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

# ── pkg-config file ───────────────────────────────────────────────────────────
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
message(STATUS "┌─ geco2 ${PROJECT_VERSION} configuration ──────────────────────────")
message(STATUS "│  Install prefix  : ${CMAKE_INSTALL_PREFIX}")
message(STATUS "│  Build type      : ${CMAKE_BUILD_TYPE}")
message(STATUS "│  Shared libs     : ${BUILD_SHARED_LIBS}")
message(STATUS "│  CUDA (optional) : ${GECO_HAS_CUDA}")
if(GECO_HAS_CUDA)
message(STATUS "│    toolkit       : ${CUDAToolkit_VERSION}")
message(STATUS "│    architectures : ${CMAKE_CUDA_ARCHITECTURES}")
endif()
message(STATUS "│  OpenMP          : ${OpenMP_CXX_VERSION}  [required]")
message(STATUS "│  OpenBLAS        : ${OpenBLAS_VERSION}  [required]")
message(STATUS "│  NLopt           : ${NLopt_VERSION}  [required]")
message(STATUS "└──────────────────────────────────────────────────────────────")
