#[[
/*
 * SpacerScope - gRNA Design and Off-Target Effect Analysis Tool
 * Copyright (C) 2026 charlesqu666
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful for CRISPR
 * research and therapeutic applications, but WITHOUT ANY WARRANTY; without
 * even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE. See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 *
 * Source code is available at: https://github.com/charlesqu666/SpacerScope
 */
]]
cmake_minimum_required(VERSION 3.10)
project(spacerscope VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# macOS deployment target
if(APPLE)
    if(NOT CMAKE_OSX_DEPLOYMENT_TARGET)
        set(CMAKE_OSX_DEPLOYMENT_TARGET "11.0" CACHE STRING "Minimum OS X deployment version" FORCE)
    endif()
endif()

# Default to Release build type
if(NOT CMAKE_CONFIGURATION_TYPES)
    set(_valid_build_types Debug Release RelWithDebInfo MinSizeRel)
    if(NOT CMAKE_BUILD_TYPE)
        set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    else()
        string(TOLOWER "${CMAKE_BUILD_TYPE}" _build_type_lower)
        if(_build_type_lower STREQUAL "release")
            set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
        elseif(_build_type_lower STREQUAL "debug")
            set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
        elseif(_build_type_lower STREQUAL "relwithdebinfo")
            set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Build type" FORCE)
        elseif(_build_type_lower STREQUAL "minsizerel")
            set(CMAKE_BUILD_TYPE MinSizeRel CACHE STRING "Build type" FORCE)
        endif()
    endif()
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS ${_valid_build_types})
endif()

option(SPACERSCOPE_ENABLE_IPO "Enable interprocedural optimization for release builds" ON)
option(SPACERSCOPE_ENABLE_NATIVE "Enable host-specific CPU tuning for release builds" OFF)

include(CheckCXXCompilerFlag)
include(CheckIPOSupported)

# MSVC: Force static CRT (/MT) instead of dynamic (/MD)
if(WIN32 AND NOT MINGW)
    foreach(_flag_var
            CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE
            CMAKE_CXX_FLAGS_RELWITHDEBINFO CMAKE_CXX_FLAGS_MINSIZEREL
            CMAKE_C_FLAGS   CMAKE_C_FLAGS_DEBUG   CMAKE_C_FLAGS_RELEASE
            CMAKE_C_FLAGS_RELWITHDEBINFO   CMAKE_C_FLAGS_MINSIZEREL)
        if(${_flag_var} MATCHES "/MD")
            string(REGEX REPLACE "/MD" "/MT" ${_flag_var} "${${_flag_var}}")
        endif()
    endforeach()
    message(STATUS "MSVC: switched to /MT (static CRT).")
endif()

# Prefer static archives for dependencies on MinGW/Linux
if(MINGW OR (UNIX AND NOT APPLE))
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a;.so" CACHE STRING "" FORCE)
endif()

find_package(OpenMP REQUIRED)
find_package(ZLIB REQUIRED)
find_package(Threads REQUIRED)

# Find static versions of OpenMP and ZLIB for static linking
if(MINGW OR (UNIX AND NOT APPLE))
    # Try to find static libgomp
    find_library(LIBGOMP_STATIC NAMES libgomp.a gomp.a
            PATHS /usr/lib /usr/local/lib /usr/lib/x86_64-linux-gnu /usr/lib/gcc
            PATH_SUFFIXES gcc x86_64-linux-gnu
            NO_DEFAULT_PATH)

    # Try to find static libz
    find_library(LIBZ_STATIC NAMES libz.a z.a
            PATHS /usr/lib /usr/local/lib /usr/lib/x86_64-linux-gnu
            NO_DEFAULT_PATH)

    if(LIBGOMP_STATIC)
        message(STATUS "Found static libgomp: ${LIBGOMP_STATIC}")
    endif()
    if(LIBZ_STATIC)
        message(STATUS "Found static libz: ${LIBZ_STATIC}")
    endif()
endif()

set(SPACERSCOPE_SOURCE_DIR  "${CMAKE_CURRENT_SOURCE_DIR}")
set(SPACERSCOPE_CORE_DIR    "${CMAKE_CURRENT_SOURCE_DIR}/core")
set(SPACERSCOPE_CLI_DIR     "${CMAKE_CURRENT_SOURCE_DIR}/cli")
set(SPACERSCOPE_WEB_DIR     "${CMAKE_CURRENT_SOURCE_DIR}/web")
set(SPACERSCOPE_RESOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/resources")

# -----------------------------------------------------------------------------
# Core Target: spacerscope
# -----------------------------------------------------------------------------
add_executable(spacerscope
        ${SPACERSCOPE_CLI_DIR}/main.cpp
        ${SPACERSCOPE_CORE_DIR}/fasta_parser.cpp
        ${SPACERSCOPE_CORE_DIR}/mismatch_search.cpp
        ${SPACERSCOPE_CORE_DIR}/indel_search.cpp
        ${SPACERSCOPE_CORE_DIR}/spacer_extraction.cpp
        ${SPACERSCOPE_CORE_DIR}/dust_filter.cpp
        ${SPACERSCOPE_CORE_DIR}/count_filter.cpp
        ${SPACERSCOPE_CORE_DIR}/anno_parser.cpp
        ${SPACERSCOPE_CORE_DIR}/on_target_filter.cpp
        ${SPACERSCOPE_CORE_DIR}/utils.cpp
        ${SPACERSCOPE_CLI_DIR}/postprocess.cpp
        ${SPACERSCOPE_CLI_DIR}/report_html.cpp
        ${SPACERSCOPE_CLI_DIR}/progress.cpp
)

if(WIN32)
    target_sources(spacerscope PRIVATE ${SPACERSCOPE_RESOURCE_DIR}/spacerscope.rc)
endif()

set(SPACERSCOPE_COMMON_INCLUDE_DIRS
        ${SPACERSCOPE_SOURCE_DIR}
        ${SPACERSCOPE_CORE_DIR}
        ${SPACERSCOPE_CLI_DIR}
        ${SPACERSCOPE_WEB_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/third_party
)

target_include_directories(spacerscope PUBLIC ${SPACERSCOPE_COMMON_INCLUDE_DIRS})

# Link libraries - use static versions if found
if(LIBGOMP_STATIC AND (MINGW OR CMAKE_SYSTEM_NAME STREQUAL "Linux"))
    target_link_libraries(spacerscope PRIVATE "${LIBGOMP_STATIC}")
else()
    target_link_libraries(spacerscope PRIVATE OpenMP::OpenMP_CXX)
endif()

if(LIBZ_STATIC AND (MINGW OR CMAKE_SYSTEM_NAME STREQUAL "Linux"))
    target_link_libraries(spacerscope PRIVATE "${LIBZ_STATIC}")
else()
    target_link_libraries(spacerscope PRIVATE ZLIB::ZLIB)
endif()

target_link_libraries(spacerscope PRIVATE Threads::Threads)
target_compile_definitions(spacerscope PRIVATE $<$<CONFIG:Release>:NDEBUG>)
target_compile_definitions(spacerscope PRIVATE SPACERSCOPE_VERSION="${PROJECT_VERSION}")

# Core Compiler Flags
if(MSVC)
    target_compile_options(spacerscope PRIVATE
            /utf-8 /Oi /Ot /Gy $<$<CONFIG:Release>:/O2>
    )

    # MSVC OpenMP Configuration
    check_cxx_compiler_flag("/openmp:llvm" MSVC_OPENMP_LLVM)
    if(MSVC_OPENMP_LLVM)
        target_compile_options(spacerscope PRIVATE /openmp:llvm)
        find_library(LIBOMP_STATIC NAMES libomp.lib omp.lib
                PATHS "$ENV{LLVM_DIR}/lib" "$ENV{VCToolsInstallDir}/lib"
                "C:/Program Files/LLVM/lib"
                "C:/Program Files (x86)/Microsoft Visual Studio/2022/BuildTools/VC/Tools/Llvm/lib"
                NO_DEFAULT_PATH)
        if(LIBOMP_STATIC)
            message(STATUS "MSVC: linking static libomp: ${LIBOMP_STATIC}")
            target_link_libraries(spacerscope PRIVATE "${LIBOMP_STATIC}")
        else()
            message(WARNING "MSVC: static libomp not found. Executable will depend on libomp.dll.")
        endif()
    else()
        message(WARNING "MSVC: /openmp:llvm not supported. Using default /openmp (vcomp*.dll).")
    endif()
    set_property(TARGET spacerscope APPEND_STRING PROPERTY LINK_FLAGS_RELEASE " /OPT:REF /OPT:ICF")
else()
    target_compile_options(spacerscope PRIVATE
            $<$<CONFIG:Release>:-O3>
            $<$<CONFIG:Release>:-ffunction-sections>
            $<$<CONFIG:Release>:-fdata-sections>
            $<$<CONFIG:Debug>:-g>
    )

    # Platform-specific static linking to ensure self-contained binaries
    if(MINGW)
        target_link_options(spacerscope PRIVATE
                -static -static-libgcc -static-libstdc++
        )
    endif()

    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        target_link_options(spacerscope PRIVATE
                -static-libgcc -static-libstdc++
        )
        set_property(TARGET spacerscope APPEND_STRING PROPERTY LINK_FLAGS_RELEASE " -Wl,--gc-sections")
    elseif(APPLE)
        # macOS: Attempt static linking for libomp via Homebrew/LLVM
        find_library(LIBOMP_STATIC NAMES libomp.a omp.a
                PATHS "/opt/homebrew/opt/libomp/lib" "/usr/local/opt/libomp/lib"
                "/opt/homebrew/opt/llvm/lib" "/usr/local/opt/llvm/lib"
                NO_DEFAULT_PATH)
        if(LIBOMP_STATIC)
            message(STATUS "macOS: linking static libomp: ${LIBOMP_STATIC}")
            set_target_properties(OpenMP::OpenMP_CXX PROPERTIES
                    IMPORTED_LOCATION "${LIBOMP_STATIC}"
                    INTERFACE_LINK_LIBRARIES "${LIBOMP_STATIC}")
        else()
            message(WARNING "macOS: static libomp.a not found. Binary will require libomp.dylib.")
        endif()
        set_property(TARGET spacerscope APPEND_STRING PROPERTY LINK_FLAGS_RELEASE " -Wl,-dead_strip")
    endif()
endif()

# -----------------------------------------------------------------------------
# SIMD Optimization
# -----------------------------------------------------------------------------
function(enable_simd_feature target feature_name gcc_flag msvc_flag)
    set(_flag "")
    if(MSVC)
        set(_flag "${msvc_flag}")
    else()
        set(_flag "${gcc_flag}")
    endif()

    if(_flag STREQUAL "")
        set(SPACERSCOPE_SIMD_ENABLED FALSE PARENT_SCOPE)
        return()
    endif()

    string(REGEX REPLACE "[^A-Za-z0-9]" "_" _feature_id "${feature_name}")
    string(TOUPPER "${_feature_id}" _feature_id)
    set(_check_var "COMPILER_SUPPORTS_${_feature_id}")
    check_cxx_compiler_flag("${_flag}" ${_check_var})

    if(${_check_var})
        message(STATUS "  + ${feature_name} enabled (${CMAKE_CXX_COMPILER_ID}: ${_flag})")
        target_compile_options(${target} PRIVATE "${_flag}")
        foreach(_def IN LISTS ARGN)
            target_compile_definitions(${target} PRIVATE "${_def}=1")
        endforeach()
        set(SPACERSCOPE_SIMD_ENABLED TRUE PARENT_SCOPE)
    else()
        set(SPACERSCOPE_SIMD_ENABLED FALSE PARENT_SCOPE)
    endif()
endfunction()

function(target_enable_best_simd target)
    message(STATUS "Detecting best available SIMD instruction set for ${target}...")

    if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64|i.86")
        enable_simd_feature(${target} "AVX2" "-mavx2" "/arch:AVX2" "__AVX2__")
        if(SPACERSCOPE_SIMD_ENABLED)
            return()
        endif()

        enable_simd_feature(${target} "SSE4.1" "-msse4.1" "/arch:AVX" "__SSE4_1__" "__SSE__")
        if(SPACERSCOPE_SIMD_ENABLED)
            return()
        endif()

        message(STATUS "  - No extra x86 SIMD flag enabled.")
    elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "arm|aarch64")
        if(APPLE)
            message(STATUS "  + NEON always-on (Apple Silicon, enabled by default)")
        else()
            enable_simd_feature(${target} "NEON" "-march=armv8-a+simd" "" "__ARM_NEON")
            if(NOT SPACERSCOPE_SIMD_ENABLED)
                message(STATUS "  - No extra ARM SIMD flag enabled.")
            endif()
        endif()
    else()
        message(STATUS "  - SIMD auto-detect skipped for processor ${CMAKE_SYSTEM_PROCESSOR}.")
    endif()
endfunction()

# Enable best SIMD support for core target
target_enable_best_simd(spacerscope)

# -----------------------------------------------------------------------------
# Web Target: spacerscope-web
# -----------------------------------------------------------------------------
add_executable(spacerscope-web
        ${SPACERSCOPE_WEB_DIR}/web_main.cpp
        ${SPACERSCOPE_WEB_DIR}/web_job.cpp
)
if(WIN32)
    target_sources(spacerscope-web PRIVATE ${SPACERSCOPE_RESOURCE_DIR}/spacerscope-web.rc)
endif()

# Enable best SIMD support for web target AFTER it is defined
target_enable_best_simd(spacerscope-web)

target_include_directories(spacerscope-web PUBLIC ${SPACERSCOPE_COMMON_INCLUDE_DIRS})
target_compile_definitions(spacerscope-web PRIVATE $<$<CONFIG:Release>:NDEBUG>)
target_compile_definitions(spacerscope-web PRIVATE SPACERSCOPE_VERSION="${PROJECT_VERSION}")

if(WIN32)
    target_link_libraries(spacerscope-web PRIVATE ws2_32 Threads::Threads)
else()
    target_link_libraries(spacerscope-web PRIVATE Threads::Threads)
endif()

# Web Compiler Flags
if(MSVC)
    target_compile_options(spacerscope-web PRIVATE
            /utf-8 /Oi /Ot /Gy $<$<CONFIG:Release>:/O2>
    )
    set_property(TARGET spacerscope-web APPEND_STRING PROPERTY LINK_FLAGS_RELEASE " /OPT:REF /OPT:ICF")
else()
    target_compile_options(spacerscope-web PRIVATE
            $<$<CONFIG:Release>:-O3>
            $<$<CONFIG:Release>:-ffunction-sections>
            $<$<CONFIG:Release>:-fdata-sections>
            $<$<CONFIG:Debug>:-g>
    )

    # Platform-specific static linking for the web target to prevent runtime DLL errors
    if(MINGW)
        target_link_options(spacerscope-web PRIVATE
                -static -static-libgcc -static-libstdc++
        )
    endif()

    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        target_link_options(spacerscope-web PRIVATE
                -static-libgcc -static-libstdc++
        )
        set_property(TARGET spacerscope-web APPEND_STRING PROPERTY LINK_FLAGS_RELEASE " -Wl,--gc-sections")
    elseif(APPLE)
        set_property(TARGET spacerscope-web APPEND_STRING PROPERTY LINK_FLAGS_RELEASE " -Wl,-dead_strip")
    endif()
endif()

# -----------------------------------------------------------------------------
# Global Build Optimizations (Native Tuning & IPO)
# -----------------------------------------------------------------------------
if(SPACERSCOPE_ENABLE_NATIVE)
    message(STATUS "Host-specific tuning enabled for release builds.")
    if(NOT MSVC)
        check_cxx_compiler_flag("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
        if(COMPILER_SUPPORTS_MARCH_NATIVE)
            target_compile_options(spacerscope PRIVATE $<$<CONFIG:Release>:-march=native>)
            target_compile_options(spacerscope-web PRIVATE $<$<CONFIG:Release>:-march=native>)
        endif()

        check_cxx_compiler_flag("-mtune=native" COMPILER_SUPPORTS_MTUNE_NATIVE)
        if(COMPILER_SUPPORTS_MTUNE_NATIVE)
            target_compile_options(spacerscope PRIVATE $<$<CONFIG:Release>:-mtune=native>)
            target_compile_options(spacerscope-web PRIVATE $<$<CONFIG:Release>:-mtune=native>)
        endif()
    else()
        message(STATUS "  - No portable MSVC equivalent for -march=native; skipping.")
    endif()
endif()

if(SPACERSCOPE_ENABLE_IPO)
    check_ipo_supported(RESULT IPO_SUPPORTED OUTPUT IPO_ERROR)
    if(IPO_SUPPORTED)
        set_property(TARGET spacerscope PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
        set_property(TARGET spacerscope PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELWITHDEBINFO TRUE)
        set_property(TARGET spacerscope-web PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
        set_property(TARGET spacerscope-web PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELWITHDEBINFO TRUE)
        message(STATUS "IPO/LTO enabled for Release and RelWithDebInfo.")
    else()
        message(STATUS "IPO/LTO not enabled: ${IPO_ERROR}")
    endif()
endif()

install(TARGETS spacerscope spacerscope-web RUNTIME DESTINATION bin)