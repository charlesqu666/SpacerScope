cmake_minimum_required(VERSION 3.10)
project(spacerscope VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(APPLE)
    if(NOT CMAKE_OSX_DEPLOYMENT_TARGET)
        set(CMAKE_OSX_DEPLOYMENT_TARGET "11.0" CACHE STRING "Minimum OS X deployment version" FORCE)
    endif()
endif()

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

find_package(OpenMP REQUIRED)
find_package(ZLIB REQUIRED)

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp")
    set(SPACERSCOPE_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")
else()
    set(SPACERSCOPE_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
endif()

if(EXISTS "${SPACERSCOPE_SOURCE_DIR}/core/fasta_parser.cpp")
    set(SPACERSCOPE_CORE_DIR "${SPACERSCOPE_SOURCE_DIR}/core")
else()
    set(SPACERSCOPE_CORE_DIR "${SPACERSCOPE_SOURCE_DIR}")
endif()

if(EXISTS "${SPACERSCOPE_SOURCE_DIR}/cli/main.cpp")
    set(SPACERSCOPE_CLI_DIR "${SPACERSCOPE_SOURCE_DIR}/cli")
else()
    set(SPACERSCOPE_CLI_DIR "${SPACERSCOPE_SOURCE_DIR}")
endif()

if(EXISTS "${SPACERSCOPE_SOURCE_DIR}/web/web_main.cpp")
    set(SPACERSCOPE_WEB_DIR "${SPACERSCOPE_SOURCE_DIR}/web")
else()
    set(SPACERSCOPE_WEB_DIR "${SPACERSCOPE_SOURCE_DIR}")
endif()

set(SPACERSCOPE_RESOURCE_DIR "${SPACERSCOPE_SOURCE_DIR}/resources")

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
target_link_libraries(spacerscope PRIVATE OpenMP::OpenMP_CXX ZLIB::ZLIB)
target_compile_definitions(spacerscope PRIVATE $<$<CONFIG:Release>:NDEBUG>)
target_compile_definitions(spacerscope PRIVATE SPACERSCOPE_VERSION="${PROJECT_VERSION}")

if(MSVC)
    target_compile_options(spacerscope PRIVATE
        /utf-8
        /Oi
        /Ot
        /Gy
        $<$<CONFIG:Release>:/O2>
    )
    set_property(TARGET spacerscope APPEND_STRING PROPERTY LINK_FLAGS_RELEASE " /OPT:REF /OPT:ICF")
else()
    target_compile_options(spacerscope PRIVATE
        $<$<CONFIG:Release>:-O3>
        $<$<CONFIG:Release>:-ffunction-sections>
        $<$<CONFIG:Release>:-fdata-sections>
        $<$<CONFIG:Debug>:-g>
    )
    if(MINGW)
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static -static-libgcc -static-libstdc++")
    endif()

    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        set_property(TARGET spacerscope APPEND_STRING PROPERTY LINK_FLAGS_RELEASE " -Wl,--gc-sections")
    elseif(APPLE)
        set_property(TARGET spacerscope APPEND_STRING PROPERTY LINK_FLAGS_RELEASE " -Wl,-dead_strip")
    endif()
endif()

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
    message(STATUS "Detecting best available SIMD instruction set...")

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
            enable_simd_feature(${target} "NEON" "" "" "__ARM_NEON")
        else()
            enable_simd_feature(${target} "NEON" "-march=armv8-a+simd" "" "__ARM_NEON")
        endif()

        message(STATUS "  - No extra ARM SIMD flag enabled.")
    else()
        message(STATUS "  - SIMD auto-detect skipped for processor ${CMAKE_SYSTEM_PROCESSOR}.")
    endif()
endfunction()

target_enable_best_simd(spacerscope)

add_executable(spacerscope-web
        ${SPACERSCOPE_WEB_DIR}/web_main.cpp
        ${SPACERSCOPE_WEB_DIR}/web_job.cpp
)
if(WIN32)
    target_sources(spacerscope-web PRIVATE ${SPACERSCOPE_RESOURCE_DIR}/spacerscope-web.rc)
endif()

target_include_directories(spacerscope-web PUBLIC ${SPACERSCOPE_COMMON_INCLUDE_DIRS})
target_compile_definitions(spacerscope-web PRIVATE $<$<CONFIG:Release>:NDEBUG>)
target_compile_definitions(spacerscope-web PRIVATE SPACERSCOPE_VERSION="${PROJECT_VERSION}")
if(WIN32)
    target_link_libraries(spacerscope-web PRIVATE ws2_32)
endif()

if(MSVC)
    target_compile_options(spacerscope-web PRIVATE
        /utf-8
        /Oi
        /Ot
        /Gy
        $<$<CONFIG:Release>:/O2>
    )
    set_property(TARGET spacerscope-web APPEND_STRING PROPERTY LINK_FLAGS_RELEASE " /OPT:REF /OPT:ICF")
else()
    target_compile_options(spacerscope-web PRIVATE
        $<$<CONFIG:Release>:-O3>
        $<$<CONFIG:Release>:-ffunction-sections>
        $<$<CONFIG:Release>:-fdata-sections>
        $<$<CONFIG:Debug>:-g>
    )
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        set_property(TARGET spacerscope-web APPEND_STRING PROPERTY LINK_FLAGS_RELEASE " -Wl,--gc-sections")
    elseif(APPLE)
        set_property(TARGET spacerscope-web APPEND_STRING PROPERTY LINK_FLAGS_RELEASE " -Wl,-dead_strip")
    endif()
endif()

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
