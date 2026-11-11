# Emscripten build of MoorhenSquab's C++ core for the moorhen-squab NPM package.
#
# Reuses Moorhen Squab's sources in ../src and ../include unchanged. Differences from
# the top-level CMakeLists.txt (which builds the standalone Vite frontend):
#   * gemmi and nlohmann_json are found or fetched, no sibling-checkout paths
#   * ENVIRONMENT=web,worker so the module can also run inside a Web Worker
#   * exceptions are enabled at compile time too, so the try/catch blocks in
#     atomExtractor.cpp actually catch gemmi parse errors
#   * SINGLE_FILE embeds the .wasm in saxs.js (no asset plumbing for consumers)
#   * cpp/squab_npm_bindings.cpp adds squab_* wrappers that return C++
#     exception text to JS (Emscripten 6 no longer exports getExceptionMessage)
#
# Configure with:  emcmake cmake -S npm -B npm/build-wasm
# Options:
#   -DSQUAB_GEMMI_DIR=/path/to/gemmi       (default: ../checkout/gemmi if populated, else fetch v0.7.0)
#   -DSQUAB_JSON_INCLUDE=/path/to/include  (dir containing nlohmann/json.hpp; default: fetch v3.11.3)
#   -DSQUAB_SINGLE_FILE=OFF                (emit a separate saxs.wasm)
#   -DSQUAB_DEBUG=ON                       (-O1 -g, ASSERTIONS=2)

cmake_minimum_required(VERSION 3.18)
project(moorhen_squab_wasm CXX)

if(NOT EMSCRIPTEN)
    message(FATAL_ERROR "Configure this with emcmake (Emscripten), e.g. 'npm run build:wasm'.")
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(SQUAB_SINGLE_FILE "Embed saxs.wasm inside saxs.js" ON)
option(SQUAB_DEBUG "Debug build with Emscripten assertions" OFF)

get_filename_component(SQUAB_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/.." ABSOLUTE)

include(FetchContent)

# ---------------------------------------------------------------- gemmi
set(SQUAB_GEMMI_DIR "" CACHE PATH "gemmi source tree (must contain src/ and include/)")
if(NOT SQUAB_GEMMI_DIR AND EXISTS "${SQUAB_ROOT}/checkout/gemmi/src/mmcif.cpp")
    set(SQUAB_GEMMI_DIR "${SQUAB_ROOT}/checkout/gemmi")
endif()
if(NOT SQUAB_GEMMI_DIR)
    # v0.7.0 matches the headers vendored in ../include/gemmi.
    # SOURCE_SUBDIR points at a non-existent dir so gemmi's own CMake is not run.
    FetchContent_Declare(gemmi_src
        GIT_REPOSITORY https://github.com/project-gemmi/gemmi.git
        GIT_TAG v0.7.0
        GIT_SHALLOW TRUE
        SOURCE_SUBDIR _no_cmake_)
    FetchContent_MakeAvailable(gemmi_src)
    set(SQUAB_GEMMI_DIR "${gemmi_src_SOURCE_DIR}")
endif()
message(STATUS "gemmi sources: ${SQUAB_GEMMI_DIR}")

# ------------------------------------------------------- nlohmann_json
set(SQUAB_JSON_INCLUDE "" CACHE PATH "directory containing nlohmann/json.hpp")
if(NOT SQUAB_JSON_INCLUDE)
    FetchContent_Declare(nlohmann_json_src
        GIT_REPOSITORY https://github.com/nlohmann/json.git
        GIT_TAG v3.11.3
        GIT_SHALLOW TRUE
        SOURCE_SUBDIR _no_cmake_)
    FetchContent_MakeAvailable(nlohmann_json_src)
    set(SQUAB_JSON_INCLUDE "${nlohmann_json_src_SOURCE_DIR}/include")
endif()
message(STATUS "nlohmann_json include: ${SQUAB_JSON_INCLUDE}")

# ---------------------------------------------------------- flags
if(SQUAB_DEBUG)
    set(SQUAB_OPT -O1 -g)
else()
    set(SQUAB_OPT -O3)
endif()
set(SQUAB_EH -fexceptions)

# ---------------------------------------------------------- gemmi lib
# Same file list as the top-level CMakeLists.txt. Static archive, so only
# the objects actually referenced end up in the wasm.
set(G "${SQUAB_GEMMI_DIR}/src")
add_library(squab_gemmi STATIC
    ${G}/align.cpp ${G}/assembly.cpp ${G}/calculate.cpp ${G}/crd.cpp
    ${G}/ddl.cpp ${G}/eig3.cpp ${G}/intensit.cpp ${G}/json.cpp
    ${G}/mmcif.cpp ${G}/monlib.cpp ${G}/mtz.cpp ${G}/mtz2cif.cpp
    ${G}/pdb.cpp ${G}/polyheur.cpp ${G}/read_cif.cpp ${G}/read_map.cpp
    ${G}/resinfo.cpp ${G}/riding_h.cpp ${G}/select.cpp ${G}/sprintf.cpp
    ${G}/to_mmcif.cpp ${G}/to_pdb.cpp ${G}/topo.cpp ${G}/xds_ascii.cpp)
target_include_directories(squab_gemmi PUBLIC
    "${SQUAB_GEMMI_DIR}/include" "${SQUAB_GEMMI_DIR}/third_party")
target_compile_options(squab_gemmi PRIVATE ${SQUAB_OPT} ${SQUAB_EH})

# ---------------------------------------------------------- saxs module
set(S "${SQUAB_ROOT}/src")
add_executable(saxs
    ${S}/main.cpp
    ${S}/atomExtractor.cpp
    ${S}/formFactor.cpp
    ${S}/DebyeCalculator.cpp
    ${S}/scatteringModel.cpp
    ${S}/hydrationShell.cpp
    ${S}/proteinAtomType.cpp
    ${S}/excludedVolumeModel.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/cpp/squab_npm_bindings.cpp)
target_include_directories(saxs PRIVATE
    "${SQUAB_ROOT}/include" "${S}" "${SQUAB_JSON_INCLUDE}")
target_link_libraries(saxs PRIVATE squab_gemmi)
target_compile_options(saxs PRIVATE ${SQUAB_OPT} ${SQUAB_EH})

set(SQUAB_LINK
    ${SQUAB_OPT} ${SQUAB_EH}
    -lembind
    -sMODULARIZE=1
    -sEXPORT_ES6=1
    -sEXPORT_NAME=createSaxsModule
    -sENVIRONMENT=web,worker
    -sALLOW_MEMORY_GROWTH=1)
if(SQUAB_SINGLE_FILE)
    list(APPEND SQUAB_LINK -sSINGLE_FILE=1)
endif()
if(SQUAB_DEBUG)
    list(APPEND SQUAB_LINK -sASSERTIONS=2)
endif()
target_link_options(saxs PRIVATE ${SQUAB_LINK})
set_target_properties(saxs PROPERTIES SUFFIX ".js")
