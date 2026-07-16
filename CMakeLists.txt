cmake_minimum_required(VERSION 3.16)

project(saxs)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(gemmi_src ${CMAKE_CURRENT_SOURCE_DIR}/checkout/gemmi)

include_directories("${gemmi_src}/include"
                    "${gemmi_src}/third_party")



add_library(gemmi_cpp
            ${gemmi_src}/src/align.cpp
            ${gemmi_src}/src/assembly.cpp
            ${gemmi_src}/src/calculate.cpp
            ${gemmi_src}/src/crd.cpp
            ${gemmi_src}/src/ddl.cpp
            ${gemmi_src}/src/eig3.cpp
            ${gemmi_src}/src/intensit.cpp
            ${gemmi_src}/src/json.cpp
            ${gemmi_src}/src/mmcif.cpp
            ${gemmi_src}/src/monlib.cpp
            ${gemmi_src}/src/mtz.cpp
            ${gemmi_src}/src/mtz2cif.cpp
            ${gemmi_src}/src/pdb.cpp
            ${gemmi_src}/src/polyheur.cpp
            ${gemmi_src}/src/read_cif.cpp
            ${gemmi_src}/src/read_map.cpp
            ${gemmi_src}/src/resinfo.cpp
            ${gemmi_src}/src/riding_h.cpp
            ${gemmi_src}/src/select.cpp
            ${gemmi_src}/src/sprintf.cpp
            ${gemmi_src}/src/to_mmcif.cpp
            ${gemmi_src}/src/to_pdb.cpp
            ${gemmi_src}/src/topo.cpp
            ${gemmi_src}/src/xds_ascii.cpp)



add_executable(saxs
    src/main.cpp
    src/atomExtractor.cpp
    src/formFactor.cpp
    src/DebyeCalculator.cpp
    src/scatteringModel.cpp
)

target_link_libraries(saxs PRIVATE gemmi_cpp)

target_compile_options(saxs PRIVATE
    -O3
    --bind
)

target_link_options(saxs PRIVATE
    -O3
    -lembind
    -sMODULARIZE=1
    -sEXPORT_ES6=1
    -sENVIRONMENT=web
    -sALLOW_MEMORY_GROWTH=1
    -sEXPORT_NAME=createSaxsModule
    -sASSERTIONS=2
    -sNO_DISABLE_EXCEPTION_CATCHING
)

target_include_directories(saxs PRIVATE
    include
    include/webgpu
    ${CMAKE_CURRENT_SOURCE_DIR}/../saxs-profile-plugin/build/_deps/nlohmann_json-src/include

)