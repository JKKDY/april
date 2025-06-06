include(FetchContent)

set(CMAKE_CXX_STANDARD 23)


# GLOB_RECURSE sources from src/
file(GLOB_RECURSE APRIL_SOURCES CONFIGURE_DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
)

# GLOB_RECURSE headers from include/
file(GLOB_RECURSE APRIL_HEADERS CONFIGURE_DEPENDS
    ${PROJECT_SOURCE_DIR}/include/*.hpp
    ${PROJECT_SOURCE_DIR}/include/*.h
)

# Show headers/sources in Visual Studio
source_group(TREE "${PROJECT_SOURCE_DIR}/include" PREFIX "Header Files" FILES ${APRIL_HEADERS})
source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${APRIL_SOURCES})

# Define the library target
add_library(April STATIC ${APRIL_SOURCES} ${APRIL_HEADERS})

# Public headers
target_include_directories(April PUBLIC
        ${PROJECT_SOURCE_DIR}/include
)

# Apply build flags
set_target_build_flags(April)

# enable sanitizer on debug builds
 enable_sanitizer_if_dbg(April)


