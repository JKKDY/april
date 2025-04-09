include(FetchContent)

# Download and add unordered_dense (Martin Ankerl's flat hash map)
FetchContent_Declare(
    unordered_dense 
    GIT_REPOSITORY https://github.com/martinus/unordered_dense.git
    GIT_TAG v4.5.0
)
FetchContent_MakeAvailable(unordered_dense)

FetchContent_Declare(
    spdlog
    GIT_REPOSITORY https://github.com/gabime/spdlog.git
    GIT_TAG v1.15.2
)
FetchContent_MakeAvailable(unordered_dense spdlog)

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

# Link unordered_dense (included via FetchContent)
target_link_libraries(
    April 
  PUBLIC 
    unordered_dense::unordered_dense
    spdlog::spdlog
)
