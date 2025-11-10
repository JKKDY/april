

function(add_example target_name source_file folder_name)
    add_executable(${target_name} ${source_file})
    target_compile_definitions(${target_name} PRIVATE PROJECT_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
    set_target_properties(${target_name} PROPERTIES FOLDER "examples/${folder_name}")

    target_link_libraries(${target_name} PRIVATE April)
    set_target_build_flags(${target_name})
    enable_sanitizer_if_dbg(${target_name})
endfunction()

add_example(halleys_comet halleys_comet.cpp halleys_comet)
add_example(two_body_collision two_body_collision.cpp two_body_collision)
add_example(falling_water_drop falling_water_drop.cpp falling_water_drop)
add_example(sandbox sandbox.cpp sandbox)



