function(add_benchmark target_name source_file folder_name)
    add_executable(${target_name} ${source_file})
    target_compile_definitions(${target_name} PRIVATE PROJECT_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
    set_target_properties(${target_name} PROPERTIES FOLDER "benchmark/${folder_name}")

    set_target_build_flags(${target_name})
    target_link_libraries(${target_name} PRIVATE April)
endfunction()

add_benchmark(cube_benchmarkLL cube_benchmark_LL.cpp cube_benchmarkLL)
add_benchmark(cube_benchmarkDS cube_benchmark_DS.cpp cube_benchmarkDS)
add_benchmark(no_overhead_DS no_overhead_DS.cpp no_overhead_DS)
add_benchmark(no_overhead_DS_AOS no_overhead_DS_AOS.cpp no_overhead_DS_AOS)

