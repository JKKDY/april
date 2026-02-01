function(add_benchmark target_name source_file folder_name)
    add_executable(${target_name} ${source_file})
    target_compile_definitions(${target_name} PRIVATE PROJECT_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
    set_target_properties(${target_name} PROPERTIES FOLDER "benchmark/${folder_name}")

    set_target_build_flags(${target_name})
    target_link_libraries(${target_name} PRIVATE April)
endfunction()

add_benchmark(cube_benchmarkLL cube_benchmark_LL.cpp cube_benchmarkLL)
add_benchmark(cube_benchmarkDS cube_benchmark_DS.cpp cube_benchmarkDS)
add_benchmark(force_kernel_benchmark force_kernel_benchmark.cpp force_kernel_benchmark)

add_benchmark(no_overhead_DS_SOA baseline_DS_SOA.cpp no_overhead_DS_SOA)
add_benchmark(baseline_DS_AOS baseline_DS_AOS.cpp baseline_DS_AOS)
add_benchmark(baseline_DS_SIMD baseline_DS_SIMD.cpp baseline_DS_SIMD)
