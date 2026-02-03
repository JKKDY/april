#include <xsimd/xsimd.hpp>
#include <iostream>

int main() {
    std::cout << "xsimd version: " << XSIMD_VERSION_MAJOR << "."
              << XSIMD_VERSION_MINOR << "." << XSIMD_VERSION_PATCH << "\n";
    std::cout << "Instruction Set: " << xsimd::default_arch::name() << "\n";
}
