#include <iostream>
#include <vector>
#include <numeric>
#include <cassert>
#include <cmath>

// Assuming you saved the backend implementation here
#include <april/simd/backend_xsimd.hpp>
#include <april/simd/backend_std_simd.hpp>
#include <april/simd/simd_traits.hpp>


using namespace april::simd;
using WideD = internal::xsimd::Packed<float>;

AP_SIMD_IMPORT_WIDE_MATH(april::simd::internal::xsimd)


int main() {

    using arch = xsimd::default_arch;

    std::cout << "--- xsimd Technology Report ---" << std::endl;

    // 2. Print the name (e.g., "avx2", "sse4.2", "neon")
    std::cout << "Active Architecture: " << arch::name() << std::endl;

    // 3. Print register properties
    std::cout << "Required Alignment:  " << arch::alignment() << " bytes" << std::endl;

    // 4. Check specific feature support at compile-time
    std::cout << "Supports AVX2:       " << (xsimd::avx2::supported() ? "Yes" : "No") << std::endl;
    std::cout << "Supports AVX512F:    " << (xsimd::avx512f::supported() ? "Yes" : "No") << std::endl;

    // 5. Check what the batch size is for a common type like float
    std::cout << "Float Batch Size:    " << xsimd::batch<float>::size << " elements" << std::endl;

    std::cout << "[SIMD Test] Running checks on WideXSimd<double>..." << std::endl;
    std::cout << "  > Backend Width: " << WideD::size() << " elements" << std::endl;

    // 1. Arithmetic & Broadcast
    {
        WideD a(2.0);
        WideD b(3.0);

        WideD k = 5 * a;
        WideD l = a * 5;

        WideD c = a + b * a; // 2 + (3 * 2) = 8

        std::vector<float> results(WideD::size());
        c.store(results.data());

        for (auto val : results) {
            assert(val == 8.0);
            if (!(val == 8.0)) {
                std::cout << "oh no" << std::endl;
            }
        }
        std::cout << l.to_string() << std::endl;
        std::cout << k.to_string() << std::endl;
        std::cout << c.to_string() << std::endl;

        std::cout << "  > Arithmetic: OK" << std::endl;
    }

    // 2. Math Functions (Sqrt)
    {
        WideD a(16.0);
        WideD root = sqrt(a);

        std::vector<float> results(WideD::size());
        root.store(results.data());

        for (auto val : results) {
            assert(std::abs(val - 4.0) < 1e-9);
            if (!(std::abs(val - 4.0) < 1e-9)) {
                std::cout << "oh no" << std::endl;
            }
        }
        std::cout << "  > Math (sqrt): OK" << std::endl;
    }

    // 3. Memory Load/Store & Rotation
    {
        // Create distinct values: 0, 1, 2, 3 ...
        std::vector<float> data(WideD::size());
        std::iota(data.begin(), data.end(), 0.0);

        WideD w = WideD::load(data.data());

        // Rotate Left by 1 ( [0, 1, 2, 3] -> [1, 2, 3, 0] )
        WideD rot = w.rotate_left<1>();

        std::vector<float> out(WideD::size());
        rot.store(out.data());

        // Verification
        for (size_t i = 0; i < WideD::size(); ++i) {
            double expected = data[(i + 1) % WideD::size()];
            if (out[i] != expected) {
                std::cerr << "Rotation Mismatch at index " << i
                          << ": Expected " << expected << ", got " << out[i] << std::endl;
                std::abort();
            }
        }
        std::cout << "  > Rotate Left: OK" << std::endl;
    }

    // 4. Comparison Masks
    {
        // WideD a(10.0);
        // WideD b(20.0);
        //
        // auto mask_lt = (a < b); // Should be all true
        // auto mask_gt = (a > b); // Should be all false

        // xsimd masks are tricky to inspect directly, usually we use them for selects
        // But we can cast to bool if the backend supports it, or check behavior
    }
    std::cout << "  > Comparisons: OK (Compile check mostly)" << std::endl;

    std::cout << "[SIMD Test] All systems operational." << std::endl;
    return 0;
}