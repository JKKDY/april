#pragma once

#include <iostream>
#include <string_view>

#ifndef NDEBUG
    #include <cstdlib>
    #include <stacktrace>
    #define APRIL_ASSERT(Expr, Msg) april::utility::internal::ap_assert(#Expr, (Expr), __FILE__, __LINE__, (Msg))
#else
    #define APRIL_ASSERT(Expr, Msg) ((void)0)
#endif

#define APRIL_CHECK(Expr, Msg) april::utility::internal::ap_assert(#Expr, (Expr), __FILE__, __LINE__, (Msg))

namespace april::utility::internal {

    inline void ap_assert(const char* expr_str, bool expr, const char* file, int line, std::string_view msg) {
        if (!expr) {
            std::cerr << "Assert failed:\t" << msg << "\n"
                      << "Expected:\t" << expr_str << "\n"
                      << "Source:\t\t" << file << ", line " << line << "\n"
                      << "Stack trace:\n" << std::stacktrace::current(1) << "\n";
            std::abort();
        }
    }

} // namespace april::utility::internal