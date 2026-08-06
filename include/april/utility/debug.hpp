#pragma once

#include <cstdlib>
#include <iostream>
#include <string_view>
#include <version>

#if defined(APRIL_ENABLE_STACKTRACE) && defined(__cpp_lib_stacktrace) && __cpp_lib_stacktrace >= 202011L
    #include <stacktrace>
    #define APRIL_HAS_STACKTRACE 1
#else
    #define APRIL_HAS_STACKTRACE 0
#endif

#ifndef NDEBUG
    #define APRIL_ASSERT(Expr, Msg) \
    april::utility::internal::ap_assert( \
    #Expr, static_cast<bool>(Expr), __FILE__, __LINE__, (Msg))
#else
    #define APRIL_ASSERT(Expr, Msg) ((void)0)
#endif

#define APRIL_CHECK(Expr, Msg) \
april::utility::internal::ap_assert( \
#Expr, static_cast<bool>(Expr), __FILE__, __LINE__, (Msg))

namespace april::utility::internal {

    inline void ap_assert(const char* expr_str, bool expr, const char* file, int line, std::string_view msg) {
        if (expr) return;

        std::cerr << "Assert failed:\t" << msg << "\n"
                  << "Expected:\t" << expr_str << "\n"
                  << "Source:\t\t" << file << ", line " << line << "\n";

#if APRIL_HAS_STACKTRACE
        std::cerr << "Stack trace:\n"
                  << std::stacktrace::current(1) << "\n";
#endif

        std::abort();
    }

} // namespace april::utility::internal


#undef APRIL_HAS_STACKTRACE