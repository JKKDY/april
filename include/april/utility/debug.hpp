#pragma once

#ifndef NDEBUG

#include <iostream>
#include <cstdlib>


#define AP_ASSERT(Expr, Msg) april::utility::internal::ap_assert(#Expr, (Expr), __FILE__, __LINE__, (Msg))

namespace april::utility::internal {

    inline void ap_assert(const char* expr_str, bool expr, const char* file, int line, const char* msg) {
        if (!expr) {
            std::cerr << "Assert failed:\t" << msg << "\n"
                << "Expected:\t" << expr_str << "\n"
                << "Source:\t\t" << file << ", line " << line << "\n";
            std::abort();
        }
    }

    inline void ap_assert(const char* expr_str, bool expr, const char* file, int line, const std::string& msg) {
        ap_assert(expr_str, expr, file, line, msg.c_str());
    }

} // namespace april::utility::debug

#else

#define AP_ASSERT(Expr, Msg)

#endif












