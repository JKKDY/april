#pragma once

#ifdef DEBUG

#define ASSERT(Expr, Msg) md_Assert_(#Expr, (Expr), __FILE__, __LINE__, (Msg))
#include <iostream>
#include <stdexcept>


static void md_Assert_(const char* expr_str, const bool expr, const char* file, const int line, const char* msg)
{
    if (!expr) {
        std::cerr << "Assert failed:\t" << msg << "\n"
                << "Expected:\t" << expr_str << "\n"
                << "Source:\t\t" << file << ", line " << line << "\n";
        throw std::runtime_error(nullptr);
    }
}

#else

#define ASSERT(Expr, Msg)
#define W_ASSERT(Expr, Msg)

#endif