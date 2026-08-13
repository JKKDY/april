#pragma once

#if __has_include(<experimental/simd>)
    #include <experimental/simd>
#else
    #error "std::experimental::simd is not available"
#endif

#ifndef __cpp_lib_experimental_parallel_simd
    #error "The standard library does not provide Parallelism TS SIMD"
#endif


#include <array>
#include <cstddef>
#include <cstdint>
#include <sstream>
#include <string>
#include <type_traits>

#include "april/simd/packed_concept.hpp"


namespace april::simd::internal::std_simd {

    namespace stdx = std::experimental;

    template<typename T, size_t Width> struct Packed;


	template<typename T, size_t Width>
	struct native_simd_selector {
		using type = stdx::simd<
			std::remove_cv_t<T>,
			stdx::simd_abi::fixed_size<Width>
		>;
	};

	template<typename T>
	struct native_simd_selector<T, 0> {
		using type = stdx::simd<std::remove_cv_t<T>>;
	};

	template<typename T, size_t Width>
	using native_simd_t =
		typename native_simd_selector<T, Width>::type;




    template<typename T, size_t Width = 0>
    struct Mask {
	    using value_type = std::remove_cv_t<T>;
	    using native_type = native_simd_t<value_type, Width>::mask_type;
	    using mask_type = Mask;

	    native_type data;

	    Mask()
		    : data(false)
	    {}

	    Mask(bool value)
		    : data(value)
	    {}

	    Mask(const native_type& value)
		    : data(value)
	    {}

	    Mask(native_type&& value)
		    : data(std::move(value))
	    {}

	    Mask(const Mask&) = default;
	    Mask(Mask&&) = default;

	    Mask& operator=(const Mask&) = default;
	    Mask& operator=(Mask&&) = default;

	    static constexpr size_t size() {
		    return native_type::size();
	    }

	    template<typename U>
	    requires (
		    !std::same_as<T, U> &&
		    (size() == Mask<U, Width>::size())
	    )
	    Mask(const Mask<U, Width>& other)
		    : data(convert_mask(other))
	    {}

	    template<typename U>
	    requires (
		    !std::same_as<T, U> &&
		    (size() == Mask<U, Width>::size())
	    )
	    Mask& operator=(const Mask<U, Width>& other) {
		    data = convert_mask(other);
		    return *this;
	    }

	    Mask& operator=(bool value) {
		    data = native_type(value);
		    return *this;
	    }

	    [[nodiscard]] const native_type& native() const noexcept {
		    return data;
	    }

	    [[nodiscard]] native_type& native() noexcept {
		    return data;
	    }



	    // -----
	    // Loads
	    // -----
    	static Mask load(const bool* ptr) {
	    	return load_unaligned(ptr);
	    }

    	static Mask load_aligned(const bool* ptr) {
	    	return load_unaligned(ptr);
	    }

    	static Mask load_unaligned(const bool* ptr) {
	    	native_type result(false);

	    	for (size_t i = 0; i < size(); ++i)
	    		result[i] = ptr[i];

	    	return Mask{std::move(result)};
	    }


	    // ------
	    // Stores
	    // ------

    	void store(bool* ptr) const {
	    	store_unaligned(ptr);
	    }

    	void store_aligned(bool* ptr) const {
	    	store_unaligned(ptr);
	    }

    	void store_unaligned(bool* ptr) const {
	    	for (size_t i = 0; i < size(); ++i)
	    		ptr[i] = data[i];
	    }


	    // ------------------
	    // Logical reductions
	    // ------------------

	    [[nodiscard]] static bool all(const Mask& mask) {
		    return stdx::all_of(mask.data);
	    }

	    [[nodiscard]] static bool any(const Mask& mask) {
		    return stdx::any_of(mask.data);
	    }

	    [[nodiscard]] static bool none(const Mask& mask) {
		    return stdx::none_of(mask.data);
	    }


	    // -----------------
	    // Logical operators
	    // -----------------

	    friend Mask operator!(const Mask& mask) {
		    return Mask{!mask.data};
	    }

	    template<typename U>
	    requires (size() == Mask<U, Width>::size())
	    friend Mask operator&&(const Mask& lhs, const Mask<U, Width>& rhs) {
		    return Mask{lhs.data && convert_mask(rhs)};
	    }

	    template<typename U>
	    requires (size() == Mask<U, Width>::size())
	    friend Mask operator||(const Mask& lhs, const Mask<U, Width>& rhs) {
		    return Mask{lhs.data || convert_mask(rhs)};
	    }


	    // -----------------
	    // Bitwise operators
	    // -----------------

	    friend Mask operator~(const Mask& mask) {
		    return Mask{!mask.data};
	    }

	    template<typename U>
	    requires (size() == Mask<U, Width>::size())
	    friend Mask operator&(const Mask& lhs, const Mask<U, Width>& rhs) {
		    return Mask{lhs.data & convert_mask(rhs)};
	    }

	    template<typename U>
	    requires (size() == Mask<U, Width>::size())
	    friend Mask operator|(const Mask& lhs, const Mask<U, Width>& rhs) {
		    return Mask{lhs.data | convert_mask(rhs)};
	    }

	    template<typename U>
	    requires (size() == Mask<U, Width>::size())
	    friend Mask operator^(const Mask& lhs, const Mask<U, Width>& rhs) {
		    return Mask{lhs.data ^ convert_mask(rhs)};
	    }


	    // ------------------
	    // Lane-wise equality
	    // ------------------

	    template<typename U>
	    requires (size() == Mask<U, Width>::size())
	    friend Mask operator==(const Mask& lhs, const Mask<U, Width>& rhs) {
		    return Mask{lhs.data == convert_mask(rhs)};
	    }

	    template<typename U>
	    requires (size() == Mask<U, Width>::size())
	    friend Mask operator!=(const Mask& lhs, const Mask<U, Width>& rhs) {
		    return Mask{lhs.data != convert_mask(rhs)};
	    }


	    // ------------------
	    // Compound operators
	    // ------------------

	    template<typename U>
	    requires (size() == Mask<U, Width>::size())
	    Mask& operator&=(const Mask<U, Width>& rhs) {
		    data = data & convert_mask(rhs);
		    return *this;
	    }

	    template<typename U>
	    requires (size() == Mask<U, Width>::size())
	    Mask& operator|=(const Mask<U, Width>& rhs) {
		    data = data | convert_mask(rhs);
		    return *this;
	    }

	    template<typename U>
	    requires (size() == Mask<U, Width>::size())
	    Mask& operator^=(const Mask<U, Width>& rhs) {
		    data = data ^ convert_mask(rhs);
		    return *this;
	    }


	    // -----------------------
	    // Mutating lane rotations
	    // -----------------------

	    template<unsigned K = 1>
	    void rotate_left() {
		    constexpr size_t Shift = K % size();

		    if constexpr (Shift != 0) {
			    native_type rotated(false);

			    for (size_t i = 0; i < size(); ++i)
				    rotated[i] = data[(i + Shift) % size()];

			    data = rotated;
		    }
	    }

	    template<unsigned K = 1>
	    void rotate_right() {
		    constexpr size_t Shift = K % size();

		    if constexpr (Shift != 0) {
			    native_type rotated(false);

			    for (size_t i = 0; i < size(); ++i)
				    rotated[i] = data[(i + size() - Shift) % size()];

			    data = rotated;
		    }
	    }


	    // ---------------------
	    // Bitmask import/export
	    // ---------------------

	    [[nodiscard]] uint64_t to_bitmask() const {
		    static_assert(
			    size() <= 64,
			    "Mask bit export supports at most 64 lanes"
		    );

		    uint64_t bits = 0;

		    for (size_t i = 0; i < size(); ++i) {
			    if (data[i])
				    bits |= uint64_t{1} << i;
		    }

		    return bits;
	    }

	    static Mask from_bitmask(const uint64_t bits) {
		    static_assert(
			    size() <= 64,
			    "Mask bit import supports at most 64 lanes"
		    );

		    native_type result(false);

		    for (size_t i = 0; i < size(); ++i)
			    result[i] = ((bits >> i) & uint64_t{1}) != 0;

		    return Mask{std::move(result)};
	    }


	    // -------------------------
	    // Debugging and inspection
	    // -------------------------

    	[[nodiscard]] std::array<bool, size()> to_array() const {
	    	std::array<bool, size()> result;
	    	store(result.data());
	    	return result;
	    }

	    [[nodiscard]] std::string to_string() const {
		    const auto arr = to_array();

		    std::stringstream ss;
		    ss << "[";

		    for (size_t i = 0; i < size(); ++i) {
			    ss << (arr[i] ? "true" : "false");
			    if (i < size() - 1)
				    ss << ", ";
		    }

		    ss << "]";
		    return ss.str();
	    }

    private:
	    template<typename U>
	    [[nodiscard]] static native_type convert_mask(
		    const Mask<U, Width>& other
	    ) {
		    static_assert(
			    size() == Mask<U, Width>::size(),
			    "Cannot convert masks with different logical widths"
		    );

		    native_type converted(false);

		    for (size_t i = 0; i < size(); ++i)
			    converted[i] = other.native()[i];

		    return converted;
	    }
    };







    // Width == 0: Use Native ABI (Best fit for hardware, e.g. 4 doubles on AVX2)
    // Width > 0:  Use Fixed Size ABI (Compiler manages register spanning, e.g. 16 doubles)
    template<typename T, size_t Width = 0>
    struct Packed {
    	using value_type = std::remove_cv_t<T>;
    	using storage_type = packed_storage_t<value_type>;
    	using native_type = native_simd_t<storage_type, Width>;
    	using abi_type = native_type::abi_type;
    	using mask_type = Mask<storage_type, Width>;
    	using packed_type = Packed;

        native_type data;

        // Construction
        Packed() = default;

        Packed(value_type scalar)
            : data(load_scalar(scalar))
        {}

        Packed(const native_type& value)
            : data(value)
        {}

        Packed(native_type&& value)
            : data(std::move(value))
        {}

        Packed(const Packed&) = default;
        Packed(Packed&&) = default;

        Packed& operator=(const Packed&) = default;
        Packed& operator=(Packed&&) = default;

        Packed& operator=(value_type scalar) {
            data = native_type(load_scalar(scalar));
            return *this;
        }

        // Lane count
        [[nodiscard]] static constexpr size_t size() {
            return native_type::size();
        }

        // Native representation access
        [[nodiscard]] const native_type& native() const noexcept {
            return data;
        }

        [[nodiscard]] native_type& native() noexcept {
            return data;
        }


    	// ----------------
		// Contiguous loads
		// ----------------
		template<typename PtrT>
		requires (
			IsPackableValue<std::remove_cv_t<PtrT>> &&
			(sizeof(PtrT) <= sizeof(storage_type))
		)
		static Packed load(const PtrT* ptr) {
			return load_unaligned(ptr);
		}

		template<typename PtrT>
		requires (
			IsPackableValue<std::remove_cv_t<PtrT>> &&
			(sizeof(PtrT) <= sizeof(storage_type))
		)
		static Packed load_aligned(const PtrT* ptr) {
			if constexpr (std::is_enum_v<std::remove_cv_t<PtrT>>) {
				alignas(stdx::memory_alignment_v<native_type, storage_type>)
				std::array<storage_type, size()> temp;

				for (size_t i = 0; i < size(); ++i)
					temp[i] = load_scalar(ptr[i]);

				native_type result;
				result.copy_from(temp.data(), stdx::vector_aligned);
				return Packed{std::move(result)};
			}
			else {
				native_type result;
				result.copy_from(ptr, stdx::vector_aligned);
				return Packed{std::move(result)};
			}
		}

		template<typename PtrT>
		requires (
			IsPackableValue<std::remove_cv_t<PtrT>> &&
			(sizeof(PtrT) <= sizeof(storage_type))
		)
		static Packed load_unaligned(const PtrT* ptr) {
			if constexpr (std::is_enum_v<std::remove_cv_t<PtrT>>) {
				alignas(stdx::memory_alignment_v<native_type, storage_type>)
				std::array<storage_type, size()> temp;

				for (size_t i = 0; i < size(); ++i)
					temp[i] = load_scalar(ptr[i]);

				native_type result;
				result.copy_from(temp.data(), stdx::vector_aligned);
				return Packed{std::move(result)};
			}
			else {
				native_type result;
				result.copy_from(ptr, stdx::element_aligned);
				return Packed{std::move(result)};
			}
		}


		// -----------------
		// Contiguous stores
		// -----------------
		template<typename PtrT>
		requires (
			!std::is_const_v<PtrT> &&
			IsPackableValue<std::remove_cv_t<PtrT>> &&
			(sizeof(PtrT) <= sizeof(storage_type))
		)
		void store(PtrT* ptr) const {
			store_unaligned(ptr);
		}

		template<typename PtrT>
		requires (
			!std::is_const_v<PtrT> &&
			IsPackableValue<std::remove_cv_t<PtrT>> &&
			(sizeof(PtrT) <= sizeof(storage_type))
		)
		void store_aligned(PtrT* ptr) const {
			if constexpr (std::is_enum_v<std::remove_cv_t<PtrT>>) {
				alignas(stdx::memory_alignment_v<native_type, storage_type>)
				std::array<storage_type, size()> temp;

				data.copy_to(temp.data(), stdx::vector_aligned);

				for (size_t i = 0; i < size(); ++i)
					ptr[i] = store_scalar<PtrT>(temp[i]);
			}
			else {
				data.copy_to(ptr, stdx::vector_aligned);
			}
		}

		template<typename PtrT>
		requires (
			!std::is_const_v<PtrT> &&
			IsPackableValue<std::remove_cv_t<PtrT>> &&
			(sizeof(PtrT) <= sizeof(storage_type))
		)
		void store_unaligned(PtrT* ptr) const {
			if constexpr (std::is_enum_v<std::remove_cv_t<PtrT>>) {
				alignas(stdx::memory_alignment_v<native_type, storage_type>)
				std::array<storage_type, size()> temp;

				data.copy_to(temp.data(), stdx::vector_aligned);

				for (size_t i = 0; i < size(); ++i)
					ptr[i] = store_scalar<PtrT>(temp[i]);
			}
			else {
				data.copy_to(ptr, stdx::element_aligned);
			}
		}


    	// ------------
    	// Gather Loads
    	// ------------
    	// Compile-time byte stride
    	template<std::ptrdiff_t ByteStride, typename PtrT>
		requires (
			IsPackableValue<std::remove_cv_t<PtrT>> &&
			(sizeof(PtrT) <= sizeof(storage_type))
		)
		static Packed gather_strided(const PtrT* ptr) {
        	static_assert(ByteStride > 0, "Byte stride must be positive");

        	if constexpr (ByteStride == sizeof(PtrT)) {
        		return load(ptr);
        	}
        	else if constexpr (ByteStride % sizeof(PtrT) == 0) {
        		constexpr std::ptrdiff_t ElementStride =
					ByteStride / sizeof(PtrT);

        		return Packed{
        			native_type([&](size_t i) {
						return load_scalar(
							ptr[static_cast<std::ptrdiff_t>(i) * ElementStride]
						);
					})
				};
        	}
        	else {
        		const auto* bytes = reinterpret_cast<const std::byte*>(ptr);

        		return Packed{
        			native_type([&](size_t i) {
						PtrT value;

						std::memcpy(
							&value,
							bytes + static_cast<std::ptrdiff_t>(i) * ByteStride,
							sizeof(PtrT)
						);

						return load_scalar(value);
					})
				};
        	}
        }


    	// Runtime byte stride
    	template<typename PtrT>
		requires (
			IsPackableValue<std::remove_cvref_t<PtrT>> &&
			(sizeof(PtrT) <= sizeof(storage_type))
		)
		static Packed gather_strided(
			const PtrT* ptr,
			const std::ptrdiff_t byte_stride
		) {
        	if (byte_stride == static_cast<std::ptrdiff_t>(sizeof(PtrT)))
        		return load(ptr);

        	if (byte_stride % static_cast<std::ptrdiff_t>(sizeof(PtrT)) == 0) {
        		const std::ptrdiff_t element_stride =
					byte_stride / static_cast<std::ptrdiff_t>(sizeof(PtrT));

        		return Packed{
        			native_type([&](size_t i) {
						return load_scalar(
							ptr[static_cast<std::ptrdiff_t>(i) * element_stride]
						);
					})
				};
        	}

        	const auto* bytes = reinterpret_cast<const std::byte*>(ptr);

        	return Packed{
        		native_type([&](size_t i) {
					PtrT value;

					std::memcpy(
						&value,
						bytes + static_cast<std::ptrdiff_t>(i) * byte_stride,
						sizeof(PtrT)
					);

					return load_scalar(value);
				})
			};
        }


    	// Arbitrary byte offsets
    	template<typename PtrT, size_t N>
		requires (
			IsPackableValue<std::remove_cvref_t<PtrT>> &&
			(sizeof(PtrT) <= sizeof(storage_type)) &&
			(N == size())
		)
		static Packed gather(
			const PtrT* ptr,
			const ByteOffsets<N>& offsets
		) {
        	bool element_aligned = true;

        	for (size_t i = 0; i < size(); ++i) {
        		if (offsets.values[i] % static_cast<std::ptrdiff_t>(sizeof(PtrT)) != 0) {
        			element_aligned = false;
        			break;
        		}
        	}

        	if (element_aligned) {
        		return Packed{
        			native_type([&](size_t i) {
						const std::ptrdiff_t index =
							offsets.values[i] / static_cast<std::ptrdiff_t>(sizeof(PtrT));

						return load_scalar(ptr[index]);
					})
				};
        	}

        	const auto* bytes =
				reinterpret_cast<const std::byte*>(ptr);

        	return Packed{
        		native_type([&](size_t i) {
					PtrT value;

					std::memcpy(
						&value,
						bytes + offsets.values[i],
						sizeof(PtrT)
					);

					return load_scalar(value);
				})
			};
        }


    	// Arbitrary pointers
    	template<typename PointerContainer>
		requires requires(const PointerContainer& pointers) {
        	pointers[size_t{}];

        	requires std::is_pointer_v<
				std::remove_cvref_t<decltype(pointers[size_t{}])>
			>;
		}
    	static Packed gather(const PointerContainer& pointers) {
        	using pointer_type =
				std::remove_cvref_t<decltype(pointers[size_t{}])>;

        	using pointed_type =
				std::remove_pointer_t<pointer_type>;

        	using source_type =
				std::remove_cv_t<pointed_type>;

        	static_assert(
				IsPackableValue<source_type>,
				"Gather pointers must point to arithmetic or enum values"
			);

        	static_assert(
				sizeof(source_type) <= sizeof(storage_type),
				"Gather source type is wider than packed storage type"
			);

        	return Packed{
        		native_type([&](size_t i) {
					return load_scalar(*pointers[i]);
				})
			};
        }

    	// Arbitrary byte offsets
    	template<typename PtrT, size_t N>
		requires (
			!std::is_const_v<PtrT> &&
			IsPackableValue<std::remove_cv_t<PtrT>> &&
			(sizeof(PtrT) <= sizeof(storage_type)) &&
			(N == size())
		)
		void scatter(
			PtrT* ptr,
			const ByteOffsets<N>& offsets
		) const {
        	bool element_aligned = true;

        	for (size_t i = 0; i < size(); ++i) {
        		if (offsets.values[i] % static_cast<std::ptrdiff_t>(sizeof(PtrT)) != 0) {
        			element_aligned = false;
        			break;
        		}
        	}

        	if (element_aligned) {
        		for (size_t i = 0; i < size(); ++i) {
        			const std::ptrdiff_t index =
						offsets.values[i] / static_cast<std::ptrdiff_t>(sizeof(PtrT));

        			ptr[index] = store_scalar<PtrT>(data[i]);
        		}

        		return;
        	}

        	auto* bytes =
				reinterpret_cast<std::byte*>(ptr);

        	for (size_t i = 0; i < size(); ++i) {
        		const PtrT value =
					store_scalar<PtrT>(data[i]);

        		std::memcpy(
					bytes + offsets.values[i],
					&value,
					sizeof(PtrT)
				);
        	}
        }


    	// Arbitrary pointers
    	template<typename PointerContainer>
		requires requires(const PointerContainer& pointers) {
        	pointers[size_t{}];

        	requires std::is_pointer_v<
				std::remove_cvref_t<decltype(pointers[size_t{}])>
			>;
		}
    	void scatter(const PointerContainer& pointers) const {
        	using pointer_type =
				std::remove_cvref_t<decltype(pointers[size_t{}])>;

        	using pointed_type =
				std::remove_pointer_t<pointer_type>;

        	using destination_type =
				std::remove_cv_t<pointed_type>;

        	static_assert(
				!std::is_const_v<pointed_type>,
				"Scatter pointers must point to writable values"
			);

        	static_assert(
				IsPackableValue<destination_type>,
				"Scatter pointers must point to arithmetic or enum values"
			);

        	static_assert(
				sizeof(destination_type) <= sizeof(storage_type),
				"Scatter destination type is wider than packed storage type"
			);

        	for (size_t i = 0; i < size(); ++i) {
        		*pointers[i] =
					store_scalar<destination_type>(data[i]);
        	}
        }



    	// --------------
    	// Scatter Stores
    	// --------------
    	// Compile-time byte stride
    	template<std::ptrdiff_t ByteStride, typename PtrT>
		requires (
			!std::is_const_v<PtrT> &&
			IsPackableValue<std::remove_cv_t<PtrT>> &&
			(sizeof(PtrT) <= sizeof(storage_type))
		)
		void scatter_strided(PtrT* ptr) const {
        	static_assert(ByteStride > 0, "Byte stride must be positive");

        	if constexpr (ByteStride == sizeof(PtrT)) {
        		store(ptr);
        	}
        	else if constexpr (ByteStride % sizeof(PtrT) == 0) {
        		constexpr std::ptrdiff_t ElementStride =
					ByteStride / sizeof(PtrT);

        		for (size_t i = 0; i < size(); ++i) {
        			ptr[static_cast<std::ptrdiff_t>(i) * ElementStride] =
						store_scalar<PtrT>(data[i]);
        		}
        	}
        	else {
        		auto* bytes =
					reinterpret_cast<std::byte*>(ptr);

        		for (size_t i = 0; i < size(); ++i) {
        			const PtrT value =
						store_scalar<PtrT>(data[i]);

        			std::memcpy(
						bytes + static_cast<std::ptrdiff_t>(i) * ByteStride,
						&value,
						sizeof(PtrT)
					);
        		}
        	}
        }


    	// Runtime byte stride
    	template<typename PtrT>
		requires (
			!std::is_const_v<PtrT> &&
			IsPackableValue<std::remove_cv_t<PtrT>> &&
			(sizeof(PtrT) <= sizeof(storage_type))
		)
		void scatter_strided(
			PtrT* ptr,
			const std::ptrdiff_t byte_stride
		) const {
        	if (byte_stride == static_cast<std::ptrdiff_t>(sizeof(PtrT))) {
        		store(ptr);
        		return;
        	}

        	if (byte_stride % static_cast<std::ptrdiff_t>(sizeof(PtrT)) == 0) {
        		const std::ptrdiff_t element_stride =
					byte_stride / static_cast<std::ptrdiff_t>(sizeof(PtrT));

        		for (size_t i = 0; i < size(); ++i) {
        			ptr[static_cast<std::ptrdiff_t>(i) * element_stride] =
						store_scalar<PtrT>(data[i]);
        		}

        		return;
        	}

        	auto* bytes =
				reinterpret_cast<std::byte*>(ptr);

        	for (size_t i = 0; i < size(); ++i) {
        		const PtrT value =
					store_scalar<PtrT>(data[i]);

        		std::memcpy(
					bytes + static_cast<std::ptrdiff_t>(i) * byte_stride,
					&value,
					sizeof(PtrT)
				);
        	}
        }



    	// ----------------------
    	// Permutes and rotations
    	// ----------------------
    	template<size_t... Indices>
		[[nodiscard]] Packed permute() const {
        	static_assert(sizeof...(Indices) > 0, "A permutation requires at least one index");
        	static_assert(((Indices < size()) && ...), "Permutation index is outside the pack");
        	static_assert(
				sizeof...(Indices) == 1 || sizeof...(Indices) == size(),
				"A permutation must provide one index or exactly one index per lane"
			);

        	constexpr std::array<size_t, sizeof...(Indices)> indices{Indices...};

        	return Packed{
        		native_type([&](size_t i) {
					if constexpr (sizeof...(Indices) == 1)
						return data[indices[0]];
					else
						return data[indices[i]];
				})
			};
        }

    	template<unsigned K = 1>
		[[nodiscard]] Packed rotate_left() const {
        	constexpr size_t Shift = K % size();

        	if constexpr (Shift == 0) {
        		return *this;
        	}
        	else {
        		return Packed{
        			native_type([&](size_t i) {
						return data[(i + Shift) % size()];
					})
				};
        	}
        }

    	template<unsigned K = 1>
		[[nodiscard]] Packed rotate_right() const {
        	constexpr size_t Shift = K % size();

        	if constexpr (Shift == 0) {
        		return *this;
        	}
        	else {
        		return Packed{
        			native_type([&](size_t i) {
						return data[(i + size() - Shift) % size()];
					})
				};
        	}
        }



    	// --------------------
		// Arithmetic operators
		// --------------------

		friend Packed operator+(const Packed& value) {
			return Packed{+value.data};
		}

		friend Packed operator-(const Packed& value) {
			return Packed{-value.data};
		}

		friend Packed operator+(const Packed& lhs, const Packed& rhs) {
			return Packed{lhs.data + rhs.data};
		}

		friend Packed operator-(const Packed& lhs, const Packed& rhs) {
			return Packed{lhs.data - rhs.data};
		}

		friend Packed operator*(const Packed& lhs, const Packed& rhs) {
			return Packed{lhs.data * rhs.data};
		}

		friend Packed operator/(const Packed& lhs, const Packed& rhs) {
			return Packed{lhs.data / rhs.data};
		}

		friend Packed operator+(const Packed& lhs, storage_type rhs) {
			return Packed{lhs.data + native_type(rhs)};
		}

		friend Packed operator+(storage_type lhs, const Packed& rhs) {
			return Packed{native_type(lhs) + rhs.data};
		}

		friend Packed operator-(const Packed& lhs, storage_type rhs) {
			return Packed{lhs.data - native_type(rhs)};
		}

		friend Packed operator-(storage_type lhs, const Packed& rhs) {
			return Packed{native_type(lhs) - rhs.data};
		}

		friend Packed operator*(const Packed& lhs, storage_type rhs) {
			return Packed{lhs.data * native_type(rhs)};
		}

		friend Packed operator*(storage_type lhs, const Packed& rhs) {
			return Packed{native_type(lhs) * rhs.data};
		}

		friend Packed operator/(const Packed& lhs, storage_type rhs) {
			return Packed{lhs.data / native_type(rhs)};
		}

		friend Packed operator/(storage_type lhs, const Packed& rhs) {
			return Packed{native_type(lhs) / rhs.data};
		}


		// ---------
		// Compounds
		// ---------

		Packed& operator+=(const Packed& rhs) {
			data += rhs.data;
			return *this;
		}

		Packed& operator-=(const Packed& rhs) {
			data -= rhs.data;
			return *this;
		}

		Packed& operator*=(const Packed& rhs) {
			data *= rhs.data;
			return *this;
		}

		Packed& operator/=(const Packed& rhs) {
			data /= rhs.data;
			return *this;
		}

		Packed& operator+=(storage_type rhs) {
			data += native_type(rhs);
			return *this;
		}

		Packed& operator-=(storage_type rhs) {
			data -= native_type(rhs);
			return *this;
		}

		Packed& operator*=(storage_type rhs) {
			data *= native_type(rhs);
			return *this;
		}

		Packed& operator/=(storage_type rhs) {
			data /= native_type(rhs);
			return *this;
		}


		// -----------
		// Comparisons
		// -----------

		friend mask_type operator==(const Packed& lhs, const Packed& rhs) {
			return compare_impl(lhs.data, rhs.data, [](const auto& a, const auto& b) {
				return a == b;
			});
		}

		friend mask_type operator!=(const Packed& lhs, const Packed& rhs) {
			return compare_impl(lhs.data, rhs.data, [](const auto& a, const auto& b) {
				return a != b;
			});
		}

		friend mask_type operator<(const Packed& lhs, const Packed& rhs) {
			return compare_impl(lhs.data, rhs.data, [](const auto& a, const auto& b) {
				return a < b;
			});
		}

		friend mask_type operator<=(const Packed& lhs, const Packed& rhs) {
			return compare_impl(lhs.data, rhs.data, [](const auto& a, const auto& b) {
				return a <= b;
			});
		}

		friend mask_type operator>(const Packed& lhs, const Packed& rhs) {
			return compare_impl(lhs.data, rhs.data, [](const auto& a, const auto& b) {
				return a > b;
			});
		}

		friend mask_type operator>=(const Packed& lhs, const Packed& rhs) {
			return compare_impl(lhs.data, rhs.data, [](const auto& a, const auto& b) {
				return a >= b;
			});
		}

		friend mask_type operator==(const Packed& lhs, storage_type rhs) {
			return compare_impl(lhs.data, native_type(rhs), [](const auto& a, const auto& b) {
				return a == b;
			});
		}

		friend mask_type operator==(storage_type lhs, const Packed& rhs) {
			return compare_impl(native_type(lhs), rhs.data, [](const auto& a, const auto& b) {
				return a == b;
			});
		}

		friend mask_type operator!=(const Packed& lhs, storage_type rhs) {
			return compare_impl(lhs.data, native_type(rhs), [](const auto& a, const auto& b) {
				return a != b;
			});
		}

		friend mask_type operator!=(storage_type lhs, const Packed& rhs) {
			return compare_impl(native_type(lhs), rhs.data, [](const auto& a, const auto& b) {
				return a != b;
			});
		}

		friend mask_type operator<(const Packed& lhs, storage_type rhs) {
			return compare_impl(lhs.data, native_type(rhs), [](const auto& a, const auto& b) {
				return a < b;
			});
		}

		friend mask_type operator<(storage_type lhs, const Packed& rhs) {
			return compare_impl(native_type(lhs), rhs.data, [](const auto& a, const auto& b) {
				return a < b;
			});
		}

		friend mask_type operator<=(const Packed& lhs, storage_type rhs) {
			return compare_impl(lhs.data, native_type(rhs), [](const auto& a, const auto& b) {
				return a <= b;
			});
		}

		friend mask_type operator<=(storage_type lhs, const Packed& rhs) {
			return compare_impl(native_type(lhs), rhs.data, [](const auto& a, const auto& b) {
				return a <= b;
			});
		}

		friend mask_type operator>(const Packed& lhs, storage_type rhs) {
			return compare_impl(lhs.data, native_type(rhs), [](const auto& a, const auto& b) {
				return a > b;
			});
		}

		friend mask_type operator>(storage_type lhs, const Packed& rhs) {
			return compare_impl(native_type(lhs), rhs.data, [](const auto& a, const auto& b) {
				return a > b;
			});
		}

		friend mask_type operator>=(const Packed& lhs, storage_type rhs) {
			return compare_impl(lhs.data, native_type(rhs), [](const auto& a, const auto& b) {
				return a >= b;
			});
		}

		friend mask_type operator>=(storage_type lhs, const Packed& rhs) {
			return compare_impl(native_type(lhs), rhs.data, [](const auto& a, const auto& b) {
				return a >= b;
			});
		}


		// ---------
		// Selection
		// ---------

		[[nodiscard]] static Packed select(
			const mask_type& mask,
			const Packed& true_value,
			const Packed& false_value
		) {
			native_type result = false_value.data;
			stdx::where(mask.native(), result) = true_value.data;
			return Packed{std::move(result)};
		}


		// ----------------
		// Basic operations
		// ----------------

		[[nodiscard]] static Packed abs(const Packed& x) {
			if constexpr (std::is_unsigned_v<storage_type>)
				return x;
			else
				return Packed{stdx::abs(x.data)};
		}

		[[nodiscard]] static Packed min(const Packed& a, const Packed& b) {
			return Packed{stdx::min(a.data, b.data)};
		}

		[[nodiscard]] static Packed max(const Packed& a, const Packed& b) {
			return Packed{stdx::max(a.data, b.data)};
		}

		[[nodiscard]] static Packed clamp(
			const Packed& x,
			const Packed& lower,
			const Packed& upper
		) {
			return Packed{stdx::clamp(x.data, lower.data, upper.data)};
		}


		// ----------------
		// Roots and powers
		// ----------------
    	[[nodiscard]] static Packed sqrt(const Packed& x) requires std::floating_point<storage_type> {
        	return Packed{stdx::sqrt(x.data)};
        }

    	[[nodiscard]] static Packed rsqrt(const Packed& x) requires std::floating_point<storage_type> {
        	return Packed{
        		native_type(storage_type{1}) / stdx::sqrt(x.data)
			 };
        }

    	[[nodiscard]] static Packed cbrt(const Packed& x) requires std::floating_point<storage_type> {
        	return Packed{stdx::cbrt(x.data)};
        }

    	[[nodiscard]] static Packed hypot(const Packed& x, const Packed& y) requires std::floating_point<storage_type> {
        	return Packed{stdx::hypot(x.data, y.data)};
        }

    	[[nodiscard]] static Packed pow(const Packed& x, const Packed& y) requires std::floating_point<storage_type> {
        	return Packed{stdx::pow(x.data, y.data)};
        }


		// ---------------------------
		// Exponential and logarithmic
		// ---------------------------
    	[[nodiscard]] static Packed exp(const Packed& x) requires std::floating_point<storage_type> {
        	return Packed{stdx::exp(x.data)};
        }

    	[[nodiscard]] static Packed exp2(const Packed& x) requires std::floating_point<storage_type> {
        	return Packed{stdx::exp2(x.data)};
        }

    	[[nodiscard]] static Packed expm1(const Packed& x) requires std::floating_point<storage_type> {
        	return Packed{stdx::expm1(x.data)};
        }

    	[[nodiscard]] static Packed log(const Packed& x) requires std::floating_point<storage_type> {
        	return Packed{stdx::log(x.data)};
        }

    	[[nodiscard]] static Packed ln(const Packed& x) requires std::floating_point<storage_type> {
        	return log(x);
        }

    	[[nodiscard]] static Packed log2(const Packed& x) requires std::floating_point<storage_type> {
        	return Packed{stdx::log2(x.data)};
        }

    	[[nodiscard]] static Packed log10(const Packed& x) requires std::floating_point<storage_type> {
        	return Packed{stdx::log10(x.data)};
        }

    	[[nodiscard]] static Packed log1p(const Packed& x) requires std::floating_point<storage_type> {
        	return Packed{stdx::log1p(x.data)};
        }


		// -------------
		// Trigonometric
		// -------------
    	[[nodiscard]] static Packed sin(const Packed& x) requires std::floating_point<storage_type> {
        	return Packed{stdx::sin(x.data)};
        }

    	[[nodiscard]] static Packed cos(const Packed& x) requires std::floating_point<storage_type> {
        	return Packed{stdx::cos(x.data)};
        }

    	[[nodiscard]] static std::pair<Packed, Packed> sincos(const Packed& x) requires std::floating_point<storage_type> {
        	return { sin(x), cos(x) };
        }

    	[[nodiscard]] static Packed tan(const Packed& x) requires std::floating_point<storage_type> {
        	return Packed{stdx::tan(x.data)};
        }

    	[[nodiscard]] static Packed asin(const Packed& x) requires std::floating_point<storage_type> {
        	return Packed{stdx::asin(x.data)};
        }

    	[[nodiscard]] static Packed acos(const Packed& x) requires std::floating_point<storage_type> {
        	return Packed{stdx::acos(x.data)};
        }

    	[[nodiscard]] static Packed atan(const Packed& x) requires std::floating_point<storage_type> {
        	return Packed{stdx::atan(x.data)};
        }

    	[[nodiscard]] static Packed atan2(const Packed& y, const Packed& x) requires std::floating_point<storage_type> {
        	return Packed{stdx::atan2(y.data, x.data)};
        }


		// ----------
		// Hyperbolic
		// ----------
    	[[nodiscard]] static Packed sinh(const Packed& x) requires std::floating_point<storage_type> {
        	return Packed{stdx::sinh(x.data)};
        }

    	[[nodiscard]] static Packed cosh(const Packed& x) requires std::floating_point<storage_type> {
        	return Packed{stdx::cosh(x.data)};
        }

    	[[nodiscard]] static Packed tanh(const Packed& x) requires std::floating_point<storage_type> {
        	return Packed{stdx::tanh(x.data)};
        }

    	[[nodiscard]] static Packed asinh(const Packed& x) requires std::floating_point<storage_type> {
        	return Packed{stdx::asinh(x.data)};
        }

    	[[nodiscard]] static Packed acosh(const Packed& x) requires std::floating_point<storage_type> {
        	return Packed{stdx::acosh(x.data)};
        }

    	[[nodiscard]] static Packed atanh(const Packed& x) requires std::floating_point<storage_type> {
        	return Packed{stdx::atanh(x.data)};
        }


		// --------
		// Rounding
		// --------

		[[nodiscard]] static Packed floor(const Packed& x) {
			if constexpr (std::integral<storage_type>)
				return x;
			else
				return Packed{stdx::floor(x.data)};
		}

		[[nodiscard]] static Packed ceil(const Packed& x) {
			if constexpr (std::integral<storage_type>)
				return x;
			else
				return Packed{stdx::ceil(x.data)};
		}

		[[nodiscard]] static Packed round(const Packed& x) {
			if constexpr (std::integral<storage_type>)
				return x;
			else
				return Packed{stdx::round(x.data)};
		}

		[[nodiscard]] static Packed trunc(const Packed& x) {
			if constexpr (std::integral<storage_type>)
				return x;
			else
				return Packed{stdx::trunc(x.data)};
		}

		[[nodiscard]] static Packed nearbyint(const Packed& x) {
			if constexpr (std::integral<storage_type>)
				return x;
			else
				return Packed{stdx::nearbyint(x.data)};
		}


		// ------------------
		// Numeric operations
		// ------------------
		[[nodiscard]] static Packed fma(
			const Packed& x,
			const Packed& y,
			const Packed& z
		) {
			if constexpr (std::integral<storage_type>)
				return Packed{x.data * y.data + z.data};
			else
				return Packed{stdx::fma(x.data, y.data, z.data)};
		}

		[[nodiscard]] static Packed fmod(const Packed& x, const Packed& y)
    	requires std::floating_point<storage_type> {
			return Packed{stdx::fmod(x.data, y.data)};
		}

		[[nodiscard]] static Packed remainder(const Packed& x, const Packed& y)
    	requires std::floating_point<storage_type> {
			return Packed{stdx::remainder(x.data, y.data)};
		}

		[[nodiscard]] static Packed copysign(const Packed& magnitude, const Packed& sign )
		requires std::floating_point<storage_type> {
			return Packed{
				stdx::copysign(magnitude.data, sign.data)
			};
		}


		// --------------
		// Classification
		// --------------

		[[nodiscard]] static mask_type isnan(const Packed& x) {
			if constexpr (std::floating_point<storage_type>)
				return mask_type{stdx::isnan(x.data)};
			else
				return mask_type{false};
		}

		[[nodiscard]] static mask_type isinf(const Packed& x) {
			if constexpr (std::floating_point<storage_type>)
				return mask_type{stdx::isinf(x.data)};
			else
				return mask_type{false};
		}

		[[nodiscard]] static mask_type isfinite(const Packed& x) {
			if constexpr (std::floating_point<storage_type>)
				return mask_type{stdx::isfinite(x.data)};
			else
				return mask_type{true};
		}

		[[nodiscard]] static mask_type signbit(const Packed& x) {
			if constexpr (std::is_unsigned_v<storage_type>)
				return mask_type{false};
			else if constexpr (std::integral<storage_type>)
				return mask_type{x.data < native_type(storage_type{0})};
			else
				return mask_type{stdx::signbit(x.data)};
		}


		// -----------------
		// Bitwise operators
		// -----------------
    	friend Packed operator~(const Packed& value) requires std::integral<storage_type> {
        	return Packed{~value.data};
        }

    	friend Packed operator&(const Packed& lhs, const Packed& rhs) requires std::integral<storage_type> {
        	return Packed{lhs.data & rhs.data};
        }

    	friend Packed operator|(const Packed& lhs, const Packed& rhs) requires std::integral<storage_type> {
        	return Packed{lhs.data | rhs.data};
        }

    	friend Packed operator^(const Packed& lhs, const Packed& rhs) requires std::integral<storage_type> {
        	return Packed{lhs.data ^ rhs.data};
        }

    	Packed& operator&=(const Packed& rhs) requires std::integral<storage_type> {
        	data &= rhs.data;
        	return *this;
        }

    	Packed& operator|=(const Packed& rhs) requires std::integral<storage_type> {
        	data |= rhs.data;
        	return *this;
        }

    	Packed& operator^=(const Packed& rhs) requires std::integral<storage_type> {
        	data ^= rhs.data;
        	return *this;
        }


    	// ----------
    	// Reductions
    	// ----------
    	[[nodiscard]] storage_type reduce_add() const {
        	return static_cast<storage_type>(stdx::reduce(data));
        }

    	[[nodiscard]] storage_type reduce_min() const noexcept {
        	return stdx::reduce(data, [](const auto& a, const auto& b) {
				if constexpr (requires { stdx::min(a, b); })
					return stdx::min(a, b);
				else
					return std::min(a, b);
			});
        }

    	[[nodiscard]] storage_type reduce_max() const noexcept {
        	return stdx::reduce(data, [](const auto& a, const auto& b) {
				if constexpr (requires { stdx::max(a, b); })
					return stdx::max(a, b);
				else
					return std::max(a, b);
			});
        }

    	// ------------------------
    	// Debugging and inspection
    	// ------------------------
    	[[nodiscard]] std::array<value_type, size()> to_array() const {
        	alignas(stdx::memory_alignment_v<native_type, storage_type>)
			std::array<storage_type, size()> storage;

        	data.copy_to(storage.data(), stdx::vector_aligned);

        	std::array<value_type, size()> result;

        	for (size_t i = 0; i < size(); ++i)
        		result[i] = store_scalar<value_type>(storage[i]);

        	return result;
        }

    	[[nodiscard]] std::string to_string() const {
        	alignas(stdx::memory_alignment_v<native_type, storage_type>)
			std::array<storage_type, size()> storage;

        	data.copy_to(storage.data(), stdx::vector_aligned);

        	std::stringstream ss;
        	ss << "[";

        	for (size_t i = 0; i < size(); ++i) {
        		ss << storage[i];

        		if (i < size() - 1)
        			ss << ", ";
        	}

        	ss << "]";
        	return ss.str();
        }


    private:
        template<typename U>
        [[nodiscard]] static constexpr storage_type load_scalar(U value) noexcept {
            if constexpr (std::is_enum_v<U>)
                return static_cast<storage_type>(std::to_underlying(value));
            else
                return static_cast<storage_type>(value);
        }

        template<typename U>
        [[nodiscard]] static constexpr U store_scalar(storage_type value) noexcept {
            return static_cast<U>(value);
        }

    	template<typename Compare>
[[nodiscard]] static mask_type compare_impl(
	const native_type& lhs,
	const native_type& rhs,
	Compare compare
) {
#if defined(__clang__) && defined(__GLIBCXX__)
        	typename mask_type::native_type result(false);

        	for (size_t i = 0; i < size(); ++i)
        		result[i] = compare(lhs[i], rhs[i]);

        	return mask_type{std::move(result)};
#else
        	return mask_type{compare(lhs, rhs)};
#endif
        }
    };


    static_assert(april::simd::IsSimdType<Packed<double>>);
    static_assert(april::simd::IsSimdType<Packed<float>>);
    static_assert(april::simd::IsSimdType<Packed<size_t>>);
    static_assert(april::simd::IsSimdType<Packed<int>>);

    static_assert(april::simd::IsSimdMask<Mask<double>>);
    static_assert(april::simd::IsSimdMask<Mask<float>>);
    static_assert(april::simd::IsSimdMask<Mask<size_t>>);
    static_assert(april::simd::IsSimdMask<Mask<int>>);

    static_assert(april::simd::IsSimdType<Packed<float, 4>>);
    static_assert(april::simd::IsSimdType<Packed<double, 2>>);
    static_assert(april::simd::IsSimdMask<Mask<float, 4>>);
    static_assert(april::simd::IsSimdMask<Mask<double, 2>>);
}
