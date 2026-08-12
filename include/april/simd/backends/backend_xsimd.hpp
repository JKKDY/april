#pragma once
#include <array>
#include <cstdint>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility>
#include <cmath>
#include <xsimd/xsimd.hpp>


#include "april/simd/packed_concept.hpp"

namespace april::simd::internal::xsimd {


	template<typename T, size_t Width>
	consteval bool has_sized_batch() {
		if constexpr (Width == 0)
			return true;
		else
			return !std::is_void_v<::xsimd::make_sized_batch_t<T, Width>>;
	}

	template<typename T, size_t Width>
	struct native_storage {
		using storage_type = std::remove_cv_t<T>;

		using int16_type = std::conditional_t<
			std::is_signed_v<storage_type>,
			std::int16_t,
			std::uint16_t
		>;

		using int32_type = std::conditional_t<
			std::is_signed_v<storage_type>,
			std::int32_t,
			std::uint32_t
		>;

		using int64_type = std::conditional_t<
			std::is_signed_v<storage_type>,
			std::int64_t,
			std::uint64_t
		>;

		using type = std::conditional_t<
			has_sized_batch<storage_type, Width>(),
			storage_type,
			std::conditional_t<
				std::integral<storage_type> &&
				sizeof(storage_type) < sizeof(int16_type) &&
				has_sized_batch<int16_type, Width>(),
				int16_type,
				std::conditional_t<
					std::integral<storage_type> &&
					sizeof(storage_type) < sizeof(int32_type) &&
					has_sized_batch<int32_type, Width>(),
					int32_type,
					std::conditional_t<
						std::integral<storage_type> &&
						sizeof(storage_type) < sizeof(int64_type) &&
						has_sized_batch<int64_type, Width>(),
						int64_type,
						void
					>
				>
			>
		>;
	};

	template<typename T, size_t Width>
	using native_storage_t = typename native_storage<T, Width>::type;


    template<typename T, size_t Width> struct Packed;


    template<typename T, size_t Width = 0>
    struct Mask {
        using value_type = std::remove_cv_t<T>;
        using native_type = typename Packed<T, Width>::native_type::batch_bool_type;

        Mask() = default;

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

        [[nodiscard]] static constexpr size_t size() {
            return native_type::size;
        }

        // Implicit conversion between masks with the same logical width.
        template<typename U>
        requires (!std::same_as<T, U> && (size() == Mask<U, Width>::size()))
        Mask(const Mask<U, Width>& other)
        : data(convert_mask(other))
            {}

        template<typename U>
        requires (!std::same_as<T, U> && (size() == Mask<U, Width>::size()))
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


        // Loads
        static Mask load(const bool* ptr) {
	        return load_unaligned(ptr);
        }

        static Mask load_aligned(const bool* ptr) {
	        return Mask{native_type::load_aligned(ptr)};
        }

        static Mask load_unaligned(const bool* ptr) {
	        return Mask{native_type::load_unaligned(ptr)};
        }

        // Stores
        void store(bool* ptr) const {
	        store_unaligned(ptr);
        }

        void store_aligned(bool* ptr) const {
	        data.store_aligned(ptr);
        }

        void store_unaligned(bool* ptr) const {
	        data.store_unaligned(ptr);
        }

        // Logical reductions
        [[nodiscard]] static bool all(const Mask& mask) {
	        return ::xsimd::all(mask.data);
        }

        [[nodiscard]] static bool any(const Mask& mask) {
	        return ::xsimd::any(mask.data);
        }

        [[nodiscard]] static bool none(const Mask& mask) {
	        return !::xsimd::any(mask.data);
        }

        // Logical operators
        friend Mask operator!(const Mask& mask) {
	        return Mask{!mask.data};
        }

        template<typename U>
        friend Mask operator&&(const Mask& lhs, const Mask<U, Width>& rhs) {
	        return Mask{lhs.data && convert_mask(rhs)};
        }

        template<typename U>
        friend Mask operator||(const Mask& lhs, const Mask<U, Width>& rhs) {
	        return Mask{lhs.data || convert_mask(rhs)};
        }

        // Bitwise operators
        friend Mask operator~(const Mask& mask) {
	        return Mask{~mask.data};
        }

        template<typename U>
        friend Mask operator&(const Mask& lhs, const Mask<U, Width>& rhs) {
	        return Mask{lhs.data & convert_mask(rhs)};
        }

        template<typename U>
        friend Mask operator|(const Mask& lhs, const Mask<U, Width>& rhs) {
	        return Mask{lhs.data | convert_mask(rhs)};
        }

        template<typename U>
        friend Mask operator^(const Mask& lhs, const Mask<U, Width>& rhs) {
	        return Mask{lhs.data ^ convert_mask(rhs)};
        }

        // Lane-wise equality
        template<typename U>
        friend Mask operator==(const Mask& lhs, const Mask<U, Width>& rhs) {
	        return Mask{lhs.data == convert_mask(rhs)};
        }

        template<typename U>
        friend Mask operator!=(const Mask& lhs, const Mask<U, Width>& rhs) {
	        return Mask{lhs.data != convert_mask(rhs)};
        }

        // Compound operators
        template<typename U>
        Mask& operator&=(const Mask<U, Width>& rhs) {
	        data = data & convert_mask(rhs);
	        return *this;
        }

        template<typename U>
        Mask& operator|=(const Mask<U, Width>& rhs) {
	        data = data | convert_mask(rhs);
	        return *this;
        }

        template<typename U>
        Mask& operator^=(const Mask<U, Width>& rhs) {
	        data = data ^ convert_mask(rhs);
	        return *this;
        }

        // Rotations
        template<unsigned K = 1>
        void rotate_right() {
            rotate<K>();
        }

        template<unsigned K = 1>
        void rotate_left() {
            rotate<(size() - (K % size())) % size()>();
        }

        // EXPORTS / DEBUGGING
        [[nodiscard]] std::uint64_t to_bitmask() const {
            static_assert(size() <= 64, "Mask bit export supports at most 64 lanes");
            return static_cast<std::uint64_t>(data.mask());
        }

        static Mask from_bitmask(const std::uint64_t bits) {
            static_assert(size() <= 64, "Mask bit import supports at most 64 lanes");
            return Mask{native_type::from_mask(bits)};
        }

        [[nodiscard]] std::array<bool, size()> to_array() const {
            alignas(alignof(native_type)) std::array<bool, size()> result;
            store_aligned(result.data());
            return result;
        }

        [[nodiscard]] std::string to_string() const {
            auto arr = to_array();
            std::stringstream ss;
            ss << "[";
            for (size_t i = 0; i < size(); ++i) {
                ss << (arr[i] ? "true" : "false");
                if (i < size() - 1) ss << ", ";
            }
            ss << "]";
            return ss.str();
        }

        native_type data;
    private:
        template<typename U>
        [[nodiscard]] static native_type convert_mask(
            const Mask<U, Width>& other
        ) {
            static_assert(
                size() == Mask<U, Width>::size(),
                "Cannot convert masks with different logical widths"
            );

            return native_type::from_mask(other.to_bitmask());
        }

        template<unsigned K>
        void rotate() {
            constexpr unsigned Shift = K % size();

            if constexpr (Shift == 0) {} else if constexpr (
                sizeof(typename native_type::register_type) <
                sizeof(typename native_type::batch_type::register_type)
            ) {
                // AVX-512 stores masks as compact predicate bits.
                constexpr uint64_t ValidBits = [] {
                    if constexpr (size() == 64) return ~uint64_t{0};
                    else return (uint64_t{1} << size()) - 1;
                }();

                const uint64_t bits = data.mask();

                // Right lane rotation corresponds to a left rotation of mask bits.
                data = native_type::from_mask(
                    ((bits << Shift) | (bits >> (size() - Shift))) & ValidBits
                );
            } else {
                // SSE/AVX masks occupy full SIMD lanes, where native rotation is cheaper.
                using Batch = native_type::batch_type;
                using BatchValue = Batch::value_type;

                data = ::xsimd::rotate_right<Shift>(Batch(data)) != Batch(BatchValue{0});
            }
        }
    };




    template<typename T, size_t Width>
	struct Packed {
    	using value_type = std::remove_cv_t<T>;
    	using storage_type = packed_storage_t<value_type>;
    	using native_storage_type = native_storage_t<storage_type, Width>;

    	static_assert(
			!std::is_void_v<native_storage_type>,
			"xsimd cannot represent this type at the requested SIMD width"
		);

    	using native_type = std::conditional_t<
			Width == 0,
			::xsimd::batch<native_storage_type>,
			::xsimd::make_sized_batch_t<native_storage_type, Width>
		>;

    	using mask_type = Mask<storage_type, Width>;
    	using packed_type = Packed;

    	static_assert(
			Width == 0 || native_type::size == Width,
			"Native SIMD width does not match logical SIMD width"
		);

        Packed() = default;

        Packed(value_type value)
            : data(native_type{load_scalar(value)})
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

        Packed& operator=(value_type value) {
            data = native_type{load_scalar(value)};
            return *this;
        }

        [[nodiscard]] static constexpr size_t size() {
            return native_type::size;
        }

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
        requires (IsPackableValue<std::remove_cv_t<PtrT>> && (sizeof(PtrT) <= sizeof(storage_type)))
        static Packed load(const PtrT* ptr) {
	        return load_unaligned(ptr);
        }

        template<typename PtrT>
        requires (IsPackableValue<std::remove_cv_t<PtrT>> && (sizeof(PtrT) <= sizeof(storage_type)))
        static Packed load_aligned(const PtrT* ptr) {
	        if constexpr (std::is_enum_v<std::remove_cv_t<PtrT>>) {
		        alignas(alignof(native_type)) storage_type temp[size()];

		        for (size_t i = 0; i < size(); ++i)
			        temp[i] = load_scalar(ptr[i]);

		        return Packed{native_type::load_aligned(temp)};
	        }
	        else {
		        return Packed{native_type::load_aligned(ptr)};
	        }
        }

        template<typename PtrT>
        requires (IsPackableValue<std::remove_cv_t<PtrT>> && (sizeof(PtrT) <= sizeof(storage_type)))
        static Packed load_unaligned(const PtrT* ptr) {
	        if constexpr (std::is_enum_v<std::remove_cv_t<PtrT>>) {
		        alignas(alignof(native_type)) storage_type temp[size()];

		        for (size_t i = 0; i < size(); ++i)
			        temp[i] = load_scalar(ptr[i]);

		        return Packed{native_type::load_aligned(temp)};
	        }
	        else {
		        return Packed{native_type::load_unaligned(ptr)};
	        }
        }


        // -----------------
        // Contiguous stores
        // -----------------
        template<typename PtrT>
        requires (
	        !std::is_const_v<PtrT> &&
	        (IsPackableValue<std::remove_cv_t<PtrT>> && (sizeof(PtrT) <= sizeof(storage_type)))
        )
        void store(PtrT* ptr) const {
	        store_unaligned(ptr);
        }

        template<typename PtrT>
        requires (
            !std::is_const_v<PtrT> &&
            (IsPackableValue<std::remove_cv_t<PtrT>> && (sizeof(PtrT) <= sizeof(storage_type)))
        )
        void store_aligned(PtrT* ptr) const {
	        if constexpr (std::is_enum_v<std::remove_cv_t<PtrT>>) {
		        alignas(alignof(native_type)) storage_type temp[size()];
		        data.store_aligned(temp);

		        for (size_t i = 0; i < size(); ++i)
			        ptr[i] = store_scalar<PtrT>(temp[i]);
	        }
	        else {
		        data.store_aligned(ptr);
	        }
        }

        template<typename PtrT>
        requires (
            !std::is_const_v<PtrT> &&
            (IsPackableValue<std::remove_cv_t<PtrT>> && (sizeof(PtrT) <= sizeof(storage_type)))
        )
        void store_unaligned(PtrT* ptr) const {
	        if constexpr (std::is_enum_v<std::remove_cv_t<PtrT>>) {
		        alignas(alignof(native_type)) storage_type temp[size()];
		        data.store_aligned(temp);

		        for (size_t i = 0; i < size(); ++i)
			        ptr[i] = store_scalar<PtrT>(temp[i]);
	        }
	        else {
		        data.store_unaligned(ptr);
	        }
        }


    	// ------------
    	// Gather Loads
    	// ------------
    	// Compile-time byte stride
    	template<std::ptrdiff_t ByteStride, typename PtrT>
		requires(IsPackableValue<std::remove_cv_t<PtrT>> && (sizeof(PtrT) <= sizeof(storage_type)))
		static Packed gather_strided(const PtrT* ptr) {
        	static_assert(ByteStride > 0, "Byte stride must be positive");

        	// Contiguous case.
        	if constexpr (ByteStride == sizeof(PtrT)) {
        		return load(ptr);
        	}
        	// xsimd gather uses element indices, so a native gather is possible
        	// whenever the byte stride is exactly representable in PtrT elements.
        	else if constexpr (
				!std::is_enum_v<std::remove_cv_t<PtrT>> &&
				ByteStride % sizeof(PtrT) == 0
			) {
        		using arch_type = native_type::arch_type;
        		using index_value_type = ::xsimd::as_integer_t<storage_type>;
        		using index_type = ::xsimd::batch<index_value_type, arch_type>;

        		static_assert(index_type::size == size());

        		constexpr std::ptrdiff_t ElementStride =
					ByteStride / sizeof(PtrT);

        		constexpr auto indices =
					::xsimd::make_batch_constant<
						index_value_type,
						StridedIndices<ElementStride>,
						arch_type
					>();

        		return Packed{
        			native_type::gather(ptr, indices.as_batch())
				};
			}
        	// Enum memory, or a stride which cannot be expressed in PtrT elements.
        	else {
        		alignas(alignof(native_type))
				std::array<storage_type, size()> values;

        		const auto* bytes =
					reinterpret_cast<const std::byte*>(ptr);

        		for (size_t i = 0; i < size(); ++i) {
        			PtrT value;

        			std::memcpy(
						&value,
						bytes + i * ByteStride,
						sizeof(PtrT)
					);

        			values[i] = load_scalar(value);
        		}

        		return Packed{
        			native_type::load_aligned(values.data())
				};
        	}
        }

    	// Runtime byte stride
		template<typename PtrT>
		requires (
			(IsPackableValue<std::remove_cvref_t<PtrT>>) &&
			(sizeof(PtrT) <= sizeof(storage_type))
		)
		static Packed gather_strided(
			const PtrT* ptr,
			const std::ptrdiff_t byte_stride
		) {
			// Contiguous case.
			if (byte_stride == static_cast<std::ptrdiff_t>(sizeof(PtrT)))
				return load(ptr);

			// Native xsimd gather when the stride is representable
			// exactly in PtrT elements.
			if constexpr (!std::is_enum_v<std::remove_cv_t<PtrT>>) {
				if (byte_stride > 0 && byte_stride % sizeof(PtrT) == 0) {
					using arch_type = native_type::arch_type;
					using index_value_type = ::xsimd::as_integer_t<storage_type>;
					using index_type = ::xsimd::batch<index_value_type, arch_type>;

					static_assert(index_type::size == size());

					const std::ptrdiff_t element_stride =
						byte_stride / sizeof(PtrT);

					const std::ptrdiff_t max_index =
						element_stride * static_cast<std::ptrdiff_t>(size() - 1);

					if (max_index <= static_cast<std::ptrdiff_t>(
						std::numeric_limits<index_value_type>::max()
					)) {
						alignas(alignof(index_type))
						std::array<index_value_type, size()> indices;

						for (size_t i = 0; i < size(); ++i) {
							indices[i] = static_cast<index_value_type>(
								static_cast<std::ptrdiff_t>(i) * element_stride
							);
						}

						const index_type index_batch =
							index_type::load_aligned(indices.data());

						return Packed{
							native_type::gather(ptr, index_batch)
						};
					}
				}
			}

			// Enum memory, non-element-aligned stride, or an index range
			// too large for xsimd's index batch.
			alignas(alignof(native_type))
			std::array<storage_type, size()> values;

			const auto* bytes =
				reinterpret_cast<const std::byte*>(ptr);

			for (size_t i = 0; i < size(); ++i) {
				PtrT value;

				std::memcpy(
					&value,
					bytes + static_cast<std::ptrdiff_t>(i) * byte_stride,
					sizeof(PtrT)
				);

				values[i] = load_scalar(value);
			}

			return Packed{
				native_type::load_aligned(values.data())
			};
		}

    	// Arbitrary byte offsets
    	template<typename PtrT, size_t N>
		requires (
			(IsPackableValue<std::remove_cvref_t<PtrT>>) &&
			(sizeof(PtrT) <= sizeof(storage_type)) &&
			(N == size())
		)
		static Packed gather(
			const PtrT* ptr,
			const ByteOffsets<N>& offsets
		) {
        	if constexpr (!std::is_enum_v<std::remove_cv_t<PtrT>>) {
        		using arch_type = native_type::arch_type;
        		using index_value_type = ::xsimd::as_integer_t<storage_type>;
        		using index_type = ::xsimd::batch<index_value_type, arch_type>;

        		static_assert(index_type::size == size());

        		alignas(alignof(index_type))
				std::array<index_value_type, size()> indices;

        		bool native_compatible = true;

        		for (size_t i = 0; i < size(); ++i) {
        			const std::ptrdiff_t offset = offsets.values[i];

        			if (offset % static_cast<std::ptrdiff_t>(sizeof(PtrT)) != 0) {
        				native_compatible = false;
        				break;
        			}

        			const std::ptrdiff_t index =
						offset / static_cast<std::ptrdiff_t>(sizeof(PtrT));

        			if (
						index < static_cast<std::ptrdiff_t>(
							std::numeric_limits<index_value_type>::min()
						) ||
						index > static_cast<std::ptrdiff_t>(
							std::numeric_limits<index_value_type>::max()
						)
					) {
        				native_compatible = false;
        				break;
					}

        			indices[i] = static_cast<index_value_type>(index);
        		}

        		if (native_compatible) {
        			const index_type index_batch =
						index_type::load_aligned(indices.data());

        			return Packed{
        				native_type::gather(ptr, index_batch)
					};
        		}
        	}

        	// Enum memory, non-element-aligned offsets, or offsets too large
        	// for the xsimd index representation.
        	alignas(alignof(native_type))
			std::array<storage_type, size()> values;

        	const auto* bytes =
				reinterpret_cast<const std::byte*>(ptr);

        	for (size_t i = 0; i < size(); ++i) {
        		PtrT value;

        		std::memcpy(
					&value,
					bytes + offsets.values[i],
					sizeof(PtrT)
				);

        		values[i] = load_scalar(value);
        	}

        	return Packed{
        		native_type::load_aligned(values.data())
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
				std::is_arithmetic_v<source_type> ||
				std::is_enum_v<source_type>,
				"Gather pointers must point to arithmetic or enum values"
			);

        	static_assert(
				sizeof(source_type) <= sizeof(storage_type),
				"Gather source type is wider than packed storage type"
			);

        	alignas(alignof(native_type))
			std::array<storage_type, size()> values;

        	for (size_t i = 0; i < size(); ++i)
        		values[i] = load_scalar(*pointers[i]);

        	return Packed{
        		native_type::load_aligned(values.data())
			};
        }


    	// --------------
    	// Scatter Stores
    	// --------------
    	// Compile-time byte stride
    	template<std::ptrdiff_t ByteStride, typename PtrT>
		requires(
			!std::is_const_v<PtrT> &&
			IsPackableValue<std::remove_cv_t<PtrT>> &&
			(sizeof(PtrT) <= sizeof(storage_type))
		)
		void scatter_strided(PtrT* ptr) const {
        	static_assert(ByteStride > 0, "Byte stride must be positive");

        	// Contiguous case.
        	if constexpr (ByteStride == sizeof(PtrT)) {
        		store(ptr);
        	}
        	// Native xsimd scatter when the byte stride can be represented
        	// exactly as a PtrT element stride.
        	else if constexpr (
				!std::is_enum_v<std::remove_cv_t<PtrT>> &&
				ByteStride % sizeof(PtrT) == 0
			) {
        		using arch_type = native_type::arch_type;
        		using index_value_type = ::xsimd::as_integer_t<storage_type>;

        		constexpr std::ptrdiff_t ElementStride =
					ByteStride / sizeof(PtrT);

        		constexpr auto indices =
					::xsimd::make_batch_constant<
						index_value_type,
						StridedIndices<ElementStride>,
						arch_type
					>();

        		data.scatter(ptr, indices.as_batch());
			}
        	// Enum memory, or a stride not representable in PtrT elements.
        	else {
        		alignas(alignof(native_type))
				std::array<storage_type, size()> values;

        		data.store_aligned(values.data());

        		auto* bytes =
					reinterpret_cast<std::byte*>(ptr);

        		for (size_t i = 0; i < size(); ++i) {
        			const PtrT value =
						store_scalar<PtrT>(values[i]);

        			std::memcpy(
						bytes + i * ByteStride,
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
			(std::is_arithmetic_v<std::remove_cv_t<PtrT>> || std::is_enum_v<std::remove_cv_t<PtrT>>) &&
			(sizeof(PtrT) <= sizeof(storage_type))
		)
		void scatter_strided(
			PtrT* ptr,
			const std::ptrdiff_t byte_stride
		) const {
			// Contiguous case.
			if (byte_stride == static_cast<std::ptrdiff_t>(sizeof(PtrT))) {
				store(ptr);
				return;
			}

			// Native xsimd scatter when the stride is representable
			// exactly in PtrT elements.
			if constexpr (!std::is_enum_v<std::remove_cv_t<PtrT>>) {
				if (byte_stride > 0 && byte_stride % sizeof(PtrT) == 0) {
					using arch_type = native_type::arch_type;
					using index_value_type = ::xsimd::as_integer_t<storage_type>;
					using index_type = ::xsimd::batch<index_value_type, arch_type>;

					static_assert(index_type::size == size());

					const std::ptrdiff_t element_stride =
						byte_stride / static_cast<std::ptrdiff_t>(sizeof(PtrT));

					const std::ptrdiff_t max_index =
						element_stride * static_cast<std::ptrdiff_t>(size() - 1);

					if (max_index <= static_cast<std::ptrdiff_t>(
						std::numeric_limits<index_value_type>::max()
					)) {
						alignas(alignof(index_type))
						std::array<index_value_type, size()> indices;

						for (size_t i = 0; i < size(); ++i) {
							indices[i] = static_cast<index_value_type>(
								static_cast<std::ptrdiff_t>(i) * element_stride
							);
						}

						const index_type index_batch =
							index_type::load_aligned(indices.data());

						data.scatter(ptr, index_batch);
						return;
					}
				}
			}

			// Enum memory, non-element-aligned stride, or an index range
			// too large for xsimd's index batch.
			alignas(alignof(native_type))
			std::array<storage_type, size()> values;

			data.store_aligned(values.data());

			auto* bytes =
				reinterpret_cast<std::byte*>(ptr);

			for (size_t i = 0; i < size(); ++i) {
				const PtrT value =
					store_scalar<PtrT>(values[i]);

				std::memcpy(
					bytes + static_cast<std::ptrdiff_t>(i) * byte_stride,
					&value,
					sizeof(PtrT)
				);
			}
		}

    	// Arbitrary byte offsets
    	template<typename PtrT, size_t N>
		requires (
			!std::is_const_v<PtrT> &&
			(std::is_arithmetic_v<std::remove_cv_t<PtrT>> || std::is_enum_v<std::remove_cv_t<PtrT>>) &&
			(sizeof(PtrT) <= sizeof(storage_type)) &&
			(N == size())
		)
		void scatter(
			PtrT* ptr,
			const ByteOffsets<N>& offsets
		) const {
        	if constexpr (!std::is_enum_v<std::remove_cv_t<PtrT>>) {
        		using arch_type = native_type::arch_type;
        		using index_value_type = ::xsimd::as_integer_t<storage_type>;
        		using index_type = ::xsimd::batch<index_value_type, arch_type>;

        		static_assert(index_type::size == size());

        		alignas(alignof(index_type))
				std::array<index_value_type, size()> indices;

        		bool native_compatible = true;

        		for (size_t i = 0; i < size(); ++i) {
        			const std::ptrdiff_t offset = offsets.values[i];

        			if (offset % static_cast<std::ptrdiff_t>(sizeof(PtrT)) != 0) {
        				native_compatible = false;
        				break;
        			}

        			const std::ptrdiff_t index =
						offset / static_cast<std::ptrdiff_t>(sizeof(PtrT));

        			if (
						index < static_cast<std::ptrdiff_t>(
							std::numeric_limits<index_value_type>::min()
						) ||
						index > static_cast<std::ptrdiff_t>(
							std::numeric_limits<index_value_type>::max()
						)
					) {
        				native_compatible = false;
        				break;
					}

        			indices[i] = static_cast<index_value_type>(index);
        		}

        		if (native_compatible) {
        			const index_type index_batch =
						index_type::load_aligned(indices.data());

        			data.scatter(ptr, index_batch);
        			return;
        		}
        	}

        	// Enum memory, non-element-aligned offsets, or offsets too large
        	// for the xsimd index representation.
        	alignas(alignof(native_type))
			std::array<storage_type, size()> values;

        	data.store_aligned(values.data());

        	auto* bytes =
				reinterpret_cast<std::byte*>(ptr);

        	for (size_t i = 0; i < size(); ++i) {
        		const PtrT value =
					store_scalar<PtrT>(values[i]);

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
				std::is_arithmetic_v<destination_type> ||
				std::is_enum_v<destination_type>,
				"Scatter pointers must point to arithmetic or enum values"
			);

        	static_assert(
				sizeof(destination_type) <= sizeof(storage_type),
				"Scatter destination type is wider than packed storage type"
			);

        	alignas(alignof(native_type))
			std::array<storage_type, size()> values;

        	data.store_aligned(values.data());

        	for (size_t i = 0; i < size(); ++i) {
        		*pointers[i] =
					store_scalar<destination_type>(values[i]);
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

	        using arch_type = native_type::arch_type;

	        constexpr auto indices =
		        ::xsimd::make_batch_constant<
			        native_storage_type,
			        PermuteIndices<Indices...>,
			        arch_type
		        >();

	        return Packed{::xsimd::swizzle(data, indices)};
        }

        template<unsigned K = 1>
        [[nodiscard]] Packed rotate_left() const {
	        constexpr unsigned Shift = K % size();
	        return Packed{::xsimd::rotate_left<Shift>(data)};
        }

        template<unsigned K = 1>
        [[nodiscard]] Packed rotate_right() const {
	        constexpr unsigned Shift = K % size();
	        return Packed{::xsimd::rotate_right<Shift>(data)};
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
	        return Packed{lhs.data + native_type{rhs}};
        }

        friend Packed operator+(storage_type lhs, const Packed& rhs) {
	        return Packed{native_type{lhs} + rhs.data};
        }

        friend Packed operator-(const Packed& lhs, storage_type rhs) {
	        return Packed{lhs.data - native_type{rhs}};
        }

        friend Packed operator-(storage_type lhs, const Packed& rhs) {
	        return Packed{native_type{lhs} - rhs.data};
        }

        friend Packed operator*(const Packed& lhs, storage_type rhs) {
	        return Packed{lhs.data * native_type{rhs}};
        }

        friend Packed operator*(storage_type lhs, const Packed& rhs) {
	        return Packed{native_type{lhs} * rhs.data};
        }

        friend Packed operator/(const Packed& lhs, storage_type rhs) {
	        return Packed{lhs.data / native_type{rhs}};
        }

        friend Packed operator/(storage_type lhs, const Packed& rhs) {
	        return Packed{native_type{lhs} / rhs.data};
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
	        data += native_type{rhs};
	        return *this;
        }

        Packed& operator-=(storage_type rhs) {
	        data -= native_type{rhs};
	        return *this;
        }

        Packed& operator*=(storage_type rhs) {
	        data *= native_type{rhs};
	        return *this;
        }

        Packed& operator/=(storage_type rhs) {
	        data /= native_type{rhs};
	        return *this;
        }


        // -----------
        // Comparisons
        // -----------
        friend mask_type operator==(const Packed& lhs, const Packed& rhs) {
	        return mask_type{lhs.data == rhs.data};
        }

        friend mask_type operator!=(const Packed& lhs, const Packed& rhs) {
	        return mask_type{lhs.data != rhs.data};
        }

        friend mask_type operator<(const Packed& lhs, const Packed& rhs) {
	        return mask_type{lhs.data < rhs.data};
        }

        friend mask_type operator<=(const Packed& lhs, const Packed& rhs) {
	        return mask_type{lhs.data <= rhs.data};
        }

        friend mask_type operator>(const Packed& lhs, const Packed& rhs) {
	        return mask_type{lhs.data > rhs.data};
        }

        friend mask_type operator>=(const Packed& lhs, const Packed& rhs) {
	        return mask_type{lhs.data >= rhs.data};
        }

        friend mask_type operator==(const Packed& lhs, storage_type rhs) {
	        return mask_type{lhs.data == native_type{rhs}};
        }

        friend mask_type operator==(storage_type lhs, const Packed& rhs) {
	        return mask_type{native_type{lhs} == rhs.data};
        }

        friend mask_type operator!=(const Packed& lhs, storage_type rhs) {
	        return mask_type{lhs.data != native_type{rhs}};
        }

        friend mask_type operator!=(storage_type lhs, const Packed& rhs) {
	        return mask_type{native_type{lhs} != rhs.data};
        }

        friend mask_type operator<(const Packed& lhs, storage_type rhs) {
	        return mask_type{lhs.data < native_type{rhs}};
        }

        friend mask_type operator<(storage_type lhs, const Packed& rhs) {
	        return mask_type{native_type{lhs} < rhs.data};
        }

        friend mask_type operator<=(const Packed& lhs, storage_type rhs) {
	        return mask_type{lhs.data <= native_type{rhs}};
        }

        friend mask_type operator<=(storage_type lhs, const Packed& rhs) {
	        return mask_type{native_type{lhs} <= rhs.data};
        }

        friend mask_type operator>(const Packed& lhs, storage_type rhs) {
	        return mask_type{lhs.data > native_type{rhs}};
        }

        friend mask_type operator>(storage_type lhs, const Packed& rhs) {
	        return mask_type{native_type{lhs} > rhs.data};
        }

        friend mask_type operator>=(const Packed& lhs, storage_type rhs) {
	        return mask_type{lhs.data >= native_type{rhs}};
        }

        friend mask_type operator>=(storage_type lhs, const Packed& rhs) {
	        return mask_type{native_type{lhs} >= rhs.data};
        }


        // ---------
        // Selection
        // ---------
        [[nodiscard]] static Packed select(
	        const mask_type& mask,
	        const Packed& true_value,
	        const Packed& false_value
        ) {
	        return Packed{
		        ::xsimd::select(mask.native(), true_value.data, false_value.data)
	        };
        }


        // ----------------
        // Basic operations
        // ----------------
        [[nodiscard]] static Packed abs(const Packed& x) {
	        if constexpr (std::is_unsigned_v<storage_type>)
		        return x;
	        else
		        return Packed{::xsimd::abs(x.data)};
        }

        [[nodiscard]] static Packed min(const Packed& a, const Packed& b) {
	        return Packed{::xsimd::min(a.data, b.data)};
        }

        [[nodiscard]] static Packed max(const Packed& a, const Packed& b) {
	        return Packed{::xsimd::max(a.data, b.data)};
        }

        [[nodiscard]] static Packed clamp(
	        const Packed& x,
	        const Packed& lower,
	        const Packed& upper
        ) {
	        return Packed{::xsimd::clip(x.data, lower.data, upper.data)};
        }


        // ----------------
        // Roots and powers
        // ----------------
        [[nodiscard]] static Packed sqrt(const Packed& x)
	        requires std::floating_point<storage_type>
        {
	        return Packed{::xsimd::sqrt(x.data)};
        }

        [[nodiscard]] static Packed rsqrt(const Packed& x)
	        requires std::floating_point<storage_type>
        {
        #if APRIL_FAST_MATH_ENABLED
	        return Packed{::xsimd::rsqrt(x.data)};
        #else
	        return Packed{
		        native_type{storage_type{1}} / ::xsimd::sqrt(x.data)
	        };
        #endif
        }

        [[nodiscard]] static Packed cbrt(const Packed& x)
	        requires std::floating_point<storage_type>
        {
	        return Packed{::xsimd::cbrt(x.data)};
        }

        [[nodiscard]] static Packed hypot(const Packed& x, const Packed& y)
	        requires std::floating_point<storage_type>
        {
	        return Packed{::xsimd::hypot(x.data, y.data)};
        }

        [[nodiscard]] static Packed pow(const Packed& x, const Packed& y)
	        requires std::floating_point<storage_type>
        {
	        return Packed{::xsimd::pow(x.data, y.data)};
        }


        // ---------------------------
        // Exponential and logarithmic
        // ---------------------------
    	[[nodiscard]] static Packed exp(const Packed& x) requires std::floating_point<storage_type> {
        	return Packed{::xsimd::exp(x.data)};
        }

    	[[nodiscard]] static Packed exp2(const Packed& x) requires std::floating_point<storage_type> {
        	return Packed{::xsimd::exp2(x.data)};
        }

    	[[nodiscard]] static Packed expm1(const Packed& x) requires std::floating_point<storage_type> {
        	return Packed{::xsimd::expm1(x.data)};
        }

    	[[nodiscard]] static Packed log(const Packed& x) requires std::floating_point<storage_type> {
        	return Packed{::xsimd::log(x.data)};
        }

    	[[nodiscard]] static Packed ln(const Packed& x) requires std::floating_point<storage_type> {
        	return log(x);
        }

    	[[nodiscard]] static Packed log2(const Packed& x) requires std::floating_point<storage_type> {
        	return Packed{::xsimd::log2(x.data)};
        }

    	[[nodiscard]] static Packed log10(const Packed& x) requires std::floating_point<storage_type> {
        	return Packed{::xsimd::log10(x.data)};
        }

    	[[nodiscard]] static Packed log1p(const Packed& x) requires std::floating_point<storage_type> {
        	return Packed{::xsimd::log1p(x.data)};
        }


        // -------------
        // Trigonometric
        // -------------
    	[[nodiscard]] static Packed sin(const Packed& x) requires std::floating_point<storage_type> {
        	return Packed{::xsimd::sin(x.data)};
        }

    	[[nodiscard]] static Packed cos(const Packed& x) requires std::floating_point<storage_type> {
        	return Packed{::xsimd::cos(x.data)};
        }

    	[[nodiscard]] static std::pair<Packed, Packed> sincos(const Packed& x) requires std::floating_point<storage_type> {
        	auto [sin_value, cos_value] = ::xsimd::sincos(x.data);
        	return {
        		Packed{sin_value},
				Packed{cos_value}
        	};
        }

    	[[nodiscard]] static Packed tan(const Packed& x) requires std::floating_point<storage_type> {
        	return Packed{::xsimd::tan(x.data)};
        }

    	[[nodiscard]] static Packed asin(const Packed& x) requires std::floating_point<storage_type> {
        	return Packed{::xsimd::asin(x.data)};
        }

    	[[nodiscard]] static Packed acos(const Packed& x) requires std::floating_point<storage_type> {
        	return Packed{::xsimd::acos(x.data)};
        }

    	[[nodiscard]] static Packed atan(const Packed& x) requires std::floating_point<storage_type> {
        	return Packed{::xsimd::atan(x.data)};
        }

    	[[nodiscard]] static Packed atan2(const Packed& y, const Packed& x) requires std::floating_point<storage_type> {
        	return Packed{::xsimd::atan2(y.data, x.data)};
        }


        // ----------
        // Hyperbolic
        // ----------
    	[[nodiscard]] static Packed sinh(const Packed& x) requires std::floating_point<storage_type> {
        	return Packed{::xsimd::sinh(x.data)};
        }

    	[[nodiscard]] static Packed cosh(const Packed& x) requires std::floating_point<storage_type> {
        	return Packed{::xsimd::cosh(x.data)};
        }

    	[[nodiscard]] static Packed tanh(const Packed& x) requires std::floating_point<storage_type> {
        	return Packed{::xsimd::tanh(x.data)};
        }

    	[[nodiscard]] static Packed asinh(const Packed& x) requires std::floating_point<storage_type> {
        	return Packed{::xsimd::asinh(x.data)};
        }

    	[[nodiscard]] static Packed acosh(const Packed& x) requires std::floating_point<storage_type> {
        	return Packed{::xsimd::acosh(x.data)};
        }

    	[[nodiscard]] static Packed atanh(const Packed& x) requires std::floating_point<storage_type> {
        	return Packed{::xsimd::atanh(x.data)};
        }


        // --------
        // Rounding
        // --------
        [[nodiscard]] static Packed floor(const Packed& x) {
	        if constexpr (std::integral<storage_type>)
		        return x;
	        else
		        return Packed{::xsimd::floor(x.data)};
        }

        [[nodiscard]] static Packed ceil(const Packed& x) {
	        if constexpr (std::integral<storage_type>)
		        return x;
	        else
		        return Packed{::xsimd::ceil(x.data)};
        }

        [[nodiscard]] static Packed round(const Packed& x) {
	        if constexpr (std::integral<storage_type>)
		        return x;
	        else
		        return Packed{::xsimd::round(x.data)};
        }

        [[nodiscard]] static Packed trunc(const Packed& x) {
	        if constexpr (std::integral<storage_type>)
		        return x;
	        else
		        return Packed{::xsimd::trunc(x.data)};
        }

        [[nodiscard]] static Packed nearbyint(const Packed& x) {
	        if constexpr (std::integral<storage_type>)
		        return x;
	        else
		        return Packed{::xsimd::nearbyint(x.data)};
        }


        // ------------------
        // Numeric operations
        // ------------------
        [[nodiscard]] static Packed fma(
	        const Packed& x,
	        const Packed& y,
	        const Packed& z
        ) {
	        return Packed{::xsimd::fma(x.data, y.data, z.data)};
        }

        [[nodiscard]] static Packed fmod(const Packed& x, const Packed& y)
	        requires std::floating_point<storage_type>
        {
	        return Packed{::xsimd::fmod(x.data, y.data)};
        }

        [[nodiscard]] static Packed remainder(const Packed& x, const Packed& y)
	        requires std::floating_point<storage_type>
        {
	        return Packed{::xsimd::remainder(x.data, y.data)};
        }

        [[nodiscard]] static Packed copysign(
	        const Packed& magnitude,
	        const Packed& sign
        )
	        requires std::floating_point<storage_type>
        {
	        return Packed{
		        ::xsimd::copysign(magnitude.data, sign.data)
	        };
        }


        // --------------
        // Classification
        // --------------
        [[nodiscard]] static mask_type isnan(const Packed& x) {
	        if constexpr (std::floating_point<storage_type>)
		        return mask_type{::xsimd::isnan(x.data)};
	        else
		        return mask_type{false};
        }

        [[nodiscard]] static mask_type isinf(const Packed& x) {
	        if constexpr (std::floating_point<storage_type>)
		        return mask_type{::xsimd::isinf(x.data)};
	        else
		        return mask_type{false};
        }

        [[nodiscard]] static mask_type isfinite(const Packed& x) {
	        if constexpr (std::floating_point<storage_type>)
		        return mask_type{::xsimd::isfinite(x.data)};
	        else
		        return mask_type{true};
        }

        [[nodiscard]] static mask_type signbit(const Packed& x) {
	        if constexpr (std::is_unsigned_v<storage_type>) {
		        return mask_type{false};
	        }
	        else if constexpr (std::integral<storage_type>) {
		        return mask_type{
			        x.data < native_type{storage_type{0}}
		        };
	        }
	        else {
		        using integer_type = std::conditional_t<
			        sizeof(storage_type) == 4,
			        std::uint32_t,
			        std::uint64_t
		        >;

		        const auto sign_bits =
			        ::xsimd::bitwise_cast<integer_type>(
				        ::xsimd::bitofsign(x.data)
			        );

		        const auto integer_mask =
			        sign_bits != decltype(sign_bits){integer_type{0}};

		        return mask_type{
			        ::xsimd::batch_bool_cast<storage_type>(integer_mask)
		        };
	        }
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
			return static_cast<storage_type>(::xsimd::reduce_add(data));
		}

		[[nodiscard]] storage_type reduce_min() const {
			return static_cast<storage_type>(::xsimd::reduce_min(data));
		}

		[[nodiscard]] storage_type reduce_max() const {
			return static_cast<storage_type>(::xsimd::reduce_max(data));
		}


		// ------------------------
		// Debugging and inspection
		// ------------------------
		[[nodiscard]] std::array<value_type, size()> to_array() const {
			alignas(alignof(native_type)) std::array<storage_type, size()> storage;
			data.store_aligned(storage.data());

			std::array<value_type, size()> result;

			for (size_t i = 0; i < size(); ++i)
				result[i] = store_scalar<value_type>(storage[i]);

			return result;
		}

		[[nodiscard]] std::string to_string() const {
			alignas(alignof(native_type)) std::array<storage_type, size()> storage;
			data.store_aligned(storage.data());

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



        native_type data;

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

    	template<std::ptrdiff_t ElementStride>
		struct StridedIndices {
        	using index_value_type =
				::xsimd::as_integer_t<native_storage_type>;

        	static constexpr index_value_type get(
				const unsigned i,
				const unsigned
			) {
        		return static_cast<index_value_type>(
					i * ElementStride
				);
        	}
        };

        template<size_t... Indices>
        struct PermuteIndices {
            static constexpr unsigned get(const unsigned i, const unsigned) {
                constexpr std::array<size_t, sizeof...(Indices)> indices{Indices...};

                if constexpr (sizeof...(Indices) == 1)
                    return static_cast<unsigned>(indices[0]);
                else
                    return static_cast<unsigned>(indices[i]);
            }
        };
    };



    static_assert(april::simd::IsSimdType<Packed<double, 0>>);
    static_assert(april::simd::IsSimdType<Packed<float, 0>>);
    static_assert(april::simd::IsSimdType<Packed<size_t, 0>>);
    static_assert(april::simd::IsSimdType<Packed<int, 0>>);

    static_assert(april::simd::IsSimdMask<Mask<double>>);
    static_assert(april::simd::IsSimdMask<Mask<float>>);
    static_assert(april::simd::IsSimdMask<Mask<size_t>>);
    static_assert(april::simd::IsSimdMask<Mask<int>>);

    static_assert(april::simd::IsSimdType<Packed<float, 4>>);
    static_assert(april::simd::IsSimdType<Packed<double, 2>>);
    static_assert(april::simd::IsSimdMask<Mask<float, 4>>);
    static_assert(april::simd::IsSimdMask<Mask<double, 2>>);
}


