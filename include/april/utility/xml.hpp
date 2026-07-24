#pragma once

#include <cstddef>
#include <ostream>
#include <string_view>
#include <type_traits>
#include <utility>

#include "april/utility/debug.hpp"


namespace april::utility {

	template<typename T>
	struct XMLAttribute {
		std::string_view name;
		T value;
	};


	template<typename T>
	[[nodiscard]] auto attribute(const std::string_view name, T && value) {
		return XMLAttribute<std::decay_t<T>> {
			name,
			std::forward<T>(value)
		};
	}


	class XMLWriter {
	public:
		explicit XMLWriter(std::ostream & output):
			out(output)
		{}


		void declaration() const {
			out << "<?xml version=\"1.0\"?>\n";
		}


		template<typename... A>
		void open(const std::string_view name, const A &... attributes) {
			write_indent();

			out << '<' << name;
			(write_attribute(attributes), ...);
			out << ">\n";

			depth++;
		}


		template<typename... A>
		void empty(const std::string_view name, const A &... attributes) {
			write_indent();

			out << '<' << name;
			(write_attribute(attributes), ...);
			out << "/>\n";
		}


		void close(const std::string_view name) {
			APRIL_ASSERT(depth > 0, "XML Writer: depth cannot be less zero or negative");

			depth--;

			write_indent();
			out << "</" << name << ">\n";
		}


		void text(const std::string_view value) const {
			write_indent();
			write_escaped(value);
			out << '\n';
		}


		template<typename T>
		void value(const T & value) {
			write_indent();
			out << value << '\n';
		}


		[[nodiscard]] std::ostream & stream() noexcept {
			return out;
		}


		[[nodiscard]] const std::ostream & stream() const noexcept {
			return out;
		}


	private:
		template<typename T>
		void write_attribute(const XMLAttribute<T> & attribute) {
			out
				<< ' '
				<< attribute.name
				<< "=\"";

			write_escaped_value(attribute.value);

			out << '"';
		}


		template<typename T>
		void write_escaped_value(const T & value) {
			if constexpr (std::is_convertible_v<T, std::string_view>) {
				write_escaped(std::string_view(value));
			}
			else {
				out << value;
			}
		}


		void write_escaped(const std::string_view value) const {
			for (const char c : value) {
				switch (c) {
					case '&':
						out << "&amp;";
						break;

					case '<':
						out << "&lt;";
						break;

					case '>':
						out << "&gt;";
						break;

					case '"':
						out << "&quot;";
						break;

					case '\'':
						out << "&apos;";
						break;

					default:
						out << c;
						break;
				}
			}
		}


		void write_indent() const {
			for (std::size_t i = 0; i < depth; i++) {
				out << '\t';
			}
		}


		std::ostream & out;
		std::size_t depth = 0;
	};

} // namespace april::io