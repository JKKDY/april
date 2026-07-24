#include <gtest/gtest.h>

#include <sstream>
#include <string>

#include "april/utility/xml.hpp"

#include "april/base/types.hpp"


using namespace april::utility;


TEST(XMLWriterTest, WritesDeclaration) {
	std::ostringstream out;
	const XMLWriter xml(out);

	xml.declaration();

	EXPECT_EQ(
		out.str(),
		"<?xml version=\"1.0\"?>\n"
	);
}


TEST(XMLWriterTest, WritesNestedElements) {
	std::ostringstream out;
	XMLWriter xml(out);

	xml.open("Root");
	xml.open("Child");
	xml.close("Child");
	xml.close("Root");

	EXPECT_EQ(
		out.str(),
		"<Root>\n"
		"\t<Child>\n"
		"\t</Child>\n"
		"</Root>\n"
	);
}


TEST(XMLWriterTest, WritesAttributes) {
	std::ostringstream out;
	XMLWriter xml(out);

	xml.empty(
		"DataArray",
		attribute("type", "Float64"),
		attribute("components", 3),
		attribute("format", "ascii")
	);

	EXPECT_EQ(
		out.str(),
		"<DataArray type=\"Float64\" components=\"3\" format=\"ascii\"/>\n"
	);
}


TEST(XMLWriterTest, EscapesAttributeValues) {
	std::ostringstream out;
	XMLWriter xml(out);

	xml.empty("Entry", attribute("name", "A&B<\"value\">'"));

	EXPECT_EQ(
		out.str(),
		"<Entry name=\"A&amp;B&lt;&quot;value&quot;&gt;&apos;\"/>\n"
	);
}


TEST(XMLWriterTest, EscapesText) {
	std::ostringstream out;
	XMLWriter xml(out);

	xml.open("Text");
	xml.text("A&B<\"value\">'");
	xml.close("Text");

	EXPECT_EQ(
		out.str(),
		"<Text>\n"
		"\tA&amp;B&lt;&quot;value&quot;&gt;&apos;\n"
		"</Text>\n"
	);
}


TEST(XMLWriterTest, WritesValues) {
	std::ostringstream out;
	XMLWriter xml(out);

	xml.open("Values");
	xml.value(42);
	xml.value(1.25);
	xml.close("Values");

	EXPECT_EQ(
		out.str(),
		"<Values>\n"
		"\t42\n"
		"\t1.25\n"
		"</Values>\n"
	);
}


TEST(XMLWriterTest, ProvidesAccessToUnderlyingStream) {
	std::ostringstream out;
	XMLWriter xml(out);

	xml.open("DataArray");

	xml.stream()
		<< '\t'
		<< "1 2 3\n";

	xml.close("DataArray");

	EXPECT_EQ(
		out.str(),
		"<DataArray>\n"
		"\t1 2 3\n"
		"</DataArray>\n"
	);
}


TEST(XMLWriterTest, WritesMinimalVTKDocument) {
	std::ostringstream out;
	XMLWriter xml(out);

	xml.declaration();

	xml.open(
		"VTKFile",
		attribute("type", "PolyData"),
		attribute("version", "0.1"),
		attribute("byte_order", "LittleEndian")
	);

	xml.open("PolyData");

	xml.empty(
		"Piece",
		attribute("NumberOfPoints", 0),
		attribute("NumberOfVerts", 0)
	);

	xml.close("PolyData");
	xml.close("VTKFile");

	EXPECT_EQ(
		out.str(),
		"<?xml version=\"1.0\"?>\n"
		"<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
		"\t<PolyData>\n"
		"\t\t<Piece NumberOfPoints=\"0\" NumberOfVerts=\"0\"/>\n"
		"\t</PolyData>\n"
		"</VTKFile>\n"
	);
}