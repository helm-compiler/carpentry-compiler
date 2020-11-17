#include "HLParser.h"

namespace HighLevelParser
{
	HLParser::HLParser()
	{

	}

	HLParser::HLParser(const std::string& file)
	{
		ParseFile(file);
	}

	bool HLParser::ParseFile(const std::string& file)
	{
		// TODO: parse a HL-HELM file to AST representation

		return true;
	}
}