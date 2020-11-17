#pragma once
#include "PreCompiled.h"
#include <QString>

#include <vector>
#include <Base/Type.h>
#include <App/DocumentObject.h>

namespace PartDesign {
	class AdditiveBox;
	class Pocket;
}

class Compiler
{
public:
	void CompileDown(std::vector<App::DocumentObject*>& sortedObjects);

	QString GetHLHELM() const;
	QString GetLLHELM() const;

	QString CompileBox(PartDesign::AdditiveBox* box, QString& varName);
	QString CompilePolycut(PartDesign::Pocket* polycut, QString& varName);

private:
	QString hlHELM, llHELM;
};

