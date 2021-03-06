/***************************************************************************
 *   Copyright (c) 2008 Jürgen Riegel (juergen.riegel@web.de)              *
 *                                                                         *
 *   This file is part of the FreeCAD CAx development system.              *
 *                                                                         *
 *   This library is free software; you can redistribute it and/or         *
 *   modify it under the terms of the GNU Library General Public           *
 *   License as published by the Free Software Foundation; either          *
 *   version 2 of the License, or (at your option) any later version.      *
 *                                                                         *
 *   This library  is distributed in the hope that it will be useful,      *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU Library General Public License for more details.                  *
 *                                                                         *
 *   You should have received a copy of the GNU Library General Public     *
 *   License along with this library; see the file COPYING.LIB. If not,    *
 *   write to the Free Software Foundation, Inc., 59 Temple Place,         *
 *   Suite 330, Boston, MA  02111-1307, USA                                *
 *                                                                         *
 ***************************************************************************/

#include "PreCompiled.h"
#ifndef _PreComp_
# include <Python.h>
#endif

#include <CXX/Extensions.hxx>
#include <CXX/Objects.hxx>

#include <Base/Console.h>
#include <Base/Interpreter.h>
#include <Gui/Application.h>
#include <Gui/Language/Translator.h>

#include "Workbench.h"
#include "ViewProviderPocket.h"
#include "ViewProviderHole.h"
#include "ViewProviderBody.h"
#include "ViewProviderSketchBased.h"
#include "ViewProviderPrimitive.h"
#include "ViewProviderBase.h"
#include "ViewProviderDatumVariable.h"

 // use a different name to CreateCommand()
void CreatePartDesignCommands(void);
void CreatePartDesignBodyCommands(void);
void CreatePartDesignPrimitiveCommands(void);

void loadPartDesignResource()
{
	// add resources and reloads the translators
	Q_INIT_RESOURCE(PartDesign);
	Gui::Translator::instance()->refresh();
}

namespace PartDesignGui {
	class Module : public Py::ExtensionModule<Module>
	{
	public:
		Module() : Py::ExtensionModule<Module>("PartDesignGui")
		{
			initialize("This module is the PartDesignGui module."); // register with Python
		}

		virtual ~Module() {}

	private:
	};

	PyObject* initModule()
	{
		return (new Module)->module().ptr();
	}
} // namespace PartDesignGui

/* Python entry */
PyMOD_INIT_FUNC(PartDesignGui)
{
	if (!Gui::Application::Instance) {
		PyErr_SetString(PyExc_ImportError, "Cannot load Gui module in console application.");
		PyMOD_Return(0);
	}

	try {
		Base::Interpreter().runString("import PartGui");
		Base::Interpreter().runString("import SketcherGui");
	}
	catch (const Base::Exception & e) {
		PyErr_SetString(PyExc_ImportError, e.what());
		PyMOD_Return(0);
	}

	PyObject* mod = PartDesignGui::initModule();
	Base::Console().Log("Loading GUI of PartDesign module... done\n");

	// instantiating the commands
	CreatePartDesignCommands();
	CreatePartDesignBodyCommands();
	CreatePartDesignPrimitiveCommands();

	PartDesignGui::Workbench::init();
	PartDesignGui::ViewProvider::init();
	PartDesignGui::ViewProviderPython::init();
	PartDesignGui::ViewProviderBody::init();
	PartDesignGui::ViewProviderSketchBased::init();
	PartDesignGui::ViewProviderPocket::init();
	PartDesignGui::ViewProviderHole::init();
	PartDesignGui::ViewProviderPrimitive::init();
	PartDesignGui::ViewProviderBase::init();
	PartDesignGui::ViewProviderDatumVariable::init();

	// add resources and reloads the translators
	loadPartDesignResource();

	PyMOD_Return(mod);
}