if(MSVC)
    add_definitions(-DHAVE_ACOSH -DHAVE_ATANH -DHAVE_ASINH)
else(MSVC)
    add_definitions(-DHAVE_LIMITS_H -DHAVE_CONFIG_H)
endif(MSVC)

add_subdirectory(sexp)

set(SEXP_INCLUDE_DIR sexp/include)

set(OPT_INCLUDE_DIR opt)
set(CLIPPER_INCLUDE_DIR clipper)
set(SOLVESPACE_INCLUDE_DIR solvespace)

include_directories(
    ${CMAKE_BINARY_DIR}
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_CURRENT_BINARY_DIR}
    ${Boost_INCLUDE_DIRS}
    ${OCC_INCLUDE_DIR}
    ${COIN3D_INCLUDE_DIRS}
    ${ZLIB_INCLUDE_DIR}
    ${PYTHON_INCLUDE_DIRS}
    ${XercesC_INCLUDE_DIRS}
    ${EIGEN3_INCLUDE_DIR}
    ${SEXP_INCLUDE_DIR}
    ${OPT_INCLUDE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CLIPPER_INCLUDE_DIR}
    ${SOLVESPACE_INCLUDE_DIR}
)

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/opt/config.hpp.in" "${CMAKE_CURRENT_SOURCE_DIR}/pagmo/config.hpp" @ONLY)

link_directories(${OCC_LIBRARY_DIR} ${SOLVESPACE_INCLUDE_DIR})

generate_from_xml(ViewProviderPy)

set(PartDesignGui_LIBS
    PartDesign
    SketcherGui
    PartGui
    SketcherGui
    FreeCADGui
)

set(PartDesignGui_MOC_HDRS
    TaskFeaturePick.h
    TaskFeatureParameters.h
    TaskSketchBasedParameters.h
    TaskPocketParameters.h
    TaskHoleParameters.h
    TaskDatumParameters.h
    TaskPrimitiveParameters.h
    CodeEditor.h
    ToolsConfiguration.h
    Highlighter.h
)
fc_wrap_cpp(PartDesignGui_MOC_SRCS ${PartDesignGui_MOC_HDRS})
SOURCE_GROUP("Moc" FILES ${PartDesignGui_MOC_SRCS})

if(BUILD_QT5)
    qt5_add_resources(PartDesignGui_SRCS Resources/PartDesign.qrc)
else()
    qt4_add_resources(PartDesignGui_SRCS Resources/PartDesign.qrc)
endif()

set(PartDesignGui_UIC_SRCS
    TaskFeaturePick.ui
    TaskPocketParameters.ui
    TaskHoleParameters.ui
    TaskPrimitiveParameters.ui
    DlgReference.ui
    ToolsConfiguration.ui
)

if(BUILD_QT5)
    qt5_wrap_ui(PartDesignGui_UIC_HDRS ${PartDesignGui_UIC_SRCS})
else()
    qt4_wrap_ui(PartDesignGui_UIC_HDRS ${PartDesignGui_UIC_SRCS})
endif()

SET(PartDesignGuiViewProvider_SRCS
    ViewProvider.cpp
    ViewProvider.h
    ViewProviderBody.cpp
    ViewProviderBody.h
    ViewProviderSketchBased.cpp
    ViewProviderSketchBased.h
    ViewProviderHole.cpp
    ViewProviderHole.h
    ViewProviderPocket.cpp
    ViewProviderPocket.h
    ViewProviderPrimitive.h
    ViewProviderPrimitive.cpp
    ViewProviderBase.h
    ViewProviderBase.cpp
    ViewProviderDatum.h
    ViewProviderDatum.cpp
    ViewProviderDatumVariable.h
    ViewProviderDatumVariable.cpp
    ViewProviderAddSub.h
    ViewProviderAddSub.cpp
)
SOURCE_GROUP("ViewProvider" FILES ${PartDesignGuiViewProvider_SRCS})

SET(PartDesignGuiTaskDlgs_SRCS
    TaskFeaturePick.ui
    TaskFeaturePick.cpp
    TaskFeaturePick.h
    ReferenceSelection.cpp
    ReferenceSelection.h
    TaskFeatureParameters.cpp
    TaskFeatureParameters.h
    TaskSketchBasedParameters.cpp
    TaskSketchBasedParameters.h
    TaskPocketParameters.ui
    TaskPocketParameters.cpp
    TaskPocketParameters.h
    TaskHoleParameters.ui
    TaskHoleParameters.cpp
    TaskHoleParameters.h
    TaskDatumParameters.cpp
    TaskDatumParameters.h
    TaskPrimitiveParameters.h
    TaskPrimitiveParameters.cpp
)
SOURCE_GROUP("TaskDialogs" FILES ${PartDesignGuiTaskDlgs_SRCS})

SET(PartDesignGuiModule_SRCS
    AppPartDesignGui.cpp
    Command.cpp
	CommandBody.cpp
    CommandPrimitive.cpp
	Geom2d_Edge.cpp
	Geom2d_Edge.h
	HLParser.cpp
	HLParser.h
    Resources/PartDesign.qrc
    PreCompiled.cpp
    PreCompiled.h
	tinyxml2.cpp
	tinyxml2.h
    Utils.cpp
    Utils.h
    Workbench.cpp
    Workbench.h
    WorkflowManager.cpp
    WorkflowManager.h
    CodeEditor.h
    CodeEditor.cpp
)
SOURCE_GROUP("Module" FILES ${PartDesignGuiModule_SRCS})

SET(PCE_SRCS
    clipper/clipper.cpp
    clipper/clipper.hpp
	PCE.cpp
	PCE.h
	PCE_arrange.cpp
	PCE_CS.cpp
	PCE_edge.cpp
	PCE_end.cpp
	PCE_functs.cpp
	PCE_Part.cpp
	PCE_position.cpp
	PCE_transform.cpp
	PCE_WoodLumber.cpp
)
SOURCE_GROUP("HELM_PCE" FILES ${PCE_SRCS})


SET(Python_SRCS
    ViewProviderPy.xml
    ViewProviderPyImp.cpp
)
SOURCE_GROUP("Python" FILES ${Python_SRCS})

# SET(SEXP_SRCS
#     error.hpp
#     float.cpp
#     float.hpp
#     io.cpp
#     io.hpp
#     lexer.cpp
#     lexer.hpp
#     parser.cpp
#     parser.hpp
#     util.cpp
#     util.hpp
#     value.cpp
#     value.hpp
# )
# SOURCE_GROUP("SEXP" FILES ${SEXP_SRCS})

SET(OPT_SRCS
    opt/optimizer.hpp 
    opt/optimizer.cpp 
    opt/egraph.hpp
    opt/egraph.cpp 
    opt/evaluator.hpp 
    opt/evaluator.cpp 
    opt/llhelm.hpp
    opt/llhelm.cpp
    opt/nsgaii.hpp 
    opt/parsers.hpp
    opt/parsers.cpp 
    opt/graphmap.hpp 
)
SOURCE_GROUP("OptCopmiler" FILES ${OPT_SRCS})

SET(HELM_COMPILER
    Compiler.cpp
    Compiler.h
    CompilerConfig.cpp
    CompilerConfig.h
    ToolsConfiguration.h
    ToolsConfiguration.cpp
    Highlighter.h
    Highlighter.cpp
    PHELM.h
    PHELM.cpp
    PHELM_Utils.cpp
    ISA.h
    ISA.cpp
)
SOURCE_GROUP("HELM" FILES ${HELM_COMPILER})

SET(PartDesignGui_SRCS
    ${PartDesignGui_SRCS}
    ${PartDesignGui_UIC_HDRS}
    ${PartDesignGuiModule_SRCS}
    ${PartDesignGuiTaskDlgs_SRCS}
    ${PartDesignGuiViewProvider_SRCS}
    ${Python_SRCS}
	${PCE_SRCS}
    ${HELM_COMPILER}
    ${OPT_SRCS}
)

SET(PartDesignGuiIcon_SVG Resources/icons/Carpentry.svg)
SET(PartDesignGuiIcon_GEOM_DLL ../Dll/carpentry_geom.dll)
SET(PartDesignGuiIcon_GMP_DLL ../Dll/libgmp-10.dll)
SET(PartDesignGuiIcon_SLVS_DLL ../Dll/slvs.dll)

add_library(PartDesignGui SHARED 
${PartDesignGui_SRCS} 
${PartDesignGuiIcon_SVG} 
${PartDesignGuiIcon_GEOM_DLL}
${PartDesignGuiIcon_GMP_DLL}
${PartDesignGuiIcon_SLVS_DLL}
)


if(WIN32)
  target_link_libraries(PartDesignGui ${PartDesignGui_LIBS} bcrypt sexp slvs)
else()
  target_link_libraries(PartDesignGui sexp)
endif()

SET_BIN_DIR(PartDesignGui PartDesignGui /Mod/PartDesign)
SET_PYTHON_PREFIX_SUFFIX(PartDesignGui)

fc_copy_sources(PartDesignGui "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_DATADIR}/Mod/PartDesign" ${PartDesignGuiIcon_SVG})

fc_copy_sources(PartDesignGui "${CMAKE_BINARY_DIR}/bin" ${PartDesignGuiIcon_GEOM_DLL})
fc_copy_sources(PartDesignGui "${CMAKE_BINARY_DIR}/bin" ${PartDesignGuiIcon_GMP_DLL})
fc_copy_sources(PartDesignGui "${CMAKE_BINARY_DIR}/bin" ${PartDesignGuiIcon_SLVS_DLL})

INSTALL(TARGETS PartDesignGui DESTINATION ${CMAKE_INSTALL_LIBDIR})
INSTALL(FILES ${PartDesignGuiIcon_SVG} DESTINATION "${CMAKE_INSTALL_DATADIR}/Mod/PartDesign/Resources/icons")
