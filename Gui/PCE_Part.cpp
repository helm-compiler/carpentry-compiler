#include <PCE.h>
#include "GeomCommonFunctions.h"
#include <Gui/Application.h>
#include <BRepBuilderAPI_GTransform.hxx>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include<tchar.h>
#include <iostream>
#include "Standard_Type.hxx"
#include <Standard.hxx>
#include <Standard_DefineAlloc.hxx>
#include <Standard_PrimitiveTypes.hxx>
#include <Geom_Circle.hxx>
#include <clipper/clipper.hpp>



using namespace PartDesignGui;


PCEShape::PCEShape()
{

}

PCEShape::PCEShape(const TopoDS_Shape& ashape_)
{
	ashape = ashape_;
	//LoadShapeGeometry(ashape_);
	GeomFunc::LoadTopoDSGeometry(ashape_, surfs, surf_normals,volume, holes, center);
}


PCEShape::PCEShape(const PCEShape& pce_shape)
{
	index = pce_shape.index;
	ashape = pce_shape.ashape;
	surfs = pce_shape.surfs;
	holes = pce_shape.holes;
	surf_normals = pce_shape.surf_normals;
	center = pce_shape.center;
	volume = pce_shape.volume;
	init_LUM = pce_shape.init_LUM;
	WLC_index = pce_shape.WLC_index;
	transforms = pce_shape.transforms;
	LGM = pce_shape.LGM;
}

void PCEShape::Clear()
{
	Vector3d2().swap(surfs);
	Vector3d1().swap(surf_normals);
	std::vector<std::pair<Vector3d1, double>>().swap(holes);
	//part_LUMs.clear();
	//std::vector<glm::dmat4>().swap(LUMs);
	for (auto& transform : transforms) transform.Clear();
	std::vector<PCETransform>().swap(transforms);
}





