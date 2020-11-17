#include "PCE.h"
#include "GeomCommonFunctions.h"
#include <Gui/Application.h>
#include <BRepBuilderAPI_GTransform.hxx>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <tchar.h>
#include <iostream>

using namespace PartDesignGui;

PCEPosition::PCEPosition()
{
	t = Vector2d(0.0, 0.0);
	style = FreeStyle;
	score = 0;
}

void PCEPosition::Clear()
{
	//transform->Clear();
	std::vector<int>().swap(edge_indexes);
	t = Vector2d(0.0, 0.0);
	style = FreeStyle;
	score = 0;
}

void PCEPosition::ComputeM()
{
	M = Math::Functs::TranslationMatrix(Vector3d(t.x, t.y, 0.));
	glm::vec4 v(Vector3d(transform->bounding_box.center_2d.x, transform->bounding_box.center_2d.y, 0.0), 1.0);
	Vector3d t_v = Vector3d(M * v);
	center = Vector2d(t_v.x, t_v.y);
}

void PCEPosition::UpdateEdges(std::vector<PCEEdge>& edges_)
{
	for (auto edge : edge_indexes)
		edges_[edge].cut.Translate(t);
}
