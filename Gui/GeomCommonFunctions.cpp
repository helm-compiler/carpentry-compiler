#include "GeomCommonFunctions.h"



Standard_Real GeomFunc::GetLength(const TopoDS_Edge& edge)
{
	Standard_Real umin, umax;
	Handle(Geom_Curve) curve = BRep_Tool::Curve(edge, umin, umax);
	GeomAdaptor_Curve adaptorCurve(curve, umin, umax);
	Standard_Real length = GCPnts_AbscissaPoint::Length(adaptorCurve, umin, umax);
	return length;
}

unsigned int GeomFunc::GetNumberOfEdges(const TopoDS_Shape& shape)
{
	TopExp_Explorer edgeExpl(shape, TopAbs_EDGE);
	unsigned int iEdges = 0;
	for (; edgeExpl.More(); edgeExpl.Next()) {
		iEdges++;
	}

	return iEdges;
}

unsigned int GeomFunc::GetNumberOfFaces(const TopoDS_Shape& shape)
{
	TopExp_Explorer faceExpl(shape, TopAbs_FACE);
	unsigned int iFaces = 0;
	for (; faceExpl.More(); faceExpl.Next()) {
		iFaces++;
	}

	return iFaces;
}


unsigned int GeomFunc::GetNumberOfSubshapes(const TopoDS_Shape& shape)
{
	if (shape.ShapeType() == TopAbs_COMPOUND) {
		unsigned int n = 0;
		for (TopoDS_Iterator anIter(shape); anIter.More(); anIter.Next()) {
			n++;
		}
		return n;
	}
	else {
		return 0;
	}
}
