#include "Compiler.h"
#include <QTextStream>
#include <unordered_map>
#include <App/Part.h>
#include <Mod/PartDesign/App/Body.h>
#include <Mod/Part/App/PartFeature.h>
#include "Mod/PartDesign/App/FeaturePrimitive.h"
#include "Mod/Sketcher/App/SketchObject.h"
#include "Mod/PartDesign/App/FeaturePocket.h"
#include "Mod/PartDesign/App/FeatureHole.h"
#include <TopoDS.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <TopExp_Explorer.hxx>
#include "GeomCommonFunctions.h"
#include "Mod/Part/App/DatumFeature.h"
#include <App/OriginFeature.h>

#ifndef sq
#define sq QString::fromUtf8
#endif

#ifndef u(x)
#define u(x) QString::fromUtf8(x)
#endif

std::string IntToString(int n, int width)
{
	const auto format = "%0" + std::to_string(width) + "d";
	char buffer[256];
	sprintf(buffer, format.c_str(), n);
	return std::string(buffer);
}

QString IntToQString(int n, int width)
{
	const auto format = "%0" + std::to_string(width) + "d";
	char buffer[256];
	sprintf(buffer, format.c_str(), n);
	return sq(buffer);
}

QString GetGeomStatement(std::unordered_map<const Part::Geometry*, QString>& mapGeomToName, const Part::Geometry* geom,
	Sketcher::PointPos pos)
{
	QString result;
	QTextStream statement(&result);
	if (mapGeomToName.find(geom) != mapGeomToName.end())
	{
		if (pos == Sketcher::end)
		{
			statement << "End(" << mapGeomToName[geom] << ")";
		}
		else if (pos == Sketcher::start)
		{
			statement << "Start(" << mapGeomToName[geom] << ")";
		}
		else
		{
			statement << mapGeomToName[geom];
		}
	}
	else
	{
		if (pos == Sketcher::start)
			statement << "Start(";
		else if (pos == Sketcher::end)
			statement << "End(";

		if (geom->Id == -1)
		{
			statement << "Horizontal";
		}
		else if (geom->Id == -2)
		{
			statement << "Vertical";
		}
		/*else if (geom->isDerivedFrom(Part::GeomLineSegment::getClassTypeId()))
		{
			const auto p = dynamic_cast<const Part::GeomLineSegment*>(geom);
			const auto l1Nd = p->getEndPoint();
			const auto l1St = p->getStartPoint();
			statement << "Line(" << l1St.x << ", " << l1St.y << ", " << l1St.z
				<< ","
				<< l1Nd.x << "," << l1Nd.y << "," << l1Nd.z << ")";
		}
		else if (geom->isDerivedFrom(Part::GeomArcOfCircle::getClassTypeId()))
		{
			std::cerr << "Other type of point on object " << geom->getTypeId().getName() << std::endl;
			system("pause");
			const auto p = dynamic_cast<const Part::GeomArcOfCircle*>(geom);
		}*/
		else
		{
			std::cerr << "Other type of point on object " << geom->getTypeId().getName() << std::endl;
			system("pause");
		}

		if (pos == Sketcher::start || pos == Sketcher::end)
			statement << ")";
	}

	return result;
}

void Compiler::CompileDown(std::vector<App::DocumentObject*>& sortedObjects)
{
	QTextStream hlCodeStream(&hlHELM);
	//hlCodeStream.setRealNumberPrecision(8);
	QTextStream llCodeStream(&llHELM);
	//llCodeStream.setRealNumberPrecision(8);
	std::unordered_map<App::DocumentObject*, bool> isFeatVisited;
	std::unordered_map<const Part::Geometry*, QString> mapGeomToName;
	int variableId = 0;
	for (auto objIt = sortedObjects.rbegin(); objIt != sortedObjects.rend(); ++objIt)
	{
		isFeatVisited[*objIt] = false;
	}
	for (auto objIt = sortedObjects.rbegin(); objIt != sortedObjects.rend(); ++objIt)
	{
		auto curObj = *objIt;

		if (!curObj->isDerivedFrom(Part::Feature::getClassTypeId()) || curObj->isDerivedFrom(
			PartDesign::Body::getClassTypeId()))
		{
			isFeatVisited[curObj] = true;
			continue;
		}

		if (!isFeatVisited[curObj])
		{
			std::cerr << curObj->getNameInDocument() << std::endl;
			isFeatVisited[curObj] = true;

			if (curObj->isDerivedFrom(PartDesign::AdditiveBox::getClassTypeId()))
			{
				auto curAddBox = dynamic_cast<PartDesign::AdditiveBox*>(curObj);
				// High-level HELM
				hlCodeStream << curObj->getNameInDocument() << " = Make_Stock(";
				hlCodeStream << curAddBox->Length.getValue() << ",";
				hlCodeStream << curAddBox->Width.getValue() << ",";
				hlCodeStream << curAddBox->Height.getValue() << ");\n";

				// Low-level HLEM
				QString nameInDoc = QString::fromUtf8(curObj->getNameInDocument());
				llCodeStream << CompileBox(curAddBox, nameInDoc);
			}
			else if (curObj->isDerivedFrom(Sketcher::SketchObject::getClassTypeId()))
			{
				auto curSkeObj = dynamic_cast<Sketcher::SketchObject*>(curObj);

				if (curSkeObj == nullptr || curSkeObj->Support.getValue() == nullptr) continue;

				auto supportNameInDoc = sq(curSkeObj->Support.getValue()->getNameInDocument());

				auto& extGeo = curSkeObj->ExternalGeo.getValues();
				auto& extGeometry = curSkeObj->ExternalGeometry.getValues();
				auto& extGeomSub = curSkeObj->ExternalGeometry.getSubValues();

				if (extGeo.size() - 2 != extGeometry.size())
				{
					std::cerr << "insame length " << std::endl;
					system("pause");
				}

				for (auto i = 0; i < extGeo.size(); ++i)
				{
					if (extGeo[i]->Id < 0) continue;

					auto tmpObj = extGeometry[i - 2];
					auto tmpSub = extGeomSub[i - 2];

					TopoDS_Shape refSubShape;
					if (tmpObj->getTypeId().isDerivedFrom(Part::Datum::getClassTypeId()))
					{
						const Part::Datum* datum = static_cast<const Part::Datum*>(tmpObj);
						refSubShape = datum->getShape();
					}
					else if (tmpObj->getTypeId().isDerivedFrom(App::Plane::getClassTypeId()))
					{
						const App::Plane* pl = static_cast<const App::Plane*>(tmpObj);
						Base::Placement plm = pl->Placement.getValue();
						Base::Vector3d base = plm.getPosition();
						Base::Rotation rot = plm.getRotation();
						Base::Vector3d normal(0, 0, 1);
						rot.multVec(normal, normal);
						gp_Pln plane(gp_Pnt(base.x, base.y, base.z), gp_Dir(normal.x, normal.y, normal.z));
						BRepBuilderAPI_MakeFace fBuilder(plane);
						if (fBuilder.IsDone())
						{
							TopoDS_Face f = TopoDS::Face(fBuilder.Shape());
							refSubShape = f;
						}
					}
					else
					{
						refSubShape = Part::Feature::getShape(
							tmpObj, tmpSub.c_str(), true, nullptr, nullptr, true, false);
					}

					switch (refSubShape.ShapeType())
					{
					case TopAbs_EDGE:
					{
						const TopoDS_Edge& edge = TopoDS::Edge(refSubShape);
						BRepAdaptor_Curve curve(edge);
						std::cerr << curve.GetType() << std::endl;
						if (curve.GetType() == GeomAbs_Line)
						{
							gp_Pnt st, nd;
							GeomFunc::GetEdgeEndPnt(edge, st, nd);
							auto mid = (st.XYZ() + nd.XYZ()) / 2.0;
							auto geomName = sq("Query_Edge_By_Closest_Point(") + supportNameInDoc +
								sq(",") + QString::number(mid.X()) + sq(",") + QString::number(mid.Y())
								+ sq(",") + QString::number(mid.Z()) + sq(")");
							mapGeomToName[extGeo[i]] = geomName;
						}
						else if (curve.GetType() == GeomAbs_Circle)
						{
							gp_Circ circle = curve.Circle();
							gp_Pnt cnt = circle.Location();
							gp_Pnt beg = curve.Value(curve.FirstParameter());
							gp_Pnt end = curve.Value(curve.LastParameter());

							auto geomName = sq("Query_Arc_By_Closest_Center_And_Radius(") +
								supportNameInDoc + sq(",") + QString::number(cnt.X()) + sq(",") +
								QString::number(cnt.Y()) + sq(",") + QString::number(cnt.Z()) + sq(",")
								+ QString::number(circle.Radius()) + sq(")");
							mapGeomToName[extGeo[i]] = geomName;

							if (beg.SquareDistance(end) < Precision::Confusion())
							{
								Part::GeomCircle* gCircle = new Part::GeomCircle();
								gCircle->setRadius(circle.Radius());
								gCircle->setCenter(Base::Vector3d(cnt.X(), cnt.Y(), cnt.Z()));

								gCircle->Construction = true;
							}
							else
							{
								Part::GeomArcOfCircle* gArc = new Part::GeomArcOfCircle();
								Handle(Geom_Curve) hCircle = new Geom_Circle(circle);
								Handle(Geom_TrimmedCurve) tCurve = new Geom_TrimmedCurve(
									hCircle, curve.FirstParameter(),
									curve.LastParameter());
								gArc->setHandle(tCurve);
								gArc->Construction = true;
							}
						}
						else
						{
							std::cerr << "Unsuppored type" << std::endl;
						}
						break;
					}
					default:
						std::cerr << "Unsuppored type" << std::endl;
						break;
					}

					//auto geomName = sq("Query_Edge_By_Closest_Point(") + sq("face,") + QString::number(mid[0]) + sq(",") + QString::number(mid[1]) + sq(")");
					//mapGeomToName[lineSeg] = geomName;
					//std::cerr << tmpGeometry << std::endl;
				}

				// First, define geom details on the skecth

				for (auto& geo : curSkeObj->Geometry.getValues())
				{
					if (geo->getTypeId() == Part::GeomLineSegment::getClassTypeId())
					{
						const auto lineSeg = dynamic_cast<const Part::GeomLineSegment*>(geo);
						const auto& st = lineSeg->getStartPoint();
						const auto& nd = lineSeg->getEndPoint();
						auto geomName = sq("MyLine") + IntToQString(variableId++, 3);
						mapGeomToName[geo] = geomName;

						hlCodeStream << geomName << " = ";
						hlCodeStream << "Line(" << st.x << "," << st.y << ","
							<< nd.x << "," << nd.y << ");\n";
					}
					else if (geo->getTypeId() == Part::GeomCircle::getClassTypeId())
					{
						const auto circ = dynamic_cast<const Part::GeomCircle*>(geo);
						const auto& center = circ->getCenter();
						const auto& radius = circ->getRadius();

						auto geomName = sq("MyCircle") + IntToQString(variableId++, 3);
						mapGeomToName[geo] = geomName;

						hlCodeStream << geomName << " = ";
						hlCodeStream << "Circle(" << center.x << "," << center.y << "," << radius << ");\n";
					}
					else if (geo->getTypeId() == Part::GeomArcOfCircle::getClassTypeId())
					{
						const Part::GeomArcOfCircle* arc = dynamic_cast<const Part::GeomArcOfCircle*>(geo);
						const auto& c = arc->getCenter();
						double u, v;
						arc->getRange(u, v, false);
						auto r = arc->getRadius();
						//const auto& st = arc->getStartPoint();
						//const auto& nd = arc->getEndPoint();

						auto geomName = sq("MyArc") + IntToQString(variableId++, 3);
						mapGeomToName[geo] = geomName;

						hlCodeStream << geomName << " = ";
						hlCodeStream << "Arc(" << c.x << "," << c.y << "," << r << "," << u << "," << v << ");\n";
					}
					else
					{
						std::cerr << "Unknown geometry: " << geo->getTypeId().getName() << std::endl;
					}
				}

				// High-level HELM
				hlCodeStream << curObj->getNameInDocument()
					<< " = Make_Sketch(Query_Face_By_Closest_Point("
					<< supportNameInDoc << ",";

				std::cerr << "Now is processing " << curObj->getNameInDocument() << std::endl;

				auto appliedFeat = curSkeObj->Support.getValue();
				const auto& sups = curSkeObj->Support.getSubValues();
				PartDesign::TopoShape shape;

				if (appliedFeat->isDerivedFrom(PartDesign::Pocket::getClassTypeId()))
				{
					auto pocketFeat = dynamic_cast<PartDesign::Pocket*>(appliedFeat);
					shape = pocketFeat->Shape.getShape();
				}
				else if (appliedFeat->isDerivedFrom(PartDesign::AdditiveBox::getClassTypeId()))
				{
					auto additiveBox = dynamic_cast<PartDesign::AdditiveBox*>(appliedFeat);
					shape = additiveBox->AddSubShape.getShape();
				}

				auto cnt = 0;
				for (auto& s : sups)
				{
					if (cnt > 0)
						hlCodeStream << ",";
					auto face = Part::Feature::getShape(appliedFeat, s.c_str(), true, nullptr, nullptr,
						true, false);
					//auto face = shape.getSubShape(s.c_str());
					auto cp = GeomFunc::GetCentralFacePoint(TopoDS::Face(face));

					hlCodeStream << cp.X() << "," << cp.Y() << "," << cp.Z();
					//hlCodeStream << s.c_str();
					++cnt;
				}

				hlCodeStream << "),Geometry(";
				cnt = 0;
				for (auto& geo : curSkeObj->Geometry.getValues())
				{
					if (geo->getTypeId() == Part::GeomLineSegment::getClassTypeId() ||
						geo->getTypeId() == Part::GeomCircle::getClassTypeId() ||
						geo->getTypeId() == Part::GeomArcOfCircle::getClassTypeId())
					{
						if (cnt > 0) hlCodeStream << ",";
						hlCodeStream << mapGeomToName[geo];

						++cnt;
					}
					else
					{
						std::cerr << "Unknown geometry: " << geo->getTypeId().getName() << std::endl;
					}
				}

				hlCodeStream << ")";

				if (curSkeObj->Constraints.getSize() > 0)
				{
					hlCodeStream << ",Constraint(";
					auto csCnt = 0;
					for (auto& cntrs : curSkeObj->Constraints.getValues())
					{
						if (cntrs->Type == Sketcher::PointOnObject)
						{
							if (csCnt > 0)
								hlCodeStream << ",";

							hlCodeStream << "PointOnObject(";

							hlCodeStream << GetGeomStatement(
								mapGeomToName, curSkeObj->getGeometry(cntrs->First), cntrs->FirstPos);

							hlCodeStream << "," << GetGeomStatement(
								mapGeomToName, curSkeObj->getGeometry(cntrs->Second),
								cntrs->SecondPos) << ")";

							++csCnt;
						}
						else if (cntrs->Type == Sketcher::Coincident)
						{
							if (csCnt > 0)
								hlCodeStream << ",";

							hlCodeStream << "Coincident(";

							hlCodeStream << GetGeomStatement(
								mapGeomToName, curSkeObj->getGeometry(cntrs->First), cntrs->FirstPos);

							hlCodeStream << "," << GetGeomStatement(
								mapGeomToName, curSkeObj->getGeometry(cntrs->Second),
								cntrs->SecondPos) << ")";

							++csCnt;
						}

						else if (cntrs->Type == Sketcher::Equal)
						{
							if (csCnt > 0)
								hlCodeStream << ",";

							hlCodeStream << "Equal(";

							hlCodeStream << GetGeomStatement(
								mapGeomToName, curSkeObj->getGeometry(cntrs->First), cntrs->FirstPos);

							hlCodeStream << "," << GetGeomStatement(
								mapGeomToName, curSkeObj->getGeometry(cntrs->Second),
								cntrs->SecondPos) << ")";

							++csCnt;
						}
						else if (cntrs->Type == Sketcher::Angle)
						{
							if (csCnt > 0)
								hlCodeStream << ",";

							hlCodeStream << "Angle(";

							if (curSkeObj->getGeometry(cntrs->First) != nullptr)
								hlCodeStream << GetGeomStatement(
									mapGeomToName, curSkeObj->getGeometry(cntrs->First), cntrs->FirstPos);

							if (curSkeObj->getGeometry(cntrs->Second) != nullptr)
								hlCodeStream << "," << GetGeomStatement(
									mapGeomToName, curSkeObj->getGeometry(cntrs->Second), cntrs->SecondPos);

							hlCodeStream << "," << cntrs->getValue() / M_PI * 180.0 << ")";
						}
						else if (cntrs->Type == Sketcher::Distance)
						{
							if (csCnt > 0)
								hlCodeStream << ",";

							hlCodeStream << "Distance(";

							if (curSkeObj->getGeometry(cntrs->First) != nullptr)
								hlCodeStream << GetGeomStatement(
									mapGeomToName, curSkeObj->getGeometry(cntrs->First), cntrs->FirstPos);

							if (curSkeObj->getGeometry(cntrs->Second) != nullptr)
								hlCodeStream << "," << GetGeomStatement(
									mapGeomToName, curSkeObj->getGeometry(cntrs->Second), cntrs->SecondPos);

							hlCodeStream << "," << cntrs->getValue() << ")";

							++csCnt;
						}
						else if (cntrs->Type == Sketcher::DistanceX)
						{
							if (csCnt > 0)
								hlCodeStream << ",";

							hlCodeStream << "DistanceX(";

							hlCodeStream << GetGeomStatement(
								mapGeomToName, curSkeObj->getGeometry(cntrs->First), cntrs->FirstPos);

							hlCodeStream << "," << GetGeomStatement(
								mapGeomToName, curSkeObj->getGeometry(cntrs->Second), cntrs->SecondPos);

							hlCodeStream << "," << cntrs->getValue() << ")";

							++csCnt;
						}
						else if (cntrs->Type == Sketcher::DistanceY)
						{
							if (csCnt > 0)
								hlCodeStream << ",";

							hlCodeStream << "DistanceY(";

							hlCodeStream << GetGeomStatement(
								mapGeomToName, curSkeObj->getGeometry(cntrs->First), cntrs->FirstPos);

							hlCodeStream << "," << GetGeomStatement(
								mapGeomToName, curSkeObj->getGeometry(cntrs->Second), cntrs->SecondPos);

							hlCodeStream << "," << cntrs->getValue() << ")";

							++csCnt;
						}
						else if (cntrs->Type == Sketcher::Perpendicular)
						{
							if (csCnt > 0)
								hlCodeStream << ",";

							hlCodeStream << "Perpendicular(";

							hlCodeStream << GetGeomStatement(
								mapGeomToName, curSkeObj->getGeometry(cntrs->First), cntrs->FirstPos);

							hlCodeStream << "," << GetGeomStatement(
								mapGeomToName, curSkeObj->getGeometry(cntrs->Second),
								cntrs->SecondPos) << ")";

							++csCnt;
						}
						else if (cntrs->Type == Sketcher::Parallel)
						{
							if (csCnt > 0)
								hlCodeStream << ",";

							hlCodeStream << "Parallel(";

							hlCodeStream << GetGeomStatement(
								mapGeomToName, curSkeObj->getGeometry(cntrs->First), cntrs->FirstPos);

							hlCodeStream << "," << GetGeomStatement(
								mapGeomToName, curSkeObj->getGeometry(cntrs->Second),
								cntrs->SecondPos) << ")";

							++csCnt;
						}
						else if (cntrs->Type == Sketcher::Radius)
						{
							// TODO: incorporate this constraint
						}
						else
						{
							std::cerr << "Unknown type " << cntrs->Type << std::endl;
							system("pause");
						}
					}
					hlCodeStream << ")";
				}

				hlCodeStream << ");\n";
			}
			else if (curObj->isDerivedFrom(PartDesign::Pocket::getClassTypeId()))
			{
				auto curCutObj = dynamic_cast<PartDesign::Pocket*>(curObj);
				if (curCutObj->Profile.getValue()->isDerivedFrom(Sketcher::SketchObject::getClassTypeId()))
				{
					auto curProfObj = dynamic_cast<Sketcher::SketchObject*>(curCutObj->Profile.getValue());
					hlCodeStream << curObj->getNameInDocument() << " = "
						<< "Make_Cut(" << curProfObj->Support.getValue()->getNameInDocument() << "," <<
						curCutObj
						->Profile.
						getValue()->
						getNameInDocument()
						<< "," << curCutObj->Part.getValue() << ");\n";
				}

				// Low-level HLEM
				QString nameInDoc = QString::fromUtf8(curObj->getNameInDocument());
				llCodeStream << CompilePolycut(curCutObj, nameInDoc);
			}
			else if (curObj->isDerivedFrom(PartDesign::Hole::getClassTypeId()))
			{
				auto curHoleObj = dynamic_cast<PartDesign::Hole*>(curObj);
				if (curHoleObj->Profile.getValue()->isDerivedFrom(Sketcher::SketchObject::getClassTypeId()))
				{
					auto curProfObj = dynamic_cast<Sketcher::SketchObject*>(curHoleObj->Profile.getValue());
					hlCodeStream << curObj->getNameInDocument() << " = "
						<< "Make_Hole(" << curProfObj->Support.getValue()->getNameInDocument() << "," <<
						curHoleObj->Profile.getValue()->getNameInDocument()
						<< ", " << curHoleObj->Diameter.getValue()
						<< ");\n";
				}

				// Low-level HLEM
				QString nameInDoc = QString::fromUtf8(curObj->getNameInDocument());
				// llCodeStream << cdWb->langCompiler->CompilePolycut(curCutObj, nameInDoc);
			}

			// visit all in-features
			for (auto inObjIt : curObj->getInList())
			{
				continue;
				if (!isFeatVisited[inObjIt])
				{
					std::cerr << inObjIt->getNameInDocument() << std::endl;
					isFeatVisited[inObjIt] = true;
				}
			}
		}
	}
}

QString Compiler::GetHLHELM() const
{
	return hlHELM;
}

QString Compiler::GetLLHELM() const
{
	return llHELM;
}

QString Compiler::CompileBox(PartDesign::AdditiveBox* box, QString& varName)
{
	QString boxHELM;
	boxHELM += varName + sq("=Box(") + sq(box->Type.getEnumVector()[box->Type.getValue()].c_str()) +
		sq(", Off=") + QString::number(box->Length.getValue()/25.4) + sq(";\n");

#ifdef OLD_VERSION
	const double MIN_JIGSAW = 2;
	const double MIN_CHOPSAW = 2;

	double side[3] = { box->Length.getValue(), 10.0, 10.0 };
	std::sort(side, side + 3);

	double maxSide = side[2];
	bool solutionFound = false;

	double objVal = -DBL_MAX;
	Material* bestMat;

	std::cout << "start compiling box" << std::endl;

	for (auto& mat : vecMaterial)
	{
		if (mat->materialType == MaterialType::LumberMaterial)
		{
			auto lumber = dynamic_cast<LumberWood*>(mat);

			if (maxSide - TOOL_EPS < lumber->height)
			{
				if (side[0] - TOOL_EPS < lumber->length && side[1] - TOOL_EPS < lumber->width)
				{
					if (std::abs(side[2] - lumber->height) > TOOL_EPS)
					{
						// Manufacturable constraint (Chopsaw)
						if (side[2] < MIN_CHOPSAW || ((lumber->height - side[2]) < MIN_CHOPSAW && (lumber->height - side[2]) > 1e-2))
						{
							continue;
						}
					}

					// Manufacturable constraint (Bandsaw 1)
					if (std::abs(side[0] - lumber->length) > TOOL_EPS)
					{
						if (side[0] < MIN_JIGSAW ||
							((lumber->length - side[0]) < MIN_JIGSAW && (lumber->length - side[0]) > 1e-2))
						{
							continue;
						}
					}

					// Manufacturable constraint (Bandsaw 2)
					if (std::abs(side[1] - lumber->width) > TOOL_EPS)
					{
						if (side[1] < MIN_JIGSAW ||
							((lumber->width - side[1]) < MIN_JIGSAW && (lumber->width - side[1]) > 1e-2))
						{
							continue;
						}
					}

					// Objective value
					double tmpVal = (side[0] * side[1]) / (lumber->length * lumber->width) * (side[2]) / (lumber->height);
					if (tmpVal > objVal)
					{
						solutionFound = true;
						objVal = tmpVal;
						bestMat = mat;
					}
				}
			}
		}
	}

	if (solutionFound)
	{
		auto lumber = dynamic_cast<LumberWood*>(bestMat);

		// exchange width and length
		auto tempVal = lumber->width;
		lumber->width = lumber->length;
		lumber->length = tempVal;

		QString strLen = QString::number(lumber->length);
		QString strWid = QString::number(lumber->width);
		QString strHei = QString::number(lumber->height);
		boxHELM += varName + sq(" = Lumber(") + strLen + sq(",") + strWid + sq(",") + strHei + sq(")\n");

		// create two faces to hold the lumber
		boxHELM += varName + sq("_s1 = Segment(Point(0, 0, 0), Point(0, ") + strWid + sq(", 0));\n");

		boxHELM += varName + sq("_s2 = Segment(Point(0, ") + strWid + sq(", 0)), Point(0, ") + strWid + sq(",") + strHei + sq("));\n");

		boxHELM += varName + sq("_s3 = Segment(Point(") + strLen + sq(", 0, 0), Point(") + strLen + sq(",") + strWid + sq(", 0));\n");

		boxHELM += varName + sq("_s4 = Segment(Point(") + strLen + sq(",") + strWid + sq(", 0)), Point(") + strLen + sq(",") + strWid + sq(",") + strHei + sq("));\n");

		gp_Pln supPlaneCS[2];
		supPlaneCS[0] = gp_Pln(gp_Pnt(0, lumber->width, 0), gp_Dir(-1, 0, 0));
		supPlaneCS[1] = gp_Pln(gp_Pnt(lumber->length, lumber->width, 0), gp_Dir(-1, 0, 0));

		boxHELM += varName + sq("_r1 = Ref(") + varName + sq("_s2, ") + QString::number(side[2]) + sq(");\n");

		boxHELM += varName + sq("_r2 = Ref(") + varName + sq("_s4, ") + QString::number(side[2]) + sq("); \n");

		boxHELM += varName + sq("_path = Path(") + varName + sq("_r1, ") + varName + sq("_r2);\n");

		if (std::abs(side[2] - lumber->height) > TOOL_EPS)
		{
			boxHELM += varName + sq("_f1 = Face(") + varName + sq("_s1, ") + varName + sq("_s2);\n");

			boxHELM += varName + sq("_f2 = Face(") + varName + sq("_s3, ") + varName + sq("_s4);\n");
			boxHELM += varName + sq(" = Chopsaw(") + varName + sq(",") + varName + sq("_f1, ") + varName + sq("_f2, ") + varName + sq("_path);\n");
		}

		if (std::abs(side[0] - lumber->width) > TOOL_EPS)
		{
			boxHELM += sq("Bandsaw(line (X ") + QString::number(side[0]) + sq(")\n");
		}

		if (std::abs(side[1] - lumber->length) > TOOL_EPS)
		{
			boxHELM += sq("Bandsaw(line (Y ") + QString::number(side[1]) + sq(")\n");
		}

		tempVal = lumber->width;
		lumber->width = lumber->length;
		lumber->length = tempVal;
	}
	else
	{
		boxHELM = sq("The object is too large to be fabricated.\n");
	}
#endif

	return boxHELM;
}

#include <Mod/Part/App/PartFeatures.h>
#include "PHELM.h"

QString Compiler::CompilePolycut(PartDesign::Pocket* polycut, QString& varName)
{
	QString polyCutHELM;

	if (polycut->Profile.getValue()->isDerivedFrom(Sketcher::SketchObject::getClassTypeId()))
	{
		auto curProfObj = dynamic_cast<Sketcher::SketchObject*>(polycut->Profile.getValue());

		std::vector<Part::GeomLineSegment*> lineSegs;
		for (auto& geo : curProfObj->Geometry.getValues())
		{
			if (geo->getTypeId() == Part::GeomLineSegment::getClassTypeId())
			{
				lineSegs.push_back(dynamic_cast<Part::GeomLineSegment*>(geo));
			}
		}

		if (lineSegs.size() == 1)
		{
			Part::GeomLineSegment* line = lineSegs[0];
			auto st = line->getStartPoint();
			auto nd = line->getEndPoint();
			Base::Vector3d dirY(0, 1, 0);

			double offset = 0;
			double theta = 0;
			if (st.y > nd.y)
				offset = st.x;
			else
				offset = nd.x;

			theta = std::acos(dirY * (st - nd).Normalize()) * 180.0 / M_PI;
			// Chopsaw
			polyCutHELM += varName + sq("=Chopsaw(") + sq(
				curProfObj->Support.getValue()->getNameInDocument()) +
				sq(", Off=") + QString::number(offset) + sq(", Theta=") + QString::number(theta) + sq(
					");\n");
			
			auto m = curProfObj->Placement.getValue();
			m.multVec(st, st);
			m.multVec(nd, nd);
	
			std::cerr << st.x << " " << st.y << " " << st.z << std::endl;
			std::cerr << nd.x << " " << nd.y << " " << nd.z << std::endl;
			gp_Dir lineDir(st.x - nd.x, st.y - nd.y, st.z - nd.z);
			TopLoc_Location invObjLoc = polycut->getLocation().Inverted();
			Base::Vector3d SketchVector = polycut->getProfileNormal();
			gp_Dir cutDir(SketchVector.x, SketchVector.y, SketchVector.z);
			cutDir.Transform(invObjLoc.Transformation());
			auto normalDir = cutDir.Crossed(lineDir);
			
			PHELM* p = new PHELM();
			p->SetKerfDistance(0.0);
			p->SetFastMode(true);

			Vector3d cutNormal(normalDir.X(), normalDir.Y(), normalDir.Z()), st3d(st.x, st.y, st.z), nd3d(nd.x, nd.y, nd.z), unk(0, 0, 0);
			PartDesignGui::PCEEdge tmpCut(st3d, nd3d, cutNormal, unk);
			std::vector<PartDesignGui::PCEEdge> cutLines = { tmpCut };
			
			std::list<TopoDS_Shape> listMaterial = {polycut->getBaseShape().getShape()};
			ListPairShapePoly listParts(listMaterial.size());

			for (auto it = listParts.begin(); it != listParts.end(); ++it)
			{
				std::get<0>(*it) = polycut->Shape.getValue();
				TopoDS_Face f;
				PHELM::GetUpDownFaces(std::get<0>(*it),true, std::get<1>(*it), f);
				Base::Vector3d cp;
				GeomFunc::GetCenterOfGravity(std::get<0>(*it), cp);
				std::get<2>(*it) = gp_Pnt(cp.x, cp.y, cp.z);
			}

			p->SimpleCompileSuccCutting(listMaterial, cutLines, listParts);
		}
		else if (lineSegs.size() > 1)
		{
			const auto& wires = polycut->getProfileWires();
			// 			if (wires.size() != 1)
			// 				;
			auto wire = TopoDS::Wire(wires[0].getShape());

			std::vector<Base::Vector3d> wireDirs;
			BRepTools_WireExplorer Ex;
			for (Ex.Init(wire); Ex.More(); Ex.Next())
			{
				Base::Vector3d endPnts[2];
				int cnt = 0;
				for (TopExp_Explorer vertexExplorer(Ex.Current(), TopAbs_VERTEX); vertexExplorer.More();
					vertexExplorer.
					Next())
				{
					const TopoDS_Vertex& vertex = TopoDS::Vertex(vertexExplorer.Current());
					gp_Pnt pt = BRep_Tool::Pnt(vertex);
					endPnts[cnt++].Set(pt.X(), pt.Y(), pt.Z());
				}
				const auto curDir = endPnts[1] - endPnts[0];
				wireDirs.push_back(curDir);
			}

			// auto& pnts = wire.points();

			// std::cerr << pnts.size() << std::endl;
			// See if the line itself is convex
			bool isConvex = true;

			const auto plNormal = polycut->getProfileNormal();
			// First dir
			const auto initDir = wireDirs[0];
			auto lastDir = wireDirs[1];
			const auto initCross = lastDir % initDir;

			bool initPos = (initCross * plNormal > 0) ? true : false;

			for (unsigned i = 2; i < lineSegs.size(); ++i)
			{
				const auto curDir = wireDirs[i];
				const auto curCross = curDir % lastDir;
				lastDir = curDir;
				bool curPos = (curCross * plNormal > 0) ? true : false;
				if (curPos != initPos)
				{
					isConvex = false;
					break;
				}
			}

			if (isConvex)
			{
				// a function to check if we can
				// convex, try to use chopsaw
				polyCutHELM += sq("convex!!!");
			}

			// non-convex, try to use bandsaw/jigsaw
			polyCutHELM += sq("non-convex!!");

			// Bandsaw or chopsaw
		}
	}

	return polyCutHELM;
}