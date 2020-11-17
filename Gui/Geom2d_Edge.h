#ifndef Geom2d_Edge_H
#define Geom2d_Edge_H Geom2d_Edge_H

#ifndef _Standard_HeaderFile
#include <Standard.hxx>
#endif

#ifndef _Standard_Macro_HeaderFile
#include <Standard_Macro.hxx>
#endif

#include <Standard_DefineHandle.hxx>
#include <gp_Pnt2d.hxx>
#include <gp_Dir2d.hxx>
#include <Geom2d_Line.hxx>
#include <ElCLib.hxx>

class gp_Pnt2d2d;
class Geom2d_Edge;
DEFINE_STANDARD_HANDLE(Geom2d_Edge, Geom2d_Line)

class Geom2d_Edge : public Geom2d_Line{
public:
	/**
	* \fn Geom2d_Edge()
	* \brief Constructs a Geom2d_Edge
	*/
	Geom2d_Edge();

	/**
	* \fn ~Geom2d_Edge()
	* \brief destructor
	*/
	~Geom2d_Edge();

	/**
	* \fn SetPoints(const gp_Pnt2d& p1,const gp_Pnt2d& p2)
	* \brief set edge  by 2 points
	* \return Standard_Boolean
	* \param p1 const gp_Pnt2d&
	* \param p2 const gp_Pnt2d&
	*/
	Standard_Boolean SetPoints(const gp_Pnt2d& p1, const gp_Pnt2d& p2);

	/**
	* \fn GetStart_Pnt()
	* \brief get start 2d point
	* \return gp_Pnt2d const
	*/
	gp_Pnt2d GetStart_Pnt() const;

	/**
	* \fn GetEnd_Pnt()
	* \brief get end 2d point
	* \return gp_Pnt2d const
	*/
	gp_Pnt2d GetEnd_Pnt() const;

	/**
	* \fn MiddlePnt()
	* \brief get middle 2d point
	* \return gp_Pnt2d const
	*/
	gp_Pnt2d MiddlePnt() const;

	/**
	* \fn StartParameter()
	* \brief get parameter of start 2d point
	* \return Standard_Real const
	*/
	Standard_Real StartParameter() const;

	/**
	* \fn EndParameter()
	* \brief get parameter of end 2d point
	* \return Standard_Real const
	*/
	Standard_Real EndParameter() const;

	DEFINE_STANDARD_RTTIEXT(Geom2d_Edge, Geom2d_Line)

private:

	gp_Pnt2d	StartPnt;
	gp_Pnt2d	EndPnt;

};

#endif