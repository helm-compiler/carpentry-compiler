/***************************************************************************
 *   Copyright (c) 2011 Juergen Riegel <FreeCAD@juergen-riegel.net>        *
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
#endif

#include "ui_TaskHoleParameters.h"
#include "TaskHoleParameters.h"
#include <Gui/Application.h>
#include <Gui/Document.h>
#include <Gui/BitmapFactory.h>
#include <Gui/ViewProvider.h>
#include <Gui/WaitCursor.h>
#include <Base/Console.h>
#include <Gui/Selection.h>
#include <Gui/Command.h>
#include <Mod/PartDesign/App/FeatureHole.h>
#include <QMessageBox>

using namespace PartDesignGui;
using namespace Gui;

/* TRANSLATOR PartDesignGui::TaskHoleParameters */

TaskHoleParameters::TaskHoleParameters(ViewProviderHole *HoleView, QWidget *parent)
    : TaskSketchBasedParameters(HoleView, parent, "PartDesign_Hole",tr("Hole parameters"))
    , observer(new Observer(this, static_cast<PartDesign::Hole*>(vp->getObject())))
    , isApplying(false)
{
    // we need a separate container widget to add all controls to
    proxy = new QWidget(this);
    ui = new Ui_TaskHoleParameters();
    ui->setupUi(proxy);
    QMetaObject::connectSlotsByName(this);

    /* Remove actual threading parameters for now */
    ui->ThreadPitch->setVisible(false);
    ui->ThreadCutOffInner->setVisible(false);
    ui->ThreadCutOffOuter->setVisible(false);
    ui->ThreadAngle->setVisible(false);

    connect(ui->ThreadPitch, SIGNAL(valueChanged(double)), this, SLOT(threadPitchChanged(double)));
    connect(ui->ThreadAngle, SIGNAL(valueChanged(double)), this, SLOT(threadAngleChanged(double)));
    connect(ui->ThreadCutOffInner, SIGNAL(valueChanged(double)), this, SLOT(threadCutOffInnerChanged(double)));
    connect(ui->ThreadCutOffOuter, SIGNAL(valueChanged(double)), this, SLOT(threadCutOffOuterChanged(double)));
    connect(ui->ThreadSize, SIGNAL(currentIndexChanged(int)), this, SLOT(threadSizeChanged(int)));
    connect(ui->ThreadClass, SIGNAL(currentIndexChanged(int)), this, SLOT(threadClassChanged(int)));
    connect(ui->ThreadFit, SIGNAL(currentIndexChanged(int)), this, SLOT(threadFitChanged(int)));
    connect(ui->Diameter, SIGNAL(valueChanged(double)), this, SLOT(threadDiameterChanged(double)));
    connect(ui->HoleCutType, SIGNAL(currentIndexChanged(int)), this, SLOT(holeCutChanged(int)));
    connect(ui->HoleCutDiameter, SIGNAL(valueChanged(double)), this, SLOT(holeCutDiameterChanged(double)));
    connect(ui->HoleCutDepth, SIGNAL(valueChanged(double)), this, SLOT(holeCutDepthChanged(double)));
    connect(ui->HoleCutCountersinkAngle, SIGNAL(valueChanged(double)), this, SLOT(holeCutCountersinkAngleChanged(double)));
    connect(ui->DepthType, SIGNAL(currentIndexChanged(int)), this, SLOT(depthChanged(int)));
    connect(ui->Depth, SIGNAL(valueChanged(double)), this, SLOT(depthValueChanged(double)));
    connect(ui->DrillPointAngle, SIGNAL(valueChanged(double)), this, SLOT(drillPointAngledValueChanged(double)));
    connect(ui->TaperedAngle, SIGNAL(valueChanged(double)), this, SLOT(taperedAngleChanged(double)));

    PartDesign::Hole* pcHole = static_cast<PartDesign::Hole*>(vp->getObject());

    pcHole->updateProps();

    ui->ThreadPitch->bind(pcHole->ThreadPitch);
    ui->ThreadAngle->bind(pcHole->ThreadAngle);
    ui->ThreadCutOffInner->bind(pcHole->ThreadCutOffInner);
    ui->ThreadCutOffOuter->bind(pcHole->ThreadCutOffOuter);
    ui->Diameter->bind(pcHole->Diameter);
    ui->HoleCutDiameter->bind(pcHole->HoleCutDiameter);
    ui->HoleCutDepth->bind(pcHole->HoleCutDepth);
    ui->HoleCutCountersinkAngle->bind(pcHole->HoleCutCountersinkAngle);
    ui->Depth->bind(pcHole->Depth);
    ui->DrillPointAngle->bind(pcHole->DrillPointAngle);
    ui->TaperedAngle->bind(pcHole->TaperedAngle);
    
    connectPropChanged = App::GetApplication().signalChangePropertyEditor.connect(
            boost::bind(&TaskHoleParameters::changedObject, this, _1, _2));
    
	this->groupLayout()->addWidget(proxy);
}

TaskHoleParameters::~TaskHoleParameters()
{
    delete ui;
}
void TaskHoleParameters::threadPitchChanged(double value)
{
    PartDesign::Hole* pcHole = static_cast<PartDesign::Hole*>(vp->getObject());

    pcHole->ThreadPitch.setValue(value);
    recomputeFeature();
}

void TaskHoleParameters::threadAngleChanged(double value)
{
    PartDesign::Hole* pcHole = static_cast<PartDesign::Hole*>(vp->getObject());

    pcHole->ThreadAngle.setValue(value);
    recomputeFeature();
}

void TaskHoleParameters::threadCutOffInnerChanged(double value)
{
    PartDesign::Hole* pcHole = static_cast<PartDesign::Hole*>(vp->getObject());

    pcHole->ThreadCutOffInner.setValue(value);
    recomputeFeature();
}

void TaskHoleParameters::threadCutOffOuterChanged(double value)
{
    PartDesign::Hole* pcHole = static_cast<PartDesign::Hole*>(vp->getObject());

    pcHole->ThreadCutOffOuter.setValue(value);
    recomputeFeature();
}

void TaskHoleParameters::holeCutChanged(int index)
{
    if (index < 0)
        return;

    PartDesign::Hole* pcHole = static_cast<PartDesign::Hole*>(vp->getObject());

    pcHole->HoleCutType.setValue(index);
    recomputeFeature();
}

void TaskHoleParameters::holeCutDiameterChanged(double value)
{
    PartDesign::Hole* pcHole = static_cast<PartDesign::Hole*>(vp->getObject());

    pcHole->HoleCutDiameter.setValue(value);
    recomputeFeature();
}

void TaskHoleParameters::holeCutDepthChanged(double value)
{
    PartDesign::Hole* pcHole = static_cast<PartDesign::Hole*>(vp->getObject());

    pcHole->HoleCutDepth.setValue(value);
    recomputeFeature();
}

void TaskHoleParameters::holeCutCountersinkAngleChanged(double value)
{
    PartDesign::Hole* pcHole = static_cast<PartDesign::Hole*>(vp->getObject());

    pcHole->HoleCutCountersinkAngle.setValue((double)value);
    recomputeFeature();
}

void TaskHoleParameters::depthChanged(int index)
{
    PartDesign::Hole* pcHole = static_cast<PartDesign::Hole*>(vp->getObject());

    pcHole->DepthType.setValue(index);
    recomputeFeature();
}

void TaskHoleParameters::depthValueChanged(double value)
{
    PartDesign::Hole* pcHole = static_cast<PartDesign::Hole*>(vp->getObject());

    pcHole->Depth.setValue(value);
    recomputeFeature();
}

void TaskHoleParameters::drillPointAngledValueChanged(double value)
{
    PartDesign::Hole* pcHole = static_cast<PartDesign::Hole*>(vp->getObject());

    pcHole->DrillPointAngle.setValue((double)value);
    recomputeFeature();
}

void TaskHoleParameters::taperedAngleChanged(double value)
{
    PartDesign::Hole* pcHole = static_cast<PartDesign::Hole*>(vp->getObject());

    pcHole->TaperedAngle.setValue(value);
    recomputeFeature();
}

void TaskHoleParameters::threadTypeChanged(int index)
{
    if (index < 0)
        return;

    PartDesign::Hole* pcHole = static_cast<PartDesign::Hole*>(vp->getObject());
    pcHole->ThreadType.setValue(index);
}

void TaskHoleParameters::threadSizeChanged(int index)
{
    if (index < 0)
        return;

    PartDesign::Hole* pcHole = static_cast<PartDesign::Hole*>(vp->getObject());

    pcHole->ThreadSize.setValue(index);
    recomputeFeature();
}

void TaskHoleParameters::threadClassChanged(int index)
{
    if (index < 0)
        return;

    PartDesign::Hole* pcHole = static_cast<PartDesign::Hole*>(vp->getObject());

    pcHole->ThreadClass.setValue(index);
    recomputeFeature();
}

void TaskHoleParameters::threadDiameterChanged(double value)
{
    PartDesign::Hole* pcHole = static_cast<PartDesign::Hole*>(vp->getObject());

    pcHole->Diameter.setValue(value);
    recomputeFeature();
}

void TaskHoleParameters::threadFitChanged(int index)
{
    PartDesign::Hole* pcHole = static_cast<PartDesign::Hole*>(vp->getObject());

    pcHole->ThreadFit.setValue(index);
    recomputeFeature();
}


void TaskHoleParameters::changeEvent(QEvent *e)
{
    TaskBox::changeEvent(e);
    if (e->type() == QEvent::LanguageChange) {
        ui->retranslateUi(proxy);
    }
}

void TaskHoleParameters::changedObject(const App::Document&, const App::Property &Prop)
{
    // happens when aborting the command
    if (vp == nullptr)
        return;
    
    PartDesign::Hole* pcHole = static_cast<PartDesign::Hole*>(vp->getObject());
    bool ro = Prop.isReadOnly();

    Base::Console().Log("Parameter %s was updated\n", Prop.getName());
    
    if (&Prop == &pcHole->ThreadPitch) {
        if (ui->ThreadPitch->value().getValue() != pcHole->ThreadPitch.getValue()) {
            ui->ThreadPitch->blockSignals(true);
            ui->ThreadPitch->setValue(pcHole->ThreadPitch.getValue());
            ui->ThreadPitch->blockSignals(false);
        }
        ui->ThreadPitch->setDisabled(ro);
    }
    else if (&Prop == &pcHole->ThreadAngle) {
        if (ui->ThreadAngle->value().getValue() != pcHole->ThreadAngle.getValue()) {
            ui->ThreadAngle->blockSignals(true);
            ui->ThreadAngle->setValue(pcHole->ThreadAngle.getValue());
            ui->ThreadAngle->blockSignals(false);
        }
        ui->ThreadAngle->setDisabled(ro);
    }
    else if (&Prop == &pcHole->ThreadCutOffInner) {
        if (ui->ThreadCutOffInner->value().getValue() != pcHole->ThreadCutOffInner.getValue()) {
            ui->ThreadCutOffInner->blockSignals(true);
            ui->ThreadCutOffInner->setValue(pcHole->ThreadCutOffInner.getValue());
            ui->ThreadCutOffInner->blockSignals(false);
        }
        ui->ThreadCutOffInner->setDisabled(ro);
    }
    else if (&Prop == &pcHole->ThreadCutOffOuter) {
        if (ui->ThreadCutOffOuter->value().getValue() != pcHole->ThreadCutOffOuter.getValue()) {
            ui->ThreadCutOffOuter->blockSignals(true);
            ui->ThreadCutOffOuter->setValue(pcHole->ThreadCutOffOuter.getValue());
            ui->ThreadCutOffOuter->blockSignals(false);
        }
        ui->ThreadCutOffOuter->setDisabled(ro);
    }
    else if (&Prop == &pcHole->ThreadSize) {
        ui->ThreadSize->setEnabled(true);
        if (ui->ThreadSize->currentIndex() != pcHole->ThreadSize.getValue()) {
            ui->ThreadSize->blockSignals(true);
            ui->ThreadSize->setCurrentIndex(pcHole->ThreadSize.getValue());
            ui->ThreadSize->blockSignals(false);
        }
        ui->ThreadSize->setDisabled(ro);
    }
    else if (&Prop == &pcHole->ThreadClass) {
        ui->ThreadClass->setEnabled(true);
        if (ui->ThreadClass->currentIndex() != pcHole->ThreadClass.getValue()) {
            ui->ThreadClass->blockSignals(true);
            ui->ThreadClass->setCurrentIndex(pcHole->ThreadClass.getValue());
            ui->ThreadClass->blockSignals(false);
        }
        ui->ThreadClass->setDisabled(ro);
    }
    else if (&Prop == &pcHole->ThreadFit) {
        ui->ThreadFit->setEnabled(true);
        if (ui->ThreadFit->currentIndex() != pcHole->ThreadFit.getValue()) {
            ui->ThreadFit->blockSignals(true);
            ui->ThreadFit->setCurrentIndex(pcHole->ThreadFit.getValue());
            ui->ThreadFit->blockSignals(false);
        }
        ui->ThreadFit->setDisabled(ro);
    }
    else if (&Prop == &pcHole->Diameter) {
        ui->Diameter->setEnabled(true);
        if (ui->Diameter->value().getValue() != pcHole->Diameter.getValue()) {
            ui->Diameter->blockSignals(true);
            ui->Diameter->setValue(pcHole->Diameter.getValue());
            ui->Diameter->blockSignals(false);
        }
        ui->Diameter->setDisabled(ro);
    }
    else if (&Prop == &pcHole->HoleCutType) {
        ui->HoleCutType->setEnabled(true);
        if (ui->HoleCutType->currentIndex() != pcHole->HoleCutType.getValue()) {
            ui->HoleCutType->blockSignals(true);
            ui->HoleCutType->setCurrentIndex(pcHole->HoleCutType.getValue());
            ui->HoleCutType->blockSignals(false);
        }
        ui->HoleCutType->setDisabled(ro);
    }
    else if (&Prop == &pcHole->HoleCutDiameter) {
        ui->HoleCutDiameter->setEnabled(true);
        if (ui->HoleCutDiameter->value().getValue() != pcHole->HoleCutDiameter.getValue()) {
            ui->HoleCutDiameter->blockSignals(true);
            ui->HoleCutDiameter->setValue(pcHole->HoleCutDiameter.getValue());
            ui->HoleCutDiameter->blockSignals(false);
        }
        ui->HoleCutDiameter->setDisabled(ro);
    }
    else if (&Prop == &pcHole->HoleCutDepth) {
        ui->HoleCutDepth->setEnabled(true);
        if (ui->HoleCutDepth->value().getValue() != pcHole->HoleCutDepth.getValue()) {
            ui->HoleCutDepth->blockSignals(true);
            ui->HoleCutDepth->setValue(pcHole->HoleCutDepth.getValue());
            ui->HoleCutDepth->blockSignals(false);
        }
        ui->HoleCutDepth->setDisabled(ro);
    }
    else if (&Prop == &pcHole->HoleCutCountersinkAngle) {
        ui->HoleCutCountersinkAngle->setEnabled(true);
        if (ui->HoleCutCountersinkAngle->value().getValue() != pcHole->HoleCutCountersinkAngle.getValue()) {
            ui->HoleCutCountersinkAngle->blockSignals(true);
            ui->HoleCutCountersinkAngle->setValue(pcHole->HoleCutCountersinkAngle.getValue());
            ui->HoleCutCountersinkAngle->blockSignals(false);
        }
        ui->HoleCutCountersinkAngle->setDisabled(ro);
    }
    else if (&Prop == &pcHole->DepthType) {
        ui->DepthType->setEnabled(true);
        if (ui->DepthType->currentIndex() != pcHole->DepthType.getValue()) {
            ui->DepthType->blockSignals(true);
            ui->DepthType->setCurrentIndex(pcHole->DepthType.getValue());
            ui->DepthType->blockSignals(false);
        }
        ui->DepthType->setDisabled(ro);
    }
    else if (&Prop == &pcHole->Depth) {
        ui->Depth->setEnabled(true);
        if (ui->Depth->value().getValue() != pcHole->Depth.getValue()) {
            ui->Depth->blockSignals(true);
            ui->Depth->setValue(pcHole->Depth.getValue());
            ui->Depth->blockSignals(false);
        }
        ui->Depth->setDisabled(ro);
    }
    else if (&Prop == &pcHole->DrillPointAngle) {
        ui->DrillPointAngle->setEnabled(true);
        if (ui->DrillPointAngle->value().getValue() != pcHole->DrillPointAngle.getValue()) {
            ui->DrillPointAngle->blockSignals(true);
            ui->DrillPointAngle->setValue(pcHole->DrillPointAngle.getValue());
            ui->DrillPointAngle->blockSignals(false);
        }
        ui->DrillPointAngle->setDisabled(ro);
    }
    else if (&Prop == &pcHole->TaperedAngle) {
        ui->TaperedAngle->setEnabled(true);
        if (ui->TaperedAngle->value().getValue() != pcHole->TaperedAngle.getValue()) {
            ui->TaperedAngle->blockSignals(true);
            ui->TaperedAngle->setValue(pcHole->TaperedAngle.getValue());
            ui->TaperedAngle->blockSignals(false);
        }
        ui->TaperedAngle->setDisabled(ro);
    }
}

void TaskHoleParameters::onSelectionChanged(const Gui::SelectionChanges& msg)
{
    Q_UNUSED(msg)
}

long   TaskHoleParameters::getThreadSize() const
{
    if ( ui->ThreadSize->currentIndex() == -1 )
        return 0;
    else
        return ui->ThreadSize->currentIndex();
}

long   TaskHoleParameters::getThreadClass() const
{
    if ( ui->ThreadSize->currentIndex() == -1 )
        return 0;
    else
        return ui->ThreadClass->currentIndex();
}


Base::Quantity TaskHoleParameters::getDiameter() const
{
    return ui->Diameter->value();
}

long   TaskHoleParameters::getHoleCutType() const
{
    if (ui->HoleCutType->currentIndex() == -1)
        return 0;
    else
        return ui->HoleCutType->currentIndex();
}

Base::Quantity TaskHoleParameters::getHoleCutDiameter() const
{
    return ui->HoleCutDiameter->value();
}

Base::Quantity TaskHoleParameters::getHoleCutDepth() const
{
    return ui->HoleCutDepth->value();
}

Base::Quantity TaskHoleParameters::getHoleCutCountersinkAngle() const
{
    return ui->HoleCutCountersinkAngle->value();
}

long   TaskHoleParameters::getDepthType() const
{
    return ui->DepthType->currentIndex();
}

Base::Quantity TaskHoleParameters::getDepth() const
{
    return ui->Depth->value();
}

Base::Quantity TaskHoleParameters::getDrillPointAngle() const
{
    return ui->DrillPointAngle->value();
}

Base::Quantity TaskHoleParameters::getTaperedAngle() const
{
    return ui->TaperedAngle->value();
}

void TaskHoleParameters::apply()
{
    auto obj = vp->getObject();
    PartDesign::Hole* pcHole = static_cast<PartDesign::Hole*>(vp->getObject());    

    isApplying = true;
            
    ui->ThreadPitch->apply();
    ui->ThreadAngle->apply();
    ui->ThreadCutOffInner->apply();
    ui->ThreadCutOffOuter->apply();
    ui->Diameter->apply();
    ui->HoleCutDiameter->apply();
    ui->HoleCutDepth->apply();
    ui->HoleCutCountersinkAngle->apply();
    ui->Depth->apply();
    ui->DrillPointAngle->apply();
    ui->TaperedAngle->apply();
   
    if (!pcHole->ThreadSize.isReadOnly())
        FCMD_OBJ_CMD(obj,"ThreadSize = " << getThreadSize());
    if (!pcHole->ThreadClass.isReadOnly())
        FCMD_OBJ_CMD(obj,"ThreadClass = " << getThreadClass());
    if (!pcHole->HoleCutType.isReadOnly())
        FCMD_OBJ_CMD(obj,"HoleCutType = " << getHoleCutType());
    if (!pcHole->DepthType.isReadOnly())
        FCMD_OBJ_CMD(obj,"DepthType = " << getDepthType());
    
    isApplying = false;
}

//**************************************************************************
//**************************************************************************
// TaskDialog
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

TaskDlgHoleParameters::TaskDlgHoleParameters(ViewProviderHole *HoleView)
    : TaskDlgSketchBasedParameters(HoleView)
{
    assert(HoleView);
    parameter  = new TaskHoleParameters(static_cast<ViewProviderHole*>(vp));

    Content.push_back(parameter);
}

TaskDlgHoleParameters::~TaskDlgHoleParameters()
{

}

#include "moc_TaskHoleParameters.cpp"

TaskHoleParameters::Observer::Observer(TaskHoleParameters *_owner, PartDesign::Hole * _hole)
    : DocumentObserver(_hole->getDocument())
    , owner(_owner)
    , hole(_hole)
{
}

void TaskHoleParameters::Observer::slotChangedObject(const App::DocumentObject &Obj, const App::Property &Prop)
{
    if (&Obj == hole) {
        Base::Console().Log("Parameter %s was updated with a new value\n", Prop.getName());
		if (Obj.getDocument())
			owner->changedObject(*Obj.getDocument(), Prop);
    }
}
