/*==============================================================================

  Program: 3D Slicer

  Portions (c) Copyright Brigham and Women's Hospital (BWH) All Rights Reserved.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

==============================================================================*/

// Qt includes
#include <QDebug>
#include <QMessageBox>

// SlicerQt includes
#include <qSlicerAbstractCoreModule.h>
#include "qSlicerApplication.h"
#include "qSlicerLayoutManager.h"
#include "qSlicerMeasureDistortionModuleWidget.h"
#include "qSlicerModuleManager.h"
#include "ui_qSlicerMeasureDistortionModuleWidget.h"


// CTK includes
#include "ctkButtonGroup.h"

// MRML includes
#include <vtkMRMLScalarVolumeNode.h>
#include <vtkMRMLSliceNode.h>
#include <vtkimagedata.h>
#include <vtkMRMLScene.h>
#include <vtkMRMLNode.h>
#include <vtkMRMLCropVolumeParametersNode.h>
#include <vtkMRMLVolumeNode.h>
#include <vtkMRMLSelectionNode.h>
#include <vtkMRMLTransformNode.h>
#include <vtkMRMLApplicationLogic.h>
#include "vtkSlicerMeasureDistortionLogic.h"




//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_Measure Distortion
class qSlicerMeasureDistortionModuleWidgetPrivate: public Ui_qSlicerMeasureDistortionModuleWidget
{
	Q_DECLARE_PUBLIC(qSlicerMeasureDistortionModuleWidget);
protected:
	qSlicerMeasureDistortionModuleWidget* const q_ptr;
public:
	qSlicerMeasureDistortionModuleWidgetPrivate(qSlicerMeasureDistortionModuleWidget& object);
	vtkSlicerMeasureDistortionLogic* logic() const;
  bool selectModule(const QString& moduleName);
};

//-----------------------------------------------------------------------------
// qSlicerMeasureDistortionModuleWidgetPrivate methods

//-----------------------------------------------------------------------------
qSlicerMeasureDistortionModuleWidgetPrivate::qSlicerMeasureDistortionModuleWidgetPrivate(qSlicerMeasureDistortionModuleWidget& object)
	: q_ptr(&object)
{
}

bool qSlicerMeasureDistortionModuleWidgetPrivate::selectModule(const QString& moduleName)
{
	Q_Q(qSlicerMeasureDistortionModuleWidget);
	qSlicerModuleManager * moduleManager = qSlicerCoreApplication::application()->moduleManager();
	if (!moduleManager)
	{
		return false;
	}
	qSlicerAbstractCoreModule * module = moduleManager->module(moduleName);
	if (!module)
	{
		QMessageBox::warning(
			q, q->tr("Raising %1 Module:").arg(moduleName),
			q->tr("Unfortunately, this requested module is not available in this Slicer session."),
			QMessageBox::Ok);
		return false;
	}
	qSlicerLayoutManager * layoutManager = qSlicerApplication::application()->layoutManager();
	if (!layoutManager)
	{
		return false;
	}
	layoutManager->setCurrentModule(moduleName);
	return true;
}

//-----------------------------------------------------------------------------
// qSlicerMeasureDistortionModuleWidget methods

//-----------------------------------------------------------------------------
qSlicerMeasureDistortionModuleWidget::qSlicerMeasureDistortionModuleWidget(QWidget* _parent)
  : Superclass( _parent )
  , d_ptr(new qSlicerMeasureDistortionModuleWidgetPrivate(*this))
{
}

//-----------------------------------------------------------------------------
qSlicerMeasureDistortionModuleWidget::~qSlicerMeasureDistortionModuleWidget()
{
}


//-----------------------------------------------------------------------------
vtkSlicerMeasureDistortionLogic* qSlicerMeasureDistortionModuleWidgetPrivate::logic() const
{
	Q_Q(const qSlicerMeasureDistortionModuleWidget);
	return vtkSlicerMeasureDistortionLogic::SafeDownCast(q->logic());
}



//-----------------------------------------------------------------------------


void qSlicerMeasureDistortionModuleWidget::setup()
{
  Q_D(qSlicerMeasureDistortionModuleWidget);
  d->setupUi(this);

 // connect(d->LoadDicomDataButton, SIGNAL(clicked()),
//	  this, SLOT(loadDicomData()));
  connect(d->CTVolumeNodeSelector, SIGNAL(currentNodeChanged(vtkMRMLNode*)),
	  this, SLOT(OnCTSelectionChanged(vtkMRMLNode*)));
  connect(d->MRVolumeNodeSelector1, SIGNAL(currentNodeChanged(vtkMRMLNode*)),
	  this, SLOT(MR1SelectionChanged(vtkMRMLNode*)));
  connect(d->MRVolumeNodeSelector2, SIGNAL(currentNodeChanged(vtkMRMLNode*)),
	  this, SLOT(MR2SelectionChanged(vtkMRMLNode*)));

  
//  QObject::connect(qSlicerApplication::application(),
//	  SIGNAL(mrmlSceneChanged(vtkMRMLScene*)),
//	  this->LayoutManager, SLOT(setMRMLScene(vtkMRMLScene*)));

#ifndef Slicer_BUILD_DICOM_SUPPORT
  d->LoadDicomDataButton->setDisabled(true);
#endif

  this->Superclass::setup();
}
//------------------------------------------------------------------------------
void qSlicerMeasureDistortionModuleWidget::OnCTSelectionChanged(vtkMRMLNode*)
{
	Q_D(const qSlicerMeasureDistortionModuleWidget);

	vtkMRMLNode* node = d->CTVolumeNodeSelector->currentNode();
	if (node)
	{
		if (this->mrmlScene() &&
			!this->mrmlScene()->IsClosing() &&
			!this->mrmlScene()->IsBatchProcessing())
		{
			// set it to be active in the slice windows
			vtkSlicerApplicationLogic *appLogic = this->module()->appLogic();
			vtkMRMLSelectionNode *selectionNode = appLogic->GetSelectionNode();
			CTvolumeNodeID = d->CTVolumeNodeSelector->currentNode()->GetID();
			selectionNode->SetReferenceActiveVolumeID(CTvolumeNodeID);
			appLogic->PropagateVolumeSelection();
		}
	}


//	qDebug() << CTnode;

//	QMessageBox msgBox;
//	msgBox.setWindowTitle("Info");
//	msgBox.setInformativeText("Test.");
//	msgBox.exec();

}

//------------------------------------------------------------------------------

void qSlicerMeasureDistortionModuleWidget::MR1SelectionChanged(vtkMRMLNode*)
{
	Q_D(const qSlicerMeasureDistortionModuleWidget);
	vtkMRMLNode* node = d->MRVolumeNodeSelector1->currentNode();
	if (node)
	{
		if (this->mrmlScene() &&
			!this->mrmlScene()->IsClosing() &&
			!this->mrmlScene()->IsBatchProcessing())
		{
			// set it to be active in the slice windows
			vtkSlicerApplicationLogic *appLogic = this->module()->appLogic();
			vtkMRMLSelectionNode *selectionNode = appLogic->GetSelectionNode();
			MR1volumeNodeID = d->MRVolumeNodeSelector1->currentNode()->GetID();
			selectionNode->SetReferenceActiveVolumeID(MR1volumeNodeID);
			appLogic->PropagateVolumeSelection();
		}
	}
}

//------------------------------------------------------------------------------
void qSlicerMeasureDistortionModuleWidget::MR2SelectionChanged(vtkMRMLNode*)
{
	Q_D(const qSlicerMeasureDistortionModuleWidget);
	vtkMRMLNode* node = d->MRVolumeNodeSelector2->currentNode();
	if (node)
	{
		if (this->mrmlScene() &&
			!this->mrmlScene()->IsClosing() &&
			!this->mrmlScene()->IsBatchProcessing())
		{
			// set it to be active in the slice windows
			vtkSlicerApplicationLogic *appLogic = this->module()->appLogic();
			vtkMRMLSelectionNode *selectionNode = appLogic->GetSelectionNode();
			MR2volumeNodeID = d->MRVolumeNodeSelector2->currentNode()->GetID();
			selectionNode->SetReferenceActiveVolumeID(MR2volumeNodeID);
			appLogic->PropagateVolumeSelection();
		}
	}
}
//-----------------------------------------------------------------------------
//bool qSlicerMeasureDistortionModuleWidget::loadDicomData()
//{
//	Q_D(qSlicerMeasureDistortionModuleWidget);
//	return d->selectModule("DICOM");
//}
//-------------------------------------------------------------
