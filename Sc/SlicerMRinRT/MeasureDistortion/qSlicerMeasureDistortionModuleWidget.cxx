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
#include <QFileDialog>


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

#include <vtkXMLPolyDataWriter.h>





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
  connect(d->CalculateButton, SIGNAL(clicked()),
	  this, SLOT(CalculateDistortionClick()));
  connect(d->GenerateReferenceButton, SIGNAL(clicked()),
	  this, SLOT(CalculateReferenceClick()));
  connect(d->LoadReferenceButton, SIGNAL(clicked()),
	  this, SLOT(LoadReferenceClick()));
  connect(d->CTVolumeNodeSelector, SIGNAL(currentNodeChanged(vtkMRMLNode*)),
	  this, SLOT(CTSelectionChanged(vtkMRMLNode*)));
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
void qSlicerMeasureDistortionModuleWidget::CTSelectionChanged(vtkMRMLNode*)
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
			CTNode = d->CTVolumeNodeSelector->currentNode();
			selectionNode->SetReferenceActiveVolumeID(CTNode->GetID());
			appLogic->PropagateVolumeSelection();
			//CTScalarVolumeNode->Copy(CTNode);
			//qDebug() << CTScalarVolumeNode;
		}
	}
//	qDebug() << CTNode;

//	qDebug() << CTvolumeNodeID;

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
			MR1Node = d->MRVolumeNodeSelector1->currentNode();
			selectionNode->SetReferenceActiveVolumeID(MR1Node->GetID());
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
			MR2Node = d->MRVolumeNodeSelector2->currentNode();
			selectionNode->SetReferenceActiveVolumeID(MR2Node->GetID());
			appLogic->PropagateVolumeSelection();
		}
	}
}
//-----------------------------------------------------------------------------
void qSlicerMeasureDistortionModuleWidget::LoadReferenceClick()
{
	Q_D(qSlicerMeasureDistortionModuleWidget); 
	QStringList path = QFileDialog::getOpenFileNames(this, tr("Open File"), "C:/S4R/Slicer-build", tr("Points Files (*.vtp)"));
	d->CurrentReferenceList->addItems(path);
}

//-------------------------------------------------------------
void qSlicerMeasureDistortionModuleWidget::CalculateReferenceClick()
{
	Q_D(qSlicerMeasureDistortionModuleWidget);
	vtkSlicerMeasureDistortionLogic* DistortionLogic;
	vtkPolyData* CTpolydata;
	ReferenceNode = DistortionLogic->CalculateReference(CTNode);

	vtkSlicerApplicationLogic *appLogic = this->module()->appLogic();
	vtkMRMLSelectionNode *selectionNode = appLogic->GetSelectionNode();
	//CTNode = d->CTVolumeNodeSelector->currentNode();
	selectionNode->SetReferenceActiveVolumeID(ReferenceNode->GetID());
	//selectionNode->SetActiveVolumeID(ReferenceNode->GetID());
	appLogic->PropagateVolumeSelection();
}
//-------------------------------------------------------------
void qSlicerMeasureDistortionModuleWidget::CalculateDistortionClick()
{
	Q_D(qSlicerMeasureDistortionModuleWidget);
	vtkSlicerMeasureDistortionLogic* DistortionLogic;
	//	qDebug() << CTNode;
	QList<QListWidgetItem *> ReferenceFile = d->CurrentReferenceList->selectedItems();
	if (ReferenceFile.size() != 1) {
		qDebug("You must select 1 reference file");
		return;
	}
	GNLDistortionNode = DistortionLogic->CalculateDistortion(MR1Node,MR2Node);
	//GNLDistortionNode = DistortionLogic->CalculateDistortion(MR1Node, MR2Node);
	//	qDebug() << ReferenceNode;
	qDebug() << "test";
	//Display Distortion Map
	vtkSlicerApplicationLogic *appLogic = this->module()->appLogic();
	vtkMRMLSelectionNode *selectionNode = appLogic->GetSelectionNode();
	selectionNode->SetReferenceActiveVolumeID(GNLDistortionNode->GetID());
	appLogic->PropagateVolumeSelection();
}
//-------------------------------------------------------------