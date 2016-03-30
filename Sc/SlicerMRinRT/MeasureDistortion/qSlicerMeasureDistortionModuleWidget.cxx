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

// MRMLWidgets includes
#include <qMRMLSliceWidget.h>
#include <qMRMLSliceControllerWidget.h>

// SlicerQt includes
#include "qSlicerApplication.h"
#include "qSlicerLayoutManager.h"
#include "qSlicerMeasureDistortionModuleWidget.h"
#include "qSlicerModuleManager.h"
#include "ui_qSlicerMeasureDistortionModuleWidget.h"

// MRML includes
#include "vtkMRMLScene.h"
#include "vtkMRMLSliceNode.h"
#include "vtkMRMLViewNode.h"

// MRMLLogic includes
#include "vtkMRMLLayoutLogic.h"

// CTK includes
#include "ctkButtonGroup.h"

// VTK includes
#include <vtkCollection.h>
#include <vtkSmartPointer.h>

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_ExtensionTemplate
class qSlicerMeasureDistortionModuleWidgetPrivate: public Ui_qSlicerMeasureDistortionModuleWidget
{
	Q_DECLARE_PUBLIC(qSlicerMeasureDistortionModuleWidget);
protected:
	qSlicerMeasureDistortionModuleWidget* const q_ptr;
public:
	qSlicerMeasureDistortionModuleWidgetPrivate(qSlicerMeasureDistortionModuleWidget& object);
	virtual ~qSlicerMeasureDistortionModuleWidgetPrivate();

  bool selectModule(const QString& moduleName);

  /// Create a Controller for a Node and pack in the widget
  void createController(vtkMRMLNode *n, qSlicerLayoutManager *lm);

  /// Remove the Controller for a Node from the widget
  void removeController(vtkMRMLNode *n);

  typedef std::map<vtkSmartPointer<vtkMRMLNode>, qMRMLViewControllerBar* > ControllerMapType;
  ControllerMapType ControllerMap;

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
void qSlicerMeasureDistortionModuleWidgetPrivate::createController(vtkMRMLNode *n, qSlicerLayoutManager *layoutManager)
{
	if (this->ControllerMap.find(n) != this->ControllerMap.end())
	{
		qDebug() << "qSlicerMeasureDistortionModuleWidgetPrivate::createController - Node already added to module";
		return;
	}

	// create the ControllerWidget and wire it to the appropriate node
	qMRMLViewControllerBar *barWidget = 0;
	vtkMRMLSliceNode *sn = vtkMRMLSliceNode::SafeDownCast(n);
	if (sn)
	{
		qMRMLSliceControllerWidget *widget =
			new qMRMLSliceControllerWidget(this->SliceControllersCollapsibleButton);
		widget->setSliceViewName(sn->GetName()); // call before setting slice node
		widget->setSliceViewLabel(sn->GetLayoutLabel());
		QColor layoutColor = QColor::fromRgbF(sn->GetLayoutColor()[0],
			sn->GetLayoutColor()[1],
			sn->GetLayoutColor()[2]);
		widget->setSliceViewColor(layoutColor);
		widget->setMRMLSliceNode(sn);
		widget->setLayoutBehavior(qMRMLViewControllerBar::Panel);

		// SliceControllerWidget needs to know the SliceLogic(s)
		qMRMLSliceWidget *sliceWidget = layoutManager->sliceWidget(sn->GetLayoutName());
		widget->setSliceLogics(layoutManager->mrmlSliceLogics());
		widget->setSliceLogic(sliceWidget->sliceController()->sliceLogic());

		// add the widget to the display
		this->SliceControllersLayout->addWidget(widget);

		barWidget = widget;
	}

	vtkMRMLViewNode *vn = vtkMRMLViewNode::SafeDownCast(n);
	if (vn)
	{
//		qMRMLThreeDViewControllerWidget *widget =
//			new qMRMLThreeDViewControllerWidget(this->ThreeDViewControllersCollapsibleButton);
//		widget->setViewLabel(vn->GetLayoutLabel());
//		widget->setMRMLViewNode(vn);
//		widget->setLayoutBehavior(qMRMLViewControllerBar::Panel);

		// ThreeDViewController needs to now the ThreeDView
//		qMRMLThreeDWidget *viewWidget = dynamic_cast<qMRMLThreeDWidget*>(layoutManager->viewWidget(vn));
//		if (viewWidget)
//		{
//			widget->setThreeDView(viewWidget->threeDView());
//		}

		// add the widget to the display
//		this->ThreeDViewControllersLayout->addWidget(widget);

//		barWidget = widget;
	}

	// cache the widget. we'll clean this up on the NodeRemovedEvent
	this->ControllerMap[n] = barWidget;
}


//-----------------------------------------------------------------------------
void qSlicerMeasureDistortionModuleWidgetPrivate::removeController(vtkMRMLNode *n)
{
	// find the widget for the SliceNode
	ControllerMapType::iterator cit = this->ControllerMap.find(n);
	if (cit == this->ControllerMap.end())
	{
		qDebug() << "qSlicerMeasureDistortionModuleWidgetPrivate::removeController - Node has no Controller managed by this module.";
		return;
	}

	// unpack the widget
	vtkMRMLSliceNode *sn = vtkMRMLSliceNode::SafeDownCast(n);
	if (sn)
	{
		SliceControllersLayout->removeWidget((*cit).second);
	}

	vtkMRMLViewNode *vn = vtkMRMLViewNode::SafeDownCast(n);
	if (vn)
	{
//		ThreeDViewControllersLayout->removeWidget((*cit).second);
	}

	// delete the widget
	delete (*cit).second;

	// remove entry from the map
	this->ControllerMap.erase(cit);
}


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
void qSlicerMeasureDistortionModuleWidget::setup()
{
  Q_D(qSlicerMeasureDistortionModuleWidget);
  d->setupUi(this);

  connect(d->LoadDicomDataButton, SIGNAL(clicked()),
	  this, SLOT(loadDicomData()));

#ifndef Slicer_BUILD_DICOM_SUPPORT
  d->LoadDicomDataButton->setDisabled(true);
#endif

  this->Superclass::setup();
}

//-----------------------------------------------------------------------------
bool qSlicerMeasureDistortionModuleWidget::loadDicomData()
{
	Q_D(qSlicerMeasureDistortionModuleWidget);
	return d->selectModule("DICOM");
}

//-----------------------------------------------------------------------------
void qSlicerMeasureDistortionModuleWidget::setMRMLScene(vtkMRMLScene *newScene)
{
	Q_D(qSlicerMeasureDistortionModuleWidget);

	vtkMRMLScene* oldScene = this->mrmlScene();

	this->Superclass::setMRMLScene(newScene);

	qSlicerApplication * app = qSlicerApplication::application();
	if (!app)
	{
		return;
	}
	qSlicerLayoutManager * layoutManager = app->layoutManager();
	if (!layoutManager)
	{
		return;
	}

	// Search the scene for the available view nodes and create a
	// Controller and connect it up
	newScene->InitTraversal();
	for (vtkMRMLNode *sn = NULL; (sn = newScene->GetNextNodeByClass("vtkMRMLSliceNode"));)
	{
		vtkMRMLSliceNode *snode = vtkMRMLSliceNode::SafeDownCast(sn);
		if (snode)
		{
			d->createController(snode, layoutManager);
		}
	}

	newScene->InitTraversal();
	for (vtkMRMLNode *sn = NULL; (sn = newScene->GetNextNodeByClass("vtkMRMLViewNode"));)
	{
		vtkMRMLViewNode *vnode = vtkMRMLViewNode::SafeDownCast(sn);
		if (vnode)
		{
			d->createController(vnode, layoutManager);
		}
	}

	// Need to listen for any new slice or view nodes being added
	this->qvtkReconnect(oldScene, newScene, vtkMRMLScene::NodeAddedEvent,
		this, SLOT(onNodeAddedEvent(vtkObject*, vtkObject*)));

	// Need to listen for any slice or view nodes being removed
	this->qvtkReconnect(oldScene, newScene, vtkMRMLScene::NodeRemovedEvent,
		this, SLOT(onNodeRemovedEvent(vtkObject*, vtkObject*)));

	// Listen to changes in the Layout so we only show controllers for
	// the visible nodes
	QObject::connect(layoutManager, SIGNAL(layoutChanged(int)), this,
		SLOT(onLayoutChanged(int)));

}

// --------------------------------------------------------------------------
void qSlicerMeasureDistortionModuleWidget::onNodeAddedEvent(vtkObject*, vtkObject* node)
{
	Q_D(qSlicerMeasureDistortionModuleWidget);

	if (!this->mrmlScene() || this->mrmlScene()->IsBatchProcessing())
	{
		return;
	}

	qSlicerApplication * app = qSlicerApplication::application();
	if (!app)
	{
		return;
	}
	qSlicerLayoutManager * layoutManager = app->layoutManager();
	if (!layoutManager)
	{
		return;
	}

	vtkMRMLSliceNode* sliceNode = vtkMRMLSliceNode::SafeDownCast(node);
	if (sliceNode)
	{

		// create the slice controller
		d->createController(sliceNode, layoutManager);
	}

	vtkMRMLViewNode* viewNode = vtkMRMLViewNode::SafeDownCast(node);
	if (viewNode)
	{
	
		// create the view controller
		d->createController(viewNode, layoutManager);
	}
}

// --------------------------------------------------------------------------
void qSlicerMeasureDistortionModuleWidget::onNodeRemovedEvent(vtkObject*, vtkObject* node)
{
	Q_D(qSlicerMeasureDistortionModuleWidget);

	if (!this->mrmlScene() || this->mrmlScene()->IsBatchProcessing())
	{
		return;
	}

	vtkMRMLSliceNode* sliceNode = vtkMRMLSliceNode::SafeDownCast(node);
	if (sliceNode)
	{
		QString layoutName = sliceNode->GetLayoutName();
		qDebug() << "qSlicerMeasureDistortionModuleModuleWidget::onNodeRemovedEvent - layoutName:" << layoutName;

		// destroy the slice controller
		d->removeController(sliceNode);
	}

	vtkMRMLViewNode* viewNode = vtkMRMLViewNode::SafeDownCast(node);
	if (sliceNode)
	{
		QString layoutName = viewNode->GetName();
		qDebug() << "qSlicerMeasureDistortionModuleModuleWidget::onNodeRemovedEvent - layoutName:" << layoutName;

		// destroy the view controller
		d->removeController(viewNode);
	}
}

// --------------------------------------------------------------------------
void qSlicerMeasureDistortionModuleWidget::onLayoutChanged(int)
{
	Q_D(qSlicerMeasureDistortionModuleWidget);

	if (!this->mrmlScene() || this->mrmlScene()->IsBatchProcessing())
	{
		return;
	}


	// add the controllers for any newly visible SliceNodes and remove
	// the controllers for any SliceNodes no longer visible

	qSlicerApplication * app = qSlicerApplication::application();
	if (!app)
	{
		return;
	}
	qSlicerLayoutManager * layoutManager = app->layoutManager();
	if (!layoutManager)
	{
		return;
	}

	vtkMRMLLayoutLogic *layoutLogic = layoutManager->layoutLogic();
	vtkCollection *visibleViews = layoutLogic->GetViewNodes();
	vtkObject *v;

	// hide Controllers for Nodes not currently visible in
	// the layout
	qSlicerMeasureDistortionModuleWidgetPrivate::ControllerMapType::iterator cit;
	for (cit = d->ControllerMap.begin(); cit != d->ControllerMap.end(); ++cit)
	{
		// is mananaged Node not currently displayed in the layout?
		if (!visibleViews->IsItemPresent((*cit).first))
		{
			// hide it
			(*cit).second->hide();
		}
	}

	// show Controllers for Nodes not currently being managed
	// by this widget
	for (visibleViews->InitTraversal(); (v = visibleViews->GetNextItemAsObject());)
	{
		vtkMRMLNode *vn = vtkMRMLNode::SafeDownCast(v);
		if (vn)
		{
			// find the controller
			qSlicerMeasureDistortionModuleWidgetPrivate::ControllerMapType::iterator cit = d->ControllerMap.find(vn);
			if (cit != d->ControllerMap.end())
			{
				// show it
				(*cit).second->show();
			}
		}
	}
}

