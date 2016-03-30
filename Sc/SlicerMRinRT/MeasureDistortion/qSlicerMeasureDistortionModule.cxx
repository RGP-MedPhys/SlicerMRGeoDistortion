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
#include <QtPlugin>
 
// MeasureDistortion Logic includes
#include <vtkSlicerMeasureDistortionLogic.h>

// MeasureDistortion includes
#include "qSlicerMeasureDistortionModule.h"
#include "qSlicerMeasureDistortionModuleWidget.h"

#include <vtkMRMLSliceNode.h>
#include <vtkMRMLViewNode.h>

//-----------------------------------------------------------------------------
Q_EXPORT_PLUGIN2(qSlicerMeasureDistortionModule, qSlicerMeasureDistortionModule);

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_ExtensionTemplate
class qSlicerMeasureDistortionModulePrivate
{
public:
  qSlicerMeasureDistortionModulePrivate();
};

//-----------------------------------------------------------------------------
// qSlicerMeasureDistortionModulePrivate methods

//-----------------------------------------------------------------------------
qSlicerMeasureDistortionModulePrivate::qSlicerMeasureDistortionModulePrivate()
{
}

//-----------------------------------------------------------------------------
// qSlicerMeasureDistortionModule methods

//-----------------------------------------------------------------------------
qSlicerMeasureDistortionModule::qSlicerMeasureDistortionModule(QObject* _parent)
  : Superclass(_parent)
  , d_ptr(new qSlicerMeasureDistortionModulePrivate) 
{
}

//-----------------------------------------------------------------------------
qSlicerMeasureDistortionModule::~qSlicerMeasureDistortionModule()
{
}

//-----------------------------------------------------------------------------
QString qSlicerMeasureDistortionModule::helpText() const
{
  return "This is a loadable module that can be bundled in an extension";
}

//-----------------------------------------------------------------------------
QString qSlicerMeasureDistortionModule::acknowledgementText() const
{
  return "Research supported in part by a grant from Philips HealthCare and technical research support from GE Healthcare.";
}

//-----------------------------------------------------------------------------
QStringList qSlicerMeasureDistortionModule::contributors() const
{
  QStringList moduleContributors;
  moduleContributors << QString("Ryan G Price (Henry Ford Health System)");
  return moduleContributors;
}

//-----------------------------------------------------------------------------
QIcon qSlicerMeasureDistortionModule::icon() const
{
  return QIcon(":/Icons/MeasureDistortion.png");
}

//-----------------------------------------------------------------------------
QStringList qSlicerMeasureDistortionModule::categories() const
{
  return QStringList() << "MRinRT";
}

//-----------------------------------------------------------------------------
QStringList qSlicerMeasureDistortionModule::dependencies() const
{
  return QStringList() << "Volumes";
}

//-----------------------------------------------------------------------------
void qSlicerMeasureDistortionModule::setup()
{

}

//-----------------------------------------------------------------------------
void qSlicerMeasureDistortionModule::setMRMLScene(vtkMRMLScene* scene)
{
	this->Superclass::setMRMLScene(scene);
	vtkSlicerMeasureDistortionLogic* logic = vtkSlicerMeasureDistortionLogic::SafeDownCast(this->logic());
	if (!logic)
	{
		qCritical() << Q_FUNC_INFO << " failed: logic is invalid";
		return;
	}
	// Update default view nodes from settings
	this->readDefaultSliceViewSettings(logic->GetDefaultSliceViewNode());
	//this->readDefaultThreeDViewSettings(logic->GetDefaultThreeDViewNode());
	this->writeDefaultSliceViewSettings(logic->GetDefaultSliceViewNode());
	//this->writeDefaultThreeDViewSettings(logic->GetDefaultThreeDViewNode());
	// Update all existing view nodes to default
	logic->ResetAllViewNodesToDefault();
}

//-----------------------------------------------------------------------------
qSlicerAbstractModuleRepresentation* qSlicerMeasureDistortionModule
::createWidgetRepresentation()
{
  return new qSlicerMeasureDistortionModuleWidget;
}

//-----------------------------------------------------------------------------
vtkMRMLAbstractLogic* qSlicerMeasureDistortionModule::createLogic()
{
  return vtkSlicerMeasureDistortionLogic::New();
}

//-----------------------------------------------------------------------------
//void qSlicerMeasureDistortionModule::readCommonViewSettings(vtkMRMLAbstractViewNode* defaultViewNode, QSettings& settings)
//{
//	if (settings.contains("OrientationMarkerType"))
//	{
//		defaultViewNode->SetOrientationMarkerType(vtkMRMLAbstractViewNode::GetOrientationMarkerTypeFromString(settings.value("OrientationMarkerType").toString().toLatin1()));
//	}
//	if (settings.contains("OrientationMarkerSize"))
//	{
//		defaultViewNode->SetOrientationMarkerSize(vtkMRMLAbstractViewNode::GetOrientationMarkerSizeFromString(settings.value("OrientationMarkerSize").toString().toLatin1()));
//	}
//	if (settings.contains("RulerType"))
//	{
//		defaultViewNode->SetRulerType(vtkMRMLAbstractViewNode::GetRulerTypeFromString(settings.value("RulerType").toString().toLatin1()));
//	}
//}

//-----------------------------------------------------------------------------
//void qSlicerMeasureDistortionModule::writeCommonViewSettings(vtkMRMLAbstractViewNode* defaultViewNode, QSettings& settings)
//{
//	settings.setValue("OrientationMarkerType", vtkMRMLAbstractViewNode::GetOrientationMarkerTypeAsString(defaultViewNode->GetOrientationMarkerType()));
//	settings.setValue("OrientationMarkerSize", vtkMRMLAbstractViewNode::GetOrientationMarkerSizeAsString(defaultViewNode->GetOrientationMarkerSize()));
//	settings.setValue("RulerType", vtkMRMLAbstractViewNode::GetRulerTypeAsString(defaultViewNode->GetRulerType()));
//}

//-----------------------------------------------------------------------------
//void qSlicerMeasureDistortionModule::readDefaultSliceViewSettings(vtkMRMLSliceNode* defaultViewNode)
//{
//	if (!defaultViewNode)
//	{
//		qCritical() << Q_FUNC_INFO << " failed: defaultViewNode is invalid";
//		return;
//	}
//	QSettings settings;
//	settings.beginGroup("DefaultSliceView");
//	readCommonViewSettings(defaultViewNode, settings);
//}

//-----------------------------------------------------------------------------
//void qSlicerMeasureDistortionModule::writeDefaultSliceViewSettings(vtkMRMLSliceNode* defaultViewNode)
//{
//	if (!defaultViewNode)
//	{
//		qCritical() << Q_FUNC_INFO << " failed: defaultViewNode is invalid";
//		return;
//	}
//	QSettings settings;
//	settings.beginGroup("DefaultSliceView");
//	writeCommonViewSettings(defaultViewNode, settings);
//}
