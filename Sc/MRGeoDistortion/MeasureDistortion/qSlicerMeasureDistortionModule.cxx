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

// Slicer includes
#include <qSlicerCoreApplication.h>
#include <qSlicerModuleManager.h>

// Qt includes
#include <QDebug>
#include <QtPlugin>
 
// MeasureDistortion Logic includes
#include <vtkSlicerMeasureDistortionLogic.h>
#include <vtkSlicerVolumesLogic.h>

// MeasureDistortion includes
#include "qSlicerMeasureDistortionModule.h"
#include "qSlicerMeasureDistortionModuleWidget.h"

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
  return QStringList() << "MRGeoDistortion";
}

//-----------------------------------------------------------------------------
QStringList qSlicerMeasureDistortionModule::dependencies() const
{
  return QStringList() << "Volumes";
}

//-----------------------------------------------------------------------------
void qSlicerMeasureDistortionModule::setup()
{
	this->Superclass::setup();

	vtkSlicerMeasureDistortionLogic* MeasureDistortionLogic =
		vtkSlicerMeasureDistortionLogic::SafeDownCast(this->logic());
	qSlicerAbstractCoreModule* volumesModule =
		qSlicerCoreApplication::application()->moduleManager()->module("Volumes");
	if (volumesModule)
	{
		vtkSlicerVolumesLogic* volumesLogic =
			vtkSlicerVolumesLogic::SafeDownCast(volumesModule->logic());
		MeasureDistortionLogic->SetVolumesLogic(volumesLogic);
	}
	else
	{
		qWarning() << "Volumes module is not found";
	}
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
