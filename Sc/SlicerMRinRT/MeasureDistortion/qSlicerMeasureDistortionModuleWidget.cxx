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
#include "qSlicerApplication.h"
#include "qSlicerLayoutManager.h"
#include "qSlicerMeasureDistortionModuleWidget.h"
#include "qSlicerModuleManager.h"
#include "ui_qSlicerMeasureDistortionModuleWidget.h"

// CTK includes
#include "ctkButtonGroup.h"

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_ExtensionTemplate
class qSlicerMeasureDistortionModuleWidgetPrivate: public Ui_qSlicerMeasureDistortionModuleWidget
{
	Q_DECLARE_PUBLIC(qSlicerMeasureDistortionModuleWidget);
protected:
	qSlicerMeasureDistortionModuleWidget* const q_ptr;
public:
	qSlicerMeasureDistortionModuleWidgetPrivate(qSlicerMeasureDistortionModuleWidget& object);

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