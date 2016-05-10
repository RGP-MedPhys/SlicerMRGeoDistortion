/*==============================================================================

  Program: 3D Slicer

  Copyright (c) Kitware Inc.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

  This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc.
  and was partially funded by NIH grant 3P41RR013218-12S1

==============================================================================*/

// FooBar Widgets includes
#include "qSlicerMeasureDistortionFooBarWidget.h"
#include "ui_qSlicerMeasureDistortionFooBarWidget.h"

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_MeasureDistortion
class qSlicerMeasureDistortionFooBarWidgetPrivate
  : public Ui_qSlicerMeasureDistortionFooBarWidget
{
  Q_DECLARE_PUBLIC(qSlicerMeasureDistortionFooBarWidget);
protected:
  qSlicerMeasureDistortionFooBarWidget* const q_ptr;

public:
  qSlicerMeasureDistortionFooBarWidgetPrivate(
    qSlicerMeasureDistortionFooBarWidget& object);
  virtual void setupUi(qSlicerMeasureDistortionFooBarWidget*);
};

// --------------------------------------------------------------------------
qSlicerMeasureDistortionFooBarWidgetPrivate
::qSlicerMeasureDistortionFooBarWidgetPrivate(
  qSlicerMeasureDistortionFooBarWidget& object)
  : q_ptr(&object)
{
}

// --------------------------------------------------------------------------
void qSlicerMeasureDistortionFooBarWidgetPrivate
::setupUi(qSlicerMeasureDistortionFooBarWidget* widget)
{
  this->Ui_qSlicerMeasureDistortionFooBarWidget::setupUi(widget);
}

//-----------------------------------------------------------------------------
// qSlicerMeasureDistortionFooBarWidget methods

//-----------------------------------------------------------------------------
qSlicerMeasureDistortionFooBarWidget
::qSlicerMeasureDistortionFooBarWidget(QWidget* parentWidget)
  : Superclass( parentWidget )
  , d_ptr( new qSlicerMeasureDistortionFooBarWidgetPrivate(*this) )
{
  Q_D(qSlicerMeasureDistortionFooBarWidget);
  d->setupUi(this);
}

//-----------------------------------------------------------------------------
qSlicerMeasureDistortionFooBarWidget
::~qSlicerMeasureDistortionFooBarWidget()
{
}
