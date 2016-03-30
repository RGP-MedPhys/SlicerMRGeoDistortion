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

#ifndef __qSlicerMeasureDistortionModule_h
#define __qSlicerMeasureDistortionModule_h

// SlicerQt includes
#include "qSlicerLoadableModule.h"

#include "qSlicerMeasureDistortionModuleExport.h"

class QSettings;

class qSlicerMeasureDistortionModulePrivate;

/// \ingroup Slicer_QtModules_ExtensionTemplate
class Q_SLICER_QTMODULES_MEASUREDISTORTION_EXPORT
qSlicerMeasureDistortionModule
  : public qSlicerLoadableModule
{
  Q_OBJECT
  Q_INTERFACES(qSlicerLoadableModule);

public:

  typedef qSlicerLoadableModule Superclass;
  explicit qSlicerMeasureDistortionModule(QObject *parent=0);
  virtual ~qSlicerMeasureDistortionModule();

  qSlicerGetTitleMacro(QTMODULE_TITLE);

  virtual QString helpText()const;
  virtual QString acknowledgementText()const;
  virtual QStringList contributors()const;

  virtual QIcon icon()const;

  virtual QStringList categories()const;
  virtual QStringList dependencies() const;

  /// Read default slice view settings from application settings (.ini file)
  /// into defaultViewNode.
  static void readDefaultSliceViewSettings(vtkMRMLSliceNode* defaultViewNode);

  /// Write default slice view settings to application settings (.ini file)
  /// from defaultViewNode.
  static void writeDefaultSliceViewSettings(vtkMRMLSliceNode* defaultViewNode);

  /// Set MRML scene for the module. Updates the default view settings based on
  /// the application settings.
  virtual void setMRMLScene(vtkMRMLScene* scene);

protected:

  /// Initialize the module. Register the volumes reader/writer
  virtual void setup();

  /// Create and return the widget representation associated to this module
  virtual qSlicerAbstractModuleRepresentation * createWidgetRepresentation();

  /// Create and return the logic associated to this module
  virtual vtkMRMLAbstractLogic* createLogic();

  /// Helper functions to read/write common view settings
  static void readCommonViewSettings(vtkMRMLAbstractViewNode* defaultViewNode, QSettings& settings);
  static void writeCommonViewSettings(vtkMRMLAbstractViewNode* defaultViewNode, QSettings& settings);

protected:
  QScopedPointer<qSlicerMeasureDistortionModulePrivate> d_ptr;

private:
  Q_DECLARE_PRIVATE(qSlicerMeasureDistortionModule);
  Q_DISABLE_COPY(qSlicerMeasureDistortionModule);

};

#endif
