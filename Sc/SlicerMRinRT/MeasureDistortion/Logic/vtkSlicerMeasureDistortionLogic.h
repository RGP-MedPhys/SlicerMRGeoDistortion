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

// .NAME vtkSlicerMeasureDistortionLogic - slicer logic class for volumes manipulation
// .SECTION Description
// This class manages the logic associated with reading, saving,
// and changing propertied of the volumes


#ifndef __vtkSlicerMeasureDistortionLogic_h
#define __vtkSlicerMeasureDistortionLogic_h

// Slicer includes
#include "vtkSlicerModuleLogic.h"

// MRML includes

// STD includes
#include <cstdlib>

#include "vtkSlicerMeasureDistortionModuleLogicExport.h"

class vtkMRMLSliceNode;
class vtkMRMLViewNode;

/// \ingroup Slicer_QtModules_ExtensionTemplate
class VTK_SLICER_MEASUREDISTORTION_MODULE_LOGIC_EXPORT vtkSlicerMeasureDistortionLogic :
  public vtkSlicerModuleLogic
{
public:

  static vtkSlicerMeasureDistortionLogic *New();
  vtkTypeMacro(vtkSlicerMeasureDistortionLogic, vtkSlicerModuleLogic);
  void PrintSelf(ostream& os, vtkIndent indent);

public:

	/// Retrieves the default slice view node from the scene.
	/// The returned node can be changed to customize the contents of
	/// new view nodes. ResetAllViewNodesToDefault() can be called to
	/// update all existing view nodes to the new default settings immediately.
	vtkMRMLSliceNode* GetDefaultSliceViewNode();

	/// Retrieves the default 3D view node from the scene.
	/// The returned node can be changed to customize the contents of
	/// new view nodes. ResetAllViewNodesToDefault() can be called to
	/// update all existing view nodes to the new default settings immediately.
	vtkMRMLViewNode* GetDefaultThreeDViewNode();

	/// Reset all existing slice and 3D view nodes to default settings.
	void ResetAllViewNodesToDefault();

protected:
  vtkSlicerMeasureDistortionLogic();
  virtual ~vtkSlicerMeasureDistortionLogic();

  virtual void SetMRMLSceneInternal(vtkMRMLScene* newScene);
  /// Register MRML Node classes to Scene. Gets called automatically when the MRMLScene is attached to this logic class.
  virtual void RegisterNodes();
  virtual void UpdateFromMRMLScene();
  virtual void OnMRMLSceneNodeAdded(vtkMRMLNode* node);
  virtual void OnMRMLSceneNodeRemoved(vtkMRMLNode* node);
private:

  vtkSlicerMeasureDistortionLogic(const vtkSlicerMeasureDistortionLogic&); // Not implemented
  void operator=(const vtkSlicerMeasureDistortionLogic&); // Not implemented
};

#endif
