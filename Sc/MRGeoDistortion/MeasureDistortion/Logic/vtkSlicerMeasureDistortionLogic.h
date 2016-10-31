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

#include "qSlicerAbstractModuleWidget.h"

#include "vtkPolyData.h"
#include <vtkDenseArray.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>

// MRML includes
#include "vtkMRML.h"
#include "vtkMRMLScene.h"
#include "vtkMRMLVolumeNode.h"

// STD includes
#include <cstdlib>
#include <list>

#include "vtkSlicerMeasureDistortionModuleLogicExport.h"

class vtkSlicerVolumesLogic;

/// \ingroup Slicer_QtModules_MeasureDistortion
class VTK_SLICER_MEASUREDISTORTION_MODULE_LOGIC_EXPORT vtkSlicerMeasureDistortionLogic :
  public vtkSlicerModuleLogic
{


public:

  static vtkSlicerMeasureDistortionLogic *New();
  vtkTypeMacro(vtkSlicerMeasureDistortionLogic, vtkSlicerModuleLogic);
  void PrintSelf(ostream& os, vtkIndent indent);
  void SetVolumesLogic(vtkSlicerVolumesLogic* logic);
  vtkSlicerVolumesLogic* GetVolumesLogic();
  vtkMRMLNode* CalculateReference(vtkMRMLNode*);
 // vtkPolyData* CalculateReference(vtkMRMLNode*);
  vtkMRMLNode* CalculateDistortion(vtkMRMLNode*, vtkMRMLNode*);
  void CalculateStats(vtkMRMLNode*, vtkSmartPointer<vtkPolyData>, vtkSmartPointer<vtkDoubleArray>, double*, int dim);
  vtkSmartPointer<vtkImageData> Distortion_polyfitSVD(vtkSmartPointer<vtkPolyData>, vtkSmartPointer<vtkDoubleArray>, int*, int);
  vnl_vector<double> Fit3DPolySVD(vtkSmartPointer<vtkPolyData>, vtkSmartPointer<vtkDoubleArray>, int);
  vtkSmartPointer<vtkImageData> Eval3DPolySVD(int*, vtkSmartPointer<vtkPolyData>, vnl_vector<double>, int);
  vnl_vector<double> vnl_vectorpow(vnl_vector<double>,int);

  /*struct stats{
	  double closest1mm;
	  double closest2mm;
	  double closest3mm;
	  double closest4mm;
	  double closest5mm;
	  double closest6mm;
	  double closest7mm;
	  double largestdist;
	  double largestdistradius;
  } Distortionx, Distortionz, Distortiony;*/


  

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
	class vtkInternal;
	vtkInternal* Internal;

//public:
	//void* pntrCTSelectionChanged = &CTSelectionChanged;
	/// The currently active mrml volume node
//	void SetActiveVolumeNode(vtkMRMLVolumeNode *ActiveVolumeNode);
//	vtkMRMLVolumeNode* GetActiveVolumeNode()const;
//	Q_INVOKABLE vtkMRMLNode* currentNode()const;
	
	
};

#endif
