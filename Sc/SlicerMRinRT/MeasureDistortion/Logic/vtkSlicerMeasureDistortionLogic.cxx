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

// MeasureDistortion Logic includes
#include "vtkSlicerMeasureDistortionLogic.h"
#include "vtkSlicerVolumesLogic.h"

// MRML includes
#include <vtkMRMLScene.h>
#include <vtkMRMLSceneViewNode.h>
#include <vtkMRMLSliceNode.h>
#include <vtkMRMLScalarVolumeNode.h>
#include <vtkMRMLScene.h>
#include <QDebug>


// VTK includes
#include <vtkIntArray.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkimagedata.h>


// STD includes
#include <cassert>

// ----------------------------------------------------------------------------
class vtkSlicerMeasureDistortionLogic::vtkInternal
{
public:
	vtkInternal();

	vtkSlicerVolumesLogic* VolumesLogic;
};


//----------------------------------------------------------------------------
vtkSlicerMeasureDistortionLogic::vtkInternal::vtkInternal()
{
	this->VolumesLogic = 0;
}
//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkSlicerMeasureDistortionLogic);

//----------------------------------------------------------------------------
vtkSlicerMeasureDistortionLogic::vtkSlicerMeasureDistortionLogic()
{
	this->Internal = new vtkInternal;
}

//----------------------------------------------------------------------------
vtkSlicerMeasureDistortionLogic::~vtkSlicerMeasureDistortionLogic()
{
}


//----------------------------------------------------------------------------
void vtkSlicerMeasureDistortionLogic::SetVolumesLogic(vtkSlicerVolumesLogic* logic)
{
	this->Internal->VolumesLogic = logic;
}

//----------------------------------------------------------------------------
vtkSlicerVolumesLogic* vtkSlicerMeasureDistortionLogic::GetVolumesLogic()
{
	return this->Internal->VolumesLogic;
}

//----------------------------------------------------------------------------
void vtkSlicerMeasureDistortionLogic::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

//---------------------------------------------------------------------------
void vtkSlicerMeasureDistortionLogic::SetMRMLSceneInternal(vtkMRMLScene * newScene)
{
  vtkNew<vtkIntArray> events;
  events->InsertNextValue(vtkMRMLScene::NodeAddedEvent);
  events->InsertNextValue(vtkMRMLScene::NodeRemovedEvent);
  events->InsertNextValue(vtkMRMLScene::EndBatchProcessEvent);
  this->SetAndObserveMRMLSceneEventsInternal(newScene, events.GetPointer());
}

//-----------------------------------------------------------------------------
void vtkSlicerMeasureDistortionLogic::RegisterNodes()
{
  assert(this->GetMRMLScene() != 0);

//  vtkMRMLSceneViewNode* viewNode = vtkMRMLSceneViewNode::New();
//  this->GetMRMLScene()->RegisterNodeClass(viewNode);

  //vtkMRMLSceneViewNode* viewNode = vtkMRMLSceneViewNode::SafeDownCast(this->GetMRMLScene()->GetNodeByID(id));
}

//---------------------------------------------------------------------------
void vtkSlicerMeasureDistortionLogic::UpdateFromMRMLScene()
{
  assert(this->GetMRMLScene() != 0);
}

//---------------------------------------------------------------------------
void vtkSlicerMeasureDistortionLogic
::OnMRMLSceneNodeAdded(vtkMRMLNode* vtkNotUsed(node))
{
}

//---------------------------------------------------------------------------
void vtkSlicerMeasureDistortionLogic
::OnMRMLSceneNodeRemoved(vtkMRMLNode* vtkNotUsed(node))
{
}
//--------------------------------------------------------------------------
//void vtkSlicerMeasureDistortionLogic::currentNode();
//{

//}
//------------------------------------------------------
void vtkSlicerMeasureDistortionLogic::CreateReference(vtkMRMLNode* CTnode)
{
	//Q_D(qSlicerMeasureDistortionModuleWidget);
	


}

//------------------------------------------------------------------------------