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

// MRML includes
#include <vtkMRMLScene.h>
#include <vtkMRMLSliceNode.h>
#include <vtkMRMLViewNode.h>

// VTK includes
#include <vtkIntArray.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>

// STD includes
#include <cassert>

//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkSlicerMeasureDistortionLogic);

//----------------------------------------------------------------------------
vtkSlicerMeasureDistortionLogic::vtkSlicerMeasureDistortionLogic()
{
}

//----------------------------------------------------------------------------
vtkSlicerMeasureDistortionLogic::~vtkSlicerMeasureDistortionLogic()
{
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

//-----------------------------------------------------------------------------
vtkMRMLSliceNode* vtkSlicerMeasureDistortionLogic::GetDefaultSliceViewNode()
{
	vtkMRMLScene *scene = this->GetMRMLScene();
	if (!scene)
	{
		vtkErrorMacro("vtkSlicerMeasureDistortionLogic::GetDefaultSliceViewNode failed: invalid scene");
		return NULL;
	}
	vtkMRMLNode* defaultNode = scene->GetDefaultNodeByClass("vtkMRMLSliceNode");
	if (!defaultNode)
	{
		defaultNode = scene->CreateNodeByClass("vtkMRMLSliceNode");
		scene->AddDefaultNode(defaultNode);
		defaultNode->Delete(); // scene owns it now
	}
	return vtkMRMLSliceNode::SafeDownCast(defaultNode);
}

//-----------------------------------------------------------------------------

vtkMRMLViewNode* vtkSlicerMeasureDistortionLogic::GetDefaultThreeDViewNode()
{
	vtkMRMLScene *scene = this->GetMRMLScene();
	if (!scene)
	{
		vtkErrorMacro("vtkSlicerViewControllersLogic::GetDefaultThreeDViewNode failed: invalid scene");
		return NULL;
	}
	vtkMRMLNode* defaultNode = scene->GetDefaultNodeByClass("vtkMRMLViewNode");
	if (!defaultNode)
	{
		defaultNode = scene->CreateNodeByClass("vtkMRMLViewNode");
		scene->AddDefaultNode(defaultNode);
		defaultNode->Delete(); // scene owns it now
	}
	return vtkMRMLViewNode::SafeDownCast(defaultNode);
}



//-----------------------------------------------------------------------------
void vtkSlicerMeasureDistortionLogic::ResetAllViewNodesToDefault()
{
	vtkMRMLScene *scene = this->GetMRMLScene();
	if (!scene)
	{
		vtkErrorMacro("vtkSlicerViewControllersLogic::ResetAllViewNodesToDefault failed: invalid scene");
		return;
	}
	scene->StartState(vtkMRMLScene::BatchProcessState);
	vtkMRMLSliceNode* defaultSliceViewNode = GetDefaultSliceViewNode();
	std::vector< vtkMRMLNode* > viewNodes;
	scene->GetNodesByClass("vtkMRMLSliceNode", viewNodes);
	for (std::vector< vtkMRMLNode* >::iterator it = viewNodes.begin(); it != viewNodes.end(); ++it)
	{
		(*it)->Reset(defaultSliceViewNode);
	}
	viewNodes.clear();
	vtkMRMLViewNode* defaultThreeDViewNode = GetDefaultThreeDViewNode();
	scene->GetNodesByClass("vtkMRMLViewNode", viewNodes);
	for (std::vector< vtkMRMLNode* >::iterator it = viewNodes.begin(); it != viewNodes.end(); ++it)
	{
		(*it)->Reset(defaultThreeDViewNode);
	}
	scene->EndState(vtkMRMLScene::BatchProcessState);
}