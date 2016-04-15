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
#include <vtkMRMLDoubleArrayNode.h>
#include <QDebug>


// VTK includes
#include <vtkSmartPointer.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkimagedata.h>
#include <vtkImageThreshold.h>
#include <vtkThreshold.h>
#include "vtkPolyDataConnectivityFilter.h"
#include "vtkConnectivityFilter.h"
#include <vtkImageGradient.h>
#include <vtkImageGradientMagnitude.h>
#include <vtkImageNonMaximumSuppression.h>
#include <vtkImageCast.h>
#include <vtkImageDataGeometryFilter.h>
#include <vtkImageExport.h>

#include <itkImage.h>
#include <itkVTKImageImport.h>
#include <vtkITKUtility.h>
#include "itkConnectedComponentImageFilter.h"





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
vtkMRMLNode* vtkSlicerMeasureDistortionLogic::CalculateReference(vtkMRMLNode* CTNode)
{
	vtkNew<vtkMRMLDoubleArrayNode> ReferencCoordsNode;
	vtkMRMLScalarVolumeNode *CTVolumeNode;
	vtkMRMLNode *ReferenceNode;
	vtkImageData* ReferenceImage;
	vtkImageData* ReferenceImage1;
	vtkImageData* CTImage;
//	double x[2]; double y[2]; double z[2];

	CTVolumeNode = vtkMRMLScalarVolumeNode::SafeDownCast(CTNode);
	CTImage = CTVolumeNode->GetImageData();
	double* PxlSpacing = CTVolumeNode->GetSpacing();
	int* Extent = CTImage->GetExtent();
	double min = CTImage->GetScalarTypeMin();
	double max = CTImage->GetScalarTypeMax();
	//vtkMatrix4x4 *inputRASToIJK = vtkMatrix4x4::New();
//	ReferencCoordsNode->SetXYValue(0, 2.5e1, 3e1, 4.5e1);
//	ReferencCoordsNode->SetXYValue(1, 2.5, 3, 4.5);

	//Thresholding
	vtkMRMLScalarVolumeNode *ReferenceVolumeNode = CTVolumeNode;
	vtkNew<vtkImageThreshold> CTThreshold;
	//vtkNew<vtkThreshold> CTThreshold;
	CTThreshold->SetInputData(CTImage);
	double lower = 56;
	CTThreshold->ThresholdByUpper(lower);
	CTThreshold->ReplaceInOn();
	CTThreshold->SetInValue(CTImage->GetScalarTypeMax());
	CTThreshold->SetOutValue(CTImage->GetScalarTypeMin());
	CTThreshold->Update();
	ReferenceImage = CTImage;
	ReferenceImage = CTThreshold->GetOutput();


	//NonMaximum Suppression
	vtkNew<vtkImageGradient> gradientFilter;
	gradientFilter->SetInputData(ReferenceImage);
	vtkNew<vtkImageCast> gradientCastFilter;
	gradientCastFilter->SetOutputScalarTypeToUnsignedShort();
	gradientCastFilter->SetInputConnection(gradientFilter->GetOutputPort());
	gradientCastFilter->Update();
	vtkNew<vtkImageGradientMagnitude> gradientMagnitudeFilter;
	gradientMagnitudeFilter->SetInputData(ReferenceImage);
	vtkNew<vtkImageCast> gradientMagnitudeCastFilter;
	gradientMagnitudeCastFilter->SetOutputScalarTypeToUnsignedShort();
	gradientMagnitudeCastFilter->SetInputConnection(gradientMagnitudeFilter->GetOutputPort());
	gradientMagnitudeCastFilter->Update();
	vtkNew<vtkImageNonMaximumSuppression> suppressionFilter;
	suppressionFilter->SetInputConnection(
		0, gradientMagnitudeFilter->GetOutputPort());
	suppressionFilter->SetInputConnection(
		1, gradientFilter->GetOutputPort());
	suppressionFilter->Update();
	vtkNew<vtkImageCast> suppressionCastFilter;
	suppressionCastFilter->SetOutputScalarTypeToUnsignedShort();
	suppressionCastFilter->SetInputConnection(suppressionFilter->GetOutputPort());
	suppressionFilter->SetDimensionality(2);
	suppressionCastFilter->Update();
	ReferenceImage1 = suppressionCastFilter->GetOutput();


	//Connectivity
	//vtkConnectivityFilter* connectivityFilter;
	//vtkImageConnectivity* connectivityFilter;
	//vtkNew<vtkImageDataGeometryFilter> imageDataGeometryFilter;
	//imageDataGeometryFilter->SetInputData(ReferenceImage);
	//imageDataGeometryFilter->Update();

//	int dims[3];
//	dims[0] = Extent[1] - Extent[0];
//	dims[1] = Extent[3] - Extent[1];
//	dims[2] = Extent[5] - Extent[4];
//	unsigned char *cImage = new unsigned char[dims[0] * dims[1] * dims[2]];
	vtkImageExport *exporter;
	exporter = vtkImageExport::New();
	exporter->SetInputData(ReferenceImage);
//	exporter->ImageLowerLeftOn();
//	exporter->Export(cImage);		
	typedef itk::Image< unsigned short, 3 >  ImageType;
	typedef itk::Image< unsigned short, 3 > OutputImageType;
	typedef itk::VTKImageImport< ImageType> ImageImportType;
	ImageImportType::Pointer importer = ImageImportType::New();
	typedef itk::ConnectedComponentImageFilter <ImageType, OutputImageType >
		ConnectedComponentImageFilterType;
//	vtkNew<vtkITKImageToImageFilter> vtkITKFilter; 
//	vtkITKFilter->SetInput(ReferenceImage);
//	ReferenceImage=vtkITKFilter->GetOutput();
	//vtkITKFilter->SetInput(ReferenceImage);
	ImageType::Pointer itkImage;
	qDebug() << "test";
	//image = vtkITKFilter->GetOutput();
	
	qDebug() << "test";
	ConnectPipelines(exporter, importer);
	qDebug() << "test";
	itkImage = importer->GetOutput();
	ConnectedComponentImageFilterType::Pointer connected = ConnectedComponentImageFilterType::New();;
	connected->SetInput(itkImage);
	connected->Update();

	qDebug() << "test";
	qDebug() << connected->GetObjectCount();


	//exporter->Delete();
	//importer->Delete();
	//connectivityFilter->SetInputConnection(CTThreshold->GetOutputPort());
	//qDebug() << "test";
	//connectivityFilter->SetExtractionModeToAllRegions();
	//connectivityFilter->SetFunctionToMeasureIsland();

	//qDebug() << "test";
	//connectivityFilter->ColorRegionsOn();
	//qDebug() << "test";
	//connectivityFilter->Update();
	//qDebug() << "test";
	//int numCtrlPnts = connectivityFilter->GetNumberOfExtractedRegions();
	//ReferenceImage = connectivityFilter->GetOutput();
	//vtkUnstructuredGrid *ControlPoints = connectivityFilter->GetOutput();
//	for (){
//		connectivityFilter->InitializeSpecifiedRegionList();
//		connectivityFilter->AddSpecifiedRegion(i);
//		connectivityFilter->Modified();
//		connectivityFilter->Update();
//	}


//	qDebug() << numCtrlPnts<<"\n";
//	qDebug() << "test";
	ReferenceVolumeNode->SetAndObserveImageData(ReferenceImage);
	ReferenceNode = vtkMRMLNode::SafeDownCast(ReferenceVolumeNode);

	exporter->Delete();
	return ReferenceNode;
}

//------------------------------------------------------------------------------
