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

#include <QFileDialog>
#include <QMessageBox>
#include <QDebug>

// MeasureDistortion Logic includes
#include "vtkSlicerMeasureDistortionLogic.h"
#include "vtkSlicerVolumesLogic.h"

// MRML includes
#include <vtkMRMLScene.h>
#include <vtkMRMLSceneViewNode.h>
#include <vtkMRMLSliceNode.h>
#include <vtkMRMLScalarVolumeNode.h>
#include <vtkMRMLDoubleArrayNode.h>



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
//#include <vtkImageGradient.h>
//#include <vtkImageGradientMagnitude.h>
//#include <vtkImageNonMaximumSuppression.h>
#include <vtkImageCast.h>
#include <vtkImageDataGeometryFilter.h>
#include <vtkImageExport.h>
//#include <vtkImageOpenClose3D.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkKdTreePointLocator.h>
#include <vtkImageMathematics.h>
#include <vtkMath.h>
#include <vtkDenseArray.h>
#include <vtkImageMedian3D.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>


// ITK includes
#include <itkImage.h>
#include <itkVTKImageImport.h>
#include <itkLabelObject.h>
#include "itkLabelMap.h"
#include "itkLabelImageToLabelMapFilter.h"
#include <vtkITKUtility.h>
#include "itkConnectedComponentImageFilter.h"
#include "itkLabelGeometryImageFilter.h"
#include "itkPoint.h"
#include "itkPointSet.h"
#include "itkRegionOfInterestImageFilter.h"






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
//vtkPolyData* vtkSlicerMeasureDistortionLogic::CalculateReference(vtkMRMLNode* CTNode)
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
	ReferenceImage = CTImage;
	vtkMRMLScalarVolumeNode *ReferenceVolumeNode = CTVolumeNode;
	double originRAS[3];
	CTVolumeNode->GetOrigin(originRAS);
	double* PxlSpacing = CTVolumeNode->GetSpacing();
	int* Extent = CTImage->GetExtent();
	double min = CTImage->GetScalarTypeMin();
	double max = CTImage->GetScalarTypeMax();
	int cpsz[3];
	cpsz[0] = round(5 / PxlSpacing[0]);
	cpsz[1] = round(5 / PxlSpacing[1]);
	cpsz[2] = round(5 / PxlSpacing[2]);
	int V = round((3.14)*(pow(cpsz[0] / 2, 2))*cpsz[2]);

	//vtkMatrix4x4 *inputRASToIJK = vtkMatrix4x4::New();
//	ReferencCoordsNode->SetXYValue(0, 2.5e1, 3e1, 4.5e1);
//	ReferencCoordsNode->SetXYValue(1, 2.5, 3, 4.5);

	//NonMaximum Suppression
//	vtkNew<vtkImageGradient> gradientFilter;
//	gradientFilter->SetInputData(ReferenceImage);
//	gradientFilter->Update();
//	vtkNew<vtkImageCast> gradientCastFilter;
//	gradientCastFilter->SetOutputScalarTypeToUnsignedShort();
//	gradientCastFilter->SetInputConnection(gradientFilter->GetOutputPort());
//	gradientCastFilter->Update();
//	vtkNew<vtkImageGradientMagnitude> gradientMagnitudeFilter;
//	gradientMagnitudeFilter->SetInputData(ReferenceImage);
//	gradientMagnitudeFilter->Update();
//	vtkNew<vtkImageCast> gradientMagnitudeCastFilter;
//	gradientMagnitudeCastFilter->SetOutputScalarTypeToUnsignedShort();
//	gradientMagnitudeCastFilter->SetInputConnection(gradientMagnitudeFilter->GetOutputPort());
//	gradientMagnitudeCastFilter->Update();
//	vtkNew<vtkImageNonMaximumSuppression> suppressionFilter;
//	suppressionFilter->SetInputConnection(
//		0, gradientMagnitudeFilter->GetOutputPort());
//	suppressionFilter->SetInputConnection(
//		1, gradientFilter->GetOutputPort());
//	suppressionFilter->Update();
//	vtkNew<vtkImageCast> suppressionCastFilter;
//	suppressionCastFilter->SetOutputScalarTypeToUnsignedShort();
//	suppressionCastFilter->SetInputConnection(suppressionFilter->GetOutputPort());
//	suppressionFilter->SetDimensionality(3);
//	suppressionCastFilter->Update();
//	ReferenceImage = suppressionCastFilter->GetOutput();

	//Median Filter

	//vtkSmartPointer<vtkImageMedian3D> medianFilter =
	//	vtkSmartPointer<vtkImageMedian3D>::New();
	//medianFilter->SetInputData(ReferenceImage);
	//medianFilter->SetKernelSize(krnsz[0], krnsz[1], krnsz[2]);
	//medianFilter->Update();
	//ReferenceImage = medianFilter->GetOutput();


	//Thresholding
	//vtkNew<vtkImageThreshold> CTThreshold;
	vtkSmartPointer<vtkImageThreshold> CTThreshold =
		vtkSmartPointer<vtkImageThreshold>::New();
	//vtkNew<vtkThreshold> CTThreshold;
	CTThreshold->SetInputData(ReferenceImage);
	double lower = -584;
	CTThreshold->ThresholdByUpper(lower);
	CTThreshold->ReplaceInOn();
	CTThreshold->SetInValue(32767);
	CTThreshold->SetOutValue(0);
	//CTThreshold->SetInValue(ReferenceImage->GetScalarTypeMax());
	//CTThreshold->SetOutValue(ReferenceImage->GetScalarTypeMin());
	CTThreshold->Update();
	ReferenceImage = CTThreshold->GetOutput();


	//Remove background with morphological Open/Close
//	vtkSmartPointer<vtkImageOpenClose3D> openClose =
//		vtkSmartPointer<vtkImageOpenClose3D>::New();
	////vtkImageOpenClose3D* openClose;
	//openClose->SetInputData(ReferenceImage);
	//openClose->SetOpenValue(0);
	//openClose->SetCloseValue(32767);
	//int krnsz[3];
	//krnsz[0] = round(5 / PxlSpacing[0]);
	//krnsz[1] = round(5 / PxlSpacing[1]);
	//krnsz[2] = round(5 / PxlSpacing[2]);
	//qDebug() << krnsz[0] << krnsz[1] << krnsz[2];
	//openClose->SetKernelSize(krnsz[0], krnsz[1], krnsz[2]);
	////openClose->SetKernelSize(5, 5, 3);
	////openClose->ReleaseDataFlagOff();
	//openClose->Update();
	//ReferenceImage = openClose->GetOutput();
	//openClose->GetCloseValue();
	//openClose->GetOpenValue();


	//Establish pipeline connection between VTK and ITK
	vtkImageExport *exporter;
	exporter = vtkImageExport::New();
	exporter->SetInputData(ReferenceImage);
//	exporter->ImageLowerLeftOn();
//	exporter->Export(cImage);		
	typedef itk::Image<short, 3>  ImageType;
	typedef itk::Image<short, 3> OutputImageType;
	typedef itk::VTKImageImport< ImageType> ImageImportType;
	ImageImportType::Pointer importer = ImageImportType::New();
	ImageType::Pointer itkImage=ImageType::New();
	ImageType::Pointer labelImage = ImageType::New();
	ConnectPipelines(exporter, importer);
	itkImage = importer->GetOutput();
	

	//ITK Connectivity Filter
	typedef itk::ConnectedComponentImageFilter <ImageType, OutputImageType>
		ConnectedComponentImageFilterType;
	ConnectedComponentImageFilterType::Pointer connected = ConnectedComponentImageFilterType::New();
	connected->SetInput(itkImage);
	connected->Update();
	//input here for future code generalizing for phantom
	labelImage=connected->GetOutput();
	//qDebug() << connected->GetObjectCount();
	vtkIdType numPoints = connected->GetObjectCount();
//	qDebug() << numPoints;
	ImageType::DirectionType DirCosinesRAS;

	//ATTENTION: direction cosines not yet generalized
	DirCosinesRAS.SetIdentity();
	DirCosinesRAS(0, 0) = -1;
	DirCosinesRAS(1, 1) = -1;


	//LabelImage Geometry processing
	typedef itk::LabelGeometryImageFilter<ImageType> LabelGeometryImageFilterType;
//	qDebug() << "test";
	LabelGeometryImageFilterType::Pointer labelGeometryImageFilter = LabelGeometryImageFilterType::New();
//	qDebug() << "test";
	labelGeometryImageFilter->SetInput(labelImage);
	//labelGeometryImageFilter->SetIntensityInput(itkImage);
	//labelGeometryImageFilter->CalculatePixelIndicesOn();
//	qDebug() << "test";
	labelGeometryImageFilter->Update();					//odd behavior here?
//	qDebug() << "test";
	LabelGeometryImageFilterType::LabelsType allLabels =
		labelGeometryImageFilter->GetLabels();

	typedef itk::PointSet< float ,3 >   PointSetType;
	typedef PointSetType::PointType PointType;
	typedef PointSetType::PointsContainerPointer PointsContainerPointer;
	//PointSetType::Pointer  PointSet = PointSetType::New();
	//PointsContainerPointer  points = PointSet->GetPoints();
	PointType pijk;
	//vtkPointSet* ReferencePoints;
	vtkSmartPointer<vtkPoints> points =
		vtkSmartPointer<vtkPoints>::New();
	//points->SetNumberOfPoints(1);
//	typedef itk::Mesh< float, 3 >   MeshType;
	//MeshType::Pointer  mesh = MeshType::New();
//	qDebug() << numPoints;
	//points->SetNumberOfPoints(numPoints);
	//qDebug() << "test1";


	vtkIdType j = 0;
	float pRAS[3];
	
	//p0[0] = 1; p0[1] = 1; p0[2] = 1;
	//points->SetPoint(j, p0);
	
	LabelGeometryImageFilterType::LabelsType::iterator allLabelsIt;
	qDebug() << labelGeometryImageFilter->GetNumberOfLabels();
	for (allLabelsIt = allLabels.begin(); allLabelsIt != allLabels.end(); allLabelsIt++)
	{	
		LabelGeometryImageFilterType::LabelPixelType labelValue = *allLabelsIt;
		pijk=labelGeometryImageFilter->GetCentroid(labelValue);
		
		
		if (((labelGeometryImageFilter->GetVolume(labelValue)) < (V-15)) || 
			((labelGeometryImageFilter->GetVolume(labelValue)) > (V+15)) ||
			(labelGeometryImageFilter->GetEccentricity(labelValue) > 0.9))
			
		{
		//	qDebug() << "Volume: " << labelGeometryImageFilter->GetVolume(labelValue);
		//	qDebug() << "Eccentricity: " << labelGeometryImageFilter->GetEccentricity(labelValue);
		}
		else
		{
			
		//	qDebug() << "Volume: " << labelGeometryImageFilter->GetVolume(labelValue);
		//	qDebug() << "Eccentricity: " << labelGeometryImageFilter->GetEccentricity(labelValue);
		//	qDebug() << "j" << j;
		//	qDebug() << "point0:" << p[0];
		//	qDebug() << "point1:" << p[1];
		//	qDebug() << "point2:" << p[2];

			pRAS[0] = pijk[0] * (DirCosinesRAS(0, 0) * PxlSpacing[0]) + originRAS[0];
			pRAS[1] = pijk[1] * (DirCosinesRAS(1, 1) * PxlSpacing[1]) + originRAS[1];
			pRAS[2] = pijk[2] * (DirCosinesRAS(2, 2) * PxlSpacing[2]) + originRAS[2];
			
			points->InsertNextPoint(pRAS);
			j++;
		}

	}

	

	//Output
	vtkSmartPointer<vtkXMLPolyDataWriter> writer =
		vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	vtkSmartPointer<vtkPolyData> CTpolydata =
		vtkSmartPointer<vtkPolyData>::New();

	CTpolydata->SetPoints(points);
	writer->SetFileName("Reference.vtp");
	writer->SetInputData(CTpolydata);
	writer->Write();
	QString filepath = QDir::cleanPath(QDir::currentPath() + QDir::separator() + "Reference.vtp");
	QMessageBox msgBox;
	msgBox.setWindowTitle("Reference Saved");
	msgBox.setInformativeText("Reference has been saved to:");
	msgBox.setDetailedText(filepath);
	msgBox.exec();
	qDebug() << "j" << j;

	//Cleanup
	exporter->Delete();
	//polydata->Delete();
	//writer->Delete();

	ReferenceVolumeNode->SetAndObserveImageData(ReferenceImage);
	ReferenceNode = vtkMRMLNode::SafeDownCast(ReferenceVolumeNode);
	return ReferenceNode;
	//return CTpolydata;
}

//------------------------------------------------------------------------------
vtkMRMLNode* vtkSlicerMeasureDistortionLogic::CalculateDistortion(vtkMRMLNode* MR1Node, vtkMRMLNode* MR2Node){
	
	vtkSmartPointer<vtkXMLPolyDataReader> reader =
		vtkSmartPointer<vtkXMLPolyDataReader>::New();
	vtkSmartPointer<vtkPolyData> CTpolydata =
		vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkPolyData> Difference1 =
		vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkPolyData> Difference2 =
		vtkSmartPointer<vtkPolyData>::New();
	vtkMRMLNode *GNLDistortionNode;
	
	

	//Read Reference Points
	//	double *ptest;
		reader->SetFileName("Reference.vtp");
		reader->Update();
		CTpolydata = reader->GetOutput();
		vtkMRMLScalarVolumeNode *MR1VolumeNode;
		vtkImageData *MR1Image;
		MR1VolumeNode = vtkMRMLScalarVolumeNode::SafeDownCast(MR1Node);
		MR1Image = MR1VolumeNode->GetImageData();
		vtkImageData *ThresholdImage=MR1Image;


	//Find MR Center Control Points (closest CP to isocenter)
		vtkSmartPointer<vtkKdTreePointLocator> KdTree =
			vtkSmartPointer<vtkKdTreePointLocator>::New();
		vtkSmartPointer<vtkImageThreshold> MRThreshold =
			vtkSmartPointer<vtkImageThreshold>::New();

		double originRAS[3];
		MR1VolumeNode->GetOrigin(originRAS);
		double* PxlSpacingRAS = MR1VolumeNode->GetSpacing();
		int* Extent = MR1Image->GetExtent();
		int cpsz[3];
		cpsz[0] = round(5 / PxlSpacingRAS[0]);
		cpsz[1] = round(5 / PxlSpacingRAS[1]);
		cpsz[2] = round(5 / PxlSpacingRAS[2]);
		int V = round((3.14)*(pow(cpsz[0] / 2, 2))*cpsz[2]);
		
		double CTbounds[6];
		CTpolydata->GetBounds(CTbounds);

		MRThreshold->SetInputData(MR1Image);
		double lower = 100;
		MRThreshold->ThresholdByUpper(lower);
		MRThreshold->ReplaceInOn();
		MRThreshold->SetInValue(MR1Image->GetScalarTypeMax());
		MRThreshold->SetOutValue(0);
		MRThreshold->Update();
		ThresholdImage = MRThreshold->GetOutput();


		//Establish pipeline connection between VTK and ITK
		vtkImageExport *exporter;
		exporter = vtkImageExport::New();
		exporter->SetInputData(ThresholdImage);
		//	exporter->ImageLowerLeftOn();
		//	exporter->Export(cImage);		
		typedef itk::Image< unsigned short, 3 >  ImageType;
		typedef itk::Image< unsigned short, 3 > OutputImageType;
		typedef itk::VTKImageImport< ImageType> ImageImportType;
		ImageImportType::Pointer importer = ImageImportType::New();
		ImageType::Pointer itkImage = ImageType::New();
		ImageType::Pointer labelImage = ImageType::New();
		ConnectPipelines(exporter, importer);
		itkImage = importer->GetOutput();
		//const ImageType::PointType & originLPS =originRAS;
		//const ImageType::SpacingType& PxlSpacingLPS = PxlSpacingRAS;
		ImageType::DirectionType DirCosinesRAS;

		//ATTENTION: direction cosines not yet generalized
		DirCosinesRAS.SetIdentity();
		DirCosinesRAS(0, 0) = -1;
		DirCosinesRAS(1, 1) = -1;
		//itkImage->SetOrigin(originLPS);
		//itkImage->SetSpacing(PxlSpacingLPS);
		//itkImage->SetDirection(DirCosinesRAS);

		
		//qDebug() << itkImage->GetOrigin()[0] << itkImage->GetOrigin()[1] << itkImage->GetOrigin()[2];
		//ITK Connectivity Filter
		typedef itk::ConnectedComponentImageFilter <ImageType, OutputImageType >
			ConnectedComponentImageFilterType;
		ConnectedComponentImageFilterType::Pointer connected = ConnectedComponentImageFilterType::New();
		

		double center[3];           //in ijk image coordinates
		center[0] = (Extent[1] + Extent[0]) / 2;
		center[1] = (Extent[3] + Extent[2]) / 2;
		center[2] = (Extent[5] + Extent[4]) / 2;

		double centerCT0[3];
		centerCT0[0] = (CTbounds[1] - (10 * 25));
		centerCT0[1] = (CTbounds[3] - (7 * 25));
		centerCT0[2] = (CTbounds[5] - (10 * 25));

		typedef itk::RegionOfInterestImageFilter< ImageType, ImageType > FilterType;
		FilterType::Pointer filter = FilterType::New();
		ImageType::IndexType start;
		start[0] = center[0] - 10;
		start[1] = center[1] - 10;
		start[2] = center[2] - 5;

		ImageType::IndexType end;
		end[0] = center[0] + 10;
		end[1] = center[1] + 10;
		end[2] = center[2] + 5;


		ImageType::RegionType region;
		region.SetIndex(start);
		region.SetUpperIndex(end);

		filter->SetInput(itkImage);
		filter->SetRegionOfInterest(region);

		//connected->SetInput(itkImage);
		connected->SetInput(filter->GetOutput());

		connected->Update();
		labelImage = connected->GetOutput();
		//qDebug() << connected->GetObjectCount();
		//vtkIdType numPoints = connected->GetObjectCount();


		//LabelImage Geometry processing
		typedef itk::LabelGeometryImageFilter< ImageType > LabelGeometryImageFilterType;
		LabelGeometryImageFilterType::Pointer labelGeometryImageFilter = LabelGeometryImageFilterType::New();
		labelGeometryImageFilter->SetInput(labelImage);
		//labelGeometryImageFilter->SetIntensityInput(itkImage);
		labelGeometryImageFilter->Update();
		LabelGeometryImageFilterType::LabelsType allLabels =
			labelGeometryImageFilter->GetLabels();
		typedef itk::PointSet< float, 3 >   PointSetType;
		typedef PointSetType::PointType PointType;
		typedef PointSetType::PointsContainerPointer PointsContainerPointer;
		//PointSetType::Pointer  PointSet = PointSetType::New();
		//PointsContainerPointer  points = PointSet->GetPoints();
		PointType p;

		//vtkPointSet* ReferencePoints;
		vtkSmartPointer<vtkPoints> points =
			vtkSmartPointer<vtkPoints>::New();
		//points->SetNumberOfPoints(1);
		//	typedef itk::Mesh< float, 3 >   MeshType;
		//MeshType::Pointer  mesh = MeshType::New();
		//qDebug() << numPoints;
		//points->SetNumberOfPoints(numPoints);
		//qDebug() << "test1";


		vtkIdType j = 0;
		float pijk[3];

		//p0[0] = 1; p0[1] = 1; p0[2] = 1;
		//points->SetPoint(j, p0);
		vtkImageData *TestImage = ThresholdImage;
		vtkSmartPointer<vtkImageMathematics> imageMath =
			vtkSmartPointer<vtkImageMathematics>::New();
		imageMath->SetOperationToMultiplyByK();
		imageMath->SetConstantK(0.0);
		imageMath->SetInputData(ThresholdImage);
		imageMath->Update();
		TestImage = imageMath->GetOutput();
	
		//find center control point position
		LabelGeometryImageFilterType::LabelPixelType labelValue = *allLabels.begin();
		p = labelGeometryImageFilter->GetCentroid(labelValue);
		pijk[0] = p[0] + start[0];
		pijk[1] = p[1] + start[1];
		pijk[2] = p[2] + start[2];
		double centerMR[3];
		centerMR[0] = pijk[0] * (DirCosinesRAS(0, 0) * PxlSpacingRAS[0]) + originRAS[0];
		centerMR[1] = pijk[1] * (DirCosinesRAS(1, 1) * PxlSpacingRAS[1]) + originRAS[1];
		centerMR[2] = pijk[2] * (DirCosinesRAS(2, 2) * PxlSpacingRAS[2]) + originRAS[2];

		KdTree->SetDataSet(CTpolydata);
		KdTree->BuildLocator();
		vtkIdType iD = KdTree->FindClosestPoint(centerCT0);
		double centerCT[3];
		KdTree->GetDataSet()->GetPoint(iD, centerCT);

		qDebug() << "centerCT0" << centerCT0[0] << centerCT0[1] << centerCT0[2];
		qDebug() << "CTbounds" << CTbounds[0] << CTbounds[1] << CTbounds[2]
			<< CTbounds[3] << CTbounds[4] << CTbounds[5];
		vtkIdType numPoints; 
		double CTp[3];
		double MRp[3];


	// iterate through Reference points and find corresponding MR CP positions
		j = 0;
		for (vtkIdType i = 0; i < points->GetNumberOfPoints(); i++)
		{
			CTpolydata->GetPoint(i, CTp);
			CTp[0] = CTp[0] - centerCT[0] + centerMR[0];
			start[0] = CTp[0] - 10;
			start[1] = CTp[1] - 10;
			start[2] = CTp[2] - 5;

			end[0] = CTp[0] + 10;
			end[1] = CTp[1] + 10;
			end[2] = CTp[2] + 5;

			region.SetIndex(start);
			region.SetUpperIndex(end);
			filter->SetInput(itkImage);
			filter->SetRegionOfInterest(region);
			connected->SetInput(filter->GetOutput());
			connected->Update();
			labelImage = connected->GetOutput();
			numPoints = connected->GetObjectCount();

			//LabelImage Geometry processing
			labelGeometryImageFilter->SetInput(labelImage);
			labelGeometryImageFilter->Update();
			allLabels = labelGeometryImageFilter->GetLabels();

			if (numPoints == 1){
				labelValue = *allLabels.begin();
				if (((labelGeometryImageFilter->GetVolume(labelValue)) < (V - 35)) ||
					((labelGeometryImageFilter->GetVolume(labelValue)) > (V + 35)) /*||   //could use more fine tuning
					(labelGeometryImageFilter->GetEccentricity(labelValue) > 0.9)*/)
				{
				}
				else
				{
					p = labelGeometryImageFilter->GetCentroid(labelValue);
					pijk[0] = p[0] + start[0];
					pijk[1] = p[1] + start[1];
					pijk[2] = p[2] + start[2];

					MRp[0] = pijk[0] * (DirCosinesRAS(0, 0) * PxlSpacingRAS[0]) + originRAS[0];
					MRp[1] = pijk[1] * (DirCosinesRAS(1, 1) * PxlSpacingRAS[1]) + originRAS[1];
					MRp[2] = pijk[2] * (DirCosinesRAS(2, 2) * PxlSpacingRAS[2]) + originRAS[2];

					points->InsertNextPoint(pijk);
					j++;
				}
			}
		}

		//----


		LabelGeometryImageFilterType::LabelsType::iterator allLabelsIt;
		for (allLabelsIt = allLabels.begin(); allLabelsIt != allLabels.end(); allLabelsIt++)
		{
			labelValue = *allLabelsIt;
			p = labelGeometryImageFilter->GetCentroid(labelValue); //in

			
			
			//ImageType::DirectionType DirCosinesRAS;
			if (((labelGeometryImageFilter->GetVolume(labelValue)) < (V - 35)) ||
				((labelGeometryImageFilter->GetVolume(labelValue)) > (V + 35)) /*||   //could use more fine tuning
				(labelGeometryImageFilter->GetEccentricity(labelValue) > 0.9)*/)
			{
			}
			else
			{

	//			qDebug() << "Volume: " << labelGeometryImageFilter->GetVolume(labelValue);
	//			qDebug() << "Eccentricity: " << labelGeometryImageFilter->GetEccentricity(labelValue);
	//			qDebug() << "j" << j;
			//	qDebug() << "point0:" << pitk[0];
			//	qDebug() << "point1:" << pitk[1];
			//	qDebug() << "point2:" << pitk[2];
				/*pijk[0] = (originRAS[0] + (DirCosinesRAS(0, 0)*pitk[0])) / PxlSpacingRAS[0];
				pijk[1] = (originRAS[1] + (DirCosinesRAS(1, 1)*pitk[1])) / PxlSpacingRAS[1];
				pijk[2] = (originRAS[2] + (DirCosinesRAS(2, 2)*pitk[2])) / PxlSpacingRAS[2];*/

				/*pijk[0] = (pitk[0] - originRAS[0]) / (DirCosinesRAS(0, 0)*PxlSpacingRAS[0]);
				pijk[1] = (pitk[1] - originRAS[1]) / (DirCosinesRAS(1, 1)*PxlSpacingRAS[1]);
				pijk[2] = (pitk[2] - originRAS[2]) / (DirCosinesRAS(2, 2)*PxlSpacingRAS[2]);*/

				pijk[0] = p[0] + start[0];
				pijk[1] = p[1] + start[1];
				pijk[2] = p[2] + start[2];

				qDebug() << "pointijk0:" << pijk[0];
				qDebug() << "pointijk1:" << pijk[1];
				qDebug() << "pointijk2:" << pijk[2];
		

			/*p0[0] = 10;
			p0[1] = 10;
			p0[2] = 10;*/
				points->InsertNextPoint(pijk);
				j++;
				TestImage->SetScalarComponentFromDouble(round(pijk[0]), round(pijk[1]), round(pijk[2]), 0, ThresholdImage->GetScalarTypeMax());
			}

		}
		//Output
		vtkSmartPointer<vtkXMLPolyDataWriter> writer =
			vtkSmartPointer<vtkXMLPolyDataWriter>::New();
		vtkSmartPointer<vtkPolyData> MR1polydata =
			vtkSmartPointer<vtkPolyData>::New();
		MR1polydata->SetPoints(points);
		// Find center of image
		//double center[3];
	//	center[0] = (Extent[1] + Extent[0]) / 2;
	//	center[1] = (Extent[3] + Extent[2]) / 2;
	//	center[2] = (Extent[5] + Extent[4]) / 2;
		// Pick a radius for the circle
	//	int radius = 5;
		

	//	double origintest[3];
	//	MR1VolumeNode->GetOrigin(origintest);
	//	double* PxlSpacingRAS = MR1VolumeNode->GetSpacing();

		//KdTree->SetDataSet(MR1polydata);
		//KdTree->BuildLocator();
		//vtkIdType iD = KdTree->FindClosestPoint(center);
		//double centerPoint[3];
		//KdTree->GetDataSet()->GetPoint(iD, centerPoint);

		

	//	for (vtkIdType i = 0; i < points->GetNumberOfPoints(); i++)
	//	{
		//	double p[3];
		//	MR1polydata->GetPoint(i, p);
		//	TestImage->SetScalarComponentFromDouble(round(p[0]), round(p[1]), round(p[2]),0, ThresholdImage->GetScalarTypeMax());


	//	}

		//ThresholdImage = TestImage;
		double bounds[6];
		MR1polydata->GetBounds(bounds);
		//KdTree->SetDataSet();



		qDebug() << "originRAS0" << originRAS[0];
		qDebug() << "originRAS1" << originRAS[1];
		qDebug() << "originRAS2" << originRAS[2];
		qDebug() << "c0" << center[0];
		qDebug() << "c1" << center[1];
		qDebug() << "c2" << center[2];
		qDebug() << "bounds0" << bounds[0];
		qDebug() << "bounds1" << bounds[1];
		qDebug() << "bounds2" << bounds[2];
		qDebug() << "bounds3" << bounds[3];
		qDebug() << "bounds4" << bounds[4];
		qDebug() << "bounds5" << bounds[5];
		qDebug() << "j" << j;
		qDebug() << "numpoints" << numPoints;

	//	ptest = polydata->GetPoint(0);
	//	qDebug() << "point0:" << ptest[0];
	//	qDebug() << "point1:" << ptest[1];
	//	qDebug() << "point2:" << ptest[2];

	//Calculate  which control point this is on Reference (Currently phantom dependent)

	//Transform Reference Coordinates to MR Coordinates

	
	//Calculate MR1 Centroid positions and Difference from Reference
		//Difference1 = CalculateMRCentroids(MR1Node, CTpolydata);

	//Calculate MR2 Centroid positions and Difference from Reference
		//Difference2 = CalculateMRCentroids(MR2Node, CTpolydata);

	//Apply Mask


	//Calculate Sequence Dependent Distortion
	//vtkPolyData *SequenceDependent = ((Difference1y - Difference2y). / 2);

	//Remove Sequence Dependent Distortion from Total Distortion
		vtkSmartPointer<vtkDenseArray<double>>  GNLDist =
			vtkSmartPointer<vtkDenseArray<double>>::New();
		GNLDist->Resize(MR1polydata->GetNumberOfPoints());
		GNLDist->Fill(1);
		//GNLDisty = Differences1y - SequenceDependent;
		//GNLDistx = Differences1x;
		//GNLDistz = Differences1z;
	//Interpolate Distortion Map
		int order = 6;
		//GNLDistortionNode = Distortion_polyfitSVD(MR1polydata, GNLDist, Extent, order);
		//Distortion_polyfitSVD(MR1polydata, GNLDist, Extent, order);
		//GNLDistortionNodex = SVD_polyfit();
		//GNLDistortionNodey = SVD_polyfit();
		//GNLDistortionNodez = SVD_polyfit();


		//Cleanup
		exporter->Delete();
		MR1VolumeNode->SetAndObserveImageData(TestImage);
		GNLDistortionNode = vtkMRMLNode::SafeDownCast(MR1VolumeNode);
		return GNLDistortionNode;
}
//-----------------------------------------------------------------------------
vtkPolyData* vtkSlicerMeasureDistortionLogic::CalculateMRCentroids(vtkMRMLNode* MRNode, vtkPolyData*  CTpolydata){
	vtkSmartPointer<vtkPolyData> MRpolydata =
		vtkSmartPointer<vtkPolyData>::New();


	return MRpolydata;
}
//------------------------------------------------------------------------------
//vtkMRMLNode* vtkSlicerMeasureDistortionLogic::Distortion_polyfitSVD(vtkPolyData* MRpolydata,
void vtkSlicerMeasureDistortionLogic::Distortion_polyfitSVD(vtkPolyData* MRpolydata,
	vtkDenseArray<double>* GNLDist, int* Extent, int order){
	vtkMRMLNode *GNLDistortionNode;

	vnl_matrix<double> coeffs = Fit3DPolySVD(MRpolydata, GNLDist, order);
	double* zbar = Eval3DPolySVD(Extent, coeffs, order);

	//return GNLDistortionNode;
}
//-----------------------------------------------------------------------------
vnl_matrix<double> vtkSlicerMeasureDistortionLogic::Fit3DPolySVD(vtkPolyData* MRpolydata,
	vtkDenseArray<double>* GNLDist, int order){
	
	//if (MRpolydata->GetNumberOfPoints() != MRpolydata->GetNumberOfPoints());

	//Scale
//	vtkSmartPointer<vtkDenseArray<double>>  data =
//		vtkSmartPointer<vtkDenseArray<double>>::New();
//	data->Resize(MRpolydata->GetNumberOfPoints(),4);
	vnl_matrix<double> data(MRpolydata->GetNumberOfPoints(), 3);
	double bounds[6];
	double maxabs[4] = { 0, 0, 0, 0 };
//	MRpolydata->GetBounds(bounds);
//	qDebug() << dim1start;
//	qDebug() << dim1end;
//	return bounds;
	for (vtkIdType i = 0; i < MRpolydata->GetNumberOfPoints(); i++){
		double p[3];	
		MRpolydata->GetPoint(i, p);
		data(i, 0) = p[0];
		data(i, 1) = p[1];
		data(i, 2) = p[2];
		data(i, 3) = GNLDist->GetValue(i);
		if (abs(data(i, 0))>maxabs[0]){ maxabs[0] = data(i, 0); }
		if (abs(data(i, 1))>maxabs[1]){ maxabs[1] = data(i, 1); }
		if (abs(data(i, 2))>maxabs[2]){ maxabs[2] = data(i, 2); }
		if (abs(data(i, 3))>maxabs[3]){ maxabs[3] = data(i, 3); }
	}
	int numCoeffs = (order + 3)*(order + 2)*(order + 1) / 6;
	vnl_matrix<double> A(MRpolydata->GetNumberOfPoints(), numCoeffs, 0);
	vnl_vector<double> B(MRpolydata->GetNumberOfPoints());
	vnl_vector<double> C(MRpolydata->GetNumberOfPoints());
	vnl_vector<double> D(MRpolydata->GetNumberOfPoints());
	int column = 0;
	for (int xpower = 0; xpower < order; xpower++){
		for (int ypower = 0; ypower < order - xpower; ypower++){
			for (int zpower = 0; zpower < order - xpower - ypower; zpower++){
				B = vnl_vectorpow(data.get_column(0).operator/= (maxabs[0]), xpower);
				C = vnl_vectorpow(data.get_column(1).operator/= (maxabs[1]), ypower);
				D = vnl_vectorpow(data.get_column(2).operator/= (maxabs[2]), zpower);
				for (int k = 0; k < MRpolydata->GetNumberOfPoints(); k++){
					A(k, column) = B(k)*C(k)*D(k);
				}

				column = column + 1;
			}
		}
	}
		
	vnl_svd<double> svd(A);
		double sigma = pow(std::numeric_limits<double>::epsilon(),1/order);
		
		vnl_matrix<double> W = svd.W();
		vnl_matrix<double> q = W;
		vnl_matrix<double> result;
		vnl_matrix<double> coeffs;
		vnl_matrix<double> comp(data.rows(),1);
		unsigned int size; 
		int sizec = W.columns();
		int sizer = W.rows();
		if (sizer >= sizec){
			size = W.columns();
		}
		else{
			size = W.rows();
		}
		for (int i = 0; i < size; i++){
			if (abs(W(i, i)) >= sigma){
				q(i, i) = 1 / W(i, i);
			}
			else{
				q(i, i) = 0;
			}
		}		
		comp.set_column(0, data.get_column(3).operator/= (maxabs[3]));
		result = svd.V().operator*(q.transpose());
		result = result.operator*(svd.U().transpose());
		result = result.operator*(comp);

		//rescale results
		int i = 0;
		for (int xpower = 0; xpower < order; xpower++){
			for (int ypower = 0; ypower < order - xpower; ypower++){
				for (int zpower = 0; zpower < order - xpower - ypower; zpower++){
					result(i,0) = result(i,0)*(pow(1 / maxabs[0], xpower))*(pow(1 / maxabs[1], ypower))*(pow(1 / maxabs[2], zpower)) / (1 / maxabs[3]);


					i = i + 1;
				}
			}
		}

		coeffs = result;
		return coeffs;
}
//-----------------------------------------------------------------------------
double* vtkSlicerMeasureDistortionLogic::Eval3DPolySVD(int* Extent, vnl_matrix<double> coeffs, int order){
	double *zbar;


	return zbar;
}
//-----------------------------------------------------------------------------
vnl_vector<double> vtkSlicerMeasureDistortionLogic::vnl_vectorpow(vnl_vector<double> v,int p){
	for (int i = 0; i < p-1; i++){
		for (int j = 0; j < v.size(); j++){
			v(j) = v(j)*v(j);
		}
	}
	return v;
}