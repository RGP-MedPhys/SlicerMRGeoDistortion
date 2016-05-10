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
#include <vtkDoubleArray.h>
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
#include <vnl/algo/vnl_svd.h>


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
	vtkMRMLNode *GNLDistortionNode =MR1Node;
	
	

	//Read Reference Points
	//	double *ptest;
		reader->SetFileName("Reference.vtp");
		reader->Update();
		CTpolydata->DeepCopy(reader->GetOutput());
		vtkMRMLScalarVolumeNode *MR1VolumeNode;
		vtkMRMLScalarVolumeNode *MR2VolumeNode;
		vtkImageData *MR1Image;
		vtkImageData *MR2Image;
		MR1VolumeNode = vtkMRMLScalarVolumeNode::SafeDownCast(MR1Node);
		MR2VolumeNode = vtkMRMLScalarVolumeNode::SafeDownCast(MR2Node);
		MR1Image = MR1VolumeNode->GetImageData();
		MR2Image = MR2VolumeNode->GetImageData();
		//vtkSmartPointer<vtkImageData> ThresholdImage1 =
		//	vtkSmartPointer<vtkImageData>::New();
		//ThresholdImage1 = MR1Image;
		//vtkSmartPointer<vtkImageData> ThresholdImage2 =
		//	vtkSmartPointer<vtkImageData>::New();
		//ThresholdImage2 = MR2Image;


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
		vtkSmartPointer<vtkImageData> ThresholdImage1 =
			vtkSmartPointer<vtkImageData>::New();
		ThresholdImage1 -> DeepCopy(MRThreshold->GetOutput());

		MRThreshold->SetInputData(MR2Image);
		MRThreshold->ThresholdByUpper(lower);
		MRThreshold->ReplaceInOn();
		MRThreshold->SetInValue(MR2Image->GetScalarTypeMax());
		MRThreshold->SetOutValue(0);
		MRThreshold->Update();
		vtkSmartPointer<vtkImageData> ThresholdImage2 =
			vtkSmartPointer<vtkImageData>::New();
		ThresholdImage2->DeepCopy(MRThreshold->GetOutput());


		//Establish pipeline connection between VTK and ITK
		vtkImageExport *exporter1;
		
		exporter1 = vtkImageExport::New();
		exporter1->SetInputData(ThresholdImage1);
		//	exporter->ImageLowerLeftOn();
		//	exporter->Export(cImage);	

		typedef itk::Image< unsigned short, 3 >  ImageType;
		typedef itk::Image< unsigned short, 3 > OutputImageType;
		typedef itk::VTKImageImport< ImageType> ImageImportType;
		ImageImportType::Pointer importer1 = ImageImportType::New();
		ImageImportType::Pointer importer2 = ImageImportType::New();
		ImageType::Pointer itkImage1 = ImageType::New();
		ImageType::Pointer labelImage1 = ImageType::New();
		
		ConnectPipelines(exporter1, importer1);
		itkImage1 = importer1->GetOutput();
		importer1->Update();
		
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
		ConnectedComponentImageFilterType::Pointer connected1 = ConnectedComponentImageFilterType::New();
		ConnectedComponentImageFilterType::Pointer connected2 = ConnectedComponentImageFilterType::New();


		double center[3];           //in ijk image coordinates
		center[0] = (Extent[1] + Extent[0]) / 2;
		center[1] = (Extent[3] + Extent[2]) / 2;
		center[2] = (Extent[5] + Extent[4]) / 2;

		double centerCT0[3];
		centerCT0[0] = (CTbounds[1] - (10 * 25));
		centerCT0[1] = (CTbounds[3] - (7 * 25));
		centerCT0[2] = (CTbounds[5] - (10 * 25));

		typedef itk::RegionOfInterestImageFilter< ImageType, ImageType > FilterType;
		//FilterType::Pointer filter = FilterType::New();
		FilterType::Pointer filter1 = FilterType::New();
		FilterType::Pointer filter2 = FilterType::New();
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
		filter1->SetInput(itkImage1);
		filter1->SetRegionOfInterest(region);
		connected1->SetInput(filter1->GetOutput());
		//connected->SetInput(itkImage1);
		connected1->Update();
		labelImage1 = connected1->GetOutput();


		//LabelImage Geometry processing
		typedef itk::LabelGeometryImageFilter< ImageType > LabelGeometryImageFilterType;
		LabelGeometryImageFilterType::Pointer labelGeometryImageFilter1 = LabelGeometryImageFilterType::New();
		labelGeometryImageFilter1->SetInput(labelImage1);
		labelGeometryImageFilter1->Update();

		LabelGeometryImageFilterType::LabelsType allLabels1;
		allLabels1 = labelGeometryImageFilter1->GetLabels();
		typedef itk::PointSet< float, 3 >   PointSetType;
		typedef PointSetType::PointType PointType;
		//typedef PointSetType::PointsContainerPointer PointsContainerPointer1;
		//typedef PointSetType::PointsContainerPointer PointsContainerPointer2;
		//PointType p;
		PointType p1;
		PointType p2;

		vtkSmartPointer<vtkPoints> MRpoints1 =
			vtkSmartPointer<vtkPoints>::New();
		vtkSmartPointer<vtkPoints> difpoints1 =
			vtkSmartPointer<vtkPoints>::New();
		/*vtkSmartPointer<vtkPoints> MRpoints2 =
			vtkSmartPointer<vtkPoints>::New();*/
		vtkSmartPointer<vtkPoints> difpoints2 =
			vtkSmartPointer<vtkPoints>::New();


		vtkIdType j = 0;
		//float pijk[3];
		double pijk1[3];
		double pijk2[3];

		vtkImageData *TestImage = ThresholdImage1;
		vtkSmartPointer<vtkImageMathematics> imageMath =
			vtkSmartPointer<vtkImageMathematics>::New();
		//imageMath->SetOperationToMultiplyByK();
		//imageMath->SetConstantK(0.0);
		imageMath->SetOperationToAdd();
		imageMath->SetInput1Data(MR1Image);
		imageMath->SetInput2Data(MR2Image);
		imageMath->Update();
		TestImage = imageMath->GetOutput();
	
		//find center control point position
		int labelValue = 1;
		//LabelGeometryImageFilterType::LabelPixelType labelValue1 = allLabels1.at(1);
		//p1 = labelGeometryImageFilter1->GetCentroid(labelValue1);
		p1 = labelGeometryImageFilter1->GetCentroid(allLabels1[labelValue]);
	
		pijk1[0] = p1[0] + start[0];
		pijk1[1] = p1[1] + start[1];
		pijk1[2] = p1[2] + start[2];
		double centerMR1[3];
		centerMR1[0] = pijk1[0] * (DirCosinesRAS(0, 0) * PxlSpacingRAS[0]) + originRAS[0];
		centerMR1[1] = pijk1[1] * (DirCosinesRAS(1, 1) * PxlSpacingRAS[1]) + originRAS[1];
		centerMR1[2] = pijk1[2] * (DirCosinesRAS(2, 2) * PxlSpacingRAS[2]) + originRAS[2];

		ImageType::Pointer itkImage2 = ImageType::New();
		ImageType::Pointer labelImage2 = ImageType::New();
		vtkImageExport *exporter2;
		exporter2 = vtkImageExport::New();
		exporter2->SetInputData(ThresholdImage2);
		ConnectPipelines(exporter2, importer2);
		itkImage2 = importer2->GetOutput();
		importer2->Update();
		filter2->SetInput(itkImage2);
		filter2->SetRegionOfInterest(region);
		connected2->SetInput(filter2->GetOutput());
		connected2->Update();
		labelImage2 = connected2->GetOutput();
		LabelGeometryImageFilterType::Pointer labelGeometryImageFilter2 = LabelGeometryImageFilterType::New();
		labelGeometryImageFilter2->SetInput(labelImage2);
		labelGeometryImageFilter2->Update();
		LabelGeometryImageFilterType::LabelsType allLabels2;
		allLabels2 = labelGeometryImageFilter2->GetLabels();

	//	LabelGeometryImageFilterType::LabelPixelType labelValue2 = allLabels2.at(1);
		//p2 = labelGeometryImageFilter2->GetCentroid(labelValue2);
		p2 = labelGeometryImageFilter2->GetCentroid(allLabels2[labelValue]);
	//	qDebug() << p1[0] << p1[1] << p1[2];
	//	qDebug() << p2[0] << p2[1] << p2[2];
		pijk2[0] = p2[0] + start[0];
		pijk2[1] = p2[1] + start[1];
		pijk2[2] = p2[2] + start[2];
		double centerMR2[3];
		centerMR2[0] = pijk2[0] * (DirCosinesRAS(0, 0) * PxlSpacingRAS[0]) + originRAS[0];
		centerMR2[1] = pijk2[1] * (DirCosinesRAS(1, 1) * PxlSpacingRAS[1]) + originRAS[1];
		centerMR2[2] = pijk2[2] * (DirCosinesRAS(2, 2) * PxlSpacingRAS[2]) + originRAS[2];

		KdTree->SetDataSet(CTpolydata);
		KdTree->BuildLocator();
		vtkIdType iD = KdTree->FindClosestPoint(centerCT0);
		double centerCT[3];
		KdTree->GetDataSet()->GetPoint(iD, centerCT);

		vtkIdType numPoints1; 
		vtkIdType numPoints2;
		double CTp[3];
		//double MRp1[3];
		//double MRp2[3];
		double difp1[3];
		double difp2[3];
		vtkSmartPointer<vtkDoubleArray> Differences1 =
			vtkSmartPointer<vtkDoubleArray>::New();
		vtkSmartPointer<vtkDoubleArray> Differences2 =
			vtkSmartPointer<vtkDoubleArray>::New();
		Differences1->SetNumberOfComponents(3); 
		Differences2->SetNumberOfComponents(3);
		//Differences1->SetNumberOfTuples(num of points);

		
	// iterate through Reference points and find corresponding MR CP positions
		j = 0;
		int check1;
		int check2;

		for (vtkIdType i = 0; i < CTpolydata->GetNumberOfPoints(); i++)
		{
			
			CTpolydata->GetPoint(i, CTp);
			CTp[0] = CTp[0] - centerCT[0] + centerMR1[0];
			CTp[1] = CTp[1] - centerCT[1] + centerMR1[1];
			CTp[2] = CTp[2] - centerCT[2] + centerMR1[2];

			/*TestImage->SetScalarComponentFromDouble(
				round(((CTp[0] - originRAS[0]) / ((DirCosinesRAS(0, 0) * PxlSpacingRAS[0])))), 
				round(((CTp[1] - originRAS[1]) / ((DirCosinesRAS(1, 1) * PxlSpacingRAS[1])))),
				round(((CTp[2] - originRAS[2]) / ((DirCosinesRAS(2, 2) * PxlSpacingRAS[2])))),
				0, ThresholdImage1->GetScalarTypeMax());*/

			//centerMR2[0] = pijk[0] * (DirCosinesRAS(0, 0) * PxlSpacingRAS[0]) + originRAS[0];
			start[0] = ((CTp[0] - originRAS[0]) / ((DirCosinesRAS(0, 0) * PxlSpacingRAS[0]))) - 10;
			start[1] = ((CTp[1] - originRAS[1]) / ((DirCosinesRAS(1, 1) * PxlSpacingRAS[1]))) - 10;
			start[2] = ((CTp[2] - originRAS[2]) / ((DirCosinesRAS(2, 2) * PxlSpacingRAS[2]))) - 5;

			end[0] = ((CTp[0] - originRAS[0]) / ((DirCosinesRAS(0, 0) * PxlSpacingRAS[0]))) + 10;
			end[1] = ((CTp[1] - originRAS[1]) / ((DirCosinesRAS(1, 1) * PxlSpacingRAS[1]))) + 10;
			end[2] = ((CTp[2] - originRAS[2]) / ((DirCosinesRAS(2, 2) * PxlSpacingRAS[2]))) + 5;

			/*qDebug() << start[0] << Extent[0];
			qDebug() << start[1] << Extent[2];
			qDebug() << start[2] << Extent[4];
			qDebug() << end[0] << Extent[1];
			qDebug() << end[1] << Extent[3];
			qDebug() << end[2] << Extent[5];*/
			if ((start[0] > Extent[0]) && (start[1] > Extent[2]) && (start[2] > Extent[4])
				&& (end[0] < Extent[1]) && (end[1] < Extent[3]) && (end[2] < Extent[5]))
			{
			
				check1 = 0;
				check2 = 0;
				region.SetIndex(start);
				region.SetUpperIndex(end);
				filter1->SetInput(itkImage1);
				filter1->SetRegionOfInterest(region);
				connected1->SetInput(filter1->GetOutput());
				connected1->Update();
				labelImage1 = connected1->GetOutput();
				numPoints1 = connected1->GetObjectCount();

				

				//LabelImage Geometry processing
				labelGeometryImageFilter1->SetInput(labelImage1);
				labelGeometryImageFilter1->Update();
				allLabels1 = labelGeometryImageFilter1->GetLabels();

				if (numPoints1 == 1){
					//LabelGeometryImageFilterType::LabelsType::iterator allLabelsIt;
					//for (allLabelsIt = allLabels1.begin(); allLabelsIt != allLabels1.end(); allLabelsIt++)
					//{
					//	labelValue1 = *allLabelsIt;
						//labelValue1 = allLabels1.at(1);
						//qDebug() << labelGeometryImageFilter->GetVolume(labelValue);
						if (((labelGeometryImageFilter1->GetVolume(allLabels1[labelValue])) < (V - 35)) ||
							((labelGeometryImageFilter1->GetVolume(allLabels1[labelValue])) > (V + 35)) ||   //could use more fine tuning
							(labelGeometryImageFilter1->GetEccentricity(allLabels1[labelValue]) > 0.9))
						{
							
							//qDebug() << labelGeometryImageFilter1->GetVolume(labelValue1);
						}
						else
						{
							check1 =+ 1;
							p1 = labelGeometryImageFilter1->GetCentroid(allLabels1[labelValue]);
							pijk1[0] = p1[0] + start[0];
							pijk1[1] = p1[1] + start[1];
							pijk1[2] = p1[2] + start[2];

							pijk1[0] = pijk1[0] * (DirCosinesRAS(0, 0) * PxlSpacingRAS[0]) + originRAS[0];
							pijk1[1] = pijk1[1] * (DirCosinesRAS(1, 1) * PxlSpacingRAS[1]) + originRAS[1];
							pijk1[2] = pijk1[2] * (DirCosinesRAS(2, 2) * PxlSpacingRAS[2]) + originRAS[2];

							difp1[0] = CTp[0] - pijk1[0];
							difp1[1] = CTp[1] - pijk1[1];
							difp1[2] = CTp[2] - pijk1[2];
							/*TestImage->SetScalarComponentFromDouble(
								round(pijk[0]),
								round(pijk[1]),
								round(pijk[2]),
								0, 100);*/
						}
					//}
				}
				else{ check1 = 0; }
				//LabelImage Geometry processing
				region.SetIndex(start);
				region.SetUpperIndex(end);
				filter2->SetInput(itkImage2);
				filter2->SetRegionOfInterest(region);
				connected2->SetInput(filter2->GetOutput());
				connected2->Update();
				labelImage2 = connected2->GetOutput();
				numPoints2 = connected2->GetObjectCount();
				labelGeometryImageFilter2->SetInput(labelImage2);
				labelGeometryImageFilter2->Update();
				allLabels2 = labelGeometryImageFilter2->GetLabels();

				if ((numPoints2 == 1) && (check1 == 1)){
				//	LabelGeometryImageFilterType::LabelsType::iterator allLabelsIt;
				//	for (allLabelsIt = allLabels2.begin(); allLabelsIt != allLabels2.end(); allLabelsIt++)
				//	{
						//labelValue2 = *allLabelsIt;
						//labelValue2 = allLabels2.at(1);
					if (((labelGeometryImageFilter2->GetVolume(allLabels2[labelValue])) < (V - 35)) ||
						((labelGeometryImageFilter2->GetVolume(allLabels2[labelValue])) > (V + 35)) ||   //could use more fine tuning
							(labelGeometryImageFilter2->GetEccentricity(allLabels2[labelValue]) > 0.9))
						{
							
						}
						else
						{
							check2 =+ 1;
							p2 = labelGeometryImageFilter2->GetCentroid(allLabels2[labelValue]);
							pijk2[0] = p2[0] + start[0];
							pijk2[1] = p2[1] + start[1];
							pijk2[2] = p2[2] + start[2];

							pijk2[0] = pijk2[0] * (DirCosinesRAS(0, 0) * PxlSpacingRAS[0]) + originRAS[0];
							pijk2[1] = pijk2[1] * (DirCosinesRAS(1, 1) * PxlSpacingRAS[1]) + originRAS[1];
							pijk2[2] = pijk2[2] * (DirCosinesRAS(2, 2) * PxlSpacingRAS[2]) + originRAS[2];

							difp2[0] = CTp[0] - pijk2[0];
							difp2[1] = CTp[1] - pijk2[1];
							difp2[2] = CTp[2] - pijk2[2];
						}
					//}
				}
				else{ check2 = 0; }
				if ((check2 == 1) && (check1 == 1))
				{
					MRpoints1->InsertNextPoint(pijk1);
					Differences1->InsertNextTupleValue(difp1);
					Differences2->InsertNextTupleValue(difp2);
					//difpoints1->InsertNextPoint(difp1);
				////	MRpoints2->InsertNextPoint(MRp2);
				//	difpoints2->InsertNextPoint(difp2);
					//qDebug() << difp1[1] << difp2[1];
					j++;
				
				}
			}
			
		}
		qDebug() << j;
		vtkSmartPointer<vtkPolyData> MR1polydata =
			vtkSmartPointer<vtkPolyData>::New();
		//vtkSmartPointer<vtkPolyData> MR2polydata =
	//		vtkSmartPointer<vtkPolyData>::New();
		MR1polydata->SetPoints(MRpoints1);
		//MR2polydata->SetPoints(MRpoints2);
		//Difference1->SetPoints(difpoints1);
		//Difference2->SetPoints(difpoints2);
	

	//Calculate MR1 Centroid positions and Difference from Reference
		//Difference1 = CalculateMRCentroids(MR1Node, CTpolydata);

	//Calculate MR2 Centroid positions and Difference from Reference
		//Difference2 = CalculateMRCentroids(MR2Node, CTpolydata);

	//Apply Mask


	//Calculate Sequence Dependent Distortion
		//and Remove Sequence Dependent Distortion from Total Distortion
		//ASSUMES READ ENCODE DIRECTIONS ARE ANTERIOR AND POSTERIOR FOR
		//MR1 AND MR2 RESPECTIVELY
		vtkSmartPointer<vtkDoubleArray> SequenceDependent =
			vtkSmartPointer<vtkDoubleArray>::New();
		vtkSmartPointer<vtkDoubleArray> GNLDist =
			vtkSmartPointer<vtkDoubleArray>::New();
		vtkSmartPointer<vtkDoubleArray>GNLDistx =
			vtkSmartPointer<vtkDoubleArray>::New();
		vtkSmartPointer<vtkDoubleArray>GNLDisty =
			vtkSmartPointer<vtkDoubleArray>::New();
		vtkSmartPointer<vtkDoubleArray>GNLDistz =
			vtkSmartPointer<vtkDoubleArray>::New();
		SequenceDependent->SetNumberOfComponents(3);
		GNLDist->SetNumberOfComponents(3);
		GNLDistx->SetNumberOfComponents(1);
		GNLDisty->SetNumberOfComponents(1);
		GNLDistz->SetNumberOfComponents(1);
		double sequencedependentp[3];
		double GNLp[3];
		
		//vtkSmartPointer<vtkPoints> seqdependentpoints =
		//	vtkSmartPointer<vtkPoints>::New();
		//vtkSmartPointer<vtkPoints> GNLpoints =
		//	vtkSmartPointer<vtkPoints>::New();
	//	qDebug() << j;
		for (vtkIdType i = 0; i < j; i++)
		{
			//double p[3];
			//Difference1->GetPoint(i, difp1);
			//Difference2->GetPoint(i, difp2);
			Differences1->GetTupleValue(i, difp1);
			Differences2->GetTupleValue(i, difp2);
			sequencedependentp[0] = 0;
			sequencedependentp[1] = ((difp1[1] - difp2[1]) / 2);
			//qDebug() << sequencedependentp[1];
			sequencedependentp[2] = 0;
			GNLp[0] = (difp1[0]);
			GNLp[1] = (difp1[1] - sequencedependentp[1]);
			GNLp[2] = (difp1[2]);
			GNLDistx->InsertComponent(i, 0, GNLp[0]);
			GNLDisty->InsertComponent(i, 0, GNLp[1]);
			GNLDistz->InsertComponent(i, 0, GNLp[2]);
			//GNLdisty[i] = GNLp[1];
			//GNLdistz[i] = GNLp[2];
			GNLDist->InsertNextTupleValue(GNLp);
			SequenceDependent->InsertNextTupleValue(sequencedependentp);
			//seqdependentpoints->InsertNextPoint(sequencedependentp);
			//GNLpoints->InsertNextPoint(GNLp);
			//qDebug() << GNLdisty[i];
		}
		
		//GNLDist->SetPoints(GNLpoints);
		//SequenceDependent->SetPoints(seqdependentpoints);
		//GNLDist->Resize(MR1polydata->GetNumberOfPoints());
		//GNLDist->Fill(1);

		//GNLDisty = Differences1y - SequenceDependent;
		//GNLDistx = Differences1x;
		//GNLDistz = Differences1z;
	//Interpolate Distortion Map
		int order = 6;
		Distortion_polyfitSVD(MR1polydata, GNLDistx, Extent, order);
		//GNLDistortionNode = Distortion_polyfitSVD(MR1polydata, GNLDisty, Extent, order);
		//GNLDistortionNode = Distortion_polyfitSVD(MR1polydata, GNLDistz, Extent, order);

	//	qDebug() << "test";
		//Output
		vtkSmartPointer<vtkXMLPolyDataWriter> writer =
			vtkSmartPointer<vtkXMLPolyDataWriter>::New();
		//writer->SetFileName("GNLDist.vtp");
		//writer->SetInputData(GNLDist);
		//writer->Write();
		writer->SetFileName("MRposition.vtp");
		writer->SetInputData(MR1polydata);
		writer->Write();
	
		MR1VolumeNode->SetAndObserveImageData(TestImage);
		GNLDistortionNode = vtkMRMLNode::SafeDownCast(MR1VolumeNode);
		//Cleanup
		exporter1->Delete();
		exporter2->Delete();

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
	vtkDoubleArray* GNLDist, int* Extent, int order){
	vtkMRMLNode *GNLDistortionNode;

	vnl_matrix<double> coeffs = Fit3DPolySVD(MRpolydata, GNLDist, order);
	double* zbar = Eval3DPolySVD(Extent, coeffs, order);

	//return GNLDistortionNode;
}
//-----------------------------------------------------------------------------
vnl_matrix<double> vtkSlicerMeasureDistortionLogic::Fit3DPolySVD(vtkPolyData* MRpolydata,
	vtkDoubleArray* GNLDist, int order){
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
		data(i, 3) = GNLDist->GetComponent(i,0);
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
	qDebug() << B.size()<<C.size()<< D.size()<< A.size();
	
	double sigma = pow(std::numeric_limits<double>::epsilon(),1/order);
	vnl_svd<double> svd(A,sigma);


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
	
	qDebug() << Extent[0] << Extent[1] << Extent[2] << Extent[3] << Extent[4] << Extent[5];
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