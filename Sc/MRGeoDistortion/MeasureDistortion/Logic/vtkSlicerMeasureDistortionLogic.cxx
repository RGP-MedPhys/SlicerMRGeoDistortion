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
#include <vtkStructuredGrid.h>
#include <vtkUnstructuredGrid.h>
#include <vtkMath.h>
#include<vtkDelaunay3D.h>
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
#include "itkOtsuThresholdImageFilter.h"






// STD includes
#include <cassert>
#include <time.h>
#include <math.h>



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
//void vtkSlicerMeasureDistortionLogic::ProcessMRMLNodesEvents(vtkObject *vtkNotUsed(caller), 
//	unsigned long event, void *callData)
//{
//	vtkDebugMacro("ProcessMRMLNodesEvents");
//
//	vtkMRMLNode* node = reinterpret_cast<vtkMRMLNode*> (callData);
//	vtkMRMLAnnotationNode* annotationNode = vtkMRMLAnnotationNode::SafeDownCast(node);
//	if (annotationNode)
//	{
//		switch (event)
//		{
//		case vtkMRMLScene::NodeAddedEvent:
//			this->OnMRMLSceneNodeAdded(annotationNode);
//			break;
//		case vtkMRMLScene::NodeRemovedEvent:
//			this->OnMRMLSceneNodeRemoved(annotationNode);
//			break;
//		}
//
//---------------------------------------------------------------------------
void vtkSlicerMeasureDistortionLogic::UpdateFromMRMLScene()
{
  assert(this->GetMRMLScene() != 0);
}

//---------------------------------------------------------------------------
void vtkSlicerMeasureDistortionLogic::OnMRMLSceneNodeAdded(vtkMRMLNode* node)
{
	
}

//---------------------------------------------------------------------------
void vtkSlicerMeasureDistortionLogic::OnMRMLSceneNodeRemoved(vtkMRMLNode* vtkNotUsed(node))
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
	ImageType::Pointer labelImage2 = ImageType::New();
	ConnectPipelines(exporter, importer);
	itkImage = importer->GetOutput();
	

	//ITK Connectivity Filter
	typedef itk::ConnectedComponentImageFilter <ImageType, OutputImageType>
		ConnectedComponentImageFilterType;
	ConnectedComponentImageFilterType::Pointer connected = ConnectedComponentImageFilterType::New();
	ConnectedComponentImageFilterType::Pointer connected2 = ConnectedComponentImageFilterType::New();
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
	LabelGeometryImageFilterType::Pointer labelGeometryImageFilter2 = LabelGeometryImageFilterType::New();
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
	PointType pijk2;
	PointType p;
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
	qDebug() << PxlSpacing[0] << PxlSpacing[1] << PxlSpacing[2];
	//p0[0] = 1; p0[1] = 1; p0[2] = 1;
	//points->SetPoint(j, p0);
	vtkImageData *TestImage = ReferenceImage;
	vtkSmartPointer<vtkImageMathematics> imageMath =
		vtkSmartPointer<vtkImageMathematics>::New();
	imageMath->SetOperationToMultiplyByK();
	imageMath->SetConstantK(0.0);
	//	imageMath->SetOperationToAdd();
	imageMath->SetInput1Data(CTImage);
	//	imageMath->SetInput2Data(MR2Image);
	imageMath->Update();
	TestImage = imageMath->GetOutput();
	
	typedef itk::RegionOfInterestImageFilter< ImageType, ImageType > FilterType;
	//FilterType::Pointer filter = FilterType::New();
	FilterType::Pointer filter = FilterType::New();
	ImageType::IndexType start;
	ImageType::IndexType end;
	double Vrange[3];
	double V1;
	Vrange[0] = 10;
	Vrange[1] = 10;
	Vrange[2] = 6;
	LabelGeometryImageFilterType::LabelPixelType labelValue2;
	ImageType::RegionType region;
	LabelGeometryImageFilterType::LabelsType::iterator allLabelsIt;
	for (allLabelsIt = allLabels.begin(); allLabelsIt != allLabels.end(); allLabelsIt++)
	{
		
		LabelGeometryImageFilterType::LabelPixelType labelValue = *allLabelsIt;
		p = labelGeometryImageFilter->GetCentroid(allLabels[labelValue]);

		start[0] = p[0] - Vrange[0];
		start[1] = p[1] - Vrange[1];
		start[2] = p[2] - Vrange[2];

		end[0] = p[0] + Vrange[0];
		end[1] = p[1] + Vrange[1];
		end[2] = p[2] + Vrange[2];

		if ((start[0] >= Extent[0]) && (start[1] >= Extent[2]) && (start[2] >= Extent[4])
			&& (end[0] <= Extent[1]) && (end[1] <= Extent[3]) && (end[2] <= Extent[5]))
		{
			int check = 0;
			region.SetIndex(start);
			region.SetUpperIndex(end);
			filter->SetInput(itkImage);
			filter->SetRegionOfInterest(region);
			filter->Update();
			connected2->SetInput(filter->GetOutput());
			connected2->Update();
			labelImage2 = connected2->GetOutput();
			labelGeometryImageFilter2->SetInput(labelImage2);
			labelGeometryImageFilter2->Update();

			LabelGeometryImageFilterType::LabelsType allLabels2 = labelGeometryImageFilter2->GetLabels();
			LabelGeometryImageFilterType::LabelsType::iterator allLabelsIt2;
			for (allLabelsIt2 = allLabels2.begin(); allLabelsIt2 != allLabels2.end(); allLabelsIt2++)
			{
				labelValue2 = *allLabelsIt2;
				pijk = labelGeometryImageFilter2->GetCentroid(allLabels2[labelValue2]);

				if (((labelGeometryImageFilter2->GetVolume(allLabels2[labelValue2])) < (80)) ||
					((labelGeometryImageFilter2->GetVolume(allLabels2[labelValue2])) > (135)) ||
					(labelGeometryImageFilter2->GetEccentricity(allLabels2[labelValue2]) > 0.92))

				{
					//	qDebug() << "Volume: " << labelGeometryImageFilter->GetVolume(labelValue);
					//	qDebug() << "Eccentricity: " << labelGeometryImageFilter->GetEccentricity(labelValue);
				}
				else
				{
					check = check + 1;
					pijk2[0] = pijk[0] + start[0];
					pijk2[1] = pijk[1] + start[1];
					pijk2[2] = pijk[2] + start[2];
					pRAS[0] = pijk2[0] * (DirCosinesRAS(0, 0) * PxlSpacing[0]) + originRAS[0];
					pRAS[1] = pijk2[1] * (DirCosinesRAS(1, 1) * PxlSpacing[1]) + originRAS[1];
					pRAS[2] = pijk2[2] * (DirCosinesRAS(2, 2) * PxlSpacing[2]) + originRAS[2];

					V1 = labelGeometryImageFilter2->GetVolume(allLabels2[labelValue2]);
				}
				

			}
			if (check == 1){
				points->InsertNextPoint(pRAS);
				j++;

				for (int tx = -Vrange[0]; tx < Vrange[0]; tx++){
					for (int ty = -Vrange[1]; ty < Vrange[1]; ty++){
						for (int tz = -Vrange[2]; tz < Vrange[2]; tz++){
							TestImage->SetScalarComponentFromDouble(
								round(pijk2[0] + tx),
								round(pijk2[1] + ty),
								round(pijk2[2] + tz),
							//	0, 100 * V1);
								0, pRAS[0]);
							/*TestImage->SetScalarComponentFromDouble(
							round(start[0] + 5 + tx),
							round(start[1] + 5 + ty),
							round(start[2] + 5),
							0, 12345);*/

						}
					}
				}
			}
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

	ReferenceVolumeNode->SetAndObserveImageData(TestImage);
	//ReferenceVolumeNode->SetAndObserveImageData(ReferenceImage);
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
	
	
	clock_t begin = clock();

	//Read Reference Points
	//	double *ptest;
		reader->SetFileName("Reference.vtp");
		reader->Update();
		CTpolydata->DeepCopy(reader->GetOutput());
		vtkMRMLScalarVolumeNode *MR1VolumeNode;
		vtkMRMLScalarVolumeNode *MR2VolumeNode;
		vtkImageData *MR1Image;
		vtkImageData *MR2Image;
		qDebug() << "test";
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


		qDebug() << "test";

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
		cpsz[0] = round(6 / PxlSpacingRAS[0]);
		cpsz[1] = round(6 / PxlSpacingRAS[1]);
		cpsz[2] = round(6 / PxlSpacingRAS[2]);
		//int V = round((3.14)*(pow(cpsz[0] / 2, 2))*cpsz[2]);
		int V = round((3.14)*(pow(cpsz[0] / 2, 2))*cpsz[2]);
		qDebug() << PxlSpacingRAS[0] << PxlSpacingRAS[1] << PxlSpacingRAS[2];
		
		double CTbounds[6];
		CTpolydata->GetBounds(CTbounds);

		MRThreshold->SetInputData(MR1Image);
	//	double lower = 120;       //1T panorama
	//	double lower = 24000;        //3T Ingenia
		double lower = 290;       //1.5T Ingenia
		MRThreshold->ThresholdByUpper(lower);
		MRThreshold->ReplaceInOn();
		MRThreshold->SetInValue(MR1Image->GetScalarTypeMax()); //1T and 1.5T
		//MRThreshold->SetInValue(1000); //3T
		MRThreshold->SetOutValue(0);
		MRThreshold->Update();
		vtkSmartPointer<vtkImageData> ThresholdImage1 =
			vtkSmartPointer<vtkImageData>::New();
		ThresholdImage1 -> DeepCopy(MRThreshold->GetOutput());

		MRThreshold->SetInputData(MR2Image);
		MRThreshold->ThresholdByUpper(lower);
		MRThreshold->ReplaceInOn();
		MRThreshold->SetInValue(MR2Image->GetScalarTypeMax());  //1T and 1.5T
		//MRThreshold->SetInValue(1000);  //3T
		MRThreshold->SetOutValue(0);
		MRThreshold->Update();
		vtkSmartPointer<vtkImageData> ThresholdImage2 =
			vtkSmartPointer<vtkImageData>::New();
		ThresholdImage2->DeepCopy(MRThreshold->GetOutput());
		qDebug() << "test";
		int t = 0;

		if (t == 1){
		
			vtkMRMLScalarVolumeNode *DistortionVolumeNodex = MR1VolumeNode;
			vtkMRMLNode *GNLDistortionNodex;

			DistortionVolumeNodex->SetAndObserveImageData(ThresholdImage1);
		
			GNLDistortionNodex = vtkMRMLNode::SafeDownCast(DistortionVolumeNodex);

		


			return GNLDistortionNodex;
		}


		//Establish pipeline connection between VTK and ITK
		vtkImageExport *exporter1;
		vtkSmartPointer<vtkImageCast> castFilter1 =
			vtkSmartPointer<vtkImageCast>::New();
		castFilter1->SetInputData(ThresholdImage1);
		castFilter1->SetOutputScalarTypeToUnsignedShort();
		castFilter1->Update();
		exporter1 = vtkImageExport::New();
		exporter1->SetInputConnection(castFilter1->GetOutputPort());
		//exporter1->SetInputData(MR1Image);
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

		qDebug() << "test";
		//qDebug() << itkImage->GetOrigin()[0] << itkImage->GetOrigin()[1] << itkImage->GetOrigin()[2];
		//ITK Connectivity Filter
		typedef itk::ConnectedComponentImageFilter <ImageType, OutputImageType >
			ConnectedComponentImageFilterType;
		ConnectedComponentImageFilterType::Pointer connected1 = ConnectedComponentImageFilterType::New();
		ConnectedComponentImageFilterType::Pointer connected2 = ConnectedComponentImageFilterType::New();


		//1T panorama
		double center[3];           //in ijk image coordinates
	//	center[0] = (Extent[1] + Extent[0]) / 2;
	//	center[1] = (Extent[3] + Extent[2]) / 2;
	//	center[2] = (Extent[5] + Extent[4]) / 2;

		double centerCT0[3];
	//	centerCT0[0] = (CTbounds[1] - (10 * 25));
	//	centerCT0[1] = (CTbounds[3] - (7 * 25));
	//	centerCT0[2] = (CTbounds[5] - (10 * 25));

		double Vrange[3];
		Vrange[0] = 11;
		Vrange[1] = 11;
		Vrange[2] = 5;
		
		//1.5T Ingenia
		center[0] = 310;
		center[1] = 352;
		center[2] = 98;

		centerCT0[0] = 1.4;
		centerCT0[1] = -30.2;
		centerCT0[2] = 341;

		//3T Ingenia
	//	center[0] = 348;
	//	center[1] = 401;
	//	center[2] = 95;

	//	centerCT0[0] = 2;
	//	centerCT0[1] = -29.5;
	//	//centerCT0[2] = 263;
	//	centerCT0[2] = -497;

	//	Vrange[0] = 13;
	//	Vrange[1] = 13;
	//	Vrange[2] = 5;

		qDebug() << "testcenterCT0" << centerCT0[0] << centerCT0[1] << centerCT0[2];
		qDebug() << "testcenter" << center[0] << center[1] << center[2];
		qDebug() << "Extent" << Extent[0] << Extent[1] << Extent[2] << Extent[3] << Extent[4] << Extent[5];

		typedef itk::RegionOfInterestImageFilter< ImageType, ImageType > FilterType;
		//FilterType::Pointer filter = FilterType::New();
		FilterType::Pointer filter1 = FilterType::New();
		FilterType::Pointer filter2 = FilterType::New();
		ImageType::IndexType start;
		ImageType::SizeType VROI;
		start[0] = center[0] - Vrange[0];
		start[1] = center[1] - Vrange[1];
		start[2] = center[2] - Vrange[2];

		ImageType::IndexType end;
		end[0] = center[0] + Vrange[0];
		end[1] = center[1] + Vrange[1];
		end[2] = center[2] + Vrange[2];

		VROI[0] = Vrange[0] * 2 + 1;
		VROI[1] = Vrange[1] * 2 + 1;
		VROI[2] = Vrange[2] * 2 + 1;


		ImageType::RegionType region;
		region.SetIndex(start);
		region.SetSize(VROI);
		//region.SetUpperIndex(end);
		
		filter1->SetRegionOfInterest(region);
		filter1->SetInput(itkImage1);	
		filter1->Update();
		typedef itk::PointSet< float, 3 >   PointSetType;
		typedef PointSetType::PointType PointType;
		PointType ROIorigin;
		ROIorigin = filter1->GetOutput()->GetOrigin();


	//	typedef itk::OtsuThresholdImageFilter <ImageType, ImageType>
	//		OFilterType;
	//	OFilterType::Pointer otsuFilter
	//		= OFilterType::New();
		//otsuFilter->SetInput(filter1->GetOutput());
		//otsuFilter->Update();
		

		/*itk::Size<3> ROISize;
		ROISize[0] = Vrange[0] * 2;
		ROISize[1] = Vrange[1] * 2;
		ROISize[2] = Vrange[2] * 2;

		itk::Index<3> ROIIndex;
		ROIIndex[0] = start[0];
		ROIIndex[1] = start[1];
		ROIIndex[2] = start[2];

		itk::ImageRegion<3> ROIRegion(ROIIndex, ROISize);*/

		connected1->SetInput(filter1->GetOutput());
		//connected1->SetInput(otsuFilter->GetOutput());
		//connected1->SetInput(itkImage1);
		//connected1->GetOutput()->SetRequestedRegion(ROIRegion);
		connected1->Update();
		labelImage1 = connected1->GetOutput();
		

		qDebug() << "test";
		//LabelImage Geometry processing
		typedef itk::LabelGeometryImageFilter< ImageType > LabelGeometryImageFilterType;
		LabelGeometryImageFilterType::Pointer labelGeometryImageFilter1 = LabelGeometryImageFilterType::New();
		labelGeometryImageFilter1->SetInput(labelImage1);
		labelGeometryImageFilter1->Update();

		LabelGeometryImageFilterType::LabelsType allLabels1;
		allLabels1 = labelGeometryImageFilter1->GetLabels();


		//typedef PointSetType::PointsContainerPointer PointsContainerPointer1;
		//typedef PointSetType::PointsContainerPointer PointsContainerPointer2;
		//PointType p;
		PointType p1;
		PointType p2;
		PointType pt;
		qDebug() << "test";
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

		//vtkImageData *TestImage = MR1Image;
		vtkSmartPointer<vtkImageData> TestImage1 =
			vtkSmartPointer<vtkImageData>::New();
		TestImage1->SetExtent(Extent[0], Extent[1], Extent[2], Extent[3], Extent[4], Extent[5]);
		TestImage1->SetOrigin(0.0, 0.0, 0.0);
		TestImage1->SetSpacing(1.0, 1.0, 1.0);
		TestImage1->AllocateScalars(VTK_DOUBLE, 1);
		vtkSmartPointer<vtkImageMathematics> imageMath =
			vtkSmartPointer<vtkImageMathematics>::New();
		imageMath->SetOperationToMultiplyByK();
		imageMath->SetConstantK(0.0);
	//	imageMath->SetOperationToAdd();
		imageMath->SetInput1Data(MR1Image);
	//	imageMath->SetInput2Data(MR2Image);
		imageMath->Update();
		//TestImage1 = imageMath->GetOutput();
		qDebug() << "test";
		//find center control point position
	//	int labelValue = 1;
		//LabelGeometryImageFilterType::LabelPixelType labelValue1 = allLabels1.at(1);
		//p1 = labelGeometryImageFilter1->GetCentroid(labelValue1);
		int check1 = 0;
		int check2 = 0;
	
		//for (int tx = -Vrange[0]; tx < Vrange[0]; tx++){
		//	for (int ty = -Vrange[1]; ty < Vrange[1]; ty++){
		//		for (int tz = -Vrange[2]; tz < Vrange[2]; tz++){
		//			//			//TestImage->SetScalarComponentFromDouble(
		//			//			//	round(MRp1ijk[0] + tx),
		//			//			//	round(MRp1ijk[1] + ty),
		//			//			//	round(MRp1ijk[2] + tz),
		//			//			//	0, 100 * labelGeometryImageFilter1->GetVolume(labelValue));
		//			TestImage1->SetScalarComponentFromDouble(
		//				round(start[0] + Vrange[0] + tx),
		//				round(start[1] + Vrange[1] + ty),
		//				round(start[2] + Vrange[2] + tz),
		//				0, 100 * labelGeometryImageFilter1->GetVolume(allLabels1[1]));
		//				//				0, 100 * difp1[0]);

		//		}
		//	}
		//}

		LabelGeometryImageFilterType::LabelsType::iterator allLabelsIt;
		for (allLabelsIt = allLabels1.begin(); allLabelsIt != allLabels1.end(); allLabelsIt++)
		{
			LabelGeometryImageFilterType::LabelPixelType labelValue = *allLabelsIt;
			//qDebug() << labelGeometryImageFilter->GetVolume(labelValue);
			if (((labelGeometryImageFilter1->GetVolume(allLabels1[labelValue])) < (20)) ||
				//((labelGeometryImageFilter1->GetVolume(allLabels1[labelValue])) > (80)) //1T panorama
				((labelGeometryImageFilter1->GetVolume(allLabels1[labelValue])) > (180)) //1.5T Ingenia
			//	((labelGeometryImageFilter1->GetVolume(allLabels1[labelValue])) > (300)) //3T Ingenia
				)
			{
				qDebug() << "testv" << labelGeometryImageFilter1->GetVolume(allLabels1[labelValue]);
			}
			else{
				p1 = labelGeometryImageFilter1->GetCentroid(allLabels1[labelValue]);
				//p1 = labelGeometryImageFilter1->GetCentroid->SetRequestedRegion(ROIRegion);
				qDebug() << "testcentervolume" << labelGeometryImageFilter1->GetVolume(allLabels1[labelValue]);

				pijk1[0] = p1[0] + start[0];
				pijk1[1] = p1[1] + start[1];
				pijk1[2] = p1[2] + start[2];
				//pijk1[0] = p1[0];
				//pijk1[1] = p1[1];
				//pijk1[2] = p1[2];
				check1 = check1 + 1;
			}
		}
		if (check1 != 1){
			qDebug() << "error: MR center not found"<< check1;
		}
		double centerMR1[3];
		centerMR1[0] = pijk1[0] * (DirCosinesRAS(0, 0) * PxlSpacingRAS[0]) + originRAS[0];
		centerMR1[1] = pijk1[1] * (DirCosinesRAS(1, 1) * PxlSpacingRAS[1]) + originRAS[1];
		centerMR1[2] = pijk1[2] * (DirCosinesRAS(2, 2) * PxlSpacingRAS[2]) + originRAS[2];
			
		
		qDebug() << "test" << pijk1[2] << DirCosinesRAS(2, 2) << PxlSpacingRAS[2] << originRAS[2]<< start[2]<< p1[2];
		qDebug() << "testcenterMR1(mm)" << centerMR1[0] << centerMR1[1] << centerMR1[2];
		int dbg = 0;
		if (dbg == 1){ 
			qDebug() << "testp1" << p1[0] << p1[1] << p1[2];
			qDebug() << "ROIorigin" << ROIorigin[0] << ROIorigin[1] << ROIorigin[2];
			vtkMRMLScalarVolumeNode *DistortionVolumeNodex = MR1VolumeNode;
			DistortionVolumeNodex->SetAndObserveImageData(ThresholdImage1);
			vtkMRMLNode* GNLDistortionNodex = vtkMRMLNode::SafeDownCast(DistortionVolumeNodex);
			return GNLDistortionNodex;
		}

		ImageType::Pointer itkImage2 = ImageType::New();
		ImageType::Pointer labelImage2 = ImageType::New();
		vtkSmartPointer<vtkImageCast> castFilter2 =
			vtkSmartPointer<vtkImageCast>::New();
		castFilter2->SetInputData(ThresholdImage2);
		castFilter2->SetOutputScalarTypeToUnsignedShort();
		castFilter2->Update();
		vtkImageExport *exporter2;
		exporter2 = vtkImageExport::New();
		exporter2->SetInputConnection(castFilter2->GetOutputPort());
		//exporter2->SetInputData(MR2Image);
		ConnectPipelines(exporter2, importer2);
		itkImage2 = importer2->GetOutput();
		importer2->Update();
		filter2->SetInput(itkImage2);
		filter2->SetRegionOfInterest(region);
		filter2->Update();
		//otsuFilter->SetInput(filter2->GetOutput());
		//otsuFilter->Update();
		connected2->SetInput(filter2->GetOutput());
		//connected2->SetInput(otsuFilter->GetOutput());
		connected2->Update();
		labelImage2 = connected2->GetOutput();
		LabelGeometryImageFilterType::Pointer labelGeometryImageFilter2 = LabelGeometryImageFilterType::New();
		labelGeometryImageFilter2->SetInput(labelImage2);
		labelGeometryImageFilter2->Update();
		LabelGeometryImageFilterType::LabelsType allLabels2;
		allLabels2 = labelGeometryImageFilter2->GetLabels();
		qDebug() << "test";
	//	LabelGeometryImageFilterType::LabelPixelType labelValue2 = allLabels2.at(1);
		//p2 = labelGeometryImageFilter2->GetCentroid(labelValue2);
		//LabelGeometryImageFilterType::LabelsType::iterator allLabelsIt;
		for (allLabelsIt = allLabels2.begin(); allLabelsIt != allLabels2.end(); allLabelsIt++)
		{
			LabelGeometryImageFilterType::LabelPixelType labelValue = *allLabelsIt;
			//qDebug() << labelGeometryImageFilter->GetVolume(labelValue);
			if (((labelGeometryImageFilter2->GetVolume(allLabels2[labelValue])) < (20)) ||
				//((labelGeometryImageFilter2->GetVolume(allLabels2[labelValue])) > (80)) //1T panorama
				((labelGeometryImageFilter1->GetVolume(allLabels1[labelValue])) > (180)) //1.5T Ingenia
			//	((labelGeometryImageFilter1->GetVolume(allLabels1[labelValue])) > (300)) //3T Ingenia
				)
			{
			}
			else{
				p2 = labelGeometryImageFilter2->GetCentroid(allLabels2[labelValue]);
				//	qDebug() << p1[0] << p1[1] << p1[2];
				//	qDebug() << p2[0] << p2[1] << p2[2];
				pijk2[0] = p2[0] + start[0];
				pijk2[1] = p2[1] + start[1];
				pijk2[2] = p2[2] + start[2];
				check2 = check2 + 1;
			}
		}
		if (check2 != 1){
			qDebug() << "error: MR center not found";
			
		}
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
		double MRp1ijk[3];
		//double MRp2[3];
		double difp1[3];
		double difp2[3];
		vtkSmartPointer<vtkDoubleArray> Differences1 =
			vtkSmartPointer<vtkDoubleArray>::New();
		vtkSmartPointer<vtkDoubleArray> Differences2 =
			vtkSmartPointer<vtkDoubleArray>::New();
		Differences1->SetNumberOfComponents(3); 
		Differences2->SetNumberOfComponents(3);
		/*vtkSmartPointer<vtkPoints> Differences1 =
			vtkSmartPointer<vtkPoints>::New();
		vtkSmartPointer<vtkPoints> Differences2 =
			vtkSmartPointer<vtkPoints>::New();*/
		//Differences1->SetNumberOfTuples(num of points);
		qDebug() << "testcenterCT(mm)" << CTpolydata->GetNumberOfPoints() << centerCT[0] << centerCT[1] << centerCT[2];
		
	// iterate through Reference points and find corresponding MR CP positions
		j = 0;
	//	int check1;
	//	int check2;
		//for (int tx = -5; tx < 5; tx++){
		//	for (int ty = -5; ty < 5; ty++){
		//		TestImage1->SetScalarComponentFromDouble(
		//			round(pijk1[0]),
		//			round(pijk1[1]),
		//			round(pijk1[2]),
		//			0, 1000);

		//		TestImage1->SetScalarComponentFromDouble(
		//			round(p1[0]),
		//			round(p1[1]),
		//			round(p1[2]),
		//			0, 1000);

		//		TestImage1->SetScalarComponentFromDouble(
		//			round(start[0]),
		//			round(start[1]),
		//			round(start[2]),
		//			0, 1000);
		//		/*TestImage->SetScalarComponentFromDouble(
		//			round(((centerCT[0] - originRAS[0]) / ((DirCosinesRAS(0, 0) * PxlSpacingRAS[0]))) + tx),
		//			round(((centerCT[1] - originRAS[1]) / ((DirCosinesRAS(1, 1) * PxlSpacingRAS[1]))) + ty),
		//			round(((centerCT[2] - originRAS[2]) / ((DirCosinesRAS(2, 2) * PxlSpacingRAS[2])))),
		//			0, 1000);

		//		TestImage->SetScalarComponentFromDouble(
		//			round(((centerMR1[0] - originRAS[0]) / ((DirCosinesRAS(0, 0) * PxlSpacingRAS[0]))) + tx),
		//			round(((centerMR1[1] - originRAS[1]) / ((DirCosinesRAS(1, 1) * PxlSpacingRAS[1]))) + ty),
		//			round(((centerMR1[2] - originRAS[2]) / ((DirCosinesRAS(2, 2) * PxlSpacingRAS[2])))),
		//			0, 2000);*/
		//	}
		//}
		double CTp1[3];
		double CTp2[3];
		qDebug() << originRAS[0] << DirCosinesRAS(0, 0) << PxlSpacingRAS[0];
		qDebug() << originRAS[1] << DirCosinesRAS(1, 1) << PxlSpacingRAS[1];
		qDebug() << originRAS[2] << DirCosinesRAS(2, 2) << PxlSpacingRAS[2];
		double V1;
		double E1;
		for (vtkIdType i = 0; i < CTpolydata->GetNumberOfPoints(); i++)
		{
			
			CTpolydata->GetPoint(i, CTp);

			CTp1[0] = CTp[0] - centerCT[0] + centerMR1[0];
			CTp1[1] = CTp[1] - centerCT[1] + centerMR1[1];
			CTp1[2] = CTp[2] - centerCT[2] + centerMR1[2];
			CTp2[0] = CTp[0] - centerCT[0] + centerMR2[0];
			CTp2[1] = CTp[1] - centerCT[1] + centerMR2[1];
			CTp2[2] = CTp[2] - centerCT[2] + centerMR2[2];

			/*TestImage->SetScalarComponentFromDouble(
				round(((CTp[0] - originRAS[0]) / ((DirCosinesRAS(0, 0) * PxlSpacingRAS[0])))), 
				round(((CTp[1] - originRAS[1]) / ((DirCosinesRAS(1, 1) * PxlSpacingRAS[1])))),
				round(((CTp[2] - originRAS[2]) / ((DirCosinesRAS(2, 2) * PxlSpacingRAS[2])))),
				0, ThresholdImage1->GetScalarTypeMax());*/

			//centerMR2[0] = pijk[0] * (DirCosinesRAS(0, 0) * PxlSpacingRAS[0]) + originRAS[0];
			start[0] = ((CTp1[0] - originRAS[0]) / ((DirCosinesRAS(0, 0) * PxlSpacingRAS[0]))) - Vrange[0];
			start[1] = ((CTp1[1] - originRAS[1]) / ((DirCosinesRAS(1, 1) * PxlSpacingRAS[1]))) - Vrange[1];
			start[2] = ((CTp1[2] - originRAS[2]) / ((DirCosinesRAS(2, 2) * PxlSpacingRAS[2]))) - Vrange[2];

			end[0] = ((CTp1[0] - originRAS[0]) / ((DirCosinesRAS(0, 0) * PxlSpacingRAS[0]))) + Vrange[0];
			end[1] = ((CTp1[1] - originRAS[1]) / ((DirCosinesRAS(1, 1) * PxlSpacingRAS[1]))) + Vrange[1];
			end[2] = ((CTp1[2] - originRAS[2]) / ((DirCosinesRAS(2, 2) * PxlSpacingRAS[2]))) + Vrange[2];
			
			

			/*qDebug() << start[0] << Extent[0];
			qDebug() << start[1] << Extent[2];
			qDebug() << start[2] << Extent[4];
			qDebug() << end[0] << Extent[1];
			qDebug() << end[1] << Extent[3];
			qDebug() << end[2] << Extent[5];*/
			if ((start[0] >= Extent[0]) && (start[1] >= Extent[2]) && (start[2] >= Extent[4])
				&& (end[0] <= Extent[1]) && (end[1] <= Extent[3]) && (end[2] <= Extent[5]))
				/*((start[0] > 72) && (start[1] > 20) && (start[2] > Extent[4])
				&& (end[0] < 414) && (end[1] < 430) && (end[2] < Extent[5]))*/
			{
				/*for (int tx = 0; tx < 10; tx++){
					for (int ty = 0; ty < 10; ty++){
						TestImage->SetScalarComponentFromDouble(
							round(start[0] + tx),
							round(start[1] + ty),
							round(start[2]),
							0, 1000);
					}
				}*/
			//	qDebug() << start[0] << start[1] << start[2];
				check1 = 0;
				check2 = 0;
				region.SetIndex(start);
				region.SetUpperIndex(end);
				filter1->SetInput(itkImage1);
				filter1->SetRegionOfInterest(region);
				//otsuFilter->SetInput(filter1->GetOutput());
				//otsuFilter->Update();
				//connected1->SetInput(otsuFilter->GetOutput());
				connected1->SetInput(filter1->GetOutput());
				connected1->Update();
				labelImage1 = connected1->GetOutput();
				numPoints1 = connected1->GetObjectCount();

			
				//LabelImage Geometry processing
				labelGeometryImageFilter1->SetInput(labelImage1);
				labelGeometryImageFilter1->Update();
				allLabels1 = labelGeometryImageFilter1->GetLabels();
			//	qDebug() << "test" << numPoints1 << labelGeometryImageFilter1->GetNumberOfObjects();
			//	if (numPoints1 == 1){
				int labelnum = 0;
				//	LabelGeometryImageFilterType::LabelsType::iterator allLabelsIt;
					for (allLabelsIt = allLabels1.begin(); allLabelsIt != allLabels1.end(); allLabelsIt++)
					{
						LabelGeometryImageFilterType::LabelPixelType labelValue = *allLabelsIt;
						if (numPoints1 > 1){
							pt = labelGeometryImageFilter1->GetCentroid(allLabels1[labelValue]);
						//	qDebug() << labelnum << labelGeometryImageFilter1->GetVolume(allLabels1[labelValue]) << pt[0] << pt[1] << pt[2];
							labelnum = labelnum + 1;
						}
						//labelValue=1;
						//labelValue1 = allLabels1.at(1);
						//qDebug() << labelGeometryImageFilter->GetVolume(labelValue);
						if (((labelGeometryImageFilter1->GetVolume(allLabels1[labelValue])) < (10)) || //(10 for 1T, 10 for 1.5, 100 for 3T)
							//	((labelGeometryImageFilter1->GetVolume(allLabels1[labelValue])) > (80)) //1T panorama
							((labelGeometryImageFilter1->GetVolume(allLabels1[labelValue])) > (180))  //1.5TIngenia
						//	((labelGeometryImageFilter1->GetVolume(allLabels1[labelValue])) > (400)) //3T Ingenia
					//		||   //could use more fine tuning
					//		(labelGeometryImageFilter1->GetEccentricity(allLabels1[labelValue]) > 0.92)
							)
						{
							
					//		qDebug() << V << labelGeometryImageFilter1->GetVolume(allLabels1[labelValue]);
					//		qDebug() << labelGeometryImageFilter1->GetEccentricity(allLabels1[labelValue]);
							/*p1 = labelGeometryImageFilter1->GetCentroid(allLabels1[labelValue]);
							MRp1ijk[0] = p1[0] + start[0];
							MRp1ijk[1] = p1[1] + start[1];
							MRp1ijk[2] = p1[2] + start[2];

							TestImage->SetScalarComponentFromDouble(
								round(MRp1ijk[0]),
								round(MRp1ijk[1]),
								round(MRp1ijk[2]),
								0, 1000);*/
						}
						else
						{
			
							check1 = check1 + 1;
							p1 = labelGeometryImageFilter1->GetCentroid(allLabels1[labelValue]);
							V1 = labelGeometryImageFilter1->GetVolume(allLabels1[labelValue]);
							E1 = labelGeometryImageFilter1->GetEccentricity(allLabels1[labelValue]);
							MRp1ijk[0] = p1[0] + start[0];
							MRp1ijk[1] = p1[1] + start[1];
							MRp1ijk[2] = p1[2] + start[2];

				

							/*MRp1ijk[0] = MRp1ijk[0] * PxlSpacingRAS[0];
							MRp1ijk[1] = MRp1ijk[1] * PxlSpacingRAS[1];
							MRp1ijk[2] = MRp1ijk[2] * PxlSpacingRAS[2];*/
							//qDebug() << "difp12" << CTp1[2] << originRAS[2] << MRp1ijk[2];


							difp1[0] = (((CTp1[0] - originRAS[0]) / ((DirCosinesRAS(0, 0))))) - (MRp1ijk[0] * PxlSpacingRAS[0]);
							difp1[1] = (((CTp1[1] - originRAS[1]) / ((DirCosinesRAS(1, 1))))) - (MRp1ijk[1] * PxlSpacingRAS[1]);
							difp1[2] = (((CTp1[2] - originRAS[2]) / ((DirCosinesRAS(2, 2))))) - (MRp1ijk[2] * PxlSpacingRAS[2]);

							/*difp1[0] = (((CTp1[0] - centerMR1[0]))) - (MRp1ijk[0] * PxlSpacingRAS[0] - centerMR1[0]);
							difp1[1] = (((CTp1[1] - centerMR1[1]))) - (MRp1ijk[1] * PxlSpacingRAS[1] - centerMR1[1]);
							difp1[2] = (((CTp1[2] - centerMR1[2]))) - (MRp1ijk[2] * PxlSpacingRAS[2] - centerMR1[2]);*/
							
							//pijk1[0] = MRp1ijk[0] * (DirCosinesRAS(0, 0) * PxlSpacingRAS[0]) + originRAS[0];
							//pijk1[1] = MRp1ijk[1] * (DirCosinesRAS(1, 1) * PxlSpacingRAS[1]) + originRAS[1];
							//pijk1[2] = MRp1ijk[2] * (DirCosinesRAS(2, 2) * PxlSpacingRAS[2]) + originRAS[2];

							//difp1[0] = CTp1[0] - pijk1[0];
							//difp1[1] = CTp1[1] - pijk1[1];
							//difp1[2] = CTp1[2] - pijk1[2];


							/*if (((start[0] + Vrange[0]) == 238) && ((start[1] + Vrange[1]) == 346) && ((start[2] + Vrange[2]) == 181)){
								qDebug() << "testdifp1" << difp1[0] << difp1[1] << difp1[2];
							}*/
							//qDebug() << labelGeometryImageFilter1->GetVolume(labelValue);
							//for (int tx = -Vrange[0]; tx < Vrange[0]; tx++){
							//	for (int ty = -Vrange[1]; ty < Vrange[1]; ty++){
							//		for (int tz = -Vrange[2]; tz < Vrange[2]; tz++){
							//			//TestImage->SetScalarComponentFromDouble(
							//			//	round(MRp1ijk[0] + tx),
							//			//	round(MRp1ijk[1] + ty),
							//			//	round(MRp1ijk[2] + tz),
							//			//	0, 100 * labelGeometryImageFilter1->GetVolume(labelValue));
							//			TestImage1->SetScalarComponentFromDouble(
							//				round(MRp1ijk[0] + tx),
							//				round(MRp1ijk[1] + ty),
							//				round(MRp1ijk[2] + tz),
							//				//0, 100 * labelGeometryImageFilter1->GetVolume(labelValue));
							//				0, (((CTp1[2] - originRAS[2]) / ((DirCosinesRAS(2, 2))))));

							//		}
							//	}
							//}
						}
					}
				//}
			//	else{ check1 = 0; }
				//LabelImage Geometry processing
					start[0] = ((CTp2[0] - originRAS[0]) / ((DirCosinesRAS(0, 0) * PxlSpacingRAS[0]))) - Vrange[0];
					start[1] = ((CTp2[1] - originRAS[1]) / ((DirCosinesRAS(1, 1) * PxlSpacingRAS[1]))) - Vrange[1];
					start[2] = ((CTp2[2] - originRAS[2]) / ((DirCosinesRAS(2, 2) * PxlSpacingRAS[2]))) - Vrange[2];

					end[0] = ((CTp2[0] - originRAS[0]) / ((DirCosinesRAS(0, 0) * PxlSpacingRAS[0]))) + Vrange[0];
					end[1] = ((CTp2[1] - originRAS[1]) / ((DirCosinesRAS(1, 1) * PxlSpacingRAS[1]))) + Vrange[1];
					end[2] = ((CTp2[2] - originRAS[2]) / ((DirCosinesRAS(2, 2) * PxlSpacingRAS[2]))) + Vrange[2];
					if ((start[0] >= Extent[0]) && (start[1] >= Extent[2]) && (start[2] >= Extent[4])
						&& (end[0] <= Extent[1]) && (end[1] <= Extent[3]) && (end[2] <= Extent[5]))
					{

						region.SetIndex(start);
						region.SetUpperIndex(end);
						filter2->SetInput(itkImage2);
						filter2->SetRegionOfInterest(region);
						//	otsuFilter->SetInput(filter2->GetOutput());
						//	otsuFilter->Update();
						//	connected2->SetInput(otsuFilter->GetOutput());
						connected2->SetInput(filter2->GetOutput());
						connected2->Update();
						labelImage2 = connected2->GetOutput();
						numPoints2 = connected2->GetObjectCount();
						labelGeometryImageFilter2->SetInput(labelImage2);
						labelGeometryImageFilter2->Update();
						allLabels2 = labelGeometryImageFilter2->GetLabels();
					}
					else{ check1 = 0; }
				if (/*(numPoints2 == 1) && */(check1 == 1)){
				//	LabelGeometryImageFilterType::LabelsType::iterator allLabelsIt;
					for (allLabelsIt = allLabels2.begin(); allLabelsIt != allLabels2.end(); allLabelsIt++)
					{
						LabelGeometryImageFilterType::LabelPixelType labelValue = *allLabelsIt;
					//	labelValue=1;
						//labelValue2 = allLabels2.at(1);
						if (((labelGeometryImageFilter2->GetVolume(allLabels2[labelValue])) < (10)) || //1T: 10 3T:100
							//	((labelGeometryImageFilter1->GetVolume(allLabels1[labelValue])) > (80)) //1T panorama
							((labelGeometryImageFilter1->GetVolume(allLabels2[labelValue])) > (180))  //1.5TIngenia
						//	((labelGeometryImageFilter2->GetVolume(allLabels2[labelValue])) > (400)) //3T Ingenia
						//	||   //could use more fine tuning
							//(labelGeometryImageFilter2->GetEccentricity(allLabels2[labelValue]) > 0.92)
							)
						{
							
						}
						else
						{
							check2 = check2 + 1;
							p2 = labelGeometryImageFilter2->GetCentroid(allLabels2[labelValue]);
							pijk2[0] = p2[0] + start[0];
							pijk2[1] = p2[1] + start[1];
							pijk2[2] = p2[2] + start[2];

							pijk2[0] = pijk2[0] * PxlSpacingRAS[0];
							pijk2[1] = pijk2[1] * PxlSpacingRAS[1];
							pijk2[2] = pijk2[2] * PxlSpacingRAS[2];


							difp2[0] = (((CTp2[0] - originRAS[0]) / ((DirCosinesRAS(0, 0))))) - (pijk2[0]);
							difp2[1] = (((CTp2[1] - originRAS[1]) / ((DirCosinesRAS(1, 1))))) - (pijk2[1]);
							difp2[2] = (((CTp2[2] - originRAS[2]) / ((DirCosinesRAS(2, 2))))) - (pijk2[2]);


							//pijk2[0] = pijk2[0] * (DirCosinesRAS(0, 0) * PxlSpacingRAS[0]) + originRAS[0];
							//pijk2[1] = pijk2[1] * (DirCosinesRAS(1, 1) * PxlSpacingRAS[1]) + originRAS[1];
							//pijk2[2] = pijk2[2] * (DirCosinesRAS(2, 2) * PxlSpacingRAS[2]) + originRAS[2];

						//	difp2[0] = CTp2[0] - pijk2[0];
						//	difp2[1] = CTp2[1] - pijk2[1];
						//	difp2[2] = CTp2[2] - pijk2[2];

						}
					}
				}
				else{ check2 = 0; }
				if ((check2 == 1) && (check1 == 1)) 
				{
				//	qDebug() << "test" << CTp[0] << pijk1[0];
					MRpoints1->InsertNextPoint(MRp1ijk);
				//	Differences1->InsertNextPoint(difp1);
				//	Differences2->InsertNextPoint(difp2);
					Differences1->InsertNextTupleValue(difp1);
					Differences2->InsertNextTupleValue(difp2);



					//qDebug() << "difp1" << difp1[2];
				for (int tx = -Vrange[0]; tx < Vrange[0]; tx++){
						for (int ty = -Vrange[1]; ty < Vrange[1]; ty++){
							for (int tz = -Vrange[2]; tz < Vrange[2]; tz++){
							//	TestImage1->SetScalarComponentFromDouble(
								//	round(MRp1ijk[0] + tx),
								//	round(MRp1ijk[1] + ty),
								//	round(MRp1ijk[2] + tz),
					//				0, 100 * pijk1[1]);
				//				//double ptest1 =difp1[0];
				//		//		double* pixel = static_cast<double*>(TestImage1->GetScalarPointer(round(start[0] + Vrange[0] + tx), round(start[1] + Vrange[1] + ty), round(start[2] + Vrange[2] + tz)));
				//			//	pixel[0] = difp2[2];
				//				//pixel = static_cast<double*>(TestImage1->GetScalarPointer(round(pijk2[0] + tx), round(pijk2[1] + ty), round(pijk2[2] + tz)));
				//				//pixel[0] = pijk2[2];
								TestImage1->SetScalarComponentFromDouble(
									round(start[0] + Vrange[0] + tx),
									round(start[1] + Vrange[1] + ty),
									round(start[2] + Vrange[2] + tz),
									0, 100 * V1);
								//	0,  100*difp1[1]);
				//					0, 100 * CTp1[2]);
				//				//pixel[0] = difp1[0];
				//			//	qDebug() << difp1[0];
				//				

							}
						}
					}
					
					
								//double ptest = TestImage->GetScalarComponentAsDouble(
								//	round(start[0] + Vrange[0]),
								//	round(start[1] + Vrange[1]),
								//	round(start[2] + Vrange[2]),
								//	//0, 100 * labelGeometryImageFilter1->GetVolume(labelValue));
								//	0);
					//double* pixel = static_cast<double*>(TestImage->GetScalarPointer(round(start[0] + Vrange[0]), round(start[1] + Vrange[1]), round(start[2] + Vrange[2])));
					//			if (((start[0] + Vrange[0]) == 238) && ((start[1] + Vrange[1]) == 346) && ((start[2] + Vrange[2]) == 181)){
					//				//qDebug() << start[0] << start[1] << start[2] << ptest;
					//			}
					//			if (abs(pixel[0]) > 15){
					//				qDebug() << start[0] << start[1] << start[2] << pixel[0]<< "testdifp2" << difp1[0] << difp1[1] << difp1[2];
					//			}

				//			}
				//		}
				//	}
				//	qDebug() << difp1[0] << CTp[0] << pijk1[0];
					//difpoints1->InsertNextPoint(difp1);
				////	MRpoints2->InsertNextPoint(MRp2);
				//	difpoints2->InsertNextPoint(difp2);
					//qDebug() << difp1[1] << difp2[1];
					j++;
				
				}
				else{
					/*if (check1 > 0 || check2 > 0){
						qDebug() << check1 << check2;
					}*/
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

	//Apply Mask

		qDebug() << "test";
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
	
		for (vtkIdType i = 0; i < j; i++)
		{
			Differences1->GetTupleValue(i, difp1);
			Differences2->GetTupleValue(i, difp2);
			//if (difp1[0] > 15){
			//	qDebug() << "difp1 is being naughty" << difp1[0];
			//}
			sequencedependentp[0] = 0;
			sequencedependentp[1] = ((difp1[1] - difp2[1]) / 2);
			sequencedependentp[2] = 0;
			GNLp[0] = ((difp1[0] ));
			GNLp[1] = (difp1[1] - sequencedependentp[1]);
			GNLp[2] = ((difp1[2]));
			GNLDistx->InsertComponent(i, 0, GNLp[0]);
			GNLDisty->InsertComponent(i, 0, GNLp[1]);
			GNLDistz->InsertComponent(i, 0, GNLp[2]);
			GNLDist->InsertNextTupleValue(GNLp);
			SequenceDependent->InsertNextTupleValue(sequencedependentp);
		}
	//Interpolate Distortion Map
		int order = 6;
		
		/*vtkSmartPointer<vtkImageData> DistortionImagex =
			vtkSmartPointer<vtkImageData>::New();
		vtkSmartPointer<vtkImageData> DistortionImagey =
			vtkSmartPointer<vtkImageData>::New();
		vtkSmartPointer<vtkImageData> DistortionImagez =
			vtkSmartPointer<vtkImageData>::New();*/
		double centerMR[3];
		centerMR[0] = ((centerMR1[0] - originRAS[0]) / ((DirCosinesRAS(0, 0) * PxlSpacingRAS[0])));
		centerMR[1] = ((centerMR1[1] - originRAS[1]) / ((DirCosinesRAS(1, 1) * PxlSpacingRAS[1])));
		centerMR[2] = ((centerMR1[2] - originRAS[2]) / ((DirCosinesRAS(2, 2) * PxlSpacingRAS[2])));
		qDebug() << "centerMR" << centerMR[0] << centerMR[1] << centerMR[2];
		qDebug() << CTp2[2] << originRAS[2];
		qDebug() << "test";

	//    vtkSmartPointer<vtkImageData> DistortionImagex = Distortion_polyfitSVD(MR1polydata, GNLDistx, Extent, order);
	//	vtkMRMLScalarVolumeNode *DistortionVolumeNodex = MR1VolumeNode;
	//	DistortionVolumeNodex->SetAndObserveImageData(DistortionImagex);
	//	vtkMRMLNode* GNLDistortionNodex = vtkMRMLNode::SafeDownCast(DistortionVolumeNodex);

	//	vtkSmartPointer<vtkImageData> DistortionImagey = Distortion_polyfitSVD(MR1polydata, GNLDisty, Extent, order);
	//	vtkMRMLScalarVolumeNode *DistortionVolumeNodey = MR1VolumeNode;
	//	DistortionVolumeNodey->SetAndObserveImageData(DistortionImagey);
	//	vtkMRMLNode* GNLDistortionNodey = vtkMRMLNode::SafeDownCast(DistortionVolumeNodey);

		vtkSmartPointer<vtkImageData> DistortionImagez = Distortion_polyfitSVD(MR1polydata, GNLDistz, Extent, order);
		vtkMRMLScalarVolumeNode *DistortionVolumeNodez = MR1VolumeNode;
		DistortionVolumeNodez->SetAndObserveImageData(DistortionImagez);
		vtkMRMLNode* GNLDistortionNodez = vtkMRMLNode::SafeDownCast(DistortionVolumeNodez);
		
		//Output
		vtkSmartPointer<vtkXMLPolyDataWriter> writer =
			vtkSmartPointer<vtkXMLPolyDataWriter>::New();
		//writer->SetFileName("GNLDist.vtp");
		//writer->SetInputData(GNLDist);
		//writer->Write();
		writer->SetFileName("MRposition.vtp");
		writer->SetInputData(MR1polydata);
		writer->Write();
	//	CalculateStats(GNLDistortionNodex, MR1polydata, GNLDistx, centerMR, 0);
	//	CalculateStats(GNLDistortionNodey, MR1polydata, GNLDisty, centerMR, 1);
		CalculateStats(GNLDistortionNodez, MR1polydata, GNLDistz, centerMR, 2);
	
		//Cleanup
		exporter1->Delete();
		exporter2->Delete();

		//clock_t end = clock();
		double elapsed_secs = (clock() - begin) / CLOCKS_PER_SEC;
		qDebug() << "time to execute" << elapsed_secs / 60;
//Debug--------------		
		//	vtkMRMLScalarVolumeNode *DistortionVolumeNodex = MR1VolumeNode;
		//	DistortionVolumeNodex->SetAndObserveImageData(TestImage1);
		//	vtkMRMLNode* GNLDistortionNodex = vtkMRMLNode::SafeDownCast(DistortionVolumeNodex);
		//	return GNLDistortionNodex;
			////-------------
	
		/*qSlicerMeasureDistortionModuleWidget *widget;
		AddMapstoScene(GNLDistortionNodex);*/
		//GNLDistortionNodex->AddToSceneOn();
		return GNLDistortionNodez;
}
//-----------------------------------------------------------------------------
void vtkSlicerMeasureDistortionLogic::CalculateStats(vtkMRMLNode* MRNode, vtkSmartPointer<vtkPolyData> MRpolydata, vtkSmartPointer<vtkDoubleArray> GNLDist, double centerMR[3], int dim){
	 
	
	//Initialization
	qDebug() << "dim:" << dim;
	 vtkMRMLScalarVolumeNode *DistVolumeNode;
	 vtkImageData *DistImage;
	 DistVolumeNode = vtkMRMLScalarVolumeNode::SafeDownCast(MRNode);
	 DistImage = DistVolumeNode->GetImageData();

	 double* PxlSpacing = DistVolumeNode->GetSpacing();
	 int* Extent = DistImage->GetExtent();
	 double *Origin = DistVolumeNode->GetOrigin();
	 double radius;

	 double largest = 0;
	 double largestloc = 0;

	 double closest1 = 3000; double closest2 = 3000; double closest3 = 3000;
	 double closest4 = 3000; double closest5 = 3000; double closest6 = 3000;
	 double closest7 = 3000;

	 double Num1 = 0; double Num2 = 0; double Num3 = 0;
	 double Num4 = 0; double Num5 = 0; double Num6 = 0;
	 double Num7 =0; double Num0 = 0;

	//Open File for editing
	 char* name; char* name1; char* name2; char* name3;
	 if (dim == 0){
		 name = "DistortionVsRadius_x.txt";
		 name1 = "innerx.txt";
		 name2 = "middlex.txt";
		 name3 = "outerx.txt";
	 }
	 else if (dim == 1){
		 name = "DistortionVsRadius_y.txt";
	 }
	 else if (dim == 2){
		 name = "DistortionVsRadius_z.txt";
	 }
	 else{
		 qCritical() << "dimension must be of value 0, 1, or 2";
		 return;
	 }
	 ofstream dvsrfile(name);
	// ofstream dvsrfile1(name1);
	 //ofstream dvsrfile2(name2);
	 //ofstream dvsrfile3(name3);
	

	 //Loop through entire distortion map and calculate summary statistics
	 qDebug() << Extent[0] << Extent[2] << Extent[4] << Origin[0] << Origin[1] << Origin[2];
	 qDebug() << "centermr" << centerMR[0] << centerMR[1] << centerMR[2];
	 qDebug() << PxlSpacing[0] << PxlSpacing[1] << PxlSpacing[2];
	 for (int k = Extent[4]; k <= Extent[5]; k++){
		 for (int j = Extent[2]; j <= Extent[3]; j++){
			 for (int i = Extent[0]; i <= Extent[1]; i++){
				 double* pixel = static_cast<double*>(DistImage->GetScalarPointer(i, j, k));
				 radius = sqrt(pow(PxlSpacing[0] * abs(centerMR[0] - i), 2) + 
					 pow(PxlSpacing[1] * abs(centerMR[1] - j), 2) +
					 pow(PxlSpacing[2] * abs(centerMR[2] - k), 2));
				 
				 if (isnan(pixel[0]) == 0){
					 if ((i % 5)==0){
						 dvsrfile << radius << "\t" << pixel[0] << "\t" << "\n";
					 }
					/* if (radius <= 100){
						 dvsrfile1 << radius << "\t" << pixel[0] << "\t" << "\n";
					 }
					 if ((radius > 100)&&(radius<=200)){
						 dvsrfile2 << radius << "\t" << pixel[0] << "\t" << "\n";
					 }
					 if (radius > 200){
						 dvsrfile3 << radius << "\t" << pixel[0] << "\t" << "\n";
					 }*/
					// if (radius<10){ qDebug() << "radius/distortion" << radius << pixel[0]; }

					 if (abs(pixel[0]) > largest){ largest = abs(pixel[0]); largestloc = radius; }


					 if (abs(pixel[0]) > 0){
						 Num0++;
					 }

					 if (abs(pixel[0]) >= 1){
						 if (radius < closest1){
							 closest1 = radius;
						 }
						 Num1++;
					 }

					 if (abs(pixel[0]) >= 2){
						 if (radius < closest2){
							 closest2 = radius;
						 }
						 Num2++;
					 }

					 if (abs(pixel[0]) >= 3){
						 if (radius < closest3){
							 closest3 = radius;
						 }
						 Num3++;
					 }

					 if (abs(pixel[0]) >= 4){
						 if (radius < closest4){
							 closest4 = radius;
						 }
						 Num4++;
					 }

					 if (abs(pixel[0]) >= 5){
						 if (radius < closest5){
							 closest5 = radius;
						 }
						 Num5++;
					 }

					 if (abs(pixel[0]) >= 6){
						 if (radius < closest6){
							 closest6 = radius;
						 }
						 Num6++;
					 }

					 if (abs(pixel[0]) >= 7){
						 if (radius < closest7){
							 closest7 = radius;
						 }
						 Num7++;
					 }
				 }
			 }
		 }
	 }
	 dvsrfile.close();
	// dvsrfile1.close();
	 //dvsrfile2.close();
	 //dvsrfile3.close();

	 ofstream statsfile("DistortionStats.txt", ios::out | ios::app);
	 if (statsfile.is_open())
	 {
		 if (dim == 0){
			 statsfile << "Left to Right Distortion\n";
		 }
		 if (dim == 1){
			 statsfile << "Anterior to Posterior Distortion\n";
		 }
		 if (dim == 2){
			 statsfile << "Superior to Inferior Distortion\n";
		 }
		 statsfile << "Largest Distortion:\t";
		 statsfile << largest << "\n";

		 statsfile << "Pct of voxels with distortion over 1mm:\t";
		 statsfile << Num1 / Num0 << "\n";
		 statsfile << "Pct of voxels with distortion over 2mm:\t";
		 statsfile << Num2 / Num0 << "\n";
		 statsfile << "Pct of voxels with distortion over 3mm:\t";
		 statsfile << Num3 / Num0 << "\n";
		 statsfile << "Pct of voxels with distortion over 4mm:\t";
		 statsfile << Num4 / Num0 << "\n";
		 statsfile << "Pct of voxels with distortion over 5mm:\t";
		 statsfile << Num5 / Num0 << "\n";
		 statsfile << "Pct of voxels with distortion over 6mm:\t";
		 statsfile << Num6 / Num0 << "\n";
		 statsfile << "Pct of voxels with distortion over 7mm:\t";
		 statsfile << Num7 / Num0 << "\n";
		
		 statsfile << "Closest 1mm Distortion:\t";
		 statsfile << closest1 << "\n";
		 statsfile << "Closest 2mm Distortion:\t";
		 statsfile << closest2 << "\n";
		 statsfile << "Closest 3mm Distortion:\t";
		 statsfile << closest3 << "\n";
		 statsfile << "Closest 4mm Distortion:\t";
		 statsfile << closest4 << "\n";
		 statsfile << "Closest 5mm Distortion:\t";
		 statsfile << closest5 << "\n";
		 statsfile << "Closest 6mm Distortion:\t";
		 statsfile << closest6 << "\n";
		 statsfile << "Closest 7mm Distortion:\t";
		 statsfile << closest7 << "\n\n";
		 
		 //Calculate goodness of fit
		 unsigned numPts = MRpolydata->GetNumberOfPoints();
		 vnl_vector<double> fitdif(numPts);
		 double fitstdev = 0;
		 double fitdifmean = 0;
		 double fitdifmax = 0;
		 int j = 0;
		 for (vtkIdType i = 0; i < numPts; i++){
			 double p[3];
			 MRpolydata->GetPoint(i, p);
			 double* pixel = static_cast<double*>(DistImage->GetScalarPointer(round(p[0]), round(p[1]), round(p[2])));
			 fitdif(i) = abs(GNLDist->GetComponent(i, 0) - pixel[0]);
			 if (isnan(fitdif(i)) == 0){
				 fitdifmean = fitdifmean + fitdif(i);
				 if (fitdif(i) > fitdifmax){ fitdifmax = fitdif(i); }
				 j++;
			 }
		 }
		 fitdifmean = fitdifmean / j;
		 for (vtkIdType i = 0; i < numPts; i++){
			 if (isnan(fitdif(i)) == 0){
				 fitstdev = fitstdev + (fitdif(i) - fitdifmean)*(fitdif(i) - fitdifmean);
			 }
		 }
		 qDebug() << "error " << fitdifmean << sqrt(fitstdev / j) << fitdifmax;
		 statsfile << "mean fit difference:\t" << fitdifmean << "+/-" << sqrt(fitstdev / j);
		 statsfile << "\n";
		 statsfile << "max fit difference:\t" << fitdifmax << "\n \n";
		 
		 statsfile.close();
	 }


	
//	return MRpolydata;
}
//------------------------------------------------------------------------------
vtkSmartPointer<vtkImageData> vtkSlicerMeasureDistortionLogic::Distortion_polyfitSVD(vtkSmartPointer<vtkPolyData> MRpolydata,
	vtkSmartPointer<vtkDoubleArray> GNLDistComp, int* Extent, int order){
	
	/*vtkSmartPointer<vtkImageData> DistortionImage =
		vtkSmartPointer<vtkImageData>::New();
	DistortionImage->SetExtent(Extent[0], Extent[1], Extent[2], Extent[3], Extent[4], Extent[5]);
	DistortionImage->SetOrigin(0.0, 0.0, 0.0);
	DistortionImage->SetSpacing(1.0, 1.0, 1.0);
	DistortionImage->AllocateScalars(VTK_DOUBLE, 1);

	for (vtkIdType i = 0; i < MRpolydata->GetNumberOfPoints(); i++){
		double p[3];
		MRpolydata->GetPoint(i, p);
		double difp = GNLDistComp->GetComponent(i, 0);
		for (int tx = -10; tx < 10; tx++){
			for (int ty = -10; ty < 10; ty++){
				for (int tz = -5; tz < 5; tz++){
					DistortionImage->SetScalarComponentFromDouble(
						round(p[0] + tx),
						round(p[1] + ty),
						round(p[2] + tz),
						0, 100 * difp);



				}
			}
		}
	}
	return DistortionImage;*/

	vnl_vector<double> coeffs((order + 3)*(order + 2)*(order + 1) / 6);
	qDebug() << coeffs.size();
	coeffs = Fit3DPolySVD(MRpolydata, GNLDistComp, order);
	qDebug() << "test_in" << coeffs[0] << coeffs[1] << coeffs[2] << coeffs[3] << coeffs[4] << coeffs[5];
	vtkSmartPointer<vtkImageData> DistortionImage = Eval3DPolySVD(Extent, MRpolydata, coeffs, order);
	qDebug() << "test_out" << coeffs[0] << coeffs[1] << coeffs[2] << coeffs[3] << coeffs[4] << coeffs[5];
	return DistortionImage;

}
//-----------------------------------------------------------------------------
vnl_vector<double> vtkSlicerMeasureDistortionLogic::Fit3DPolySVD(vtkSmartPointer<vtkPolyData> MRpolydata,
	vtkSmartPointer<vtkDoubleArray> GNLDist, int order){

	unsigned numCoeffs = (order + 3)*(order + 2)*(order + 1) / 6;
	unsigned numPts = MRpolydata->GetNumberOfPoints();

	//Scale
	vnl_matrix<double> data(numPts, 4);
	double bounds[6];
	double maxabs[4] = { 0, 0, 0, 0 };

	for (vtkIdType i = 0; i < numPts; i++){
		double p[3];	
		MRpolydata->GetPoint(i, p);
		data.put(i, 0, p[0]);
		data.put(i, 1, p[1]);
		data.put(i, 2, p[2]);
		data.put(i, 3, GNLDist->GetComponent(i, 0));

	//	qDebug() << GNLDist->GetComponent(i, 0) << data.get(i, 3);

		if (abs(data.get(i, 0))>maxabs[0]){ maxabs[0] = data.get(i, 0); }
		if (abs(data.get(i, 1))>maxabs[1]){ maxabs[1] = data.get(i, 1); }
		if (abs(data.get(i, 2))>maxabs[2]){ maxabs[2] = data.get(i, 2); }
		if (abs(data.get(i, 3))>maxabs[3]){ maxabs[3] = data.get(i, 3); }
		//qDebug() << "test";
	}
	vnl_matrix<double> A(numPts, numCoeffs,0.0);
	vnl_vector<double> B(numPts);
	vnl_vector<double> C(numPts);
	vnl_vector<double> D(numPts);
	int column = 0;
	for (int xpower = 0; xpower < order + 1; xpower++){
		for (int ypower = 0; ypower < order - xpower + 1; ypower++){
			for (int zpower = 0; zpower < order - xpower - ypower + 1; zpower++){
				B = vnl_vectorpow((data.get_column(0))/(maxabs[0]), xpower);
				C = vnl_vectorpow((data.get_column(1))/(maxabs[1]), ypower);
				D = vnl_vectorpow((data.get_column(2))/(maxabs[2]), zpower);
				for (int k = 0; k < numPts; k++){
					A.put(k, column, B.get(k)*C.get(k)*D.get(k));
				}
				column = column + 1;
			}
		}
	}

	double sigma = pow(std::numeric_limits<double>::epsilon(),1/order);
	vnl_svd<double> svd(A);
	vnl_matrix<double> W = svd.W();
	//vnl_matrix<double> q = W;	
	vnl_matrix<double> q(numCoeffs, numCoeffs, 0);
	vnl_matrix<double> result;
	vnl_vector<double> coeffs(numCoeffs);
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
	//	if (abs(W.get(i, i)) >= sigma){
			q.put(i, i, 1 / W.get(i, i));
//		}
	//	else{
	//		q.put(i, i, 0);
	//	}
	}		
	//vnl_matrix<double> qn(numPts, numCoeffs,0);
	if (numPts > numCoeffs){
		qDebug() << "numPts > numCoeffs" << numPts << numCoeffs << W.rows() << W.columns();
		/*for (int i = numCoeffs; i < numPts; i++){
			for (int j = 0; j < numCoeffs; j++){
				q.put(i, j, 0);
			}
		}*/
	}
	/*qDebug() << svd.V().rows() << svd.V().columns();
	qDebug() << q.transpose().rows() << q.transpose().columns();
	qDebug() << svd.U().transpose().rows() << svd.U().transpose().columns();*/
	
	comp.set_column(0, (data.get_column(3))/(maxabs[3]));
	//result = svd.V()*(q.transpose());
	//result = result*(svd.U().transpose());
	//qDebug() << result.rows() << result.columns();
	//result = result*(comp);                //double check this algorithm
	result = (svd.V()*(q.transpose())*(svd.U().transpose())*(comp));
	qDebug() << result.rows() << result.columns();
	coeffs = result.get_column(0);
	
	

	//rescale results
	int row = 0;
	for (int xpower = 0; xpower < order + 1; xpower++){
		for (int ypower = 0; ypower < order - xpower + 1; ypower++){
			for (int zpower = 0; zpower < order - xpower - ypower + 1; zpower++){
				
				coeffs.put(row, coeffs.get(row)*(pow(1/maxabs[0], xpower))*(pow(1/maxabs[1], ypower))*(pow(1/maxabs[2], zpower)) / (1/maxabs[3]));
				row = row + 1;
			}
		}
	}
		return coeffs;
}
//-----------------------------------------------------------------------------
vtkSmartPointer<vtkImageData> vtkSlicerMeasureDistortionLogic::Eval3DPolySVD(int* Extent, vtkSmartPointer<vtkPolyData> MRpolydata, vnl_vector<double> coeffs, int order){

	
	vtkSmartPointer<vtkImageData> DistortionImage =
		vtkSmartPointer<vtkImageData>::New();
	DistortionImage->SetExtent(Extent[0], Extent[1], Extent[2], Extent[3], Extent[4], Extent[5]);
	DistortionImage->SetOrigin(0.0, 0.0, 0.0);
	DistortionImage->SetSpacing(1.0, 1.0, 1.0);
	//DistortionImage->SetSpacing(0.955882, 0.955882, 2.0);
	DistortionImage->AllocateScalars(VTK_DOUBLE, 1);

	vtkSmartPointer<vtkDelaunay3D> delaunay3D = 
		vtkSmartPointer<vtkDelaunay3D>::New();
	delaunay3D->SetInputData(MRpolydata);
	delaunay3D->Update();

	double point[3];
	vtkDoubleArray* test;

	qDebug() << "test_inside" << coeffs[0] << coeffs[1] << coeffs[2] << coeffs[3] << coeffs[4] << coeffs[5];
	int ti=0;
	vtkIdType cellId;
	//vtkSmartPointer<vtkDoubleArray>celltest =
	//	vtkSmartPointer<vtkDoubleArray>::New();
	//celltest->SetNumberOfComponents(1);
		
		
	int row = 0;	
	for (int xpower = 0; xpower < order + 1; xpower++){
		for (int ypower = 0; ypower < order - xpower + 1; ypower++){
			for (int zpower = 0; zpower < order - xpower - ypower + 1; zpower++){
				
		
				for (unsigned int k = Extent[4]; k <= Extent[5]; k++){					
					for (unsigned int j = Extent[2]; j <= Extent[3]; j++){					
						for (unsigned int i = Extent[0]; i <= Extent[1]; i++){		
							
						//	double* pixel = static_cast<double*>(DistortionImage->GetScalarPointer(i, j, k));
							
								if ((row == 0)){
									point[0] = i; point[1] = j; point[2] = k;
									double pcoords[3];
									double weights[4];
									int subId;
									cellId = delaunay3D->GetOutput()->FindCell(point, NULL, 0, .1, subId, pcoords, weights);
									//celltest->InsertComponent(ti, 0, cellId);
									
									if (cellId >= 0)
									{
										DistortionImage->SetScalarComponentFromDouble(i,j,k,0, 
											((coeffs.get(row) * (pow(i, xpower))*(pow(j, ypower))*(pow(k, zpower)))));
										//pixel[0] = ((coeffs[row]*(pow(i, xpower))*(pow(j, ypower))*(pow(k, zpower))));
							
									}
									else{
										//pixel[0] = 0.0;
										DistortionImage->SetScalarComponentFromDouble(i,j,k,0, NAN);
									//
									}
								}
								else{
									double pixel0 = DistortionImage->GetScalarComponentAsDouble(i, j, k, 0);
									//if (celltest->GetComponent(ti, 0) >= 0){
									if (isnan(pixel0)==0 ){
										double old;
											
										double pixel = (pixel0 + (coeffs.get(row)*(pow(i, xpower))*(pow(j, ypower))*(pow(k, zpower))));
										//old = pixel[0];
										//pixel[0] = (old + (coeffs[row] * (pow(i, xpower))*(pow(j, ypower))*(pow(k, zpower))));
											DistortionImage->SetScalarComponentFromDouble(i, j, k, 0, pixel);
									}
									//else{
									//	//pixel[0] = 0.0;
									//	DistortionImage->SetScalarComponentFromDouble(i, j, k, 0, NAN);
									//}
								}
								//ti = ti + 1;
								
							}
				
						}
					}
				row = row + 1;
				
				}
			}
		}
	


	return DistortionImage;
}
//-----------------------------------------------------------------------------
vnl_vector<double> vtkSlicerMeasureDistortionLogic::vnl_vectorpow(vnl_vector<double> v,int p){
	if (p < 0){
		qCritical() << "Error at vnl_vectorpow: p must be >=0";
		return v;
	}
	//if (p == 0){
//		for (int j = 0; j < v.size(); j++){
	//	v.fill(1);
		//}
	//	return v;
	//}

	if (p == 1){
			return v;
	}

	vnl_vector<double> vn = v;
	//for (int i = 0; i < p-1; i++){
		for (int j = 0; j < v.size(); j++){
			//v(j) = v(j)*vn(j);
			v(j) = pow(v(j), p);
		}
	//}
	return v;
}
//---------------------------------------------
