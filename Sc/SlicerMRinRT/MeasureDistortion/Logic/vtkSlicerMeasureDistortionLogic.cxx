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
#include <vtkImageGradient.h>
#include <vtkImageGradientMagnitude.h>
#include <vtkImageNonMaximumSuppression.h>
#include <vtkImageCast.h>
#include <vtkImageDataGeometryFilter.h>
#include <vtkImageExport.h>
#include <vtkImageOpenClose3D.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkKdTreePointLocator.h>
#include <vtkImageMathematics.h>
#include <vtkMath.h>
#include <vtkDenseArray.h>
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
	double* PxlSpacing = CTVolumeNode->GetSpacing();
	int* Extent = CTImage->GetExtent();
	double min = CTImage->GetScalarTypeMin();
	double max = CTImage->GetScalarTypeMax();
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

	
	//Thresholding
	vtkMRMLScalarVolumeNode *ReferenceVolumeNode = CTVolumeNode;
	//vtkNew<vtkImageThreshold> CTThreshold;
	vtkSmartPointer<vtkImageThreshold> CTThreshold =
		vtkSmartPointer<vtkImageThreshold>::New();
	//vtkNew<vtkThreshold> CTThreshold;
	CTThreshold->SetInputData(ReferenceImage);
	double lower = 56;
	CTThreshold->ThresholdByUpper(lower);
	CTThreshold->ReplaceInOn();
	CTThreshold->SetInValue(ReferenceImage->GetScalarTypeMax());
	CTThreshold->SetOutValue(ReferenceImage->GetScalarTypeMin());
	CTThreshold->Update();
	ReferenceImage = CTThreshold->GetOutput();


	//Remove background with morphological Open/Close
	vtkSmartPointer<vtkImageOpenClose3D> openClose =
		vtkSmartPointer<vtkImageOpenClose3D>::New();
	//vtkImageOpenClose3D* openClose;
	openClose->SetInputData(ReferenceImage);
	openClose->SetOpenValue(ReferenceImage->GetScalarTypeMin());
	openClose->SetCloseValue(ReferenceImage->GetScalarTypeMax());
	openClose->SetKernelSize(5, 5, 3);
	//openClose->ReleaseDataFlagOff();
	openClose->Update();
	ReferenceImage = openClose->GetOutput();
	//openClose->GetCloseValue();
	//openClose->GetOpenValue();
	qDebug() << "test";


	//Establish pipeline connection between VTK and ITK
	vtkImageExport *exporter;
	exporter = vtkImageExport::New();
	exporter->SetInputData(ReferenceImage);
//	exporter->ImageLowerLeftOn();
//	exporter->Export(cImage);		
	typedef itk::Image< unsigned short, 3 >  ImageType;
	typedef itk::Image< unsigned short, 3 > OutputImageType;
	typedef itk::VTKImageImport< ImageType> ImageImportType;
	ImageImportType::Pointer importer = ImageImportType::New();
	ImageType::Pointer itkImage=ImageType::New();
	ImageType::Pointer labelImage = ImageType::New();;
	ConnectPipelines(exporter, importer);
	itkImage = importer->GetOutput();


	//ITK Connectivity Filter
	typedef itk::ConnectedComponentImageFilter <ImageType, OutputImageType >
		ConnectedComponentImageFilterType;
	ConnectedComponentImageFilterType::Pointer connected = ConnectedComponentImageFilterType::New();;
	connected->SetInput(itkImage);
	connected->Update();
	labelImage=connected->GetOutput();
	//qDebug() << connected->GetObjectCount();
	vtkIdType numPoints = connected->GetObjectCount();


	//LabelImage Geometry processing
	typedef itk::LabelGeometryImageFilter< ImageType > LabelGeometryImageFilterType;
	LabelGeometryImageFilterType::Pointer labelGeometryImageFilter = LabelGeometryImageFilterType::New();
	labelGeometryImageFilter->SetInput(labelImage);
	//labelGeometryImageFilter->SetIntensityInput(itkImage);
	labelGeometryImageFilter->Update();
	LabelGeometryImageFilterType::LabelsType allLabels =
		labelGeometryImageFilter->GetLabels();
	typedef itk::PointSet< float ,3 >   PointSetType;
	typedef PointSetType::PointType PointType;
	typedef PointSetType::PointsContainerPointer PointsContainerPointer;
	//PointSetType::Pointer  PointSet = PointSetType::New();;
	//PointsContainerPointer  points = PointSet->GetPoints();
	PointType p;
	
	//vtkPointSet* ReferencePoints;
	vtkSmartPointer<vtkPoints> points =
		vtkSmartPointer<vtkPoints>::New();
	//points->SetNumberOfPoints(1);
//	typedef itk::Mesh< float, 3 >   MeshType;
	//MeshType::Pointer  mesh = MeshType::New();
	qDebug() << numPoints;
	//points->SetNumberOfPoints(numPoints);
	//qDebug() << "test1";


	vtkIdType j = 0;
	float p0[3];
	
	//p0[0] = 1; p0[1] = 1; p0[2] = 1;
	//points->SetPoint(j, p0);
	
	LabelGeometryImageFilterType::LabelsType::iterator allLabelsIt;
	for (allLabelsIt = allLabels.begin(); allLabelsIt != allLabels.end(); allLabelsIt++)
	{	
		LabelGeometryImageFilterType::LabelPixelType labelValue = *allLabelsIt;
		p=labelGeometryImageFilter->GetCentroid(labelValue);
		
		
		if (((labelGeometryImageFilter->GetVolume(labelValue)) < 26) || (labelGeometryImageFilter->GetEccentricity(labelValue) > 0.8))
		{
		}
		else
		{
			
			qDebug() << "Volume: " << labelGeometryImageFilter->GetVolume(labelValue);
			qDebug() << "Eccentricity: " << labelGeometryImageFilter->GetEccentricity(labelValue);
			qDebug() << "j" << j;
			qDebug() << "point0:" << p[0];
			qDebug() << "point1:" << p[1];
			qDebug() << "point2:" << p[2];
			p0[0] = p[0];
			p0[1] = p[1];
			p0[2] = p[2];
			points->InsertNextPoint(p0);
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
		vtkImageData *ThresholdImage;
		MR1VolumeNode = vtkMRMLScalarVolumeNode::SafeDownCast(MR1Node);
		MR1Image = MR1VolumeNode->GetImageData();


	//Find MR Center Control Points (closest CP to isocenter)
		vtkSmartPointer<vtkKdTreePointLocator> KdTree =
			vtkSmartPointer<vtkKdTreePointLocator>::New();
		vtkSmartPointer<vtkImageThreshold> MRThreshold =
			vtkSmartPointer<vtkImageThreshold>::New();

		MRThreshold->SetInputData(MR1Image);
		double lower = 100;
		MRThreshold->ThresholdByUpper(lower);
		MRThreshold->ReplaceInOn();
		MRThreshold->SetInValue(MR1Image->GetScalarTypeMax());
		MRThreshold->SetOutValue(MR1Image->GetScalarTypeMin());
		MRThreshold->Update();
		ThresholdImage = MRThreshold->GetOutput();

		
		
		double iso[3];
		MR1VolumeNode->GetOrigin(iso);
		double* PxlSpacing = MR1VolumeNode->GetSpacing();
		int* Extent = MR1Image->GetExtent();

		//Remove background with morphological Open/Close
		vtkSmartPointer<vtkImageOpenClose3D> openClose =
			vtkSmartPointer<vtkImageOpenClose3D>::New();
		//vtkImageOpenClose3D* openClose;
		openClose->SetInputData(ThresholdImage);
		openClose->SetOpenValue(ThresholdImage->GetScalarTypeMin());
		openClose->SetCloseValue(ThresholdImage->GetScalarTypeMax());
		openClose->SetKernelSize(5, 5, 3);
		//openClose->ReleaseDataFlagOff();
		openClose->Update();
		ThresholdImage = openClose->GetOutput();


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
		ImageType::Pointer labelImage = ImageType::New();;
		ConnectPipelines(exporter, importer);
		itkImage = importer->GetOutput();


		//ITK Connectivity Filter
		typedef itk::ConnectedComponentImageFilter <ImageType, OutputImageType >
			ConnectedComponentImageFilterType;
		ConnectedComponentImageFilterType::Pointer connected = ConnectedComponentImageFilterType::New();;
		connected->SetInput(itkImage);
		connected->Update();
		labelImage = connected->GetOutput();
		//qDebug() << connected->GetObjectCount();
		vtkIdType numPoints = connected->GetObjectCount();


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
		//PointSetType::Pointer  PointSet = PointSetType::New();;
		//PointsContainerPointer  points = PointSet->GetPoints();
		PointType p;

		//vtkPointSet* ReferencePoints;
		vtkSmartPointer<vtkPoints> points =
			vtkSmartPointer<vtkPoints>::New();
		//points->SetNumberOfPoints(1);
		//	typedef itk::Mesh< float, 3 >   MeshType;
		//MeshType::Pointer  mesh = MeshType::New();
		qDebug() << numPoints;
		//points->SetNumberOfPoints(numPoints);
		//qDebug() << "test1";


		vtkIdType j = 0;
		float p0[3];

		//p0[0] = 1; p0[1] = 1; p0[2] = 1;
		//points->SetPoint(j, p0);

		LabelGeometryImageFilterType::LabelsType::iterator allLabelsIt;
		for (allLabelsIt = allLabels.begin(); allLabelsIt != allLabels.end(); allLabelsIt++)
		{
			LabelGeometryImageFilterType::LabelPixelType labelValue = *allLabelsIt;
			p = labelGeometryImageFilter->GetCentroid(labelValue);


			if (((labelGeometryImageFilter->GetVolume(labelValue)) < 26) || (labelGeometryImageFilter->GetEccentricity(labelValue) > 0.8))
			{
			}
			else
			{

				qDebug() << "Volume: " << labelGeometryImageFilter->GetVolume(labelValue);
				qDebug() << "Eccentricity: " << labelGeometryImageFilter->GetEccentricity(labelValue);
				qDebug() << "j" << j;
				qDebug() << "point0:" << p[0];
				qDebug() << "point1:" << p[1];
				qDebug() << "point2:" << p[2];
				p0[0] = p[0];
				p0[1] = p[1];
				p0[2] = p[2];
				points->InsertNextPoint(p0);
				j++;
			}

		}

		//Output
		vtkSmartPointer<vtkXMLPolyDataWriter> writer =
			vtkSmartPointer<vtkXMLPolyDataWriter>::New();
		vtkSmartPointer<vtkPolyData> MR1polydata =
			vtkSmartPointer<vtkPolyData>::New();

		// Find center of image
		int center[2];
		center[0] = (Extent[1] + Extent[0]) / 2;
		center[1] = (Extent[3] + Extent[2]) / 2;
		//center[2] = (Extent[5] + Extent[4]) / 2;
		// Pick a radius for the circle
		int radius = 5;
		vtkImageData *TestImage=ThresholdImage;
		vtkSmartPointer<vtkImageMathematics> imageMath =
			vtkSmartPointer<vtkImageMathematics>::New();
		imageMath->SetOperationToMultiplyByK();
		imageMath->SetConstantK(0.0);
		imageMath->SetInputData(ThresholdImage);
		imageMath->Update();
		TestImage = imageMath->GetOutput();
		MR1polydata->SetPoints(points);

		//for (vtkIdType i = 0; i < points->GetNumberOfPoints(); i++)
		//{
		//	double p[3];
		//	MR1polydata->GetPoint(i, p);
		//	TestImage->SetScalarComponentFromDouble(round(p[0]), round(p[1]), round(p[2]),0, ThresholdImage->GetScalarTypeMax());


		//}

		//ThresholdImage = TestImage;
		double bounds[6];
		MR1polydata->GetBounds(bounds);
		//KdTree->SetDataSet();



		qDebug() << "iso0" << iso[0];
		qDebug() << "iso1" << iso[1];
		qDebug() << "iso2" << iso[2];
		qDebug() << "bounds0" << bounds[0];
		qDebug() << "bounds1" << bounds[1];
		qDebug() << "bounds2" << bounds[2];
		qDebug() << "bounds3" << bounds[3];
		qDebug() << "bounds4" << bounds[4];
		qDebug() << "bounds5" << bounds[5];

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
		Distortion_polyfitSVD(MR1polydata, GNLDist, Extent, order);
		//GNLDistortionNodex = SVD_polyfit();
		//GNLDistortionNodey = SVD_polyfit();
		//GNLDistortionNodez = SVD_polyfit();


		//Cleanup
		exporter->Delete();
	
		MR1VolumeNode->SetAndObserveImageData(ThresholdImage);
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