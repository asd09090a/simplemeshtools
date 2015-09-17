#include "VTKMeshClipPipeline.h"
#include <map>
#include "vtkPolyDataWriter.h"
#include <vtkSmartPointer.h>
#include <vtkImplicitVolume.h>
#include <vtkImageDataToPointSet.h>
#include <vtkClipPolyData.h>

using namespace std;

template<class ImageType>
VTKMeshClipPipeline<ImageType>
::VTKMeshClipPipeline()
{
  // Initialize all the filters involved in the transaction, but do not
  // pipe the inputs and outputs between these filters. The piping is quite
  // complicated and depends on the set of options that the user wishes to 
  // apply
  
  // Initialize the VTK Exporter
  m_VTKExporter = VTKExportType::New();
  //m_VTKExporter->ReleaseDataFlagOn();
  
  // Initialize the VTK Importer
  m_VTKImporter = vtkSmartPointer<vtkImageImport>::New();
  //m_VTKImporter->ReleaseDataFlagOn();

  // Pipe the importer into the exporter (that's a lot of code)
  m_VTKImporter->SetUpdateInformationCallback(
    m_VTKExporter->GetUpdateInformationCallback());
  m_VTKImporter->SetPipelineModifiedCallback(
    m_VTKExporter->GetPipelineModifiedCallback());
  m_VTKImporter->SetWholeExtentCallback(
    m_VTKExporter->GetWholeExtentCallback());
  m_VTKImporter->SetSpacingCallback(
    m_VTKExporter->GetSpacingCallback());
  m_VTKImporter->SetOriginCallback(
    m_VTKExporter->GetOriginCallback());
  m_VTKImporter->SetScalarTypeCallback(
    m_VTKExporter->GetScalarTypeCallback());
  m_VTKImporter->SetNumberOfComponentsCallback(
    m_VTKExporter->GetNumberOfComponentsCallback());
  m_VTKImporter->SetPropagateUpdateExtentCallback(
    m_VTKExporter->GetPropagateUpdateExtentCallback());
  m_VTKImporter->SetUpdateDataCallback(
    m_VTKExporter->GetUpdateDataCallback());
  m_VTKImporter->SetDataExtentCallback(
    m_VTKExporter->GetDataExtentCallback());
  m_VTKImporter->SetBufferPointerCallback(
    m_VTKExporter->GetBufferPointerCallback());  
  m_VTKImporter->SetCallbackUserData(
    m_VTKExporter->GetCallbackUserData());

}
template<class ImageType>
VTKMeshClipPipeline<ImageType>
::~VTKMeshClipPipeline()
{
  // Destroy the filters
  // m_VTKImporter->Delete();
  // m_VTKGaussianFilter->Delete();
  // m_PolygonSmoothingFilter->Delete();
  // m_StripperFilter->Delete();

  // m_MarchingCubesFilter->Delete();
  // m_TransformFilter->Delete();
  // m_Transform->Delete();
  // m_DecimateFilter->Delete();
}


#include <ctime>

template<class ImageType>
void
VTKMeshClipPipeline<ImageType>
::ComputeMesh(vtkPolyData *inMesh, vtkPolyData *outMesh)
{
  // Graft the polydata to the last filter in the pipeline
  //m_StripperFilter->SetOutput(outMesh);
  //m_MarchingCubesFilter->SetOutput(outMesh);
  
  // Connect importer and exporter
  m_VTKImporter->SetCallbackUserData(
    m_VTKExporter->GetCallbackUserData());

  // Update the ITK portion of the pipeline
  m_VTKExporter->SetInput(m_InputImage);
  m_VTKImporter->Modified();

  // Update the importer
  m_VTKImporter->Update();

  m_VTKImporter->Modified();


  vtkSmartPointer<vtkImplicitVolume> vol = vtkSmartPointer<vtkImplicitVolume>::New();
  vol->SetVolume(m_VTKImporter->GetOutput());

  vtkSmartPointer<vtkClipPolyData> clipper = vtkSmartPointer<vtkClipPolyData>::New();
  clipper->SetClipFunction(vol);
  clipper->SetInputData(inMesh);
  clipper->Update();
  outMesh->ShallowCopy(clipper->GetOutput());

}

template<class ImageType>
void
VTKMeshClipPipeline<ImageType>
::SetImage(ImagePointer image)
{
  // Store the image 
  m_InputImage = image;

  m_VTKExporter->SetInput(m_InputImage);
}


