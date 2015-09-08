#include "ioutils.h"
#include "VTKMeshPipeline.h"
#include <vtkSmartPointer.h>
#include "vtkPolyDataWriter.h"
#include "vtkImageData.h"
#include <itkVTKImageExport.h>
#include <vtkImageImport.h>
#include <vtkMarchingCubes.h>


int main(int argc, char * argv[])
{
  typedef itk::Image<float, 3> ImageType;
  ImageType::Pointer input = readIm<ImageType>(argv[1]);
  VTKMeshPipeline<ImageType> MeshPipeline;
  MeshPipeline.SetImage(input);
  MeshPipeline.SetUseGaussianSmoothing(true);
  MeshPipeline.SetMeshOptions();
  vtkSmartPointer<vtkPolyData> m = vtkSmartPointer<vtkPolyData>::New();
  MeshPipeline.ComputeMesh(m);
  vtkPolyDataWriter* writer1 = vtkPolyDataWriter::New();
  writer1->SetFileTypeToBinary(); //sets vtk polydata to binary
  writer1->SetInputData(m);
  writer1->SetFileName(argv[2]);
  writer1->Write();
  
}
