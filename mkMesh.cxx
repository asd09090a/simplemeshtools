#include "tclap/CmdLine.h"
#include "ioutils.h"
#include "VTKMeshPipeline.h"
#include <vtkSmartPointer.h>
#include "vtkPolyDataWriter.h"
#include "vtkImageData.h"
#include <itkVTKImageExport.h>
#include <vtkImageImport.h>
#include <vtkMarchingCubes.h>

typedef class CmdLineType
{
public:
  std::string InputIm, OutputMesh;
  float  GaussSmooth, MarchThresh;
} CmdLineType;

void ParseCmdLine(int argc, char* argv[],
                  CmdLineType &CmdLineObj
  )
{
  using namespace TCLAP;
  try
    {
    // Define the command line object.
    CmdLine cmd("mkSurf", ' ', "0.9");

    ValueArg<std::string> inArg("i","input","T1 input image",true,"result","string");
    cmd.add( inArg );

    ValueArg<std::string> outArg("o","output","output mesg",true,"result","string");
    cmd.add( outArg );


    ValueArg<float> smoothArg("", "smoothing", "size of Gaussian smoothing on image (mm)", false, 2, "float");
    cmd.add(smoothArg);

    ValueArg<float> marchArg("", "marchthresh", "value of threshold used in marching cubes", 
                              false, 0.5, "float");
    cmd.add(marchArg);

    // Parse the args.
    cmd.parse( argc, argv );

    CmdLineObj.InputIm = inArg.getValue();
    CmdLineObj.OutputMesh = outArg.getValue();
    CmdLineObj.GaussSmooth = smoothArg.getValue();
    CmdLineObj.MarchThresh = marchArg.getValue();
    }
  catch (ArgException &e)  // catch any exceptions
    {
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    }
}



int main(int argc, char * argv[])
{
  CmdLineType CmdLineObj;
  ParseCmdLine(argc, argv, CmdLineObj);

  typedef itk::Image<float, 3> ImageType;
  ImageType::Pointer input = readIm<ImageType>(CmdLineObj.InputIm);
  VTKMeshPipeline<ImageType> MeshPipeline;
  MeshPipeline.SetImage(input);
  MeshPipeline.SetUseGaussianSmoothing(true);
  MeshPipeline.SetGaussianStandardDeviation(CmdLineObj.GaussSmooth);
  MeshPipeline.SetMarchThresh(CmdLineObj.MarchThresh);
  MeshPipeline.SetMeshOptions();
  vtkSmartPointer<vtkPolyData> m = vtkSmartPointer<vtkPolyData>::New();
  MeshPipeline.ComputeMesh(m);
  vtkPolyDataWriter* writer1 = vtkPolyDataWriter::New();
  writer1->SetFileTypeToBinary(); //sets vtk polydata to binary
  writer1->SetInputData(m);
  writer1->SetFileName(CmdLineObj.OutputMesh.c_str());
  writer1->Write();
  
}
