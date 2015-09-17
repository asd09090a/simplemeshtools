#include "tclap/CmdLine.h"
#include "ioutils.h"
#include "itkBinaryThresholdImageFilter.h"

#include "VTKMeshClipPipeline.h"
#include <vtkSmartPointer.h>
#include "vtkPolyDataWriter.h"
#include "vtkPolyDataReader.h"
#include "vtkImageData.h"
#include <itkVTKImageExport.h>
#include <vtkImageImport.h>

#include <itkSmartPointer.h>
namespace itk
{
template <typename T>
class Instance : public T::Pointer {
public:
  Instance() : SmartPointer<T>( T::New() ) {}
};
}


typedef class CmdLineType
{
public:
  std::string InputIm, InputMesh, OutputMesh;
  float  thresh;
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

    ValueArg<std::string> inMArg("m","meshinput","The mesh to be clipped",true,"result","string");
    cmd.add( inMArg );

    ValueArg<std::string> outArg("o","output","output mesh",true,"result","string");
    cmd.add( outArg );

    ValueArg<float> tArg("", "threshold", "T1 threshold used to clip surgical bits", false, 100, "float");
    cmd.add(tArg);


    // Parse the args.
    cmd.parse( argc, argv );

    CmdLineObj.InputIm = inArg.getValue();
    CmdLineObj.InputMesh = inMArg.getValue();
    CmdLineObj.OutputMesh = outArg.getValue();
    CmdLineObj.thresh = tArg.getValue();
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
  typedef itk::Image<unsigned char, 3> MaskImType;

  // Image loading and thresholding
  ImageType::Pointer input = readIm<ImageType>(CmdLineObj.InputIm);
  itk::Instance<itk::BinaryThresholdImageFilter<ImageType,MaskImType> > Thresh;
  Thresh->SetInput(input);
  Thresh->SetUpperThreshold(CmdLineObj.thresh);
  Thresh->SetInsideValue(0);
  Thresh->SetOutsideValue(1);

  vtkSmartPointer<vtkPolyDataReader> vtkReader =   vtkSmartPointer<vtkPolyDataReader>::New();
  vtkReader->SetFileName(CmdLineObj.InputMesh.c_str());
  vtkReader->Update();

  VTKMeshClipPipeline<MaskImType> MeshPipeline;
  MeshPipeline.SetImage(Thresh->GetOutput());
  vtkSmartPointer<vtkPolyData> m = vtkSmartPointer<vtkPolyData>::New();
  MeshPipeline.ComputeMesh(vtkReader->GetOutput(), m);
  vtkPolyDataWriter* writer1 = vtkPolyDataWriter::New();
  writer1->SetFileTypeToBinary(); //sets vtk polydata to binary
  writer1->SetInputData(m);
  writer1->SetFileName(CmdLineObj.OutputMesh.c_str());
  writer1->Write();
  
}
