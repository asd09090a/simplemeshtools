#include "tclap/CmdLine.h"

#include <vtkSmartPointer.h>
#include "vtkPolyDataWriter.h"
#include "vtkPolyDataReader.h"
#include <vtkDistancePolyDataFilter.h>
#include <itkSmartPointer.h>
#include <vtkTriangleFilter.h>

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
  std::string InputMesh1, InputMesh2, OutputMesh1, OutputMesh2;
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

    ValueArg<std::string> inArg1("","inmesh1","first input mesh",true,"result","string");
    cmd.add( inArg1 );

    ValueArg<std::string> inArg2("","inmesh2","second input mesh",true,"result","string");
    cmd.add( inArg2 );

    ValueArg<std::string> outArg1("","output1","first output mesh",true,"result","string");
    cmd.add( outArg1 );

    ValueArg<std::string> outArg2("","output2","second output mesh",true,"result","string");
    cmd.add( outArg2 );




    // Parse the args.
    cmd.parse( argc, argv );

    CmdLineObj.InputMesh1 = inArg1.getValue();
    CmdLineObj.InputMesh2 = inArg2.getValue();
    CmdLineObj.OutputMesh1 = outArg1.getValue();
    CmdLineObj.OutputMesh2 = outArg2.getValue();
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

  vtkSmartPointer<vtkPolyDataReader> vtkReader1 =   vtkSmartPointer<vtkPolyDataReader>::New();
  vtkReader1->SetFileName(CmdLineObj.InputMesh1.c_str());
  vtkReader1->Update();

  vtkSmartPointer<vtkPolyDataReader> vtkReader2 =   vtkSmartPointer<vtkPolyDataReader>::New();
  vtkReader2->SetFileName(CmdLineObj.InputMesh2.c_str());
  vtkReader2->Update();

  vtkSmartPointer<vtkTriangleFilter> tf1 = vtkSmartPointer<vtkTriangleFilter>::New();
  vtkSmartPointer<vtkTriangleFilter> tf2 = vtkSmartPointer<vtkTriangleFilter>::New();

  tf1->SetInputConnection(vtkReader1->GetOutputPort());
  tf2->SetInputConnection(vtkReader2->GetOutputPort());

  vtkSmartPointer<vtkDistancePolyDataFilter> distfilt = vtkSmartPointer<vtkDistancePolyDataFilter>::New();

  distfilt->ComputeSecondDistanceOn();
  distfilt->SetInputConnection(0, tf1->GetOutputPort());
  distfilt->SetInputConnection(1, tf2->GetOutputPort());

  vtkPolyDataWriter* writer1 = vtkPolyDataWriter::New();
  writer1->SetFileTypeToBinary(); //sets vtk polydata to binary
  writer1->SetInputConnection(distfilt->GetOutputPort(0));
  writer1->SetFileName(CmdLineObj.OutputMesh1.c_str());
  writer1->Write();
  
  vtkPolyDataWriter* writer2 = vtkPolyDataWriter::New();
  writer2->SetFileTypeToBinary(); //sets vtk polydata to binary
  writer2->SetInputConnection(distfilt->GetOutputPort(1));
  writer2->SetFileName(CmdLineObj.OutputMesh2.c_str());
  writer2->Write();

}

