// some simple binary filtering to create fiducials from the holes in
// the scalp used to mount the surgical frame.
// Threshold - small closing. Largeclosing - small closing, mask with
// manual indicator.
#include <iostream>
#include <cstdio>
#include <vector>

#include "tclap/CmdLine.h"

#include "ioutils.h"
#include <itkBinaryShapeOpeningImageFilter.h>
#include <itkGradientMagnitudeImageFilter.h>
#include <itkOtsuThresholdImageFilter.h>
#include <itkMaskImageFilter.h>
#include "itkBinaryDilateParaImageFilter.h"
#include "itkBinaryCloseParaImageFilter.h"
#include "itkBinaryOpenParaImageFilter.h"
#include "itkBinaryFillholeImageFilter.h"
#include <itkSubtractImageFilter.h>
#include <itkBinaryShapeKeepNObjectsImageFilter.h>
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
  std::string InputIm, OutputIm;
  float  smallclose, largeopen;
} CmdLineType;

bool debug=false;
std::string debugprefix="/tmp/align";
std::string debugsuffix=".nii.gz";

template <class TImage>
void writeImDbg(const typename TImage::Pointer Im, std::string filename)
{
  if (debug)
    {
    writeIm<TImage>(Im, debugprefix + "_" + filename + debugsuffix);
    }
}

void ParseCmdLine(int argc, char* argv[],
                  CmdLineType &CmdLineObj
  )
{
  using namespace TCLAP;
  try
    {
    // Define the command line object.
    CmdLine cmd("mkFiducials", ' ', "0.9");

    ValueArg<std::string> inArg("i","input","T1 input image",true,"result","string");
    cmd.add( inArg );

    ValueArg<std::string> outArg("o","output","output image",true,"result","string");
    cmd.add( outArg );



    ValueArg<float> scloseArg("", "smallclose", "size of small closing (mm)", false, 15, "float");
    cmd.add(scloseArg);

    ValueArg<float> lcloseArg("", "largeopen", "size of large opening (mm)", 
                              false, 10, "float");
    cmd.add(lcloseArg);

    SwitchArg debugArg("d", "debug", "save debug images", debug);
    cmd.add( debugArg );

    // Parse the args.
    cmd.parse( argc, argv );

    CmdLineObj.InputIm = inArg.getValue();
    CmdLineObj.OutputIm = outArg.getValue();
    debug = debugArg.getValue();
    CmdLineObj.smallclose = scloseArg.getValue();
    CmdLineObj.largeopen = lcloseArg.getValue();

    }
  catch (ArgException &e)  // catch any exceptions
    {
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    }
}

/////////////////////////////////////////////////////////////////

template <class PixType, int dimension>
void doSeg(const CmdLineType &CmdLineObj)
{
  // floating point image for spm prob maps
  typedef itk::Image<PixType, dimension> ImageType;
  typedef typename ImageType::Pointer IPtr;
  typedef typename itk::Image<unsigned char, ImageType::ImageDimension> MaskImType;
  typedef typename MaskImType::Pointer MIPtr;

  IPtr T1 = readIm<ImageType>(CmdLineObj.InputIm);

  // threshold first
  itk::Instance<itk::OtsuThresholdImageFilter <ImageType, MaskImType> > Thresh;
  Thresh->SetInput(T1);
  Thresh->SetInsideValue(0);
  Thresh->SetOutsideValue(1);

  // small closing
  itk::Instance<itk::BinaryCloseParaImageFilter<MaskImType> > smClose;
  itk::Instance<itk::BinaryOpenParaImageFilter<MaskImType> > lgOpen;

  smClose->SetInput(Thresh->GetOutput());
  smClose->SetRadius(CmdLineObj.smallclose);
  smClose->SetUseImageSpacing(true);
  // writeIm<MaskImType>(smClose->GetOutput(), CmdLineObj.OutputIm);
  // return;


  lgOpen->SetInput(smClose->GetOutput());
  lgOpen->SetRadius(CmdLineObj.largeopen);
  lgOpen->SetUseImageSpacing(true);


  itk::Instance<itk::BinaryShapeKeepNObjectsImageFilter <MaskImType> > Biggest;
  Biggest->SetInput(lgOpen->GetOutput());
  Biggest->SetForegroundValue(1);
  Biggest->SetNumberOfObjects(1);

  itk::Instance<itk::BinaryDilateParaImageFilter<MaskImType> > Dilate;
  Dilate->SetInput(Biggest->GetOutput());
  Dilate->SetRadius(5);

  writeIm<MaskImType>(Dilate->GetOutput(), CmdLineObj.OutputIm);

}

int main(int argc, char * argv[])
{
  CmdLineType CmdLineObj;
  ParseCmdLine(argc, argv, CmdLineObj);

  const int dimension = 3;
  int dim1 = 0;
  itk::ImageIOBase::IOComponentType ComponentType;
  if (!readImageInfo(CmdLineObj.InputIm, &ComponentType, &dim1)) 
    {
    std::cerr << "Failed to open " << CmdLineObj.InputIm << std::endl;
    return(EXIT_FAILURE);
    }
  if (dim1 != dimension) 
    {
      std::cerr << CmdLineObj.InputIm << "isn't 3D" << std::endl;
      return(EXIT_FAILURE);
    }
  // These tolerances are being set high because we rely on spm to pass in
  // appropriate data. Problems arise because spm uses the sform when
  // the sform and qform are different while ITK seems to use the
  // qform.
  // This code doesn't use orientation info, so we'll ignore the
  // headers.
  // We'll use spm code to copy the headers in spm style.

  itk::ImageToImageFilterCommon::SetGlobalDefaultCoordinateTolerance(1000.0);
  itk::ImageToImageFilterCommon::SetGlobalDefaultDirectionTolerance(1000.0);

  switch (ComponentType) 
    {
    case (itk::ImageIOBase::SHORT):
      doSeg<short, dimension>(CmdLineObj);
      break;
    case (itk::ImageIOBase::USHORT):
      doSeg<unsigned short, dimension>(CmdLineObj);
      break;
    case (itk::ImageIOBase::INT):
      doSeg<int, dimension>(CmdLineObj);
      break;
    default:
      doSeg<float, dimension>(CmdLineObj);
      break;
    }

  return(EXIT_SUCCESS);
}
