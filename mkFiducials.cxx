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
  std::string InputIm, OutputIm, InputMaskIm;
  float markerdilate, smallclose, largeclose;
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

    ValueArg<std::string> maskArg("m","mask","input fiducial mask (points in the middle of the hole)",false,"","string");
    cmd.add( maskArg );

    ValueArg<float> dilateArg("", "fdilate", "size of dilation applied to fiducal mask (mm)", false, 20, "float");
    cmd.add(dilateArg);

    ValueArg<float> scloseArg("", "smallclose", "size of small closing (mm)", false, 3, "float");
    cmd.add(scloseArg);

    ValueArg<float> lcloseArg("", "largeclose", "size of large closing (mm)", 
                              false, 20, "float");
    cmd.add(lcloseArg);

    SwitchArg debugArg("d", "debug", "save debug images", debug);
    cmd.add( debugArg );

    // Parse the args.
    cmd.parse( argc, argv );

    CmdLineObj.InputIm = inArg.getValue();
    CmdLineObj.OutputIm = outArg.getValue();
    CmdLineObj.InputMaskIm = maskArg.getValue();
    debug = debugArg.getValue();
    CmdLineObj.markerdilate = dilateArg.getValue();
    CmdLineObj.smallclose = scloseArg.getValue();
    CmdLineObj.largeclose = lcloseArg.getValue();

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

  MIPtr fmarker = readIm<MaskImType>(CmdLineObj.InputMaskIm);

  // threshold first
  itk::Instance<itk::OtsuThresholdImageFilter <ImageType, MaskImType> > Thresh;
  Thresh->SetInput(T1);
  Thresh->SetInsideValue(0);
  Thresh->SetOutsideValue(1);

  // small closing
  itk::Instance<itk::BinaryCloseParaImageFilter<MaskImType> > smClose;
  itk::Instance<itk::BinaryCloseParaImageFilter<MaskImType> > lgClose;

  smClose->SetInput(Thresh->GetOutput());
  smClose->SetRadius(CmdLineObj.smallclose);
  smClose->SetUseImageSpacing(true);
  // writeIm<MaskImType>(smClose->GetOutput(), CmdLineObj.OutputIm);
  // return;

  //itk::Instance<itk::BinaryOpenParaImageFilter<MaskImType> > firstsmOpen;

  // firstsmOpen->SetInput(smClose->GetOutput());
  // firstsmOpen->SetRadius(CmdLineObj.smallclose);
  // firstsmOpen->SetUseImageSpacing(true);

  lgClose->SetInput(smClose->GetOutput());
  lgClose->SetRadius(CmdLineObj.largeclose);
  lgClose->SetUseImageSpacing(true);

  itk::Instance<itk::SubtractImageFilter<MaskImType, MaskImType, MaskImType> > Diff;
  Diff->SetInput(lgClose->GetOutput());
  Diff->SetInput2(smClose->GetOutput());

  // mask with the rough fiducials
  //itk::Instance<itk::BinaryDilateParaImageFilter<MaskImType> > Dilater;
  //Dilater->SetInput(fmarker);
  //Dilater->SetRadius(CmdLineObj.markerdilate);
  
  itk::Instance<itk::MaskImageFilter<MaskImType, MaskImType> > Masker;
  Masker->SetInput(Diff->GetOutput());
  //Masker->SetInput2(Dilater->GetOutput());
  Masker->SetInput2(fmarker);

  // // little opening on the result
  // itk::Instance<itk::BinaryOpenParaImageFilter<MaskImType> > smOpen;
  // smOpen->SetInput(Masker->GetOutput());
  // smOpen->SetRadius(1);
  // smOpen->SetUseImageSpacing(true);
  // smOpen->SafeBorderOff();

  itk::Instance<itk::BinaryShapeKeepNObjectsImageFilter <MaskImType> > Biggest;
  Biggest->SetInput(Masker->GetOutput());
  Biggest->SetForegroundValue(1);
  Biggest->SetNumberOfObjects(3);

  itk::Instance<itk::BinaryShapeOpeningImageFilter<MaskImType> > SizeFilt;
  SizeFilt->SetInput(Biggest->GetOutput());
  SizeFilt->SetLambda(10);
  SizeFilt->SetForegroundValue(1);

  // write out the centroids
  typedef typename itk::BinaryImageToShapeLabelMapFilter<MaskImType> LabellerType;

  typename LabellerType::Pointer Labeller = LabellerType::New();
  Labeller->SetInput(SizeFilt->GetOutput());
  Labeller->SetInputForegroundValue(1);
  Labeller->SetFullyConnected(true);

  typedef typename LabellerType::OutputImageType LabelMapType;
  typedef typename LabellerType::OutputImagePointer LabelMapPointerType;
  typedef typename LabelMapType::LabelObjectType LabelObjectType;
  LabelMapPointerType labmap = Labeller->GetOutput();
  labmap->Update();
  //std::cout << "Objects: " << labmap->GetNumberOfLabelObjects() << std::endl;

  for(unsigned int i = 0; i < labmap->GetNumberOfLabelObjects(); i++)
    {
    // Get the ith region
    LabelObjectType* labelObject = labmap->GetNthLabelObject(i);
    std::cout << labelObject->GetCentroid() << std::endl;
    }
  writeIm<MaskImType>(Biggest->GetOutput(), CmdLineObj.OutputIm);

  writeImDbg<MaskImType>(Thresh->GetOutput(), "pinthresh");
  writeImDbg<MaskImType>(smClose->GetOutput(), "pinsmclose");
  writeImDbg<MaskImType>(lgClose->GetOutput(), "pinlgclose");
  writeImDbg<MaskImType>(Diff->GetOutput(), "pindiff");
  writeImDbg<MaskImType>(Masker->GetOutput(), "pinmask");
  writeImDbg<MaskImType>(Biggest->GetOutput(), "pinbiggest");

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
