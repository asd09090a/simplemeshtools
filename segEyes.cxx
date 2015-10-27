// segmentation of eyeballs/lenses give a marker image
// with approximate centres marked
#include <iostream>
#include <cstdio>
#include <vector>

#include "tclap/CmdLine.h"

#include "ioutils.h"
#include <itkBinaryShapeOpeningImageFilter.h>
#include <itkGradientMagnitudeImageFilter.h>
#include <itkOtsuThresholdImageFilter.h>
#include <itkTriangleThresholdImageFilter.h>
#include <itkMaskImageFilter.h>
#include "itkBinaryDilateParaImageFilter.h"
#include "itkBinaryErodeParaImageFilter.h"
#include "itkBinaryCloseParaImageFilter.h"
#include "itkBinaryOpenParaImageFilter.h"
#include "itkBinaryFillholeImageFilter.h"
#include <itkSubtractImageFilter.h>
#include <itkBinaryShapeKeepNObjectsImageFilter.h>
#include <itkBinaryReconstructionByDilationImageFilter.h>

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
  std::string StructuralIm, OutputIm, InputEyeIm;
  float markerdilate;
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
    CmdLine cmd("segEyes", ' ', "0.9");

    ValueArg<std::string> inArg("i","input","T1 input image",true,"result","string");
    cmd.add( inArg );

    ValueArg<std::string> outArg("o","output","output image",true,"result","string");
    cmd.add( outArg );

    ValueArg<std::string> maskArg("m","eye","input eye mask (points in the middle of the eye)",false,"","string");
    cmd.add( maskArg );

    ValueArg<float> dilateArg("", "fdilate", "rough size of eyeball - should be an over estimate", false, 20, "float");
    cmd.add(dilateArg);

    SwitchArg debugArg("d", "debug", "save debug images", debug);
    cmd.add( debugArg );

    // Parse the args.
    cmd.parse( argc, argv );

    CmdLineObj.StructuralIm = inArg.getValue();
    CmdLineObj.OutputIm = outArg.getValue();
    CmdLineObj.InputEyeIm = maskArg.getValue();
    debug = debugArg.getValue();
    CmdLineObj.markerdilate = dilateArg.getValue();

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

  IPtr T1 = readIm<ImageType>(CmdLineObj.StructuralIm);

  MIPtr emarker = readIm<MaskImType>(CmdLineObj.InputEyeIm);

  // dilate the markers
  itk::Instance<itk::BinaryDilateParaImageFilter <MaskImType> > DilEye ;
  DilEye->SetInput(emarker);
  DilEye->SetRadius(CmdLineObj.markerdilate);
  DilEye->SetUseImageSpacing(true);

  itk::Instance<itk::TriangleThresholdImageFilter<ImageType, MaskImType, MaskImType> > Thresh;
  Thresh->SetInput(T1);
  Thresh->SetMaskImage(DilEye->GetOutput());
  Thresh->SetMaskValue(1);
  Thresh->SetInsideValue(1);
  Thresh->SetOutsideValue(0);

  // erode
  itk::Instance<itk::BinaryErodeParaImageFilter <MaskImType> > EroEye ;
  EroEye->SetInput(Thresh->GetOutput());
  EroEye->SetRadius(3);
  
  // reconstruct from the seed points to clean up
  itk::Instance<itk::BinaryReconstructionByDilationImageFilter<MaskImType> > BRecon;
  BRecon->SetMaskImage(EroEye->GetOutput());
  BRecon->SetMarkerImage(emarker);
  BRecon->SetBackgroundValue(0);
  BRecon->SetForegroundValue(1);

  itk::Instance<itk::BinaryDilateParaImageFilter <MaskImType> > DilEye2 ;
  DilEye2->SetInput(BRecon->GetOutput());
  DilEye2->SetRadius(3);

  writeIm<MaskImType>(DilEye2->GetOutput(), CmdLineObj.OutputIm);
  typedef typename itk::BinaryImageToShapeLabelMapFilter<MaskImType> LabellerType;

  typename LabellerType::Pointer Labeller = LabellerType::New();
  Labeller->SetInput(DilEye2->GetOutput());
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

  // Saving images for paper
  writeImDbg<MaskImType>(DilEye->GetOutput(), "dilatedeye");
  writeImDbg<MaskImType>(Thresh->GetOutput(), "trianglethresh");
  writeImDbg<MaskImType>(EroEye->GetOutput(), "erodedeye");
  writeImDbg<MaskImType>(BRecon->GetOutput(), "recon");


}

int main(int argc, char * argv[])
{
  CmdLineType CmdLineObj;
  ParseCmdLine(argc, argv, CmdLineObj);

  const int dimension = 3;
  int dim1 = 0;
  itk::ImageIOBase::IOComponentType ComponentType;
  if (!readImageInfo(CmdLineObj.StructuralIm, &ComponentType, &dim1)) 
    {
    std::cerr << "Failed to open " << CmdLineObj.StructuralIm << std::endl;
    return(EXIT_FAILURE);
    }
  if (dim1 != dimension) 
    {
      std::cerr << CmdLineObj.StructuralIm << "isn't 3D" << std::endl;
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
