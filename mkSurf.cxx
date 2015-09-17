// Try to locate the highest gradient near the surface of the scalp.
// Uses a marker derived from a presurgical scan
#include <iostream>
#include <cstdio>
#include <vector>

#define DIRGRAD

#include "tclap/CmdLine.h"

#include "ioutils.h"
#include <itkBinaryShapeOpeningImageFilter.h>
#include <itkMaskImageFilter.h>
#include "itkBinaryDilateParaImageFilter.h"
#include "itkBinaryErodeParaImageFilter.h"
#include <itkMaximumImageFilter.h>
#include "itkBinaryFillholeImageFilter.h"
#include <itkSubtractImageFilter.h>
#include <itkBinaryShapeKeepNObjectsImageFilter.h>
#include <itkMorphologicalWatershedFromMarkersImageFilter.h>

#ifdef DIRGRAD

#include "itkDirectionalGradientImageFilter.h"
#include <itkSmoothingRecursiveGaussianImageFilter.h>

#else

#include <itkGradientMagnitudeRecursiveGaussianImageFilter.h>

#endif
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
  std::string InputIm, OutputIm, MaskIm;
  float  erodesize, dilatesize, smoothsize;
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
    CmdLine cmd("mkSurf", ' ', "0.9");

    ValueArg<std::string> inArg("i","input","T1 input image",true,"result","string");
    cmd.add( inArg );

    ValueArg<std::string> outArg("o","output","output image",true,"result","string");
    cmd.add( outArg );

    ValueArg<std::string> maskArg("m","mask","mask image - used to generate markers",true,"result","string");
    cmd.add( maskArg );

    ValueArg<float> smoothArg("", "smoothing", "size of small gradient smoothing (mm)", false, 2, "float");
    cmd.add(smoothArg);

    ValueArg<float> eroArg("", "erode", "size of erosion to create marker(mm)", 
                              false, 3, "float");
    cmd.add(eroArg);

    ValueArg<float> dilArg("", "dilate", "size of dilation background marker (mm)", 
                              false, 3, "float");
    cmd.add(dilArg);

    SwitchArg debugArg("d", "debug", "save debug images", debug);
    cmd.add( debugArg );

    // Parse the args.
    cmd.parse( argc, argv );

    CmdLineObj.InputIm = inArg.getValue();
    CmdLineObj.OutputIm = outArg.getValue();
    CmdLineObj.MaskIm = maskArg.getValue();
    debug = debugArg.getValue();
    CmdLineObj.erodesize = eroArg.getValue();
    CmdLineObj.smoothsize = smoothArg.getValue();
    CmdLineObj.dilatesize = dilArg.getValue();
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
  MIPtr mask = readIm<MaskImType>(CmdLineObj.MaskIm);

  // fill holes
  itk::Instance<itk::BinaryFillholeImageFilter<MaskImType > > Fill;
  Fill->SetInput(mask);
  // set up the marker
  itk::Instance<itk::BinaryErodeParaImageFilter<MaskImType> > Erode;
  itk::Instance<itk::BinaryDilateParaImageFilter<MaskImType> > Dilate;

  Erode->SetInput(Fill->GetOutput());
  Erode->SetRadius(CmdLineObj.erodesize);
  Erode->SetUseImageSpacing(true);

  Dilate->SetInput(Fill->GetOutput());
  Dilate->SetRadius(CmdLineObj.dilatesize);
  Dilate->SetUseImageSpacing(true);
  
  itk::Instance<itk::BinaryThresholdImageFilter<MaskImType,MaskImType> > Invert;
  Invert->SetInput(Dilate->GetOutput());
  Invert->SetUpperThreshold(0);
  Invert->SetLowerThreshold(0);
  Invert->SetInsideValue(2);
  Invert->SetOutsideValue(0);

  itk::Instance<itk::MaximumImageFilter<MaskImType, MaskImType, MaskImType> >Comb;
  Comb->SetInput(Erode->GetOutput());
  Comb->SetInput2(Invert->GetOutput());

  // Directional gradient
#ifdef DIRGRAD
  itk::Instance< itk::DirectionalGradientImageFilter<ImageType, MaskImType, ImageType> > GradD;
  GradD->SetInput(T1);
  GradD->SetMaskImage(mask);
  GradD->SetOutsideValue(0);
  // No need to threshold (marker will cover most of the negative
  // bits), but do smooth a little - name the smoother Grad to keep
  // following code simple.
  itk::Instance< itk::SmoothingRecursiveGaussianImageFilter <ImageType, ImageType> > Grad;
  Grad->SetInput(GradD->GetOutput());
  Grad->SetSigma(CmdLineObj.smoothsize);
#else  
  itk::Instance< itk::GradientMagnitudeRecursiveGaussianImageFilter<ImageType, ImageType> > Grad;
  Grad->SetInput(T1);
  Grad->SetSigma(CmdLineObj.smoothsize);
#endif
  itk::Instance<itk::MorphologicalWatershedFromMarkersImageFilter<ImageType, MaskImType> > WS;
  WS->SetInput(Grad->GetOutput());
  WS->SetMarkerImage(Comb->GetOutput());
  WS->SetMarkWatershedLine(false);

  itk::Instance<itk::BinaryThresholdImageFilter<MaskImType,MaskImType> > Select;
  Select->SetInput(WS->GetOutput());
  Select->SetUpperThreshold(1);
  Select->SetLowerThreshold(1);
  Select->SetInsideValue(1);
  Select->SetOutsideValue(0);

  writeIm<MaskImType>(Select->GetOutput(), CmdLineObj.OutputIm);
  writeIm<ImageType>(Grad->GetOutput(), "/tmp/grad.nii.gz");
  writeIm<MaskImType>(Comb->GetOutput(), "/tmp/marker.nii.gz");

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
