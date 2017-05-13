#include <iomanip>
 
#include "itkNumericTraits.h"
#include "itkImage.h"
#include "itkRandomImageSource.h"
#include "itkImageFileWriter.h"
#include "itkPlaneSpatialObject.h"
#include "itkEllipseSpatialObject.h"
#include "itkSpatialObjectToImageStatisticsCalculator.h"
#include "itkSpatialObjectToImageFilter.h"
 
 
typedef float           PixelType;  //  Operations
const   unsigned int     Dimension = 2;
 
typedef itk::Image<PixelType, Dimension>    ImageType;
typedef itk::ImageFileWriter<ImageType> WriterType;
typedef itk::EllipseSpatialObject<Dimension> EllipseSpatialObjectType;
typedef itk::PlaneSpatialObject<Dimension> PlaneSpatialObjectType;
typedef itk::SpatialObjectToImageStatisticsCalculator<ImageType, EllipseSpatialObjectType > EllipseCalculatorType;
typedef itk::SpatialObjectToImageFilter<EllipseSpatialObjectType, ImageType> EllipseSpatialObjectToImageFilterType;
typedef itk::SpatialObjectToImageStatisticsCalculator<ImageType, PlaneSpatialObjectType > PlaneCalculatorType;
typedef itk::SpatialObjectToImageFilter<PlaneSpatialObjectType, ImageType> PlaneSpatialObjectToImageFilterType;
 
typedef ImageType::IndexType ImageIndexType;
typedef ImageType::SpacingType ImageSpacingType;
typedef ImageType::RegionType ImageRegionType;
typedef ImageType::IndexType ImageIndexType;
typedef ImageType::SizeType  ImageSizeType;
typedef ImageType::PointType PointType;
 
int main(int argc, char* argv[])
{
 
  if( argc < 3)
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << "   spacing outputFileNamePrefix length"
 	      << std::endl;
    return EXIT_FAILURE;
    }
 
  int ArgumentNumber = 0;
  const char* outputFileNamePrefix = argv[++ArgumentNumber];
  const float spacing = atof (argv[++ArgumentNumber]);
  const float length = atof (argv[++ArgumentNumber]);
 
  std::cout << "Image spacing : " << spacing << std::endl;
  std::cout << "edge of the square, or diameter : " << length << std::endl;
   
  ImageSizeType imageSize;
  imageSize[0] = 100;
  imageSize[1] = 100;
  ImageSpacingType imageSpacing;
  imageSpacing[0] = spacing;
  imageSpacing[1] = spacing;
  ImageType::IndexType start;
  start[0] =  0;
  start[1] =  0;
 
  ImageType::Pointer image = ImageType::New();
  image->SetSpacing (imageSpacing);
  ImageType::RegionType region;
  region.SetSize (imageSize);
  region.SetIndex (start);
  image->SetRegions (region);
 
  image->Allocate ();
  image->FillBuffer (itk::NumericTraits<PixelType>::One);
 
 
  // Find the center of the image
  ImageIndexType imageCenterIndex;
  for (int i = 0 ; i < Dimension; ++i)
    {
    imageCenterIndex[i] = imageSize[i]/2;
    }
  PointType imageCenterPosition;
  image->TransformIndexToPhysicalPoint (imageCenterIndex, 
 					imageCenterPosition);
 
  PointType LLCornerPosition, UUCornerPosition, zeroPosition;
  for (int i = 0; i < Dimension; ++i)
    {
    LLCornerPosition[i] = imageCenterPosition[i] - length/2 - imageSpacing[i]/2;
    UUCornerPosition[i] = imageCenterPosition[i] + length/2 + imageSpacing[i]/2;
    zeroPosition[i] = 0.0;
    }
 
  std::cout << "imageCenterPosition : " << imageCenterPosition << std::endl;
  std::cout << "1st square corner : " << LLCornerPosition << std::endl;
  std::cout << "2nd square corner : " << UUCornerPosition << std::endl;
 
 
  // Define a plane spatial object
  PlaneSpatialObjectType::Pointer plane = PlaneSpatialObjectType::New();
  plane->SetLowerPoint (LLCornerPosition);
  plane->SetUpperPoint (UUCornerPosition);
 
  PlaneCalculatorType::Pointer planeCalculator = PlaneCalculatorType::New();
  planeCalculator->SetSpatialObject (plane);
  planeCalculator->SetImage (image);
  planeCalculator->Update ();
 
  // Define an ellipse spatial object
  EllipseSpatialObjectType::Pointer ellipse = EllipseSpatialObjectType::New();
  ellipse->SetRadius (length/2);
  EllipseSpatialObjectType::VectorType offset;
  offset[0] = imageCenterPosition[0];  
  offset[1] = imageCenterPosition[1];
  ellipse->GetIndexToObjectTransform()->SetOffset (offset);
  ellipse->ComputeObjectToParentTransform ();
 
  EllipseCalculatorType::Pointer ellipseCalculator = EllipseCalculatorType::New();
  ellipseCalculator->SetSpatialObject (ellipse);
  ellipseCalculator->SetImage (image);
  ellipseCalculator->Update ();
 
  std::cout << "Is the center of the image inside  " << std::endl
 	    << "   - the PlaneSpatialObject : " 
 	    << plane->IsInside (imageCenterPosition)  << std::endl
 	    << "   - the EllipseSpatialObject : " 
 	    << ellipse->IsInside (imageCenterPosition) << std::endl;
  std::cout << "Is (0,0) inside  " << std::endl
 	    << "   - the PlaneSpatialObject : " 
 	    <<  plane->IsInside (zeroPosition) << std::endl
 	    << "   - the EllipseSpatialObject : " 
 	    << ellipse->IsInside (zeroPosition) << std::endl;
  std::cout << "Number of pixels in the region " << std::endl
 	    << "   - PlaneCalculator : " << planeCalculator->GetNumberOfPixels()
	    << std::endl
 	    << "   - EllipseCalculator : " << ellipseCalculator->GetNumberOfPixels()
 	    << std::endl;
  std::cout << "Sample mean  " << std::endl
 	    << "   - the PlaneCalculator : " 
 	    << planeCalculator->GetMean() << std::endl
 	    << "   - the EllipseCalculator : "
 	    << ellipseCalculator->GetMean() << std::endl;
 
  PlaneSpatialObjectToImageFilterType::Pointer planeSpatialObjectToImageFilter = PlaneSpatialObjectToImageFilterType::New();
  planeSpatialObjectToImageFilter->SetInput (plane);
  planeSpatialObjectToImageFilter->SetOrigin (image->GetOrigin());
  planeSpatialObjectToImageFilter->SetSpacing (image->GetSpacing());
  planeSpatialObjectToImageFilter->SetSize (image->GetLargestPossibleRegion().GetSize());
   
  WriterType::Pointer writer1 = WriterType::New();
  std::ostringstream outputFileName1;
  outputFileName1 << outputFileNamePrefix
		  << "PlaneSpatialImageFilter.vtk";
  writer1->SetFileName (outputFileName1.str().c_str());
  writer1->SetInput (planeSpatialObjectToImageFilter->GetOutput());
  writer1->Update ();
 
  EllipseSpatialObjectToImageFilterType::Pointer ellipseSpatialObjectToImageFilter = EllipseSpatialObjectToImageFilterType::New();
  ellipseSpatialObjectToImageFilter->SetInput (ellipse);
  ellipseSpatialObjectToImageFilter->SetOrigin (image->GetOrigin());
  ellipseSpatialObjectToImageFilter->SetSpacing (image->GetSpacing());
  ellipseSpatialObjectToImageFilter->SetSize (image->GetLargestPossibleRegion().GetSize());
   
  WriterType::Pointer writer2 = WriterType::New();
  std::ostringstream outputFileName2;
  outputFileName2 << outputFileNamePrefix
		  << "EllipseSpatialImageFilter.vtk";
  writer2->SetFileName (outputFileName2.str().c_str());
  writer2->SetInput (ellipseSpatialObjectToImageFilter->GetOutput());
  writer2->Update ();
 
  WriterType::Pointer writer3 = WriterType::New();
  std::ostringstream outputFileName3;
  outputFileName3 << outputFileNamePrefix
		  << "Image.vtk";
  writer3->SetFileName (outputFileName3.str().c_str());
  writer3->SetInput (image);
  writer3->Update ();
 
}
