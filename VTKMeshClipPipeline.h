#ifndef __VTKMeshClipPipeline_h_
#define __VTKMeshClipPipeline_h_

#include <itkVTKImageExport.h>
#include <vtkSmartPointer.h>

// VTK includes
#include <vtkCellArray.h>
#include <vtkImageData.h>
#include <vtkImageImport.h>
#include <vtkImageToStructuredPoints.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>

#ifndef vtkFloatingPointType
# define vtkFloatingPointType vtkFloatingPointType
typedef float vtkFloatingPointType;
#endif

class VTKProgressAccumulator;

/**
 * \class VTKMeshPipeline
 * \brief A small pipeline used to convert an ITK image with a level set into
 * a VTK contour, with optional blurring
 */
template <class ImageType>
class VTKMeshClipPipeline
{
public:
  /** Input image type */
  typedef typename ImageType::Pointer ImagePointer;
  
  /** Set the input segmentation image */
  void SetImage(ImagePointer image);

  /** Compute a mesh for a particular color label */
  void ComputeMesh(vtkPolyData *inData, vtkPolyData *outData);


  /** Constructor, which builds the pipeline */
  VTKMeshClipPipeline();

  /** Deallocate the pipeline filters */
  ~VTKMeshClipPipeline();

private:

  // VTK-ITK Connection typedefs
  typedef itk::VTKImageExport<ImageType> VTKExportType;
  typedef itk::SmartPointer<VTKExportType> VTKExportPointer;
  
  // The input image
  ImagePointer m_InputImage;

  // The VTK exporter for the data
  VTKExportPointer m_VTKExporter;

  // The VTK importer for the data
  vtkSmartPointer<vtkImageImport> m_VTKImporter;


};
#ifndef ITK_MANUAL_INSTANTIATION
#include "VTKMeshClipPipeline.cxx"
#endif


#endif // __VTKMeshPipeline_h_
