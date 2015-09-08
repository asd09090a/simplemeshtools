/*=========================================================================

  Program:   ITK-SNAP
  Module:    $RCSfile: VTKMeshPipeline.h,v $
  Language:  C++
  Date:      $Date: 2009/01/23 20:09:38 $
  Version:   $Revision: 1.4 $
  Copyright (c) 2007 Paul A. Yushkevich
  
  This file is part of ITK-SNAP 

  ITK-SNAP is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
 
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

  -----

  Copyright (c) 2003 Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notices for more information. 

  Modified for vtk 6 by Richard Beare


=========================================================================*/
#ifndef __VTKMeshPipeline_h_
#define __VTKMeshPipeline_h_

#include <itkVTKImageExport.h>
#include <vtkSmartPointer.h>

// VTK includes
#include <vtkCellArray.h>
#include <vtkImageData.h>
#include <vtkImageImport.h>
#include <vtkImageGaussianSmooth.h>
#include <vtkImageThreshold.h>
#include <vtkImageToStructuredPoints.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkStripper.h>
#include <vtkCallbackCommand.h>
#include <vtkMarchingCubes.h>
#include <vtkDecimatePro.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkTransform.h>

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
class VTKMeshPipeline
{
public:
  /** Input image type */
  typedef typename ImageType::Pointer ImagePointer;
  
  /** Set the input segmentation image */
  void SetImage(ImagePointer image);


  /** Compute a mesh for a particular color label */
  void ComputeMesh(vtkPolyData *outData);

  void SetMeshOptions();

  /** Constructor, which builds the pipeline */
  VTKMeshPipeline();

  /** Deallocate the pipeline filters */
  ~VTKMeshPipeline();

  itkSetMacro(MarchThresh, float);

  itkSetMacro(UseGaussianSmoothing, bool);
  itkSetMacro(UseDecimation, bool);
  itkSetMacro(UseMeshSmoothing, bool);

  itkSetMacro(GaussianStandardDeviation, float);
  itkSetMacro(GaussianError, float);
  itkSetMacro(DecimateTargetReduction, float);
  itkSetMacro(DecimateInitialError, float);
  itkSetMacro(DecimateAspectRatio, float);
  itkSetMacro(DecimateFeatureAngle, float);
  itkSetMacro(DecimateErrorIncrement, float);
  itkSetMacro(DecimateMaximumIterations, unsigned int);
  itkSetMacro(DecimatePreserveTopology, bool);
  
  itkSetMacro(MeshSmoothingRelaxationFactor, float);
  itkSetMacro(MeshSmoothingIterations, unsigned int );
  itkSetMacro(MeshSmoothingConvergence, float);
  itkSetMacro(MeshSmoothingFeatureAngle, float);
  itkSetMacro(MeshSmoothingFeatureEdgeSmoothing, bool);
  itkSetMacro(MeshSmoothingBoundarySmoothing, bool);

private:
  bool m_UseGaussianSmoothing;
  bool m_UseDecimation;
  bool m_UseMeshSmoothing;
  
  // threshold for marching cubes
  float m_MarchThresh;
  // Begin gsmooth params
  float m_GaussianStandardDeviation;
  float m_GaussianError;

  // Begin decimate parameters
  float m_DecimateTargetReduction;
  float m_DecimateInitialError;
  float m_DecimateAspectRatio;
  float m_DecimateFeatureAngle;
  float m_DecimateErrorIncrement;
  unsigned int m_DecimateMaximumIterations;
  bool m_DecimatePreserveTopology;
  
  // Begin msmooth params
  float m_MeshSmoothingRelaxationFactor;
  unsigned int m_MeshSmoothingIterations;
  float m_MeshSmoothingConvergence;
  float m_MeshSmoothingFeatureAngle;
  bool m_MeshSmoothingFeatureEdgeSmoothing;
  bool m_MeshSmoothingBoundarySmoothing;

  // VTK-ITK Connection typedefs
  typedef itk::VTKImageExport<ImageType> VTKExportType;
  typedef itk::SmartPointer<VTKExportType> VTKExportPointer;
  
  // to fool the itk macro generator
  void Modified(){};

  // The input image
  ImagePointer m_InputImage;

  // The VTK exporter for the data
  VTKExportPointer m_VTKExporter;

  // The VTK importer for the data
  vtkSmartPointer<vtkImageImport> m_VTKImporter;

  // VTK Gaussian (because we care about the speed and not so much about
  // precision)
  vtkSmartPointer<vtkImageGaussianSmooth> m_VTKGaussianFilter;

  // The polygon smoothing filter
  vtkSmartPointer<vtkSmoothPolyDataFilter> m_PolygonSmoothingFilter;

  // Triangle stripper
  vtkSmartPointer<vtkStripper> m_StripperFilter;  
  
  // Marching cubes filter
  vtkSmartPointer<vtkMarchingCubes>     m_MarchingCubesFilter;

  // Transform filter used to map to RAS space
  vtkSmartPointer<vtkTransformPolyDataFilter> m_TransformFilter;

  // The transform used
  vtkSmartPointer<vtkTransform> m_Transform;
  
  // The triangle decimation driver
  vtkSmartPointer<vtkDecimatePro>   m_DecimateFilter;
  /** 
   * This static function constructs a NIFTI matrix from the ITK direction
   * cosines matrix and Spacing and Origin vectors
   */

  static vnl_matrix_fixed<double,4,4> ConstructNiftiSform(
    vnl_matrix<double> m_dir, 
    vnl_vector<double> v_origin,
    vnl_vector<double> v_spacing)
  {
    // Set the NIFTI/RAS transform
    vnl_matrix<double> m_ras_matrix;
    vnl_diag_matrix<double> m_scale, m_lps_to_ras;
    vnl_vector<double> v_ras_offset;

    // Compute the matrix
    m_scale.set(v_spacing);
    m_lps_to_ras.set(vnl_vector<double>(3, 1.0));
    m_lps_to_ras[0] = -1;
    m_lps_to_ras[1] = -1;
    m_ras_matrix = m_lps_to_ras * m_dir * m_scale;

    // Compute the vector
    v_ras_offset = m_lps_to_ras * v_origin;

    // Create the larger matrix
    vnl_vector<double> vcol(4, 1.0);
    vcol.update(v_ras_offset);
    
    vnl_matrix_fixed<double,4,4> m_sform;
    m_sform.set_identity();
    m_sform.update(m_ras_matrix);
    m_sform.set_column(3, vcol);
    return m_sform;
  }

  static vnl_matrix_fixed<double,4,4> ConstructVTKtoNiftiTransform(
    vnl_matrix<double> m_dir, 
    vnl_vector<double> v_origin,
    vnl_vector<double> v_spacing)
    {
    vnl_matrix_fixed<double,4,4> vox2nii = ConstructNiftiSform(m_dir, v_origin, v_spacing);
    vnl_matrix_fixed<double,4,4> vtk2vox; 
    vtk2vox.set_identity();
    for(size_t i = 0; i < 3; i++)
      {
      vtk2vox(i,i) = 1.0 / v_spacing[i];
      vtk2vox(i,3) = - v_origin[i] / v_spacing[i];
      }
    return vox2nii * vtk2vox;
    }


};
#ifndef ITK_MANUAL_INSTANTIATION
#include "VTKMeshPipeline.cxx"
#endif


#endif // __VTKMeshPipeline_h_
