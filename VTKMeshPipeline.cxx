/*=========================================================================

  Program:   ITK-SNAP
  Module:    $RCSfile: VTKMeshPipeline.cxx,v $
  Language:  C++
  Date:      $Date: 2010/06/28 18:45:08 $
  Version:   $Revision: 1.7 $
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

=========================================================================*/
#include "VTKMeshPipeline.h"
#include <map>
#include "vtkPolyDataWriter.h"
#include <vtkSmartPointer.h>

using namespace std;

template<class ImageType>
VTKMeshPipeline<ImageType>
::VTKMeshPipeline()
{
  // Initialize all the filters involved in the transaction, but do not
  // pipe the inputs and outputs between these filters. The piping is quite
  // complicated and depends on the set of options that the user wishes to 
  // apply
  
  // Initialize the VTK Exporter
  m_VTKExporter = VTKExportType::New();
  //m_VTKExporter->ReleaseDataFlagOn();
  
  // Initialize the VTK Importer
  m_VTKImporter = vtkSmartPointer<vtkImageImport>::New();
  //m_VTKImporter->ReleaseDataFlagOn();

  this->SetMarchThresh(0.5);

  // Pipe the importer into the exporter (that's a lot of code)
  m_VTKImporter->SetUpdateInformationCallback(
    m_VTKExporter->GetUpdateInformationCallback());
  m_VTKImporter->SetPipelineModifiedCallback(
    m_VTKExporter->GetPipelineModifiedCallback());
  m_VTKImporter->SetWholeExtentCallback(
    m_VTKExporter->GetWholeExtentCallback());
  m_VTKImporter->SetSpacingCallback(
    m_VTKExporter->GetSpacingCallback());
  m_VTKImporter->SetOriginCallback(
    m_VTKExporter->GetOriginCallback());
  m_VTKImporter->SetScalarTypeCallback(
    m_VTKExporter->GetScalarTypeCallback());
  m_VTKImporter->SetNumberOfComponentsCallback(
    m_VTKExporter->GetNumberOfComponentsCallback());
  m_VTKImporter->SetPropagateUpdateExtentCallback(
    m_VTKExporter->GetPropagateUpdateExtentCallback());
  m_VTKImporter->SetUpdateDataCallback(
    m_VTKExporter->GetUpdateDataCallback());
  m_VTKImporter->SetDataExtentCallback(
    m_VTKExporter->GetDataExtentCallback());
  m_VTKImporter->SetBufferPointerCallback(
    m_VTKExporter->GetBufferPointerCallback());  
  m_VTKImporter->SetCallbackUserData(
    m_VTKExporter->GetCallbackUserData());

  // Initialize the Gaussian filter
  m_VTKGaussianFilter = vtkSmartPointer<vtkImageGaussianSmooth>::New();
  //m_VTKGaussianFilter->ReleaseDataFlagOn();
  
  // Create and configure a filter for polygon smoothing
  m_PolygonSmoothingFilter = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
  //m_PolygonSmoothingFilter->ReleaseDataFlagOn();

  // Create and configure a filter for triangle strip generation
  m_StripperFilter = vtkSmartPointer<vtkStripper>::New();
  //m_StripperFilter->ReleaseDataFlagOn();

  // Create and configure the marching cubes filter
  m_MarchingCubesFilter = vtkSmartPointer<vtkMarchingCubes>::New();
  //m_MarchingCubesFilter->ReleaseDataFlagOn();
  m_MarchingCubesFilter->ComputeScalarsOff();
  m_MarchingCubesFilter->ComputeGradientsOff();
  m_MarchingCubesFilter->SetNumberOfContours(1);
  m_MarchingCubesFilter->SetValue(0,m_MarchThresh);

  // Create the transform filter
  m_TransformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
  //m_TransformFilter->ReleaseDataFlagOn();

  // Create the transform
  m_Transform = vtkSmartPointer<vtkTransform>::New();

  // Create and configure a filter for triangle decimation
  m_DecimateFilter = vtkSmartPointer<vtkDecimatePro>::New();
  //m_DecimateFilter->ReleaseDataFlagOn();  

  // default settings
  this->SetUseGaussianSmoothing(false);
  this->SetUseDecimation(false);
  this->SetUseMeshSmoothing(false);

  this->SetGaussianStandardDeviation(0.8);
  this->SetGaussianError(0.03);

  this->SetDecimateTargetReduction(0.95);
  this->SetDecimateInitialError(0.002);
  this->SetDecimateAspectRatio(20);
  this->SetDecimateFeatureAngle(45);
  this->SetDecimateErrorIncrement(0.002);
  this->SetDecimatePreserveTopology(true);
  this->SetDecimateMaximumIterations(1);

  this->SetMeshSmoothingRelaxationFactor(0.01);
  this->SetMeshSmoothingConvergence(0);
  this->SetMeshSmoothingIterations(1);
  this->SetMeshSmoothingFeatureAngle(45);
  this->SetMeshSmoothingFeatureEdgeSmoothing(false);
  this->SetMeshSmoothingBoundarySmoothing(false);
}
template<class ImageType>
VTKMeshPipeline<ImageType>
::~VTKMeshPipeline()
{
  // Destroy the filters
  // m_VTKImporter->Delete();
  // m_VTKGaussianFilter->Delete();
  // m_PolygonSmoothingFilter->Delete();
  // m_StripperFilter->Delete();

  // m_MarchingCubesFilter->Delete();
  // m_TransformFilter->Delete();
  // m_Transform->Delete();
  // m_DecimateFilter->Delete();
}

template<class ImageType>
void
VTKMeshPipeline<ImageType>
::SetMeshOptions()
{
  // assumes the options are already set


  // Define the current pipeline end-point
  // vtkSmartPointer<vtkImageData> pipeImageTail = m_VTKImporter->GetOutputPort();
  // vtkSmartPointer<vtkPolyData> pipePolyTail;

  vtkAlgorithmOutput *pipeImageTail = m_VTKImporter->GetOutputPort();
  vtkAlgorithmOutput *pipePolyTail;

  // Route the pipeline according to the settings
  // 1. Check if Gaussian smoothing will be used

  if(m_UseGaussianSmoothing) 
    {    
    // The Gaussian filter is enabled
    m_VTKGaussianFilter->SetInputConnection(pipeImageTail);
    pipeImageTail = m_VTKGaussianFilter->GetOutputPort();

    // Apply parameters to the Gaussian filter
    float sigma = m_GaussianStandardDeviation;

    // Sigma is in millimeters
    const double *spacing = m_InputImage->GetSpacing().GetDataPointer();
    m_VTKGaussianFilter->SetStandardDeviation(
      sigma / spacing[0], sigma / spacing[1], sigma / spacing[2]);
    m_VTKGaussianFilter->SetRadiusFactors(
      3 * sigma / spacing[0], 3 * sigma / spacing[1], 3 * sigma / spacing[2]);
    }

  // 2. Set input to the appropriate contour filter

  // Marching cubes gets the tail
  m_MarchingCubesFilter->SetValue(0,m_MarchThresh);
  m_MarchingCubesFilter->SetInputConnection(pipeImageTail);
  pipePolyTail = m_MarchingCubesFilter->GetOutputPort();

  // 2.5 Pipe marching cubes output to the transform
  m_TransformFilter->SetInputConnection(m_MarchingCubesFilter->GetOutputPort());
  pipePolyTail = m_TransformFilter->GetOutputPort();

  // 3. Check if decimation is required
  if(m_UseDecimation)
    {
    // Decimate filter gets the pipe tail
    m_DecimateFilter->SetInputConnection(pipePolyTail);
    pipePolyTail = m_DecimateFilter->GetOutputPort();

    // Apply parameters to the decimation filter
    m_DecimateFilter->SetTargetReduction(
      m_DecimateTargetReduction);

    m_DecimateFilter->SetMaximumError(
      m_DecimateInitialError);

    m_DecimateFilter->SetFeatureAngle(
      m_DecimateFeatureAngle);

    m_DecimateFilter->SetPreserveTopology(
      m_DecimatePreserveTopology);

    } // If decimate enabled

  // 4. Compute the normals (non-patented only)

  // 5. Include/exclude mesh smoothing filter
  if(m_UseMeshSmoothing)
    {
    // Pipe smoothed output into the pipeline
    m_PolygonSmoothingFilter->SetInputConnection(pipePolyTail);
    pipePolyTail = m_PolygonSmoothingFilter->GetOutputPort();

    // Apply parameters to the mesh smoothing filter
    m_PolygonSmoothingFilter->SetNumberOfIterations(
      m_MeshSmoothingIterations);

    m_PolygonSmoothingFilter->SetRelaxationFactor(
      m_MeshSmoothingRelaxationFactor); 

    m_PolygonSmoothingFilter->SetFeatureAngle(
      m_MeshSmoothingFeatureAngle);

    m_PolygonSmoothingFilter->SetFeatureEdgeSmoothing(
      m_MeshSmoothingFeatureEdgeSmoothing);

    m_PolygonSmoothingFilter->SetBoundarySmoothing(
      m_MeshSmoothingBoundarySmoothing);

    m_PolygonSmoothingFilter->SetConvergence(
      m_MeshSmoothingConvergence);
    }

  // 6. Pipe in the final output into the stripper
  m_StripperFilter->SetInputConnection(pipePolyTail);
}

#include <ctime>

template<class ImageType>
void
VTKMeshPipeline<ImageType>
::ComputeMesh(vtkPolyData *outMesh)
{
  // Graft the polydata to the last filter in the pipeline
  //m_StripperFilter->SetOutput(outMesh);
  //m_MarchingCubesFilter->SetOutput(outMesh);
  
  // Connect importer and exporter
  m_VTKImporter->SetCallbackUserData(
    m_VTKExporter->GetCallbackUserData());

  // Update the ITK portion of the pipeline
  m_VTKExporter->SetInput(m_InputImage);
  m_VTKImporter->Modified();

  // Update the importer
  m_VTKImporter->Update();

  m_VTKImporter->Modified();
  // Update the pipeline
  // m_StripperFilter->Update();

  // Set the source of outMesh to null
  // outMesh->UpdateInformation();
  // outMesh->Update();
  // m_StripperFilter->UnRegister(m_StripperFilter->GetOutput());
  //m_StripperFilter->Update();
  //
  // m_MarchingCubesFilter->Update();
  // m_TransformFilter->Update();
  m_StripperFilter->Update();

  outMesh->ShallowCopy(m_StripperFilter->GetOutput());
  // In the case that the jacobian of the transform is negative,
  // flip the normals around
  if(m_Transform->GetMatrix()->Determinant() < 0)
    {
    vtkPointData *pd = outMesh->GetPointData();
    vtkDataArray *nrm = pd->GetNormals();
    for(size_t i = 0; i < (size_t)nrm->GetNumberOfTuples(); i++)
      for(size_t j = 0; j < (size_t)nrm->GetNumberOfComponents(); j++)
        nrm->SetComponent(i,j,-nrm->GetComponent(i,j));
    nrm->Modified();
    }
}

template<class ImageType>
void
VTKMeshPipeline<ImageType>
::SetImage(ImagePointer image)
{
  // Store the image 
  m_InputImage = image;

  // Compute the transform from VTK coordinates to NIFTI/RAS coordinates
  vnl_matrix_fixed<double, 4, 4> vtk2nii = 
    ConstructVTKtoNiftiTransform(
      image->GetDirection().GetVnlMatrix(),
      image->GetOrigin().GetVnlVector(),
      image->GetSpacing().GetVnlVector());


  // Update the VTK transform to match
  m_Transform->SetMatrix(vtk2nii.data_block());
  //std::cout << vtk2nii;
  // Pass the transform to the transform filter
  m_TransformFilter->SetTransform(m_Transform);
  m_VTKExporter->SetInput(m_InputImage);
}


