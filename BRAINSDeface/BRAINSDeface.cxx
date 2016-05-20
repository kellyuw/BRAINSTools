//
// Created by Leinoff, Alexander on 5/20/16.
//


#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <BRAINSDefaceCLP.h>
#include <Slicer3LandmarkIO.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkMinimumMaximumImageCalculator.h>


int main(int argc, char **argv)
{
  PARSE_ARGS;

  //Read in the imagefile. Right now assumed to be t1 weighted image
  typedef itk::Image<double,3> ImageType;
  typedef itk::ImageFileReader<ImageType> ImageReaderType;
  ImageReaderType::Pointer imageReader = ImageReaderType::New();
  imageReader->SetFileName(inputImage);

  //Read in the atlas label file
  typedef itk::Image<double,3> LabelAtlasType;
  typedef itk::ImageFileReader<LabelAtlasType> LabelAtlasReaderType;
  LabelAtlasReaderType::Pointer labelAtlasReader = LabelAtlasReaderType::New();
  labelAtlasReader->SetFileName(labelmap);

  //Read in the landmarks file
  LandmarksMapType myLandmarks = ReadSlicer3toITKLmk(landmarks);


  //Turn Label map into binary image. Use a threshold Image filter?? or brainscut?
  typedef itk::Image<unsigned char, 3> MaskAtlasType;
  typedef itk::BinaryThresholdImageFilter< LabelAtlasType, MaskAtlasType>  MaskFilterType;
  MaskFilterType::Pointer maskFilter = MaskFilterType::New();

  maskFilter->SetInput( labelAtlasReader->GetOutput() );
  maskFilter->SetOutsideValue(1);
  maskFilter->SetInsideValue(0);
  maskFilter->SetLowerThreshold(0);
  maskFilter->SetUpperThreshold(1);

  //Write to a file
  typedef itk::ImageFileWriter<MaskAtlasType> MaskAtlasWriterType;
  MaskAtlasWriterType::Pointer maskAtlasWriter = MaskAtlasWriterType::New();

  maskAtlasWriter->SetInput(maskFilter->GetOutput());
  maskAtlasWriter->SetFileName("/scratch/aleinoff/DefaceOutput.nii.gz");
  maskAtlasWriter->Update();

  return EXIT_SUCCESS;
}