//
// Created by Leinoff, Alexander on 5/20/16.
//


#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <BRAINSDefaceCLP.h>
#include <Slicer3LandmarkIO.h>
#include <itkBinaryThresholdImageFilter.h>


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


  return EXIT_SUCCESS;
}