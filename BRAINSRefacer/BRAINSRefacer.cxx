//
// Created by Leinoff, Alexander on 5/20/16.
//

#include <itkTransformToDisplacementFieldFilter.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <BRAINSRefacerCLP.h>
#include <Slicer3LandmarkIO.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkBSplineTransform.h>
#include <itkResampleImageFilter.h>
#include <itkDanielssonDistanceMapImageFilter.h>
#include <itkDivideImageFilter.h>
#include <itkVectorIndexSelectionCastImageFilter.h>
#include <itkTransformFileWriter.h>
#include <itkMultiplyImageFilter.h>
#include <itkBSplineTransformInitializer.h>
#include <itkComposeImageFilter.h>
#include <itkDisplacementFieldTransform.h>
#include <itkSubtractImageFilter.h>
#include <map>
#include <string>

#include "CreateRandomBSpline.h"
#include "CombineBSplineWithDisplacement.h"

//Convienience function to write images
template< typename TImageType >
void WriteImage(std::string filename, TImageType *image)
{
  std::cout << "Writing Image: " << filename << std::endl;
  typedef itk::ImageFileWriter<TImageType> FileWriterType;
  typename FileWriterType::Pointer fileWriter = FileWriterType::New();

  fileWriter->SetInput(image);
  fileWriter->SetFileName(filename);
  fileWriter->Update();
  std::cout << "\tdone writing Image: " << filename << std::endl;
}

//Convienience function to write images
template< typename TImageType >
void WriteSmartImage(std::string filename, typename TImageType::Pointer image)
{
  std::cout << "Writing Image: " << filename << std::endl;
  typedef itk::ImageFileWriter<TImageType> FileWriterType;
  typename FileWriterType::Pointer fileWriter = FileWriterType::New();

  fileWriter->SetInput(image);
  fileWriter->SetFileName(filename);
  fileWriter->Update();
  std::cout << "\tdone writing Image: " << filename << std::endl;
}


//Convienience function to write transforms
template< typename TTransformType >
void WriteTransform(std::string transformFileName, TTransformType transform )
{
  std::cout << "Writing Transform: " << transformFileName << std::endl;
  typedef itk::TransformFileWriter TransformWriterType;
  TransformWriterType::Pointer transformWriter = TransformWriterType::New();

  transformWriter->SetInput(transform);
  transformWriter->SetFileName(transformFileName);
  transformWriter->Update();
  std::cout << "\t done writing Transform: " << transformFileName << std::endl;
}

int main(int argc, char **argv)
{
  PARSE_ARGS;

  //Basic typedef's
  typedef  double PixelType;
  const unsigned int Dimension = 3;
  typedef itk::Image<PixelType, Dimension> ImageType;

  //Read in subject image
  typedef itk::ImageFileReader<ImageType> ImageReaderType;
  ImageReaderType::Pointer imageReader = ImageReaderType::New();
  imageReader->SetFileName(inputImage);
  ImageType::Pointer subject = imageReader->GetOutput();
  imageReader->Update();

  //Read in the atlas label file
  typedef itk::Image<PixelType, Dimension> LabelAtlasType;
  typedef itk::ImageFileReader<LabelAtlasType> LabelAtlasReaderType;
  LabelAtlasReaderType::Pointer labelAtlasReader = LabelAtlasReaderType::New();
  labelAtlasReader->SetFileName(labelmap);

  typedef itk::Image<unsigned char, Dimension> MaskAtlasType;

  //Read in the landmarks file
  LandmarksMapType myLandmarks = ReadSlicer3toITKLmk(landmarks);
  //print out landmarks
  typedef LandmarksMapType::const_iterator LandmarksIteratorType;
  LandmarksIteratorType lmIter;
/*


  for(lmIter = myLandmarks.begin(); lmIter != myLandmarks.end(); ++lmIter )
    {
    std::cout <<"Landmark: " << lmIter->first << "\t\t";
    std::cout << "Point: " << lmIter->second << std::endl;
    }
  */

  // Get the points for the right eye, the left eye and the dens_axis
  typedef itk::Point<double, 3> PointType;
#if 1
  PointType rightEye   = myLandmarks.find("RE")->second;
  PointType leftEye    = myLandmarks.find("LE")->second;
  PointType dens_axis  = myLandmarks.find("dens_axis")->second;
#endif

#if 0
  PointType rightEye;
  rightEye[0] = 2;
  rightEye[1] = 1;
  rightEye[2] = -1;

  PointType leftEye;
  leftEye[0] = 0;
  leftEye[1] = -2;
  leftEye[2] = 0;

  PointType dens_axis;
  dens_axis[0] = 1;
  dens_axis[1] = -1;
  dens_axis[2] = 2;
#endif

  std::cout << "rightEye:\t" << rightEye << std::endl;
  std::cout << "leftEye:\t" << leftEye << std::endl;
  std::cout << "dens_axis:\t" << dens_axis << std::endl;


  // find the a,b,c for the plane equation ax + by + c = 0
  // first get two vectors in the plane u and v
  typedef itk::Vector<double,3> VectorType;
  VectorType u = rightEye - leftEye;
  VectorType v = dens_axis - leftEye;

  std::cout << "vector u: \t" << u << std::endl;
  std::cout << "vector v: \t" << v << std::endl;

  VectorType cross = itk::CrossProduct(u, v);
  std::cout << "vector cross: \t" << cross <<std::endl;

  //for ax + by + cz = d plug in cross for abc and one of the points (dens_axis) for xyz
  VectorType leftEyeVector;
  leftEyeVector[0] = leftEye[0];
  leftEyeVector[1] = leftEye[1];
  leftEyeVector[2] = leftEye[2];

  VectorType::ComponentType d = cross * leftEyeVector;
  std::cout << "d:\t" << d << std::endl;

  std::cout << "equation of plane is:" << std::endl;
  std::cout << cross[0] << "x + " << cross[1] << "y + " << cross[2] << "z = " << d << std::endl;

  // make a mask based on the plane:

  // first create mask image
  MaskAtlasType::Pointer maskImageLM = MaskAtlasType::New();
  maskImageLM->SetOrigin(subject->GetOrigin());
  maskImageLM->SetSpacing(subject->GetSpacing());
  maskImageLM->SetDirection(subject->GetDirection());
  maskImageLM->SetRegions(subject->GetLargestPossibleRegion());
  maskImageLM->Allocate();

  typedef itk::ImageRegionConstIterator<ImageType> SubjectIteratorType;
  SubjectIteratorType subjectIt(subject, subject->GetLargestPossibleRegion());

  typedef itk::ImageRegionIterator<MaskAtlasType> MaskIteratorType;
  MaskIteratorType maskIt(maskImageLM, maskImageLM->GetLargestPossibleRegion());

  subjectIt.GoToBegin();
  maskIt.GoToBegin();

  // go through image to see if pixel is in mask or not??

  while(!subjectIt.IsAtEnd())
    {
    PointType currentPoint;
    subject->TransformIndexToPhysicalPoint(subjectIt.GetIndex(), currentPoint);
  //  std::cout << currentPoint <<std::endl;

    VectorType currentVector;
    currentVector[0] = currentPoint[0];
    currentVector[1] = currentPoint[1];
    currentVector[2] = currentPoint[2];

    if( d > cross * currentVector )
      {
      maskIt.Set(0);
      }
    else
      {
      maskIt.Set(1);
      }

    ++maskIt;
    ++subjectIt;
    }



  WriteSmartImage<MaskAtlasType>("/scratch/aleinoff/temp/maskTest.nii.gz", maskImageLM);

  return 0;


  //Turn Label map into binary image. Use a threshold Image filter?? or brainscut?
  //Write a new filter for this??
  //typedef itk::Image<unsigned char, Dimension> MaskAtlasType;
  typedef itk::BinaryThresholdImageFilter< LabelAtlasType, MaskAtlasType>  MaskFilterType;
  MaskFilterType::Pointer maskFilter = MaskFilterType::New();

  maskFilter->SetInput( labelAtlasReader->GetOutput() );
  maskFilter->SetOutsideValue(1);
  maskFilter->SetInsideValue(0);
  maskFilter->SetLowerThreshold(0);
  maskFilter->SetUpperThreshold(0);

  //Write to a file
  WriteImage(outputMask, maskFilter->GetOutput());
  //Get a distance map to the Brain region:
  typedef itk::DanielssonDistanceMapImageFilter<MaskAtlasType, ImageType, ImageType> DistanceMapFilter;
  DistanceMapFilter::Pointer distanceMapFilter = DistanceMapFilter::New();
  distanceMapFilter->SetInput(maskFilter->GetOutput());
  distanceMapFilter->InputIsBinaryOn();
  distanceMapFilter->SetSquaredDistance(false);

  //Write the distance map to a file so we can see what it did:
  WriteImage(distanceMapFileName, distanceMapFilter->GetOutput());

  //Try to scale distance map
  typedef itk::MultiplyImageFilter<ImageType, ImageType, ImageType> ScalingFilterType;
  ScalingFilterType::Pointer distanceMapScaler = ScalingFilterType::New();
  distanceMapScaler->SetInput(distanceMapFilter->GetOutput());
  distanceMapScaler->SetConstant(scaleDistanceMap);

  //Perform some kind of BSpline on Image
  const int BSplineOrder = 3;

  typedef CreateRandomBSpline<ImageType, PixelType, Dimension, BSplineOrder> BSplineCreator; //, BSTransformType> Test;
  BSplineCreator::Pointer bSplineCreator = BSplineCreator::New();
  bSplineCreator->SetInput(subject);
  bSplineCreator->SetBSplineControlPoints(bsplineControlPoints);
  bSplineCreator->SetRandMax(maxRandom);
  bSplineCreator->SetRandMin(minRandom);
  bSplineCreator->SetRandScale(scaleRandom);
  bSplineCreator->Update();

  typedef itk::BSplineTransform<PixelType, Dimension, BSplineOrder> BSTransformType;
  BSTransformType::Pointer bSpline = bSplineCreator->GetBSplineOutput();

  WriteTransform(bSplineFileName, bSpline);
//return 0;
  typedef itk::Vector<PixelType, Dimension > VectorPixelType;
  typedef itk::Image< VectorPixelType, Dimension> DisplacementFieldImageType;

  typedef CombineBSplineWithDisplacement<ImageType, DisplacementFieldImageType, PixelType, 3,3> CombinerType;

  CombinerType::Pointer combiner = CombinerType::New();
  combiner->SetBSplineInput(bSpline);
  combiner->SetInput(subject);
  combiner->SetDistanceMap(distanceMapScaler->GetOutput());
  combiner->Update();

  //write the new displacement image
  DisplacementFieldImageType* composedDisplacementField_rawPtr = combiner->GetComposedImage();
  WriteImage(smoothDisplacementName, composedDisplacementField_rawPtr);


#if 0
  //write composed displacement field into a displacement transform
  std::cout<<"printing composed image info" <<std::endl;
  combiner->GetComposedImage()->Print(std::cerr,5);
  std::cout<<"done printing composed image" <<std::endl;

  std::cout<<"printing distance map info" <<std::endl;
  combiner->GetDistanceMap()->Print(std::cerr, 0);
  std::cout<<"done printing distance map" <<std::endl;
#endif

  typedef itk::DisplacementFieldTransform<PixelType, Dimension> FinalTransformType;
  FinalTransformType::Pointer finalTransform = FinalTransformType::New();
  finalTransform->SetDisplacementField(composedDisplacementField_rawPtr);
  finalTransform->Print(std::cerr,5);

  WriteTransform(finalTransformFileName, finalTransform);

  // Apply transform to image with resampler:
  typedef itk::ResampleImageFilter<ImageType, ImageType> ResampleFilterType;
  ResampleFilterType::Pointer resampler = ResampleFilterType::New();

  typedef itk::LinearInterpolateImageFunction<ImageType, PixelType > InterpolatorType;
  InterpolatorType::Pointer interpolater = InterpolatorType::New();

  ImageType::RegionType subjectRegion = subject->GetBufferedRegion();

  resampler->SetInterpolator(interpolater);
  resampler->SetOutputSpacing(subject->GetSpacing());
  resampler->SetOutputOrigin(subject->GetOrigin());
  resampler->SetOutputDirection(subject->GetDirection());
  resampler->SetSize(subjectRegion.GetSize());
  resampler->SetOutputStartIndex(subjectRegion.GetIndex());

  resampler->SetInput(imageReader->GetOutput());
  resampler->SetTransform(finalTransform);

  WriteImage(deformedImageName, resampler->GetOutput());

  //Get the difference image
  typedef itk::SubtractImageFilter<ImageType, ImageType> SubtractFilter;
  SubtractFilter::Pointer subtractFilter = SubtractFilter::New();
  subtractFilter->SetInput1(subject);
  subtractFilter->SetInput2(resampler->GetOutput());

  //write the difference Image
  WriteImage( diffImageName, subtractFilter->GetOutput());

  std::cout << "done" << std::endl;

  return EXIT_SUCCESS;
}
