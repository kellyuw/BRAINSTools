//
// Created by Leinoff, Alexander on 5/20/16.
//


#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <BRAINSDefaceCLP.h>
#include <Slicer3LandmarkIO.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkBSplineTransform.h>
#include <stdlib.h>
#include <time.h>
#include <itkResampleImageFilter.h>

int main(int argc, char **argv)
{
  PARSE_ARGS;

  //Read in the imagefile. Right now assumed to be t1 weighted image

  typedef double PixelType;
  //TODO: how to activate cpp11 in BRAINSTools so I can use constexpr???
  const unsigned int Dimension = 3;

  typedef itk::Image<PixelType, Dimension> ImageType;
  typedef itk::ImageFileReader<ImageType> ImageReaderType;
  ImageReaderType::Pointer imageReader = ImageReaderType::New();
  imageReader->SetFileName(inputImage);

  //

  //Read in the atlas label file
  typedef itk::Image<PixelType, Dimension> LabelAtlasType;
  typedef itk::ImageFileReader<LabelAtlasType> LabelAtlasReaderType;
  LabelAtlasReaderType::Pointer labelAtlasReader = LabelAtlasReaderType::New();
  labelAtlasReader->SetFileName(labelmap);

  //Read in the landmarks file
  LandmarksMapType myLandmarks = ReadSlicer3toITKLmk(landmarks);


  //Turn Label map into binary image. Use a threshold Image filter?? or brainscut?
  typedef itk::Image<unsigned char, Dimension> MaskAtlasType;
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
  maskAtlasWriter->SetFileName("/scratch/aleinoff/defaceOutput/binaryBrainMask.nii.gz");
  maskAtlasWriter->Update();


  //Perform some kind of BSpline on Image
  //Extract subject image from reader
  ImageType::Pointer subject = imageReader->GetOutput();
  imageReader->Update();

  //Setup BSpline
  const unsigned int BSplineOrder = 3;
  typedef itk::BSplineTransform<PixelType, Dimension, BSplineOrder> BSplineTransform;
  BSplineTransform::Pointer bSpline = BSplineTransform::New();
  //Set BSpline basic parameters
  const unsigned int BSplineControlPoints = 8;

  typedef ImageType::RegionType ImageRegionType;
//  ImageRegionType subjectRegion = subject->GetLargestPossibleRegion();
  ImageRegionType subjectRegion = subject->GetBufferedRegion();

  bSpline->SetTransformDomainOrigin(subject->GetOrigin());           //Origin
  bSpline->SetTransformDomainDirection(subject->GetDirection());     //Direction
  bSpline->SetTransformDomainPhysicalDimensions((                    //PhysicalDimensions
    subject->GetSpacing()[0]*(subjectRegion.GetSize()[0]-1),         //Should all be set to the same as the subject Image
    subject->GetSpacing()[1]*(subjectRegion.GetSize()[1]-1),
    subject->GetSpacing()[2]*(subjectRegion.GetSize()[2]-1)
    ));

  BSplineTransform::MeshSizeType meshSize;                  //Setup a mesh that contains the number of controlpoints
  meshSize.Fill(BSplineControlPoints-BSplineOrder);         //Inspired from itk example "BSplineWarping2.cxx"

  bSpline->SetTransformDomainMeshSize(meshSize);            //TODO: ask if it is possible to have "uneven" control points.
                                                            // EG. 8 on LR axis 7 on SI 6 on AP. We did this in simple itk
                                                            //so it should be possible in cpp itk

  //Get the number of paramaters/nodes required for this BSpline
  const unsigned int numberOfParameters = bSpline->GetNumberOfParameters();
  const unsigned int numberOfNodes = numberOfParameters/Dimension;

  //print out
  std::cout << "Number of params: " << numberOfParameters << std::endl;
  std::cout << "Number of nodes:  " << numberOfNodes << std::endl;

  //Setup a paramaters variable for the bspline
  BSplineTransform::ParametersType bSplineParams( numberOfParameters );

  //  From ITK Example "BSplineWarping2"
  //  The B-spline grid should now be fed with coeficients at each node. Since
  //  this is a two dimensional grid, each node should receive two coefficients.
  //  Each coefficient pair is representing a displacement vector at this node.
  //  The coefficients can be passed to the B-spline in the form of an array where
  //  the first set of elements are the first component of the displacements for
  //  all the nodes, and the second set of elemets is formed by the second
  //  component of the displacements for all the nodes.

  //  In the ITK Example, the read the points in from a file. Here, they will be
  //  generated randomly. This should put the xyz coordinates in the correct space
  //  a better way would probably be to use the image coefficient array of the bspline,
  //  but for now I will use the method from the ITK example


  //initialize random number
  std::srand(time(nullptr));

  for( unsigned int n = 0; n < numberOfNodes; ++ n)
    {
    bSplineParams[n] = rand() % 5 + 1;                           // "x" coord;   rand number between 1 and 10
    bSplineParams[n + numberOfNodes] = rand() % 5 + 1;      // "y" coord;
    bSplineParams[n + numberOfNodes * 2] = rand() % 5 + 1;  // "z" coord;
                                                                  // TODO: x,y,z seem like they are the wrong coordinate system. Get a better model
    }

  bSpline->SetParameters(bSplineParams);

  // Apply transform to image with resampler:
  typedef itk::ResampleImageFilter<ImageType, ImageType> ResampleFilterType;
  ResampleFilterType::Pointer resampler = ResampleFilterType::New();

  typedef itk::LinearInterpolateImageFunction<ImageType, PixelType > InterpolatorType;
  InterpolatorType::Pointer interpolater = InterpolatorType::New();

  resampler->SetInterpolator(interpolater);
  resampler->SetOutputSpacing(subject->GetSpacing());
  resampler->SetOutputOrigin(subject->GetOrigin());
  resampler->SetOutputDirection(subject->GetDirection());
  resampler->SetSize(subjectRegion.GetSize());
  resampler->SetOutputStartIndex(subjectRegion.GetIndex());

  resampler->SetInput(imageReader->GetOutput());
  resampler->SetTransform(bSpline);


  //Setup ImageWriter
  typedef itk::ImageFileWriter<ImageType> DeformedSubjectFileWriterType;
  DeformedSubjectFileWriterType::Pointer deformedWriter = DeformedSubjectFileWriterType::New();
  deformedWriter->SetFileName(deformedImageName);
  deformedWriter->SetInput(resampler->GetOutput());
  deformedWriter->Update();
  std::cout << "Finished writing file " << std::endl;
  std::cout << "done" << std::endl;


  return EXIT_SUCCESS;
}