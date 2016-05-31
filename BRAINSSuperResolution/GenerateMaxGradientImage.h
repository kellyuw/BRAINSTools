/*=========================================================================
 *
 *  Copyright SINAPSE: Scalable Informatics for Neuroscience, Processing and Software Engineering
 *            The University of Iowa
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
/*
 * Author: Ali Ghayoor
 * at SINAPSE Lab,
 * The University of Iowa 2016
 */

#ifndef __GenerateMaxGradientImage_h
#define __GenerateMaxGradientImage_h

//#include "itkImage.h"
#include "itkGradientMagnitudeImageFilter.h"
#include "itkScalarImageToHistogramGenerator.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkIntensityWindowingImageFilter.h"
#include "itkMaximumImageFilter.h"
#include "itkTimeProbe.h"
#include <vector>

// Class to print the proper exception message
class EmptyVectorException
{
public:
  EmptyVectorException(const char* pStr = "The list of input images was empty.  Nothing to find maximum!") :
    pMessage(pStr)
    {
    }

  const char * what() const
  {
  return pMessage;
  }

private:
  const char * pMessage;
};

// Auxiliary function to find the maximum of the rescaled gradient image list
template <typename TImage>
typename TImage::Pointer
MaxOfImageList(const std::vector<typename TImage::Pointer> & inputImageList) // inputImageList = rescaledGradientImageList
{
  typedef itk::MaximumImageFilter<TImage, TImage, TImage>   MaximumFilterType;

  if( inputImageList.empty() )
    {
    // No images, something went wrong.
    throw EmptyVectorException();
    }
  if( inputImageList.size() == 1 )
    {
    // Only one image, no need to find maximum image.
    return inputImageList[0];
    }

  // Initialize the maximum image with the first image in the list
  typename TImage::Pointer maxImage = inputImageList[0];

  for(unsigned int i = 1; i < inputImageList.size(); ++i)
     {
     typename MaximumFilterType::Pointer myMax = MaximumFilterType::New();
     myMax->SetInput1( maxImage );
     myMax->SetInput2( inputImageList[i] );
     try
       {
       myMax->Update();
       }
     catch( itk::ExceptionObject & exp )
       {
       std::cerr << "ExceptionObject with Iterator" << std::endl;
       std::cerr << exp << std::endl;
       }
     maxImage = myMax->GetOutput();
     }

  return maxImage;
}

/*
 * Main function to generate the maximum gradient image
 */
template <class InputImageType, class OutputImageType>
typename OutputImageType::Pointer
GenerateMaxGradientImage(const std::vector<typename InputImageType::Pointer> & inputImages)
{

  typedef itk::GradientMagnitudeImageFilter<InputImageType,
                                            InputImageType>                   GradientFilterType;
  typedef itk::MinimumMaximumImageCalculator<InputImageType>                  MinMaxCalculatorType;
  typedef itk::Statistics::ScalarImageToHistogramGenerator<InputImageType>    HistogramGeneratorType;
  typedef typename HistogramGeneratorType::HistogramType                      HistogramType;
  typedef typename itk::IntensityWindowingImageFilter<InputImageType,
                                                      OutputImageType>        RescalerType;
  typedef std::vector<typename OutputImageType::Pointer>                      RescaledImageGradientVectorType;

  itk::TimeProbe MaxGradientImageTimer;
  MaxGradientImageTimer.Start();


  const float UpperPercentileMatching = 0.95F; // Map 95th Quantile and above to 127
  const float LowerPercentileMatching = 0.50F; // Map 25th Quantile and below to 0

  const unsigned int numberOfImageModalities =
    inputImages.size(); // number of modality images

  RescaledImageGradientVectorType rescaledGradientImageList( numberOfImageModalities );

  for( size_t i = 0; i < numberOfImageModalities; i++ )
     {
     //
     // TODO: Find quantiles within brain mask
     //
     typename GradientFilterType::Pointer gradientFilter = GradientFilterType::New();
     gradientFilter->SetInput( inputImages[i] );
     gradientFilter->Update();

     typename MinMaxCalculatorType::Pointer myMinMax = MinMaxCalculatorType::New();
     myMinMax->SetImage( gradientFilter->GetOutput() );
     myMinMax->Compute();

     typename HistogramGeneratorType::Pointer histogramGenerator =
       HistogramGeneratorType::New();
     histogramGenerator->SetInput( gradientFilter->GetOutput() );
     histogramGenerator->SetNumberOfBins(1024); // 4x oversampling to put into an unsigned char image.
     histogramGenerator->SetMarginalScale(10);
     histogramGenerator->SetHistogramMin( myMinMax->GetMinimum() );
     histogramGenerator->SetHistogramMax( myMinMax->GetMaximum() );
     histogramGenerator->Compute();
     typename HistogramType::ConstPointer ith_histogram = histogramGenerator->GetOutput();

     typename RescalerType::Pointer rescaler = RescalerType::New();
     rescaler->SetInput( gradientFilter->GetOutput() );
     rescaler->SetOutputMinimum(0);
     rescaler->SetOutputMaximum(127);
     const float slope =
       ( ith_histogram->Quantile(0, UpperPercentileMatching)
         - ith_histogram->Quantile(0, LowerPercentileMatching) ) / 80.0F;
     rescaler->SetWindowMinimum( ith_histogram->Quantile(0, LowerPercentileMatching) );
     rescaler->SetWindowMaximum( ith_histogram->Quantile(0, UpperPercentileMatching) + slope * 100.0
                                   * ( 1.0 - UpperPercentileMatching ) );

     rescaledGradientImageList[i] = rescaler->GetOutput();
     }

  typename OutputImageType::Pointer maxGI = MaxOfImageList<OutputImageType>(rescaledGradientImageList);

  typedef itk::IntensityWindowingImageFilter<OutputImageType,
                                             OutputImageType>           RescaleFilterType;
  typename RescaleFilterType::Pointer outputRescaler = RescaleFilterType::New();
  outputRescaler->SetOutputMinimum(0);
  outputRescaler->SetOutputMaximum(255);
  outputRescaler->SetInput( maxGI );
  outputRescaler->Update();

  MaxGradientImageTimer.Stop();
  itk::RealTimeClock::TimeStampType elapsedTime = MaxGradientImageTimer.GetTotal();
  std::cout << "Generating maximum gradient edgemap took " << elapsedTime
            << " " << MaxGradientImageTimer.GetUnit() << "." << std::endl;

  return outputRescaler->GetOutput();
}

#endif // __GenerateMaxGradientImage_h
