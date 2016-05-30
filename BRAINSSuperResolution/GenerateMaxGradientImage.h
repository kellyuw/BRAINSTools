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

template <class InputImageType, class OutputImageType>
typename OutputImageType::Pointer
GenerateMaxGradientImage(const std::vector<typename InputImageType::Pointer> & inputImages)
{

  typedef itk::GradientMagnitudeImageFilter<InputImageType,
                                            InputImageType>     GradientFilterType;
  typedef itk::MinimumMaximumImageCalculator<InputImageType>    MinMaxCalculatorType;
  typedef itk::Statistics::ScalarImageToHistogramGenerator<
    InputImageType>                                             HistogramGeneratorType;
  typedef typename HistogramGeneratorType::HistogramType        HistogramType;
  typedef typename itk::IntensityWindowingImageFilter<
    InputImageType,
    OutputImageType>                                            RescalerType;

  typedef std::vector<typename OutputImageType::Pointer>        RescaledImageGradientVector;

  typedef itk::MaximumImageFilter<OutputImageType, OutputImageType, OutputImageType>
    MaximumFilterType;
  MaximumFilterType::Pointer myMax = MaximumFilterType::New();

  itk::TimeProbe MaxGradientImageTimer;
  MaxGradientImageTimer.Start();


  const float UpperPercentileMatching = 0.95F; // Map 95th Quantile and above to 127
  const float LowerPercentileMatching = 0.50F; // Map 25th Quantile and below to 0

  const unsigned int numberOfImageModalities =
    inputImages.size(); // number of modality images

  RescaledImageGradientVector rescaledGradientImages( numberOfImageModalities )

  for( size_t i = 0; i < numberOfImageModalities; i++ )
     {
     typename GradientFilterType::Pointer gradientFilter = GradientFilterType::New();
     gradientFilter->SetInput( inputImages[i] );
     gradientFilter->Update();

     typename MinMaxCalculatorType::Pointer myMinMax = MinMaxCalculatorType::New();
     myMinMax->SetImage( gradientFilter->GetOutput() );
     myMinMax->Compute();

     typename HistogramGeneratorType::Pointer histogramGenerator =
       HistogramGeneratorType::New();
     histogramGenerator->SetInput( gradientFilter->GetOutput() );
     histogramGenerator->SetNumberOfBins(1024);
     histogramGenerator->SetMarginalScale(10);
     histogramGenerator->SetHistogramMin( myMinMax->GetMinimum() );
     histogramGenerator->SetHistogramMax( myMinMax->GetMaximum() );
     histogramGenerator->Compute();
     HistogramType::ConstPointer ith_histogram = histogramGenerator->GetOutput();

     typename RescalerType::Pointer rescaler = RescalerType::New();
     rescaler->SetInput( gradientFilter->GetOutput() );
     rescaler->SetOutputMinimum(0);
     rescaler->SetOutputMaximum(127);
     const float ithSlope =
       ( ith_histogram->Quantile(0, UpperPercentileMatching)
         - ith_histogram->Quantile(0, LowerPercentileMatching) ) / 80.0F;
     rescaler->SetWindowMinimum( ith_histogram->Quantile(0, LowerPercentileMatching) );
     rescaler->SetWindowMaximum( ith_histogram->Quantile(0, UpperPercentileMatching) + ithSlope * 100.0
                                   * ( 1.0 - UpperPercentileMatching ) );

     rescaledGradientImages[i] = rescaler->GetOutput();

     // TODO: find maximum of all rescaled gradient images

     //myMax->SetInput1( rescaler->GetOutput() );
     }

  try
    {
    myMax->Update();
    }
  catch( itk::ExceptionObject & exp )
    {
    std::cerr << "ExceptionObject with Iterator" << std::endl;
    std::cerr << exp << std::endl;
    return EXIT_FAILURE;
    }

  typename RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
  rescaler->SetOutputMinimum(0);
  rescaler->SetOutputMaximum(255);
  rescaler->SetInput( myMax->GetOutput() );
  rescaler->Update();

  typename OutputImageType::Pointer edgemap = rescaler->GetOutput();

  MaxGradientImageTimer.Stop();
  itk::RealTimeClock::TimeStampType elapsedTime = MaxGradientImageTimer.GetTotal();
  std::cout << "Generating maximum gradient edgemap took " << elapsedTime
            << " " << MaxGradientImageTimer.GetUnit() << "." << std::endl);

  return edgemap;
}

#endif // __GenerateMaxGradientImage_h
