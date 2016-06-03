
#add_executable(BRAINSDeface_tester BRAINSDeface_FirstTest.cxx)
#add_test(BRAINSDeface_test BRAINSDeface_tester)

MakeTestDriverFromSEMTool(BRAINSDeface BRAINSDefaceTest.cxx)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME BRAINSDefaceTest1
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSDefaceTestDriver>
  --compare
  DATA{${TestData_DIR}/BRAINSFitTest_AffineRotationMasks.result.nii.gz}
  ${CMAKE_CURRENT_BINARY_DIR}/applyWarp_test4.nii.gz
  --compareIntensityTolerance 30
  --compareRadiusTolerance 5
  --compareNumberOfPixelsTolerance 1
  BRAINSDefaceTest
  --inputVolume DATA{${TestData_DIR}/rotation.test.nii.gz}
  --referenceVolume DATA{${TestData_DIR}/test.nii.gz}
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/applyWarp_test4.nii.gz
  --pixelType uchar
  --warpTransform DATA{${TestData_DIR}/Transforms_h5/BRAINSFitTest_AffineRotationMasks.${XFRM_EXT}}
  )