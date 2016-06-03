MakeTestDriverFromSEMTool(BRAINSDeface BRAINSDefaceTest.cxx)
message(WARNING "TestData should be placed into directory /scratch/testData/BRAINSDefaceData")

ExternalData_add_test( ${PROJECT_NAME}FetchData NAME BRAINSDefaceTest_RunTest
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSDefaceTestDriver>
  BRAINSDefaceTest
  --inputImage
    ${DefaceTestData_DIR}/t1_average_BRAINSABC.nii.gz
  --labelmap
    ${DefaceTestData_DIR}/JointFusion_HDAtlas20_2015_dustCleaned_label.nii.gz
  --landmarks
    ${DefaceTestData_DIR}/BCD_ACPC_Landmarks.fcsv
  --deformedImageName
    ${DefaceTestData_DIR}/defacedSubject.nii.gz
  --outputMask
    ${DefaceTestData_DIR}/outputMask.nii.gz
  --distanceMapFileName
    ${DefaceTestData_DIR}/distanceMap.nii.gz
  --bSplineFileName
    ${DefaceTestData_DIR}/bSpline.h5
  --smoothDisplacementName
    ${DefaceTestData_DIR}/bSplineTimesDistanceMap.nii.gz
  --diffImageName
    ${DefaceTestData_DIR}/diffImage.nii.gz
  --finalTransformFileName
    ${DefaceTestData_DIR}/finalTransform.h5
  )
