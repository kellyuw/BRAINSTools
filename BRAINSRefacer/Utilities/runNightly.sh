#!/bin/bash
# Script to run the Test Suite for BRAINSRefacer nightly
# Change the paths of the directories below
# This needs to be copied to and run from a location outside of the source directory.

BASE_DIR=~/code/brains/
SOURCE_DIR=${BASE_DIR}/BRAINSTools
BUILD_DIR=${BASE_DIR}/BRAINSTools-Refacer-Nightly-Release
EXTERNAL_SOURCES_DIR=${BASE_DIR}/BRAINSTools-SuperBuild-RelWithDebugInfo

if [[ ! -d ${BASE_DIR} ]]; then
  echo "base directory doesn't exist";
  exit -1;
fi

echo "updating source repo"
if [[ ! -d ${SOURCE_DIR} ]]; then
  echo "source directory doesn't exist. Cloning.";
  git clone http://github.com/BRAINSIa/BRAINSTools.git;
  cd ${SOURCE_DIR}
  git checkout BRAINSRefacer
else
  cd ${SOURCE_DIR};
  git fetch origin;
  git checkout master;
  git branch -D BRAINSRefacer
  git checkout BRAINSRefacer
fi

echo "changing to build directory"
if [[ ! -d ${BUILD_DIR} ]]; then
  echo "build directory doesnt exist. Creating ..."
  mkdir -p ${BUILD_DIR};
fi

cd ${BUILD_DIR};

echo "Configuring with cmake"
cmake ${SOURCE_DIR} -DBRAINSTools_SUPERBUILD:BOOL=OFF -DSlicerExecutionModel_DIR:PATH=${EXTERNAL_SOURCES_DIR}/SlicerExecutionModel-build -DITK_DIR:PATH=${EXTERNAL_SOURCES_DIR}/ITKv4-build/ -DCMAKE_CXX_STANDARD=11 -DUSE_BRAINSRefacer=ON -DUSE_ANTS=OFF -DUSE_AutoWorkup=OFF -DUSE_BRAINSABC=OFF -DUSE_BRAINSConstellationDetector=OFF -DUSE_BRAINSDWICleanup=OFF -DUSE_BRAINSFit=OFF -DUSE_BRAINSInitializedControlPoints=OFF -DUSE_BRAINSLabelStats=OFF -DUSE_BRAINSLandmarkInitializer=OFF -DUSE_BRAINSROIAuto=OFF -DUSE_BRAINSResample=OFF -DUSE_BRAINSSnapShotWriter=OFF -DUSE_BRAINSStripRotation=OFF -DUSE_BRAINSTransformConvert=OFF -DUSE_DWIConvert=OFF -DUSE_ConvertBetweenFileFormats=OFF -DUSE_ImageCalculator=OFF -DUSE_ReferenceAtlas=OFF

echo "Building";
make -j32 -k;
make -j32 -k;
make -j32 -k;

echo "Testing"
ctest -j16 -D NightlyStart;
ctest -j16 -D NightlyConfigure;
ctest -j16 -D NightlyBuild;
ctest -j16 -D NightlyTest;
ctest -j16 -D NightlySubmit;
