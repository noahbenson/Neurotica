(* MRImage.m
 *
 * Utility functions for dealing with (mostly FreeSurfer) cortical volume image data in Mathematica.
 *
 * Copyright (C) 2014-2015 by Noah C. Benson.
 * This file is part of the Neurotica library.
 *
 * This program is free software: you can redistribute it and/or modify it under the terms of the
 * GNU General Public License as published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See
 * the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program.
 * If not, see <http://www.gnu.org/licenses/>.
 *)

(**************************************************************************************************)
BeginPackage["Neurotica`MRImage`", {"Neurotica`Global`", "Neurotica`Util`"}];
Unprotect["Neurotica`MRImage`*", "Neurotica`MRImage`Private`*"];
ClearAll["Neurotica`MRImage`*", "Neurotica`MRImage`Private`*"];

MRImage3D::usage = "MRImage3D[data] yields an Image3D-like form that can be used with Image3D functions but which stores additional relevant data regarding MR images. MRImage3D accepts all the options that can be passed to Image3D as well as the following:
  * Indeterminate (defailt: Min) specifies what value should be filled in for Indeterminate voxel values when drawing the Image3D representation; any key of MRImageStatistics may be given.
  * None (default: Min) specifies the value that should be filled in for None voxel values.
  * RightDirectionVector, AnteriorDirectionVector, and SuperiorDirectionVector (defaults: Indeterminate) are necessary for acurrate MRISlice and MRIOrient function calls; these each specify a vector telling which direction in the image is right, anterior, or superior, respectively. The vectors are interpreted to be in a space such that the vectors {1,0,0}, {0,1,0}, and {0,0,1} point from the center of voxel [[i,j,k]] to the centers of voxels [[i+1,j,k]], [[i,j+1,k]], and [[i,j,k+1]], respectively.
  * Center (default: Automatic) specifies the center of the brain or item of interest in the MRImage; if Automatic is given, then ImageDimensions/2 is used, otherwise a 3D vector in terms of voxel indices should be given (same coordinate orientation as for the direction vectors, e.g. RightDirectionVector).
  * VoxelDimensions (default: {1,1,1}) specifies the spatial dimensions of the voxels. Quantity units may be given.";

MRImage3D::badarg = "Bad argument given to MRImage3D: `1`";

MRImageStatistics::usage = "MRImageStatistics[mrimg] yields an association of a few critical statistics of all the valid numerical values in the image; these statistics are Min, Max, Mean, Median, Variance, Count, and Missing. For more information, see MRImageMax, MRImageCount, etc.";
MRImageMax::usage = "MRImageMax[mrimg] yields the maximum value of all voxels in the given MRImage object mrimg. This value ignores Ideterminate and None values.";
MRImageMin::usage = "MRImageMin[mrimg] yields the minimum value of all voxels in the given MRImage object mrimg. This value ignores Ideterminate and None values.";
MRImageMean::usage = "MRImageMean[mrimg] yields the mean value of all voxels in the given MRImage object mrimg. This value ignores Ideterminate and None values.";
MRImageMedian::usage = "MRImageMedian[mrimg] yields the median value of all voxels in the given MRImage object mrimg. This value ignores Ideterminate and None values.";
MRImageVariance::usage = "MRImageVariance[mrimg] yields the variance of all voxels in the given MRImage object mrimg. This value ignores Ideterminate and None values.";
MRImageCount::usage = "MRImageCount[mrimg] yields the number of valid (numerical) values, including separate frames of the image, in the given MRImage object mrimg. This value ignores Ideterminate and None values.";
MRImageMissing::usage = "MRImageMissing[mrimg] yields the number of invalid (Indeterminate or None) values, including separate frames of the image, in the given MRImage object mrimg.";

MRImageQ::usage = "MRImageQ[img] yields True if img is a valid MRImage3D object and yields False otherwise.";
MRImageSliceQ::usage = "MRImageSliceQ[img] yields True if img is a valid MRImage object and yields False otherwise.";
MRImageObjectQ::usage = "MRImageObjectQ[img] is equivalent to Or[MRImageQ[img], MRImageSliceQ[img]].";

RightDirectionVector::usage = "RightDirectionVector[img] yields the vector that points to the right for the given MRImage3D object img or yields Indeterminate if none can be determined. RightDirectionVector is also an optional argument to MRImage3D that specifies the image's rightward direction. In both cases, the vector should assume that the basis vecotrs i, j, and k point in the direction of increasing index in the ImageData[img] matrix; i.e. i, j, and k are the vectors from the center of voxel ImageData[img][[1,1,1]] to the centers of ImageData[img][[2,1,1]], ImageData[img][[1,2,1]], and ImageData[img][[1,1,2]], respectively.";
AnteriorDirectionVector::usage = "AnteriorDirectionVector[img] yields the vector that points to the anterior for the given MRImage3D object img or yields Indeterminate if none can be determined. AnteriorDirectionVector is also an optional argument to MRImage3D that specifies the image's forward direction (in terms of the brain). In both cases, the vector should assume that the basis vecotrs i, j, and k point in the direction of increasing index in the ImageData[img] matrix; i.e. i, j, and k are the vectors from the center of voxel ImageData[img][[1,1,1]] to the centers of ImageData[img][[2,1,1]], ImageData[img][[1,2,1]], and ImageData[img][[1,1,2]], respectively.";
SuperiorDirectionVector::usage = "SuperiorDirectionVector[img] yields the vector that points to the superior for the given MRImage3D object img or yields Indeterminate if none can be determined. SuperiorDirectionVector is also an optional argument to MRImage3D that specifies the image's upward direction (in terms of the brain). In both cases, the vector should assume that the basis vecotrs i, j, and k point in the direction of increasing index in the ImageData[img] matrix; i.e. i, j, and k are the vectors from the center of voxel ImageData[img][[1,1,1]] to the centers of ImageData[img][[2,1,1]], ImageData[img][[1,2,1]], and ImageData[img][[1,1,2]], respectively.";

VoxelDimensions::usage = "VoxelDimensions[img] yields the {i,j,k} voxel image dimensions for the given MRImage3D img or the {i,j} voxel image dimensions for the MRImage slice img. VoxelDimensions is also an argument that can be passed to both MRImage3D or MRImage to specify the size of the voxels.";
PixelDimensions::usage = "PixelDimensions[img] is equivalent to VoxelDimensions but for MRImage slices.";

VoxelIndexToRASMatrix::usage = "VoxelIndexToRASMatrix[img] yields a 4x4 matrix that will transform the voxel index {i,j,k,1} into the coordinate {x,y,z,1} for the given MRImage3D, img. The resulting matrix will orient the right, anterior, and superior directions in the positive x, y, and z axes respectively.
VoxelIndexToRASMatrix[img, spec] yields a transformation matrix that orients x, y, and z in the directions specified in spec, which is identical to the first argument of MRIOrient.";
VoxelIndexToRASTransform::usage = "VoxelIndexToRASTransform[img] yields a TransformationFunction that will transform the voxel index {i,j,k} into the coordinate {x,y,z} for the given MRImage3D, img.
VoxelIndexToRASTransform[img, spec] yields a transformation function that orients x, y, and z in the directions specified in spec, which is identical to the first argument of MRIOrient.";
VoxelIndexToRAS::usage = "VoxelIndexToRAS[img, {i,j,k}] yields the coordinate at which the voxel (i, j, k) is centered in the given MRImage3D, img.
VoxelIndexToRAS[img, {{i0,j0,k0}, {i1,j1,k1}, ...}] yields the coordinates for all of the indices.
VoxelIndexToRAS[img, coords, spec] performes the transformation but orients the positive x, y, and z axes in the directions given by spec, which is identical in format to the second argument of MRIOrient.";
VoxelIndexToRASTr::usage = "VoxelIndexToRASTr[img, Q, spec] is equivalent to Transpose @ VoxelIndexToRAS[img, Transpose @ Q, spec].";
VoxelIndexToRAS::badimg = "Given MRImage3D object does not define a full orientation";

RASToVoxelIndexMatrix::usage = "RASToVoxelIndexMatrix[img] yields a matrix that will transform the voxel index {x,y,z,1} into the voxel coordinate {i,j,k} for the given MRImage3D, img. Note that the coordinate may not be even integer numbers; use Round to find indices. The resulting indices will assume that the x, y, and z values are specified as coordinates in which the positive x, y, and z axes are pointed in the right, anterior, and superior directions, respectively.
RASToVoxelIndexMatrix[img, spec] yields a matrix that performs the same conversion but assumes that the x, y, and z values are specified in coordinates in which the positive x, y, and z axes are pointed in the directions indicated by spec, which is identical in format to the second argument of MRIOrient.";
RASToVoxelIndexTransform::usage = "RASToVoxelIndexTransform[img] yields a TransformationFunction that will transform the voxel index {x,y,z} into the voxel coordinate {i,j,k} for the given MRImage3D, img. Note that the coordinate may not be even integer numbers; use Round to find indices. The resulting indices will assume that the x, y, and z values are specified as coordinates in which the positive x, y, and z axes are pointed in the right, anterior, and superior directions, respectively.
RASToVoxelIndexTransform[img, spec] yields a transformation function that performs the same conversion but assumes that the x, y, and z values are specified in coordinates in which the positive x, y, and z axes are pointed in the directions indicated by spec, which is identical in format to the second argument of MRIOrient.";
RASToVoxelIndex::usage = "RASToVoxelIndex[img, {x,y,z}] yields the voxel index (i, j, k) in the given MRImage3D, img, that corresponds to the given (x, y, z) values. Note that the coordinate may not be even integer numbers; use Round to find indices. The resulting indices will assume that the x, y, and z values are specified as coordinates in which the positive x, y, and z axes are pointed in the right, anterior, and superior directions, respectively.
RASToVoxelIndex[img, {{x0,y0,z0}, {x1,y1,z1}, ...}] yields the voxel indices for all of the coordinates.
RASToVoxelIndex[img, coords, spec] performs the same conversion but assumes that the x, y, and z values are specified in coordinates in which the positive x, y, and z axes are pointed in the directions indicated by spec, which is identical in format to the second argument of MRIOrient.";
RASToVoxelIndexTr::usage = "RASToVoxelIndexTr[img, Q, spec] is equivalent to Transpose @ RASToVoxelIndex[img, Transpose @ Q, spec].";
RASToVoxelIndex::badimg = "Given MRImage3D object does not define a full orientation";

MRITransformation::usage = "MRITransformation[img, tx] is identical to ImageTransformation[img, tx] for MRImage3D or MRImage slice object img and transformation tx, but performs the transformation in a manner similar to how FreeSurfer and other neuroimaging software do, resulting in a slightly modified transformation than the typical Mathematica mode. Critically, this maintains orientation information through the transformation, so should be used with MRImages.";

MRIOrient::usage = "MRIOrient[img, {rows, columns, slices}] yields an MRImage3D identical to the given MRImage3D, img, but with the rows, columns, and slices occurring according to the given directives. The three arguments may take the form Left/LH, Right/RH, Inferior/Bottom, Superior/Top, Posterior/Back, or Anterior/Front. Arguments may only appear once in the list and may not appear with their antonyms.
MRIOrient[img, {rows, columns, slices}, dimensions] orients the given img as instructed and additionally centers it in an image with the number of rows, columns, and slices given in the list dimensions.
Examples:
  * MRIOrient[img, {Right, Anterior, Superior}] yields img in a RAS orientation.
  * MRIOrient[img, {RH, Front, Top}, {256, 256, 256}] is equivalent to the example above except that it also centers the image in 256 x 256 x 256 voxel space.
  * MRIOrient[img, {Left, Posterior, Superior}] yields img in a scanner orientation.
  * MRIOrient[img, {Anterior, Superior, Back}] is invalid because Anterior and Back are antonyms.
MRIOrient may also be called with string arguments such as \"RAS\" or \"LIA\", which are automatically converted into the appropriate lists.  
MRIOrient[img, {rows, columns}] orients the rows and columns of an MRImage slice as with an MRImage3D object; note that the rows and columns must correspond to axes that the given image represents; e.g., you cannot slice parallel the Saggital plane then orient the resulting image slices in left or right directions.";
MRIOrient::badarg = "Invalid instructions given to MRIOrient: `1`";
MRIOrientMatrix::usage = "MRIOrientMatrix[img, spec] and MRIOrientMatrix[img, spec, dims] both yield the matrix used in the transformation undertaken by MRIOrient[img, spec, dims]. Note that these yield matrices that are appropriate for use with ImageForwardTransformation or MRITransformation, but not ImageTransformation.";
MRIOrientTransform::usage = "MRIOrientTransform[img, spec] and MRIOrientMatrix[img, spec, dims] both yield the transformation function used in the transformation undertaken by MRIOrient[img, spec, dims]. Note that these yield transform functions that are appropriate for use with ImageForwardTransformation or MRITransformation, but not ImageTransformation.";

MRISlices::usage = "MRISlices[img, plane] yields a list of MRImage objects (2D MRImage slices) corresponding to the given plane. Plane may be specified in any of the following ways:
  * \"Sagittal\" or \"Lateral\" indicate the slices parallel to the Sagittal plane, which divides left from right, ordering the slices from left to right.
  * \"Transverse\", \"Axial\", or \"Horizontal\" indicate the slices parallel to the Horizontal plane, which divides top superior from inferior, ordering the slices from inferior to superior.
  * \"Coronal\" or \"Frontal\" indicate the slices parallel to the Coronal plane, which divides posterior from anterior, ordering the slices from posterior to anterior.
  * Left -> Right, Anterior -> Back, LH -> RH, Top -> Inferior, etc. may also be used to indicate that the slices should be ordered in a particular direction.
The option MRIOrient may be specified as well; if not specified, the default value, Automatic, indicates that the orientation provided in the image should be preserved in the resulting slices; otherwise, the value must be a 2-element list which is a valid argument to the MRIOrient[] function with any of the 2D image slices.
MRISlices[img, plane, slices] is equivalent to MRISlices[img, plane][[slices]].";
MRISlices::badplane = "Could not recognized plane given to MRISlices: `1`";

(**************************************************************************************************)
Begin["`Private`"];

(*  Notes:
 *  Mathematica plots Image3D objects in an odd fashion; the real question is now x,y, and z line
 *  up with the voxel indices i, j, k.
 *  Assuming that, for an image with the data array d, voxel i, j, k is d[[i,j,k]], then the 
 *  transformation between (x,y,z) and (i,j,k) is:
 *  {x -> k, y -> -j, z -> -i}
 *  and the reverse is:
 *  {i -> -z, j -> -y, k -> x}
 *  Where the negatives just indicate a direction reversal (really they should be max - <x>).
 *)
MMAImageToXYZTransform[img_Image3D] := Fold[
  ImageReflect,
  img,
  {Top -> Bottom, Front -> Back}];
XYZToMMAImageTransform[img_Image3D] := Fold[
  ImageReflect,
  img,
  {Back -> Front, Bottom -> Top}];
Protect[MMAImageToXYZTransform, XYZToMMAImageTransform];

(* #MRImage3D *************************************************************************************)

(* Because the MRImage3D and MRImage objects are so similar, we have a lot of overlapping functions,
 * which we store in this large Hold form for use in DefineImmutable[...].
 *)
$MRImageSharedMethods = Hold[
  (* Simple accessors for Image3D computed data *)
  ImageAspectRatio[img, opts___] := ImageAspectRatio[Image3D[img], opts],
  ImageChannels[img, opts___] := ImageChannels[Image3D[img], opts],
  ImageColorSpace[img, opts___] := ImageColorSpace[Image3D[img], opts],
  ImageDimensions[img, opts___] := ImageDimensions[Image3D[img], opts],
  ImageLevels[img, opts___] := ImageLevels[Image3D[img], opts],
  ImageMeasurements[img, opts___] := ImageMeasurements[Image3D[img], opts],

  (* Directly access values *)
  ImageValue[img, opts___] := ImageValue[Image3D[img], opts],
  PixelValue[img, opts___] := PixelValue[Image3D[img], opts],
  ImageValuePosition[img, opts___] := ImageValuePosition[Image3D[img], opts],
  PixelValuePosition[img, opts___] := PixelValuePosition[Image3D[img], opts],

  (* More complex query operations over the image *)
  ImageAccumulate[img, opts___] := ImageAccumulate[Image3D[img], opts],
  ImageHistogram[img, opts___] := ImageHistogram[Image3D[img], opts],
  DominantColors[img, opts___] := DominantColors[Image3D[img], opts],
  ImageCooccurance[img, opts___] := ImageCooccurance[Image3D[img], opts],
  BorderDimensions[img, opts___] := BorderDimensions[Image3D[img], opts],
  BinaryImageQ[img, opts___] := BinaryImageQ[Image3D[img], opts],
  AlphaChannel[img, opts___] := AlphaChannel[Image3D[img], opts],
  FindThreshold[img, opts___] := Plus[
    MRImageMin[img],
    (MRImageMax[img] - MRImageMin[img]) * FindThreshold[Image3D[img], opts]],

  (* Image Transformations *)
  Blur[img, opts___] := With[
    {res = Check[Blur[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img]},
        MRImage3D[img, ImageData -> newdat]]]],
  ImageAdd[img, imgs__] := With[
    {res = Check[
       ImageData[img] + Total[If[NumberQ[#] || ArrayQ[#], #, ImageData[#]]& /@ {imgs}],
       $Failed]},
    If[res === $Failed, res, MRImage3D[res, Sequence @@ Options[img]]]],
  ImageConvolve[img, opts___] := With[
    {res = Check[ImageConvolve[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, Evaluate[Sequence@@theOpts]]]]],
  ImageCorrelate[img, opts___] := ImageCorrelate[Image3D[img], opts],
  ImageCorrespondingPoints[img, opts___] := ImageCorrespondingPoints[Image3D[img], opts],
  ImageCrop[img, opts___] := With[
    {res = Check[ImageCrop[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, Evaluate[Sequence@@theOpts]]]]],
  ImageDeconvolve[img, opts___] := With[
    {res = Check[ImageDeconvolve[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, Evaluate[Sequence@@theOpts]]]]],
  ImageDistance[img, opts___] := ImageDistance[Image3D[img], opts],
  ImageForwardTransformation[img, opts___] := With[
    {res = Check[ImageForwardTransformation[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, Evaluate[Sequence@@theOpts]]]]],
  ImagePad[img, opts___] := With[
    {res = Check[ImagePad[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, Evaluate[Sequence@@theOpts]]]]],
  ImagePartition[img, opts___] := With[
    {res = Check[ImagePartition[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, Evaluate[Sequence@@theOpts]]]]],
  ImagePerspectiveTransformation[img, opts___] := With[
    {res = Check[ImagePerspectiveTransformation[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, Evaluate[Sequence@@theOpts]]]]],
  ImageReflect[img, opts___] := With[
    {res = Check[ImageReflect[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, Evaluate[Sequence@@theOpts]]]]],
  ImageResize[img, opts___] := With[
    {res = Check[ImageResize[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, Evaluate[Sequence@@theOpts]]]]],
  ImageRotate[img, opts___] := With[
    {res = Check[ImageRotate[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, Evaluate[Sequence@@theOpts]]]]],
  ImageTake[img, opts___] := With[
    {res = Check[ImageTake[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, Evaluate[Sequence@@theOpts]]]]],
  ImageTransformation[img, opts___] := With[
    {res = Check[ImageTransformation[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img]},
        MRImage3D[img, ImageData -> newdat]]]],
  ImageTrim[img, opts___] := With[
    {res = Check[ImageTrim[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, Evaluate[Sequence@@theOpts]]]]],
  Sharpen[img, opts___] := With[
    {res = Check[Sharpen[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, Evaluate[Sequence@@theOpts]]]]],

  (* Image Filters *)
  BandpassFilter[img, opts___] := With[
    {res = Check[BandpassFilter[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, Evaluate[Sequence@@theOpts]]]]],
  BandstopFilter[img, opts___] := With[
    {res = Check[BandstopFilter[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, Evaluate[Sequence@@theOpts]]]]],
  BilateralFilter[img, opts___] := With[
    {res = Check[BilateralFilter[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, Evaluate[Sequence@@theOpts]]]]],
  CommonestFilter[img, opts___] := With[
    {res = Check[CommonestFilter[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, Evaluate[Sequence@@theOpts]]]]],
  CurvatureFlowFilter[img, opts___] := With[
    {res = Check[CurvatureFlowFilter[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, Evaluate[Sequence@@theOpts]]]]],
  DerivativeFilter[img, opts___] := With[
    {res = Check[DerivativeFilter[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, Evaluate[Sequence@@theOpts]]]]],
  DifferentiatorFilter[img, opts___] := With[
    {res = Check[DifferentiatorFilter[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, Evaluate[Sequence@@theOpts]]]]],
  EntropyFilter[img, opts___] := With[
    {res = Check[EntropyFilter[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, Evaluate[Sequence@@theOpts]]]]],
  GaborFilter[img, opts___] := With[
    {res = Check[GaborFilter[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, Evaluate[Sequence@@theOpts]]]]],
  GaussianFilter[img, opts___] := With[
    {res = Check[GaussianFilter[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, Evaluate[Sequence@@theOpts]]]]],
  GeometricMeanFilter[img, opts___] := With[
    {res = Check[GeometricMeanFilter[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, Evaluate[Sequence@@theOpts]]]]],
  GradientFilter[img, opts___] := With[
    {res = Check[GradientFilter[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, Evaluate[Sequence@@theOpts]]]]],
  GradientOrientedFilter[img, opts___] := With[
    {res = Check[GradientOrientedFilter[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, Evaluate[Sequence@@theOpts]]]]],
  HarmonicMeanFilter[img, opts___] := With[
    {res = Check[HarmonicMeanFilter[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, Evaluate[Sequence@@theOpts]]]]],
  HighpassFilter[img, opts___] := With[
    {res = Check[HighpassFilter[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, Evaluate[Sequence@@theOpts]]]]],
  HilbertFilter[img, opts___] := With[
    {res = Check[HilbertFilter[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, Evaluate[Sequence@@theOpts]]]]],
  KuwaharaFilter[img, opts___] := With[
    {res = Check[KuwaharaFilter[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, Evaluate[Sequence@@theOpts]]]]],
  LaplacianFilter[img, opts___] := With[
    {res = Check[LaplacianFilter[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, Evaluate[Sequence@@theOpts]]]]],
  LaplacianGaussianFilter[img, opts___] := With[
    {res = Check[LaplacianGaussianFilter[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, Evaluate[Sequence@@theOpts]]]]],
  LowpassFilter[img, opts___] := With[
    {res = Check[LowpassFilter[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, Evaluate[Sequence@@theOpts]]]]],
  MaxFilter[img, opts___] := With[
    {res = Check[MaxFilter[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, Evaluate[Sequence@@theOpts]]]]],
  MeanFilter[img, opts___] := With[
    {res = Check[MeanFilter[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, Evaluate[Sequence@@theOpts]]]]],
  MeanShiftFilter[img, opts___] := With[
    {res = Check[MeanShiftFilter[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, Evaluate[Sequence@@theOpts]]]]],
  MedianFilter[img, opts___] := With[
    {res = Check[MedianFilter[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, Evaluate[Sequence@@theOpts]]]]],
  MinFilter[img, opts___] := With[
    {res = Check[MinFilter[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, Evaluate[Sequence@@theOpts]]]]],
  PeronaMalikFilter[img, opts___] := With[
    {res = Check[PeronaMalikFilter[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, Evaluate[Sequence@@theOpts]]]]],
  RangeFilter[img, opts___] := With[
    {res = Check[RangeFilter[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, Evaluate[Sequence@@theOpts]]]]],
  RidgeFilter[img, opts___] := With[
    {res = Check[RidgeFilter[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, Evaluate[Sequence@@theOpts]]]]],
  StandardDeviationFilter[img, opts___] := With[
    {res = Check[StandardDeviationFilter[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, Evaluate[Sequence@@theOpts]]]]],
  TotalVariationFilter[img, opts___] := With[
    {res = Check[TotalVariationFilter[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, Evaluate[Sequence@@theOpts]]]]],
  WienerFilter[img, opts___] := With[
    {res = Check[WienerFilter[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, Evaluate[Sequence@@theOpts]]]]],

  (* Morphological Operations *)
  Binarize[img, opts___] := Binarize[Image3D[img], opts],
  BottomHatTransform[img, opts___] := BottomHatTransform[Image3D[img], opts],
  Closing[img, opts___] := With[
    {res = Check[Closing[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, Evaluate[Sequence@@theOpts]]]]],
  ColorNegate[img, opts___] := ColorNegate[Image3D[img], opts],
  ComponentMeasurements[img, opts___] := ComponentMeasurements[Image3D[img], opts],
  Colorize[img, opts___] := Colorize[Image3D[img], opts],
  DeleteBorderComponents[img, opts___] := With[
    {res = Check[DeleteBorderComponents[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, Evaluate[Sequence@@theOpts]]]]],
  DeleteSmallComponents[img, opts___] := With[
    {res = Check[DeleteSmallComponents[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, Evaluate[Sequence@@theOpts]]]]],
  Dilation[img, opts___] := With[
    {res = Check[Dilation[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, Evaluate[Sequence@@theOpts]]]]],
  DistanceTranform[img, opts___] := DistanceTransgorm[Image3D[img], opts],
  Erosion[img, opts___] := With[
    {res = Check[Erosion[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, Evaluate[Sequence@@theOpts]]]]],
  FillingTransform[img, opts___] := With[
    {res = Check[FillingTransform[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, Evaluate[Sequence@@theOpts]]]]],
  HighlightImage[img, opts___] := HighlightImage[Image3D[img], opts],
  HitMissTransform[img, opts___] := With[
    {res = Check[HitMissTransform[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, Evaluate[Sequence@@theOpts]]]]],
  InverseDistanceTransform[img, opts___] := With[
    {res = Check[InverseDistanceTransform[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, Evaluate[Sequence@@theOpts]]]]],
  MaxDetect[img, opts___] := MaxDetect[Image3D[img], opts],
  MinDetect[img, opts___] := MinDetect[Image3D[img], opts],
  MorphologicalBinarize[img, opts___] := MorphologicalBinarize[Image3D[img], opts],
  MorphologicalBranchPoints[img, opts___] := MorphologicalBranchPoints[Image3D[img], opts],
  MorphologicalComponents[img, opts___] := MorphologicalComponents[Image3D[img], opts],
  MorphologicalEulerNumber[img, opts___] := MorphologicalEulerNumber[Image3D[img], opts],
  MorphologicalGraph[img, opts___] := MorphologicalGraph[Image3D[img], opts],
  MorphologicalPerimeter[img, opts___] := MorphologicalPerimeter[Image3D[img], opts],
  MorphologicalTransform[img, opts___] := MorphologicalTransform[Image3D[img], opts],
  Opening[img, opts___] := With[
    {res = Check[Opening[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, Evaluate[Sequence@@theOpts]]]]],
  Pruning[img, opts___] := With[
    {res = Check[Pruning[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, Evaluate[Sequence@@theOpts]]]]],
  SelectComponents[img, opts___] := SelectComponents[Image3D[img], opts],
  SkeletonTranform[img, opts___] := SkeletonTransgorm[Image3D[img], opts],
  Thinning[img, opts___] := Thinning[Image3D[img], opts],
  TopHatTransform[img, opts___] := TopHatTransform[Image3D[img], opts]];


(* #MRImage3D *************************************************************************************)
Options[MRImage3D] = Join[
  {RightDirectionVector -> Indeterminate,
   AnteriorDirectionVector -> Indeterminate,
   SuperiorDirectionVector -> Indeterminate,
   Center -> Automatic,
   VoxelDimensions -> {1,1,1},
   Indeterminate -> Min,
   None -> Min},
  Replace[
    Options[Image3D],
    {(ColorFunction -> _) -> (ColorFunction -> "XRay"),
     (Boxed -> _) -> (Boxed -> False),
     (Background -> _) -> (Background -> Black)},
    {1}]];
(* For redefining options of MRImage3D *)
MRImage3D[img_Image3D, opts___Rule] := MRImage3D[
  ImageData[img],
  opts,
  Sequence@@Options[img]];
MRImage3D[img_MRImage3D, opts___Rule] := With[
  {settables = FilterRules[{opts}, {ImageData}],
   options = DeleteCases[{opts}, HoldPattern[ImageData -> _], {1}]},
  Clone[
    img,
    Sequence @@ settables,
    Options -> Join[options, Options[img]]]];
DefineImmutable[
  MRImage3D[data_List, OptionsPattern[]] :> img,
  Evaluate[
    Join[
      $MRImageSharedMethods,
      Hold[
        (* Retreive (or edit) the raw image data *)
        ImageData[img] = N @ If[Length@Dimensions@data == 4 && Last@Dimensions@data == 1,
          data[[All,All,All,1]],
          data],
        (* Access the scale factors of min and max for the data *)
        MRImageStatistics[img] -> Fold[
          Function @ If[AssociationQ[#1],
            #1,
            With[
              {min = Min[#1]},
              Which[
                NumericQ[min], Association[
                  {Min -> min, Max -> Max[#1],
                   Mean -> Mean[#1], Median -> Median[#1],
                   Variance -> Variance[#1], Count -> Length[#1],
                   Missing -> ((Times @@ Dimensions[ImageData[img]]) - Length[#1])}],
                MatchQ[max, #2[[1]]], DeleteCases[#1, #2[[2]], {3}],
                True, #1]]],
          Flatten @ ImageData[img],
          {Indeterminate -> Indeterminate, Max[_, None] -> None}],
        MRImageMax[img] := MRImageStatistics[img][Max],
        MRImageMin[img] := MRImageStatistics[img][Min],
        MRImageMean[img] := MRImageStatistics[img][Mean],
        MRImageMedian[img] := MRImageStatistics[img][Median],
        MRImageVariance[img] := MRImageStatistics[img][Variance],
        MRImageCount[img] := MRImageStatistics[img][Count],
        MRImageMissing[img] := MRImageStatistics[img][Missing],
        (* Retreive (or edit) the raw image options *)
        Options[img] = Map[(# -> OptionValue[#])&, Options[MRImage3D][[All, 1]]],
        Options[img, opt_] := Replace[opt, Options[img]],
        (* Obtain a few options specifically! *)
        Center[img] -> Replace[Center /. Options[img], Automatic :> (ImageDimensions[img]/2)],
        RightDirectionVector[img] -> (RightDirectionVector /. Options[img]),
        AnteriorDirectionVector[img] -> (AnteriorDirectionVector /. Options[img]),
        SuperiorDirectionVector[img] -> (SuperiorDirectionVector /. Options[img]),
        VoxelDimensions[img] -> (VoxelDimensions /. Options[img]),
        
        (* Matrices and Transformations... *)
        VoxelIndexToRASMatrix[img] :> With[
          {mtx = {RightDirectionVector[img],
                  AnteriorDirectionVector[img],
                  SuperiorDirectionVector[img]}},
          With[
            {x0 = -mtx.(Center[img] + 0.5*ImageDimensions[img])},
            If[!MatrixQ[Append[mtx,x0], NumericQ],
              Indeterminate,
              Append[MapThread[Append, {mtx, x0}], {0,0,0,1}]]]],
        RASToVoxelIndexMatrix[img] :> With[
          {mtx = VoxelIndexToRASMatrix[img]},
          If[mtx === Indeterminate || mtx === $Failed, mtx, Inverse[mtx]]],
        VoxelIndexToRASTransform[img] := With[
          {mtx = VoxelIndexToRASMatrix[img]},
          If[mtx === Indeterminate || mtx === $Failed,
            mtx,
            AffineTransform[{mtx[[1;;3, 1;;3]], mtx[[1;;3, 4]]}]]],
        RASToVoxelIndexTransform[img] := With[
          {mtx = RASToVoxelIndexMatrix[img]},
          If[mtx === Indeterminate || mtx === $Failed,
            mtx,
            AffineTransform[{mtx[[1;;3, 1;;3]], mtx[[1;;3, 4]]}]]],
                
        (* This is an MRImage... Checks of data quality should go here as well. *)
        MRImageQ[img] -> With[
          {dat = ImageData[img],
           opts = Options[img],
           stats = MRImageStatistics[img]},
          With[
            {rv = RightDirectionVector /. opts,
             av = AnteriorDirectionVector /. opts,
             sv = SuperiorDirectionVector /. opts,
             c = Center /. opts,
             vdims = VoxelDimensions /. opts},
            Which[
              !ArrayQ[dat, 3|4, NumericQ], Message[
                MRImage3D::badarg,
                "ImageData must be a 3 or 4D numeric array"],
              rv =!= Indeterminate && (!ArrayQ[rv, 1, NumericQ] || Length[rv] != 3), Message[
                MRImage3D::badarg,
                "RightDirectionVector must be Indeterminate or a 3-element numeric vector"],
              av =!= Indeterminate && (!ArrayQ[av, 1, NumericQ] || Length[av] != 3), Message[
                MRImage3D::badarg,
                "AnteriorDirectionVector must be Indeterminate or a 3-element numeric vector"],
              sv =!= Indeterminate && (!ArrayQ[sv, 1, NumericQ] || Length[sv] != 3), Message[
                MRImage3D::badarg,
                "SuperiorDirectionVector must be Indeterminate or a 3-element numeric vector"],
              c =!= Automatic && (!ArrayQ[c, 1, NumericQ] || Length[c] != 3), Message[
                MRImage3D::badarg,
                "Center must be Automatic or a 3-element numeric vector"],
              !ArrayQ[vdims, 1, NumericQ] || Length[vdims] != 3, Message[
                MRImage3D::badarg,
                "VodelDimensions must be Indeterminate or a 3-element numeric vector"],
              !NumericQ[stats[Min]] || !NumericQ[stats[Max]], Message[
                MRImage3D::badarg,
                "Image data may only contain numeric, Indeterminate, or None values"],
              True, True]]],

        (* The image as it is represented internally *)
        Image3D[img] :> With[
          {stats = MRImageStatistics[img],
           opts = Options[img]},
          Image3D[
            Rescale[
              If[stats[Missing] == 0,
                ImageData[img], 
                ImageData[img] /. {
                  Indeterminate -> (If[KeyExistsQ[stats, #], stats[#], #]&@(Indeterminate /. opts)),
                  None -> (If[KeyExistsQ[stats, #], stats[#], #]& @ (None /. opts))}],
              {stats[Min], stats[Max]}],
            Sequence@@FilterRules[opts, Options[Image3D][[All,1]]]]],
        
        (* This produces slices, but the preferred method is to use MRISlices[], below. *)
        Image3DSlices[img, opts___] := With[
          {max = MRImageMax[img],
           min = MRImageMin[img]},
          Map[
            Function[
              MRImage[
                ImageData[#] * (max - min) + min,
                MRImage3D -> img,
                Sequence @@ FilterRules[
                  Options[img],
                  Cases[Options[MRImage][[All, 1]], Except[MRImage3D|PixelDimensions]]],
                PixelDimensions -> Rest[VoxelDimensions[img]]]],
            Image3DSlices[Image3D[img], opts]]]]]],
  SetSafe ->True,
  Symbol -> MRImage3D];


(* #MRImage ***************************************************************************************)
Options[MRImage] = Join[
  {MRImage3D -> None,
   PixelDimensions -> {1,1}},
  Replace[
    Options[Image],
    {(ColorFunction -> _) -> (ColorFunction -> "XRay"),
     (Boxed -> _) -> (Boxed -> False),
     (Background -> _) -> (Background -> Black)},
    {1}]];
(* For redefining options of MRImage *)
MRImage[img_Image, opts___Rule] := MRImage[
  ImageData[img],
  opts,
  Sequence@@Options[img]];
MRImage[img_MRImage, opts___Rule] := Clone[
  img,
  Options -> Join[{opts}, Options[img]]];
DefineImmutable[
  MRImage[data_, OptionsPattern[]] :> img,
  Evaluate[
    Join[
      ReplaceAll[$MRImageSharedMethods, {Image3D -> Image, MRImage3D -> MRImage}],
      Hold[
        (* Retreive (or edit) the raw image data *)
        ImageData[img] = data,
        (* Retreive (or edit) the raw image options *)
        Options[img] = Map[(# -> OptionValue[#])&, Options[MRImage][[All, 1]]],
        (* Access the scale factors of min and max for the data *)
        MRImageMax[img] -> With[
          {source = MRImage3D /. Options[img]},
          If[MRImageQ[source], MRImageMax[source], Max[Flatten@ImageData[img]]]],
        MRImageMin[img] -> With[
          {source = MRImage3D /. Options[img]},
          If[MRImageQ[source], MRImageMin[source], Min[Flatten@ImageData[img]]]],
        (* The specific options... *)
        MRImage3D[img] := With[
          {source = MRImage3D /. Options[img]},
          If[!MRImageQ[source], None, source]],
        PixelDimensions[img] := With[
          {vdims = PixelDimensions /. Options[img]},
          If[vdims === Automatic, Rest@VoxelDimensions[MRImage3D[img]], vdims]],
        
        (* The image as it is represented internally *)
        Image[img] -> With[
          {dat = ImageData[img],
           opts = Options[img]},
          Which[
            !ArrayQ[dat, 2|3, NumericQ],  Message[
              MRImage::badarg,
              "ImageData must be a 2 or 3D numeric array"],
            (# =!= None && !MRImageQ[#])&[MRImage3D /. opts], Message[
               MRImage::badarg,
               "MRImage3D option must provide an MRImage3D object or None"],
            (!ArrayQ[#,1,NumericQ] || Length[#] != 2)&[PixelDimensions /. opts], Message[
              MRImage::badarg,
              "PixelDimensions must be a 2D array"],
            True, With[
              {min = MRImageMin[img],
               max = MRImageMax[img]},
              Image[
                (dat - min) / Replace[(max - min), 0|0. -> 1],
                Sequence@@FilterRules[opts, Options[Image][[All,1]]]]]]],
        
        (* This is an MRImage... Checks of data quality should go here as well. *)
        MRImageSliceQ[img] -> True]]],
  SetSafe ->True,
  Symbol -> MRImage];

(* ImageApply uses mechanics that have to be overridden separately *)
Unprotect[MRImage3D, MRImage];
MRImage3D /: ImageApply[f_, img_MRImage3D] := MRImage3D[
    img,
    ImageData -> Map[f, ImageData[img], {3}]];
MRImage /: ImageApply[f_, img_MRImage] := MRImage[
  img,
  ImageData -> Map[f, ImageData[img], {2}]];
Protect[MRImage3D, MRImage];

(* #MakeBoxes allows us to display the image as a 3d image *)
MakeBoxes[img_MRImage3D, form_] := MakeBoxes[#, form]&[Image3D[img]];
MakeBoxes[img_MRImage, form_] := MakeBoxes[#, form]&[Image[img]];

(* We want to define the #MRImageObjectQ and related functions now... *)
MRImageQ[x_] := False;
MRImageSliceQ[x_] := False;
MRImageObjectQ[img_] := Or[MRImageQ[img], MRImageSliceQ[img]];
Protect[MRImageObjectQ, MRImageQ, MRImageSliceQ];

(* and some others *)
Protect[
  RightDirectionVector, AnteriorDirectionVector, SuperiorDirectionVector, 
  VoxelDimensions, PixelDimensions,
  MRImageStatistics, MRImageMin, MRImageMax,
  MRImageMean, MRImageMedian, MRImageVariance,
  MRImageCount, MRImageMissing];

(* #MRIORientMatrix *******************************************************************************)
MRIOrientMatrix[img_?MRImageQ, spec:{_,_,_}] := Catch @ With[
  {dims = ImageDimensions[img],
   mtx0 = {
     RightDirectionVector[img],
     AnteriorDirectionVector[img],
     SuperiorDirectionVector[img]},
   dir = Replace[
     spec,
     {Right|RH|Rule[LH|Left, RH|Right] :> {1,0,0},
      Left|LH|Rule[RH|Right, LH|Left] :> {-1,0,0},
      Anterior|Front|Rule[Posterior|Back, Anterior|Front] :> {0,1,0},
      Posterior|Back|Rule[Anterior|Front, Posterior|Back] :> {0,-1,0},
      Superior|Top|Rule[Inferior|Bottom, Superior|Top] :> {0,0,1},
      Inferior|Bottom|Rule[Superior|Top, Inferior|Bottom] :> {0,0,-1},
      x_ :> Throw[
        Message[MRIOrient::badarg, "Unrecognized spec: " <> ToString[x]];
        $Failed]},
     {1}],
   center = 0.5*(ImageDimensions[img] + 1)},
  If[!MatrixQ[mtx0, NumericQ],
    Message[MRIOrient::badarg, "MRImage object does not include orientation information"]];
  If[Abs[Total[dir]] != {1,1,1},
    Message[MRIOrient::badarg, "repeated or opposite specs given"]];
  (* mtx0 tells us which direction right/left etc are in... get the transform and the center out *)
  With[
    {tx = mtx0 . dir,
     newDims = {Floor@Min[#], Ceiling@Max[#]}& /@ Transpose[
       Dot[mtx0, dir, #]& /@ Tuples@Tuples[{{1}, dims}]],
     newCenter = Dot[mtx0, dir, center]},
    With[
      {endCenter = (#[[2]] - #[[1]] + 2)/2.0& /@ newDims,
       endDims = {1, #[[2]] - #[[1]] + 1}& /@ newDims},
      (* We don't want to change the center; just the direction transform *)
      Append[MapThread[Append, {tx, endCenter - newCenter}], {0,0,0,1}]]]];
MRIOrientMatrix[img_?MRImageQ, spec_] := MRIOrientMatrix[img, spec, ImageDimensions@img];
Protect[MRIOrientMatrix];

(* #MRIORientMatrix *******************************************************************************)
MRIOrientTransform[img_?MRImageQ, spec:{_,_,_}, dims:{_,_,_}] := With[
  {mtx = MRIOrientMatrix[img, spec, dims]},
  If[ListQ[mtx], AffineTransform[{mtx[[1;;3, 1;;3]], mtx[[4, 1;;3]]}], mtx]];
MRIOrientTransform[img_?MRImageQ, spec_] := MRIOrientTransform[img, spec, ImageDimensions@img];
Protect[MRIOrientTransform];

(* Here, we want to do all the coordinate transformations... *)

(* #$VoxelIndexToRASMatrix *****************************************************************)
$VoxelIndexToRASMatrix = {
  { 0,  0,  1},
  { 0, -1,  0},
  {-1,  0,  0}};
Protect[$VoxelIndexToRASMatrix];

(* #VoxelIndexToRASMatrix ******************************************************************)
VoxelIndexToRASMatrix[img_?MRImageQ, spec_List] := Check[
  With[
    {mtx0 = MRIOrientMatrix[img, spec, ImageDimensions[img]]},
    If[!ListQ[mtx0],
      mtx0,
      (* This matrix needs to be modified only in that it needs to center things at 0 *)
      With[
        {mtx = mtx0[[1;;3, 1;;3]]},
        Append[
          MapThread[
            Append,
            {mtx, ImageDimensions[img]/2}],
          {0,0,0,1}]]]],
  $Failed];
Protect[VoxelIndexToRASMatrix];
(* #VoxelIndexToRASTransform ***************************************************************)
VoxelIndexToRASTransform[img_?MRImageQ, spec_List] := Check[
  With[
    {mtx = VoxelIndexToRASMatrix[img, spec]},
    If[!ListQ[mtx],
      mtx,
      AffineTransform[{mtx[[1;;3,1;;3]], mtx[[4,1;;3]]}]]],
  $Failed];
Protect[VoxelIndexToRASTransform];  

(* #VoxelIndexToRASTr **********************************************************************)
VoxelIndexToRASTr[
  img_?MRImageQ,
  coords_ /; ArrayQ[coords,2] && Length[coords] == 3,
  spec:{_,_,_}
 ] := Check[
  With[
    {mtx = Most@VoxelIndexToRASMatrix[img, spec],
     x0 = ImageDimensions[img] / 2},
    (*If[!ListQ[mtx], mtx, Dot[mtx, Append[coords, ConstantArray[1, Length[coords[[1]]]]]]]],*)
    If[!ListQ[mtx], mtx, Dot[mtx[[1;;3, 1;;3]], (coords - x0)]]],
  $Failed];
VoxelIndexToRASTr[img_?MRImageQ, coords:{_List,_List,_List}] := VoxelIndexToRASTr[
  img,
  coords,
  {Right, Anterior, Superior}];
Protect[VoxelIndexToRASTr];

(* #VoxelIndexToRAS ************************************************************************)
VoxelIndexToRAS[img_?MRImageQ, coords_ /; ArrayQ[coords,2], spec:{_,_,_}] := With[
  {tmp = VoxelIndexToRASTr[img, Transpose @ coords, spec]},
  If[ListQ[tmp], Transpose @ tmp, tmp]];
VoxelIndexToRAS[img_?MRImageQ, coord:{_,_,_}, spec:{_,_,_}] := With[
  {tmp = VoxelIndexToRASTr[img, List /@ coord, spec]},
  If[ListQ[tmp], First[tmp], tmp]];
VoxelIndexToRAS[img_?MRImageQ, coords_] := VoxelIndexToRAS[
  img,
  coords,
  {Right, Anterior, Superior}];
Protect[VoxelIndexToRAS];


(* #RASToVoxelIndexMatrix ******************************************************************)
RASToVoxelIndexMatrix[img_?MRImageQ, spec_List] := With[
  {mtx = VoxelIndexToRASMatrix[img, spec]},
  If[!ListQ[mtx], mtx, Inverse[mtx]]];
Protect[RASToVoxelIndexMatrix];

(* #RASToVoxelIndexTransform ***************************************************************)
RASToVoxelIndexTransform[img_?MRImageQ, spec_List] := With[
  {mtx = VoxelIndexToRASTransform[img, spec]},
  If[mtx === Indeterminate || mtx === $Failed, mtx, InverseFunction[mtx]]];
Protect[RASToVoxelIndexTransform];  

(* #RASToVoxelIndexTr **********************************************************************)
RASToVoxelIndexTr[
  img_?MRImageQ,
  coords_ /; ArrayQ[coords,2] && Length[coords] == 3,
  spec:{_,_,_}
 ] := Check[
  With[
    {mtx = RASToVoxelIndexMatrix[img, spec]},
    If[!ListQ[mtx], mtx, Dot[mtx, Append[coords, ConstantArray[1, Length[coords[[1]]]]]]]],
  $Failed];
RASToVoxelIndexTr[img_?MRImageQ, coords:{_List,_List,_List}] := RASToVoxelIndexTr[
  img,
  coords,
  {Right, Anterior, Superior}];
Protect[RASToVoxelIndexTr];

(* #RASToVoxelIndex ************************************************************************)
RASToVoxelIndex[img_?MRImageQ, coords_ /; ArrayQ[coords,2], spec:{_,_,_}] := With[
  {tmp = RASToVoxelIndexTr[img, Transpose @ coords, spec]},
  If[ListQ[tmp], Transpose @ tmp, tmp]];
RASToVoxelIndex[img_?MRImageQ, coord:{_,_,_}, spec:{_,_,_}] := With[
  {tmp = RASToVoxelIndexTr[img, List /@ coord, spec]},
  If[ListQ[tmp], First[tmp], tmp]];
RASToVoxelIndex[img_?MRImageQ, coords_] := RASToVoxelIndex[
  img,
  coords,
  {Right, Anterior, Superior}];
Protect[RASToVoxelIndex];

(* #MRIInterpolate ********************************************************************************)
MRIInterpolate[img_?MRImageQ,
               pt_ /; VectorQ[pt, NumericQ] && Length[pt] == 3] := First@MRIInterpolate[img, {pt}];
MRIInterpolate[img_?MRImageQ, pts_ /; MatrixQ[pts, NumericQ] && Length@First[pts] == 3] := With[
  {voxelDims = VoxelDimensions[img],
   dat = ImageData[img]},
  With[
    {X = MapThread[
       Clip[#1, {1,#2}]&,
       {If[N[voxelDims] == {1.0,1.0,1.0}, #, #/voxelDims]&@Transpose[pts],
        ImageDimensions[img]}]},
    With[
      {fl = Floor[X], cl = Ceiling[X]},
      With[
        {dc = cl - X, df = X - fl},
      With[
        {volumes = Map[(Times@@#)&, Tuples@Transpose[{df, dc}]],
         values  = Map[Extract[data, #]&, Tuples@Transpose[{cl, fl}]]},
        Total[volumes * values]]]]]];      

(* #MRITransformation *****************************************************************************)
MRITransformation[img_?MRImageQ, tx0_, opts___Rule] := Check[
  With[
    {txs = Switch[tx0,
       _?MatrixQ, Switch[Dimensions[tx0],
         {3,3}, {AffineTransform[tx0], AffineTransform@Inverse[tx0]},
         {3,4}, With[
            {m = Inverse@Append[tx0, {0,0,0,1}]},
            {AffineTransform[{tx0[[1;;3, 1;;3]], tx0[[1;;3, 4]]}],
             AffineTransform[{  m[[1;;3, 1;;3]],   m[[1;;3, 4]]}]}],
         {4,4}, With[
            {m = Inverse[tx0]},
            {AffineTransform[{tx0[[1;;3, 1;;3]], tx0[[1;;3, 4]]}],
             AffineTransform[{  m[[1;;3, 1;;3]],   m[[1;;3, 4]]}]}],
         _, Message[MRITransformation::badtx]],
       {_?MatrixQ, _List}, {AffineTransform[tx0], InverseFunction@AffineTransform[tx0]},
       _, {tx0, InverseFunction[tx0]}],
     dims = ImageDimensions[img]},
    With[
      {newDims = Ceiling@Max[#]& /@ Transpose[tx@Tuples@Tuples[{{1}, dims}]],
       tx = txs[[1]],
       itx = txs[[2]],
       rv = RightDirectionVector[img],
       av = AnteriorDirectionVector[img],
       sv = SuperiorDirectionVector[img],
       c = Center[img],
       mid = 0.5*(dims + 1),
       origin = {0,0,0}},
      MRImage3D[
        ImageTransformation[img, itx, DataRange->Thread[{1, dims}]],
        RightDirectionVector    -> If[ListQ[rv], tx[rv] - tx[origin], rv],
        AnteriorDirectionVector -> If[ListQ[av], tx[av] - tx[origin], av],
        SuperiorDirectionVector -> If[ListQ[sv], tx[sv] - tx[origin], sv],
        Center                  -> If[ListQ[c],  tx[c + mid] - tx[mid], c]]]],
  $Failed];
Protect[MRITransformation];

(* #MRIOrient *************************************************************************************)
MRIOrient[img_?MRImageQ, spec:{_,_,_}, dims:{_,_,_}] := Check[
  MRITransformation[img, MRIOrientTransform[img, spec, dims]],
  $Failed];
MRIOrient[img_?MRImageQ, spec:{_,_,_}] := MRIOrient[img, spec, ImageDimensions[img]];
MRIOrient[img_?MRImageQ, type_String /; StringLength[type] == 3] := Check[
  MRIOrient[
    img,
    Replace[
      Characters @ ToUpperCase[type],
      {"R" -> Right,    "L" -> Left,
       "A" -> Anterior, "P" -> Posterior,
       "S" -> Superior, "I" -> Inferior,
       c_ :> Message[MRIOrient::badarg, "Unrecognized orientation character: " <> c]},
      {1}]],
    $Failed];
MRIOrient[img_?MRImageQ, s_String] := MRIOrient[img, s, ImageDimensions[img]];
Protect[MRIOrient];

(* #MRISlices *************************************************************************************)
Options[MRISlices] = {MRIOrient -> Automatic};
MRISlices[img_?MRImageQ, plane_, which_, OptionsPattern[]] := Check[
  With[
    {orient = OptionValue[MRIOrient],
     dir = Switch[plane,
       "Sagittal"|"Lateral", Left,
       Right|RH|(LH|Left -> RH|Right), Left,
       Left|LH|(RH|Right -> LH|Left), Right,
       "Coronal"|"Frontal", Posterior,
       Anterior|Front|(Back|Posterior -> Front|Anterior), Posterior,
       Posterior|Back|(Front|Anterior -> Back|Posterior), Anterior,
       "Transverse"|"Axial"|"Horizontal", Inferior,
       Superior|Top|(Bottom|Inferior -> Top|Superior), Inferior,
       Inferior|Bottom|(Top|Superior -> Bottom|Inferior), Superior,
       _, Message[MRISlices::badplane, plane]]},
    With[
      {tx = MRIOrient[
         img,
         Append[
           Which[
             orient === Automatic, Switch[dir,
               Right, {Anterior, Superior},
               Left, {Posterior, Superior},
               Superior, {Right, Anterior},
               Inferior, {Right, Anterior},
               Anterior, {Left, Superior},
               Posterior, {Right, Superior}],
             ListQ[orient] && Length[orient] == 2, orient,
             True, Message[MRISlices::badarg, "MRIOrient argument must be a 2-element list"]],
           dir]]},
      Image3DSlices[tx, which]]],
  $Failed];
MRISlices[img_?MRImageQ, plane_, opts:OptionsPattern[]] := MRISlices[img, plane, All, opts];
Protect[MRISlices];



(* For all of the functions below, we want to write them only once, for the 3D case, but to define
 * them twice, also for the 2D case. Their code is almost identical, however, so a hack is employed
 * to prevent code duplication.
 *)
(ReleaseHold[#]; ReleaseHold[# //. {Image3D -> Image, MRImage3D -> MRImage}];)& @ Hold[

   (* #ImageAlign is a slightly special case *)
   Unprotect[ImageAlign];
   ImageAlign[a_, b_MRImage3D, opts___] := With[
     {res = Check[ImageAlign[a, Image3D[b], opts], $Failed]},
     If[res === $Failed, 
       res,
       With[
         {res = Check[ImageCrop[Image3D[img], opts], $Failed]},
         If[res === $Failed, 
           res,
           With[
             {newdat = ImageData[res] * (MRImageMax[b] - MRImageMin[b]) + MRImageMin[b],
              theOpts = Options[res]},
             MRImage3D[newdat, Evaluate[Sequence@@theOpts]]]]]]];
   Protect[ImageAlign];

   (* #ImageCorrespondingPoints also is a special case *)
   Unprotect[ImageCorrespondingPoints];
   ImageCorrespondingPoints[a_, b_MRImage3D, opts___] := ImageCorrespondingPoints[
     a,
     Image3D[b],
     opts];
   Protect[ImageCorrespondingPoints];
   
   (* #ImageDistance is also a special case *)
   Unprotect[ImageDistance];
   ImageDistance[a_, b_MRImage3D, opts___] := ImageDistance[a, Image3D[b], opts];
   Protect[ImageDistance];

   (* As is #ImageFilter *)
   Unprotect[ImageFilter];
   (* I have no idea why this is necessary, but it seems ImageFilter reprotects itself *)
   ImageFilter;
   Unprotect[ImageFilter];
   ImageFilter[f_, img_MRImage3D, opts___] := With[
     {res = Check[ImageFilter[f, Image3D[img], opts], $Failed]},
     If[res === $Failed, 
       res,
       With[
         {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
          theOpts = Options[res]},
         MRImage3D[newdat, Evaluate[Sequence@@theOpts]]]]];
   Protect[ImageFilter];];

End[];
EndPackage[];

