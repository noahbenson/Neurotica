(* NifTI.m
 *
 * The Neurotica`NifTI namespace contains functions for reading and interpreting the NIH standard
 * neuroimaging file formats, including NifTI, GifTI, and SifTI.
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
BeginPackage[
  "Neurotica`NifTI`",
  {"Neurotica`Global`","Neurotica`Util`","Neurotica`Mesh`", "Neurotica`MRImage`"}];
Unprotect["Neurotica`NifTI`*", "Neurotica`NifTI`Private`*"];
ClearAll[ "Neurotica`NifTI`*", "Neurotica`NifTI`Private`*"];

ImportNifTI::usage = "ImportNifTI[source, options...] is equivalent to Import[source, \"NIFTI\", optios...].";
ImportNifTI::badarg = "Bad argument given to ImportNifTI: `1`";
ImportNifTI::badfmt = "Bad NifTI file format: `1`";

ReadBinaryStructure::usage = "ReadBinaryStructure[stream, instructionsList] yields the data structure result of importing the given instructions via BinaryReadList. Instructions may take the following format:
  * type (e.g., \"Integer32\" or \"Real64\") are read as is;
  * {type, n} reads in a list of n of the given type (if n is 1 then a 1-element list is returned);
  * {type, n, postprocessFn} yields postprocessFn[data] where data is the list of n of the given type;
  * name -> spec yields the rule name -> result where result is the result of reading in the given instruction spec;
  * {spec1, spec2...} yields a list of the given instruction specs, each read in order.";
ReadBinaryStructure::badinstr = "Bad instruction given to ReadBinaryStructure: `1`";

Begin["`Private`"];

(* #ReadBinaryStructure ***************************************************************************)
$ReadBinaryStructureInstructionKey = {
  type_String :> Function[BinaryRead[#, type]],
  {type_String} :> Function[BinaryReadList[#, type, 1]],
  {type_String, n_Integer} :> Function[BinaryReadList[#, type, n]],
  {type_String, n_Integer, f_} :> Function[f[BinaryReadList[#, type, n]]],
  l_List :> With[
    {els = Replace[l, $ReadBinaryStructureInstructionKey, 1]},
    Function[With[{stream = #}, #[stream]& /@ els]]],
  Rule[name_, instr_] :> With[
    {r = (instr /. $ReadBinaryStructureInstructionKey)}, 
    Function[name -> r[#]]],
  q_ :> Message[ReadBinaryStructure::badinstr, ToString[q]]};
Protect[$ReadBinaryStructureInstructionKey];
  
ReadBinaryStructure[stream_, instructions_List] := Check[
  With[
    {instFns = Replace[instructions, $ReadBinaryStructureInstructionKey, {1}]},
    Map[
      #[stream]&,
      instFns]],
  $Failed];
Protect[ReadBinaryStructure];

(* ============================================================================================== *)
(* Here we include the NifTI file spec header:
   int   sizeof_hdr;    /*!< MUST be 348           */  /* int sizeof_hdr;      */
   char  data_type[10]; /*!< ++UNUSED++            */  /* char data_type[10];  */
   char  db_name[18];   /*!< ++UNUSED++            */  /* char db_name[18];    */
   int   extents;       /*!< ++UNUSED++            */  /* int extents;         */
   short session_error; /*!< ++UNUSED++            */  /* short session_error; */
   char  regular;       /*!< ++UNUSED++            */  /* char regular;        */
   char  dim_info;      /*!< MRI slice ordering.   */  /* char hkey_un0;       */
 
   short dim[8];        /*!< Data array dimensions.*/  /* short dim[8];        */
   float intent_p1 ;    /*!< 1st intent parameter. */  /* short unused8;       */
                                                       /* short unused9;       */
   float intent_p2 ;    /*!< 2nd intent parameter. */  /* short unused10;      */
                                                       /* short unused11;      */
   float intent_p3 ;    /*!< 3rd intent parameter. */  /* short unused12;      */
                                                       /* short unused13;      */
   short intent_code ;  /*!< NIFTI_INTENT_* code.  */  /* short unused14;      */
   short datatype;      /*!< Defines data type!    */  /* short datatype;      */
   short bitpix;        /*!< Number bits/voxel.    */  /* short bitpix;        */
   short slice_start;   /*!< First slice index.    */  /* short dim_un0;       */
   float pixdim[8];     /*!< Grid spacings.        */  /* float pixdim[8];     */
   float vox_offset;    /*!< Offset into .nii file */  /* float vox_offset;    */
   float scl_slope ;    /*!< Data scaling: slope.  */  /* float funused1;      */
   float scl_inter ;    /*!< Data scaling: offset. */  /* float funused2;      */
   short slice_end;     /*!< Last slice index.     */  /* float funused3;      */
   char  slice_code ;   /*!< Slice timing order.   */
   char  xyzt_units ;   /*!< Units of pixdim[1..4] */
   float cal_max;       /*!< Max display intensity */  /* float cal_max;       */
   float cal_min;       /*!< Min display intensity */  /* float cal_min;       */
   float slice_duration;/*!< Time for 1 slice.     */  /* float compressed;    */
   float toffset;       /*!< Time axis shift.      */  /* float verified;      */
   int   glmax;         /*!< ++UNUSED++            */  /* int glmax;           */
   int   glmin;         /*!< ++UNUSED++            */  /* int glmin;           */
  
   char  descrip[80];   /*!< any text you like.    */  /* char descrip[80];    */
   char  aux_file[24];  /*!< auxiliary filename.   */  /* char aux_file[24];   */
  
   short qform_code ;   /*!< NIFTI_XFORM_* code.   */  /*-- all ANALYZE 7.5 ---*/
   short sform_code ;   /*!< NIFTI_XFORM_* code.   */  /*   fields below here  */
                                                       /*   are replaced       */
   float quatern_b ;    /*!< Quaternion b param.   */
   float quatern_c ;    /*!< Quaternion c param.   */
   float quatern_d ;    /*!< Quaternion d param.   */
   float qoffset_x ;    /*!< Quaternion x shift.   */
   float qoffset_y ;    /*!< Quaternion y shift.   */
   float qoffset_z ;    /*!< Quaternion z shift.   */
  
   float srow_x[4] ;    /*!< 1st row affine transform.   */
   float srow_y[4] ;    /*!< 2nd row affine transform.   */
   float srow_z[4] ;    /*!< 3rd row affine transform.   */

   char intent_name[16];/*!< 'name' or meaning of data.  */

   char magic[4] ;      /*!< MUST be "ni1\0" or "n+1\0". */                                       *)
(* ============================================================================================== *)

$NIFTIHeaderSize = 348;

(* #ImportNifTIHeader *****************************************************************************)
ImportNifTIHeader[stream_, opts___Rule] := Check[
  SetStreamPosition[stream, 0];
  With[
    {raw = ReadBinaryStructure[
       stream,
       {"HeaderSize" -> {"Integer32", 1, Function[
          If[#[[1]] != $NIFTIHeaderSize, 
            Message[ImportNifTI::badfmt, "Header does not start with 348"],
            #[[1]]]]},
        "DataType" -> {"Integer8", 10},
        "DBName" -> {"Character8", 18, StringJoin},
        "Extents" -> "Integer32",
        "SessionError" -> "Integer16",
        "Regular" -> "Integer8",
        "DimensionsInformation" -> "Integer8",
        "Dimensions" -> {"Integer16", 8},
        "IntentParameters" -> {"Real32", 3},
        "IntentCode" -> "Integer16",
        "Datatype" -> "Integer16",
        "BitsPerVoxel" -> "Integer16",
        "SliceStart" -> "Integer16",
        "GridSpacings" -> {"Real32", 8},
        "VoxelOffset" -> "Real32",
        "ScaleSlope" -> "Real32",
        "ScaleIntercept" -> "Real32",
        "SliceEnd" -> "Integer16",
        "SliceCode" -> "Integer8",
        "DimensionsUnits" -> "Integer8",
        "DisplayIntensityRange" -> {"Real32", 2},
        "SliceDuration" -> "Real32",
        "TimeAxisShift" -> "Real32",
        "GLRange" -> {"Integer32", 2, Reverse},
        "Description" -> {"Character8", 80, StringJoin[TakeWhile[#, ToCharacterCode[#] != {0}&]]&},
        "AuxiliaryFilename" -> {"Character8", 24, 
                                StringJoin[TakeWhile[#, ToCharacterCode[#] != {0}&]]&},
        "QFormCode" -> "Integer16",
        "SFormCode" -> "Integer16",
        "Quaternions" -> {"Real32", 6},
        "AffineTransform" -> {"Real32", 12, Partition[#, 4]&},
        "IntentName" -> {"Character8", 16, StringJoin[TakeWhile[#, ToCharacterCode[#] != {0}&]]&},
        "MagicTerminus" -> {
          "Character8", 4,
          Function[
            With[
              {str = StringJoin@TakeWhile[#, ToCharacterCode[#] != {0}&]},
              If[str == "nil" || str == "n+1", 
                str,
                Message[ImportNifTI::badfmt, "header does not end in nil or n+1 string"]]]]}}]},
    raw],
  $Failed];
    
End[];
EndPackage[];

