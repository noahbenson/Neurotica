(* Util.m
 *
 * Basic utility functions, e.g. for automatically caching data in quick-to-load mx files.
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
BeginPackage["Neurotica`Util`", {"Neurotica`Global`"}];
Unprotect["Neurotica`Util`*", "Neurotica`Util`Private`*"];
ClearAll ["Neurotica`Util`*", "Neurotica`Util`Private`*"];

$CacheDirectory::usage = "$CacheDirectory is the directory in which cached data for this user is placed. If the directory does not exist at cache time and $AutoCreateCacheDirectory is True, then this directory is automatically created. By default it is the directory FileNameJoin[{$UserBaseDirectory, \"AutoCache\"}]. If set to Temporary, then the next use of AutoCache will create a temporary directory and set the $CacheDirectory to this path.";

$AutoCreateCacheDirectory::usage = "$AutoCreateCacheDirectory is True if and only if the $CacheDirectory should be automatically created when AutoCache requires it.";

RefreshOnChange::usage = "RefreshOnChange is an option to AutoCache that, if true, instructs AutoCache that it should refresh the cache if the cache file is older than the last time of modification of the evaluating cell. If false, then the cache is never refreshed due to modification of the cell.";

AutoCache::usage = "AutoCache[name, body] checks first to see if a cache file exist for the file named by the given name and yields its contents if so. Otherwise, yields the result of evaluating body and caches it in the given filename. If the cache file is out of date relative to the evaluating cell, the contents are erased and recalculated. If the given filename is an absolute path, it is used as such, otherwise it is localized to the cache directory. Note that the mx extension is added automatically if not included.
The following options may be used:
 * RefreshOnChange (default: True) prevents the function from overwriting the cache file when the cell is changed if False.
 * Quiet (default: False) prevents the function from producing messages when updating or overwriting the cache file if False.
 * Check (default: True) if True instructs AutoCache to yield $Failed and not generate the cache file whenever messages are produced during the execution of body.
 * Directory (default: Automatic) names the directory in which the cache file should go; if Automatic, then uses $CacheDirectory.
 * CreateDirectory (default: Automatic) determines whether AutoCache can automatically create the $CacheDirectory if it is needed. If True, then AutoCache will always create the directory; if False it will never create the directory; if Automatically, it will defer to $AutoCreateCacheDirectory.";
AutoCache::expired = "AutoCache file `1` has expired; recalculating";
AutoCache::nodir = "AutoCache directory (`1`) does not exist and cannot be created";
AutoCache::nomkdir = "AutoCache directory (`1`) does not exist and AutoCache is not allowed to create it";
AutoCache::badopt = "Bad option to AutoCache: `1`"

SetSafe::usage = "SetSafe[a, b] is equivalent t0 the expression (a = b) except that if any messages are generated during the evaluation of b, $Failed is yielded and a is not set.";

Memoize::usage = "Memoize[f := expr] constructs a version of the pattern f that will auto-memoize its values and remember them.
Memoize[\!\(\*SubscriptBox[\(f\),\(1\)]\) := \!\(\*SubscriptBox[\(expr\),\(1\)]\), \!\(\*SubscriptBox[\(f\),\(2\)]\) := \!\(\*SubscriptBox[\(expr\),\(2\)]\), ...] memoizes each of the given functions.
The function Forget may be used to forget the memoized values of functions defined within a Memoize form.";
MemoizeSafe::usage = "MemoizeSafe[f := expr, ...] is identical to calling Memoize except that it uses SetSafe instead of Set to memoize values, so of the evaluation of the expression fails or generates a message, the result is not memoized and $Failed is yielded.";
Forget::usage = "Forget[f] forces all memoized values associated with the symbol f to be forgotten.
Forget[\!\(\*SubscriptBox[\(f\),\(1\)]\), \!\(\*SubscriptBox[\(f\),\(2\)]\), ...] forgets each of the symbols.";

Indices::usage = "Indices[list, patt] yields a list of the indices of elements that match patt in list; this is equivalent to Flatten[Position[list, patt, {1}, Heads->False]], but is considerably optimized.";

TemporarySymbol::usage = "TemporarySymbol[] yields a unique symbol that is marked as temporary.
TemporarySymbol[string] yields a unique symbol that begins with the given string and is marked temporary.";

DefineImmutable::usage = "DefineImmutable[constructor[args...] :> symbol, {members...}] defines an immutable data type with the given constructor and the accessor functions given by the list of members. Each member must be of the form (f[symbol] ** form) where f is the function that is used to access the data, the form yields the value and may refer to the arguments of the constructor or the other member functions, and ** must be either =, ->, :>, or :=, each of which is interpreted differently:
 * any value assigned with a := will become a generic function whose value is calculated each time the function is called;
 * any value assigned with a :> rule becomes a function that is memoized for each object once it has been calculated; note, however, that memoized values in an immutable object cannot be forgotted by the Forget[] function;
 * any value assigned with a -> rule will be calculated and saved to the object when it is constructed;
 * any value assigned with a = will also be assigned when the object is constructed, just as with the -> rule forms; however, values assigned with an = form may be modified in the Clone[] function.
In addition, the functions assigned using the := form may have additional arguments, as long as the first argument is the symbol.
DefineImmutable[constructor[args...] :> symbol, Hold[members...]] is identical to DefineImmutable[constructor[args...] :> symbol, {members...}].
The following options may be given:
 * Symbol (default: Automatic) specifies the box form used to hold the type; if no box type is given, then a unique symbol is generated to encapsulate the type.
 * SetSafe (default: True) specifies whether memoization performed on :> rule forms should be done using SetSafe (True) or Set (False).
See also Clone.";
DefineImmutable::badarg = "Bad argument given to DefineImmutable: `1`";
Clone::usage = "Clone[immutable, edits...] yields an immutable object (created via a constructor defined by the DefineImmutable interface) identical to the given object immutable with the changes given by the sequence of edits. Each edit must take the form accessor -> value where accessor is the function name for a member assigned using an = form in the DefineImmutable function and value is the new value it should take.";
Clone::badarg = "Bad argument given to Clone: `1`";

MimicAssociation::usage = "MimicAssociation[...] is identical to Association[...] if the Mathematica version is at least 10.0; otherwise, it yields a symbol that imitates an association for most basic intents and purposes.";
MimicAssociation::badarg = "Bad argument given to MimicAssociation: `1`";

NormalizeRows::usage = "NormalizeRows[X] yields a transformation of the matrix X in which each row of X has been normalized; this is the equivalent of (Normalize /@ X) but is significantly optimized.";
NormalizeColumns::usage = "NormalizeColumns[X] yields a transformation of the matrix X in which each column of X has been normalized. This is equivalent to Transpose[Normalze /@ Transpose[X]], but has been significantly optimized.";
RowNorms::usage = "RowNorms[X] yields the equivalent of Norm /@ X but has been optimized for speed.";
ColumnNorms::usage = "ColumnNorms[X] yields the equivalent of Norm /@ Transpose[X] but has been optimized for speed.";

QuaternionToRotationMatrix::usage = "QuaternionToRotationMatrix[{x,y,z,w}] yields the rotation matrix associated with the quaternion {x,y,z,w}.
QuaternionToRotationMatrix[{y, z, q}] is equivalent to QuaternionToRotationMatrix[{x,y,z,q}] where x is Sqrt[1 - (x^2 + y^2 + z^2)].";

ReadBinaryStructure::usage = "ReadBinaryStructure[stream, instructionsList] yields the data structure result of importing the given instructions via BinaryReadList. Instructions may take the following format:
  * type (e.g., \"Integer32\" or \"Real64\") are read as is;
  * {type, n} reads in a list of n of the given type (if n is 1 then a 1-element list is returned);
  * {type, n, postprocessFn} yields postprocessFn[data] where data is the list of n of the given type;
  * name -> spec yields the rule name -> result where result is the result of reading in the given instruction spec;
  * {spec1, spec2...} yields a list of the given instruction specs, each read in order.";
ReadBinaryStructure::badinstr = "Bad instruction given to ReadBinaryStructure: `1`";

BinaryStringFix::usage = "BinaryStringFix[{c1, c2, ...}] yields the string formed from the given list of characters, truncating the string upon encountering the first null character (character code 0).";

NeuroticaPermanentDatum::usage = "NeuroticaPermanentDatum[name, value] yields Null but adds the given name/value pair to the current user's Neurotica auto-initialization file. The data contained in this file are loaded at startup and must consist only of string data. The file itself is stored as a \"Table\" style file at FileNameJoin[{$UserBasePath, \".NeuroticaPermanentData.dat\"}].
If the value given in the second argument is None, then the datum is removed from the permanent cache.
NeuroticaPermanentDatum[name] yields the value of the permanent datum with the given name; if no such datum is found, then None is returned.";

$NeuroticaPermanentData::usage = "$NeuroticaPermanentData is an Association of the permanent Neurotica data, as saved using the NeuroticaPermanentDatum function.";

Begin["`Private`"];

(* #$CacheDirectory *******************************************************************************)
$CacheDirectory = FileNameJoin[{$UserBaseDirectory, "AutoCache"}];

(* #$AutoCreateCacheDirectory *********************************************************************)
$AutoCreateCacheDirectory = True;
$AutoCreateCacheDirectory /: Set[$AutoCreateCacheDirectory, b:(True|False)] := (
  Unprotect[$AutoCreateCacheDirectory];
  $AutoCreateCacheDirectory = b;
  Protect[$AutoCreateCacheDirectory];
  b);
Protect[$AutoCreateCacheDirectory];

(* #RefreshOnChange *******************************************************************************)
Protect[RefreshOnChange];

(* #AutoCache *************************************************************************************)
Options[AutoCache] = {
  RefreshOnChange -> True,
  Quiet -> False,
  Check -> True,
  Directory -> Automatic,
  CreateDirectory -> Automatic};
Attributes[AutoCache] = {HoldRest};
AutoCache[name_String, body_, OptionsPattern[]] := Catch[
  With[
    {splitName = FileNameSplit[name],
     canCreate = Replace[
       OptionValue[CreateDirectory],
       {Automatic :> $AutoCreateCacheDirectory,
        Except[True|False] :> (
          Message[AutoCache::badopt, "CreateDirectory must be True, False, or Automatic"];
          Throw[$Failed])}],
     quiet = TrueQ[OptionValue[Quiet]],
     check = TrueQ[OptionValue[Check]],
     refresh = TrueQ[OptionValue[RefreshOnChange]]},
    With[
      {dir = Replace[
         OptionValue[Directory],
         {Automatic :> Which[
            First[splitName] == "", FileNameJoin[Most[splitName]],
            First[splitName] == "~", FileNameJoin[Most[splitName]],
            $CacheDirectory === Temporary, ($CacheDirectory = CreateDirectory[]),
            True, $CacheDirectory],
          Except[_String] :> (
            Message[AutoCache::badopt, "Directory must be Automatic or a string"];
            Throw[$Failed])}],
       file = StringReplace[Last[splitName], {s:(__ ~~ ".mx") :> s, s__ :> (s <> ".mx")}, 1]},
      With[
        {fileName = Catch[
           FileNameJoin[
             {Which[
                DirectoryQ[dir], dir,
                canCreate, Replace[
                  CreateDirectory[dir],
                  $Failed :> (Message[AutoCache::nodir, dir]; Throw[None])],
                True, (Message[AutoCache::nomkdir, dir]; Throw[None])],
              file}]]},
        With[
          {fileDate = Which[
             !refresh, 0,
             fileName === None, 0,
             !FileExistsQ[FileNameJoin[{dir, file}]], 0,
             True, AbsoluteTime[FileNameJoin[{dir, file}], TimeZone -> 0]],
           cellDate = If[refresh,
             Last@Replace[
               CellChangeTimes /. AbsoluteOptions[EvaluationCell[]],
               {t0_, t1_} :> t1,
               {1}],
             0]},
          With[
            {cachedRes = Which[
               refresh && fileDate == 0, None,
               refresh && cellDate >= fileDate, (Message[AutoCache::expired]; None),
               !refresh && !FileExistsQ[FileNameJoin[{dir, file}]], None,
               True, Block[{Global`data = None}, Get[fileName]; Global`data]]},
            With[
              {theRes = If[cachedRes === None || cachedRes === $Failed,
                  If[check,
                    Check[body, $Failed],
                    body],
                  cachedRes]},
              If[cachedRes === None || cachedRes == $Failed,
                Block[{Global`data = theRes}, DumpSave[fileName, Global`data]]];
              theRes]]]]]]];
Protect[AutoCache];

(* #SetSafe ***************************************************************************************)
SetAttributes[SetSafe, HoldAll];
SetSafe[f_, expr_] := With[
  {res = Check[expr, $Failed]},
  If[res === $Failed, $Failed, Set[f, res]]];
Protect[SetSafe];

(* #Memoize ***************************************************************************************)
$MemoizePatternName = $MemoizePatternName;
Protect[$MemoizePatternName];
SetAttributes[Memoize, HoldAll];
Memoize[forms:(SetDelayed[_,_]..)] := CompoundExpression@@Replace[
  Hold[forms],
  HoldPattern[SetDelayed[f_, expr_]] :> ReplacePart[
    Hold[
      Evaluate[Pattern[Evaluate[$MemoizePatternName], f]],
      Evaluate[Hold[Evaluate[$MemoizePatternName], expr]]],
    {{2,0} -> Set, 0 -> SetDelayed}],
  {1}];
Protect[Memoize];

(* #MemoizeSafe ***********************************************************************************)
SetAttributes[MemoizeSafe, HoldAll];
MemoizeSafe[forms:(SetDelayed[_,_]..)] := CompoundExpression@@Replace[
  Hold[forms],
  HoldPattern[SetDelayed[f_, expr_]] :> ReplacePart[
    Hold[
      Evaluate[Pattern[Evaluate[$MemoizePatternName], f]],
      Evaluate[Hold[Evaluate[$MemoizePatternName], expr]]],
    {{2,0} -> SetSafe, 0 -> SetDelayed}],
  {1}];
Protect[MemoizeSafe];

(* #Forget ****************************************************************************************)
SetAttributes[Forget, HoldAll];
Forget[symbol_] := With[
  {sym1 = Unique[],
   sym2 = Unique[],
   down = DownValues[symbol]},
  With[
    {patterns = Cases[
       down,
       ReplacePart[
         Hold[
           Evaluate[
             Pattern @@ {
               sym1,
               Hold[
                 Evaluate[HoldPattern@@{Pattern[Evaluate[$MemoizePatternName], _]}],
                 Evaluate[Hold[Evaluate[$MemoizePatternName], _]]]}],
           Evaluate[Hold[Evaluate[{sym1}], Evaluate[{1,1} -> sym2]]]],
         {{1,2,2,0} -> (Set|SetSafe),
          {1,2,0} -> RuleDelayed,
          {2,1,0} -> First,
          {2,0} -> ReplacePart,
          0 -> RuleDelayed}],
       {1}]},
    With[
      {memoized = Cases[
         down[[All, 1]], 
         Evaluate[(Alternatives @@ patterns) :> Evaluate[Hold[Evaluate[sym2]]]],
         {2}]},
      Scan[Apply[Unset, #]&, memoized]]]];
Protect[Forget];      

(* #Indices ***************************************************************************************)
Indices[list_List, patt_] := Pick[Range[Length@list], list, patt];
Protect[Indices];

(* #TemporarySymbol *******************************************************************************)
TemporarySymbol[] := With[{sym = Unique[]}, SetAttributes[Evaluate[sym], Temporary]; sym];
TemporarySymbol[s_String] := With[{sym = Unique[s]}, SetAttributes[Evaluate[sym], Temporary]; sym];
Protect[TemporarySymbol];

(* #DefineImmutable *******************************************************************************)
SetAttributes[DefineImmutable, HoldAll];
Options[DefineImmutable] = {
  SetSafe -> True,
  Symbol -> Automatic};
DefineImmutable[RuleDelayed[pattern_, sym_Symbol], args_, OptionsPattern[]] := Check[
  Block[
    {sym},
    With[
      {instructions0 = If[#[[1,0]] === Hold, ReplacePart[#, {1,0} -> List], #]&[Hold[args]],
       type = Unique[SymbolName[Head[pattern]]],
       box = Replace[
         OptionValue[Symbol],
         {Automatic :> Unique[ToString[sym]],
          Except[_Symbol] :> Message[DefineImmutable::badarg, "Symbol option must be a symbol"]}],
       setFn = Replace[OptionValue[SetSafe], {True -> SetSafe, Except[True] -> Set}]},
      With[
        {instructions = Thread @ NestWhile[#[[{1},2]]&, instructions0, #[[1,0]] === With&],
         outerWithPos = Last @ NestWhile[
           {(#[[1]])[[{1}, 2]], Append[#[[2]], 2]}&,
           {instructions0, {}},
           (#[[1]])[[1, 0]] === With&]},
        With[
          {patterns = Apply[HoldPattern, #]& /@ instructions[[All, {1}, 1]],
           heads = Replace[
             instructions[[All, {1}, 0]],
             {Hold[Set] -> Set,
              Hold[Rule] -> Rule,
              Hold[RuleDelayed] -> RuleDelayed,
              Hold[SetDelayed] -> SetDelayed,
              x_ :> Message[
                DefineImmutable::badarg, 
                "Instructions bust be ->, =, :=, or :> forms: " <> ToString[x]]},
             {1}],
           bodies = instructions[[All, {1}, 2]]},
          If[Count[patterns, form_ /; form[[1,1]] === sym, {1}] < Length[patterns],
            Message[
              DefineImmutable::badarg, 
              "All field patterns must use the object as their first parameter"]];
          MapThread[
            Function[
              If[#1 =!= SetDelayed && #1 !MatchQ[Hold @@ #2, Hold[_[sym]]],
                Message[DefineImmutable::badarg, "rule and set functions must be unary"]]],
            {heads, patterns}];
          (* Setup the graph of dependencies *)
          With[
            {deps = Flatten@Last@Reap[
               (* We sort instructions into 4 groups:
                *   Init: instructions that must be calculated at startup (a -> b or a = b)
                *   Delay: instructions that needn't be calculated/saved (a := b)
                *   Settable: instructions that can be changed via Clone (a = b)
                *   Lazy: instructions that are delayed then memoized (a :> b)
                * We also collect Edges:
                *   An edge occurs when the RHS of an instruction references another instruction
                *)
               Sow[None, "Edge"];
               Sow[None, "Init"];
               Sow[None, "Delay"];
               Sow[None, "Settable"];
               Sow[None, "Lazy"];
               MapThread[
                 Function[{head, body, id},
                   If[head === SetDelayed,
                     Sow[id, "Delay"],
                     With[
                       {matches = Select[
                          Range[Length@patterns],
                          (heads[[#]] =!= SetDelayed && Count[body,patterns[[#]],Infinity] > 0)&]},
                       If[head =!= RuleDelayed, Sow[id, "Init"], Sow[id, "Lazy"]];
                       If[head === Set, Sow[id, "Settable"]];
                       Scan[
                         Function[
                           If[id == #,
                             Message[DefineImmutable::badarg, "Self-loop detected in dependencies"],
                             Sow[DirectedEdge[#, id], "Edge"]]],
                         matches]]]],
                 {heads, bodies, Range[Length@heads]}],
               {"Edge", "Init", "Delay", "Settable", "Lazy"},
               Rule[#1, Rest[#2]]&]},
            With[
              {delay = "Delay" /. deps,
               members = Complement[Range[Length@heads], "Delay" /. deps],
               init = "Init" /. deps,
               lazy = "Lazy" /. deps,
               lazyIdx = Part[
                 SortBy[
                   Join[
                     MapIndexed[#1 -> {2, #2[[1]]}&, "Lazy" /. deps],
                     MapIndexed[#1 -> {1, #2[[1]]}&, "Init" /. deps]],
                   First],
                 All, 2]},
              (* The easiest part to setup is the delayed items; do them now *)
              Do[
                TagSetDelayed @@ Join[
                  Hold[box],
                  ReplacePart[Hold @@ {patterns[[id]]}, {1,1,1} -> Pattern@@{sym, Blank[box]}],
                  bodies[[id]]],
                {id, "Delay" /. deps}];
              (* The direct accessors are also pretty easy to create... 
                 We create these slightly differently for lazy and non-lazy members. *)
              Do[
                TagSetDelayed @@ Join[
                  Hold[box],
                  ReplacePart[
                    Hold @@ {patterns[[members[[id]]]]}, 
                    {1,1,1} -> Pattern@@{sym, Blank[box]}],
                  ReplacePart[
                    (* For lazy: sym[[2, k]]; for init: sym[[1, k]] *)
                    Hold @@ {Hold[sym, Evaluate[lazyIdx[[id, 1]]], Evaluate[lazyIdx[[id, 2]]]]},
                    {1,0} -> Part]],
                {id, Range[Length@members]}];
              (* The constructor requires that we evaluate things in certain orders... 
                 we start by making some symols to use for each value *)
              With[
                {syms = Table[Unique[], {Length@members}],
                 (* We also make a graph of the dependencies *)
                 depGraph = With[
                   {tmp = Graph[members, "Edge" /. deps]},
                   If[!AcyclicGraphQ[tmp],
                     Message[DefineImmutable::badarg, "Cyclical dependency detected"],
                     tmp]]},
                (* with the transitive closure of the dependency graph, we know all the downstream
                   dependencies of any given element *)
                With[
                  {depStream = With[
                     {tc = TransitiveClosureGraph[depGraph]},
                     Map[
                       Function[
                         {EdgeList[tc, DirectedEdge[_,#]][[All, 1]],
                          EdgeList[tc, DirectedEdge[#,_]][[All, 2]]}],
                       members]],
                   (* Given the dependency graph, we can also determine the order non-delayed values
                      need to be evaluated in to prevent dependency errors *)
                   initOrder = Reverse@Part[
                     NestWhileList[
                       Function@With[
                         {goal = #[[1]], known = #[[3]]},
                         With[
                           {knowable = Select[
                              goal,
                              Function[
                                0 == Length@Complement[
                                  EdgeList[depGraph, DirectedEdge[_,#]][[All,1]], 
                                  known]]]},
                           {Complement[goal, knowable], knowable, Union[Join[known, knowable]]}]],
                       With[
                         {easy = Select[
                            "Init" /. deps,
                            EdgeCount[depGraph, DirectedEdge[_,#]] == 0&]},
                         {Complement["Init" /. deps, easy], easy, easy}],
                       Length[#[[1]]] > 0&],
                     All, 2]},
                  (* okay, now we build the constructor;
                     this is a big ugly code chunk due to the macro-nature of it... *)
                  SetDelayed @@ Join[
                    Hold[pattern],
                    With[
                      {constructorBody = ReplaceAll[
                         (* To construct the body of the custructor:
                          * (1) make a core: a body that sets the lazy symbols and makes the form
                          * (2) wrap the core in layers of With[{init...}, core]
                          *)
                         Fold[
                           (* (2): wrap the core-so-far in the next outer layer of init vars *)
                           Function[
                             ReplacePart[
                               Hold@Evaluate@Join[
                                 Hold@Evaluate@Map[
                                   Join[Hold@@{syms[[#]]}, bodies[[members[[#]]]]]&,
                                   #2],
                                 #1],
                               {{1,1,_,0} -> Set,
                                {1,0} -> With}]],
                           (* (1): make a core that sets the lazy symbols and yields the form *)
                           With[
                             {idcs = Map[
                                Position[members,#][[1,1]]&,
                                Complement[members, Flatten[initOrder]]],
                              formSym = Unique[],
                              lazyRevIdx = Part[
                                Map[
                                  Rest[SortBy[#[[All, 2 ;; 3]], Last]]&,
                                  SortBy[
                                    GatherBy[
                                      Join[
                                        MapIndexed[{#1[[1]], #2[[1]], #1[[2]]}&, lazyIdx],
                                        {{1, 0, 0}, {2, 0, 0}}],
                                      First],
                                    #[[1,1]]&]],
                                All, All, 1]},
                             ReplacePart[
                               Hold@Evaluate@Join[
                                 ReplacePart[
                                   Hold@Evaluate[Hold[#1,Unique[]]& /@ syms[[idcs]]],
                                   {1,_,0} -> Set],
                                 (* Here we make the actual form, prior to setting the symbols so
                                    that they are easier to hold: *)
                                 With[
                                   {form = Hold@Evaluate[
                                      box @@ ReplacePart[
                                        syms[[#]]& /@ lazyRevIdx,
                                        {2,0} -> Hold]]},
                                   ReplacePart[
                                     Hold@Evaluate@Join[
                                       ReplacePart[
                                         Hold @@ {{Flatten[Hold[Hold@Evaluate@formSym, form]]}},
                                         {1,1,0} -> Set],
                                       ReplacePart[
                                         Hold@Evaluate@Join[
                                           Flatten[
                                             Hold @@ MapThread[
                                               Function[
                                                 ReplacePart[
                                                   ReplaceAll[
                                                     Hold@Evaluate@Hold[
                                                       {#1}, 
                                                       Evaluate[Join[Hold[#1], #2]]],
                                                     {sym -> formSym}],
                                                   {{1,2,0} -> setFn,
                                                    {1,0} -> SetDelayed,
                                                    {1,1,0} -> Evaluate}]],
                                               {syms[[idcs]], bodies[[members[[idcs]]]]}]],
                                           Hold@Evaluate@formSym],
                                         {1,0} -> CompoundExpression]],
                                     {1,0} -> With]]],
                               {1,0} -> With]],
                           Map[Position[members, #][[1,1]]&, initOrder, {2}]],
                         MapThread[Rule, {patterns[[members]], syms}]]},
                      (* At this point, we need to return a Hold[] pattern, and constructorBody is
                         all that we really need; however, we may need to wrap some extra data
                         around our definitions... *)
                      If[outerWithPos == {}, 
                        constructorBody,
                        ReplacePart[
                          ReplacePart[instructions0, Prepend[outerWithPos, 1] -> constructorBody],
                          Join[{1}, outerWithPos, {0}] -> Identity]]]];
                  (* we also want to build a clone function *)
                  With[
                    {settables = "Settable" /. deps,
                     patternNames = #[[1,0]]& /@ patterns},
                    Evaluate[box] /: Clone[b_box, changes___Rule] := Check[
                      With[
                        {rules = {changes}},
                        Which[
                          (* No changes to be made; just return the current object *)
                          rules == {}, b,
                          (* multiple changes: do them one at a time *)
                          Length[rules] > 1, Fold[Clone, b, rules],
                          (* otherwise, we have one change; just handle it *)
                          True, With[
                            {id = Replace[
                               Position[patternNames, rules[[1,1]]], 
                               {{{i_}} /; MemberQ[settables, i] :> i,
                                _ :> Message[
                                  Clone::badarg,
                                  "unrecognized or unsettable rule: " <> ToString[rules[[1,1]]]]}],
                             val = rules[[1,2]]},
                            (* Iterate through the values that depend on this one *)
                            With[
                              {harvest = Reap[
                                 Fold[
                                   Function@With[
                                     {bCur = #1, iids = #2},
                                     ReplacePart[
                                       bCur,
                                       Map[
                                         Function[
                                           lazyIdx[[ Position[members, #1][[1,1]] ]] -> If[
                                             heads[[#1]] === RuleDelayed,
                                             With[
                                               {tmp = TemporarySymbol[]},
                                               Sow[tmp -> bodies[[#1]]];
                                               tmp],
                                             ReleaseHold@ReplaceAll[
                                               bodies[[#1]],
                                               sym -> bCur]]],
                                         iids]]],
                                   ReplacePart[b, lazyIdx[[ Position[members, id][[1,1]] ]] -> val],
                                   Rest@First@Last@Reap[
                                     NestWhile[
                                       Function[{G},
                                         With[
                                           {vs = Select[VertexList[G], VertexInDegree[G,#] == 0&]},
                                           Sow[vs];
                                           VertexDelete[G, vs]]],
                                       Subgraph[
                                         depGraph,
                                         VertexOutComponent[depGraph, id, Infinity]],
                                       VertexCount[#] > 0&]]]]},
                              (* harvest[[1]] is the new box'ed object; we need to setup the 
                                 temporary symbols, though. *)
                              With[
                                {finalSym = Unique[]},
                                Scan[
                                  Function@With[
                                    {memSym = #[[1]], memBody = #[[2]]},
                                    ReplacePart[
                                      {memSym, Join[
                                         Hold@Evaluate@memSym,
                                         ReplaceAll[memBody, sym -> finalSym]]},
                                      {{2,0} -> Set,
                                       0 -> SetDelayed}]],
                                  Flatten[harvest[[2]]]];
                                Set @@ Join[Hold[Evaluate@finalSym], Hold[harvest][[{1},1]]];
                                finalSym]]]]],
                      $Failed]];
                SetAttributes[Evaluate[box], Protected];
                True]]]]]]]],
  $Failed];

(* #MimicAssociation ******************************************************************************)
If[$MathematicaVersion >= 10.0,
  (MimicAssociation[args___] := Association[args]),
  (MimicAssociation[args___] := With[
     {sym = Unique["association"],
      assocs = Flatten[{args}]},
     If[!MatchQ[assocs, {_Rule...}],
       Message[
         MimicAssociation::badarg,
         "MimicAssociation must be given a list or sequence of rules"],
       (Scan[
          Function[sym /: sym[#[[1]]] = #[[2]]];
          assocs];
        sym /: sym[k_] := Missing["KeyAbsent", k];
        sym /: Lookup[sym, k_] := sym[k];
        sym /: Lookup[sym, k_, dflt_] := With[
          {tmp = sym[k]},
          If[tmp == Missing["KeyAbsent", k],
            dflt,
            tmp]];
        sym /: Keys[sym] = assocs[[All,1]];
        sym /: Values[sym] = assocs[[All,2]];
        sym /: KeyExistsQ[sym, k] = sym[k] != Missing["KeyAbsent", k];
        sym)]])];
Protect[MimicAssociation];

(* #NormalizeRows *********************************************************************************)
NormalizeRows[X_] := With[
  {tr = Transpose[X]},
  Transpose[tr / Table[#, {Length[tr]}]&@Sqrt@Total[tr^2]]];
NormalizeRows[{}] = {};

(* #NormalizeColumns ******************************************************************************)
NormalizeColumns[X_] := X / Table[#, {Length[X]}]&@Sqrt@Total[X^2];
NormalizeColumns[{}] = {};
Protect[NormalizeRows, NormalizeColumns];

(* #RowNorms **************************************************************************************)
RowNorms[X_] := With[{tr = Transpose[X]}, Sqrt @ Total[tr^2]];
RowNorms[{}] = {};
Protext[RowNorms];

(* #ColumnNorms ***********************************************************************************)
ColumnNorms[Xt_] := Sqrt @ Total[Xt^2];
ColumnNorms[{}] = {};
Protext[ColumnNorms];

(* #QuaternionToRotationMatrix *******************************************************************)
QuaternionToRotationMatrix[{a_, b_, c_, d_}] := {
  {a^2 + b^2 - c^2 - d^2, 2*b*c - 2*a*d, 2*b*d + 2*a*c},
  {2*b*c + 2*a*d, a^2 + c^2 - b^2 - d^2,   2*c*d - 2*a*b},
  {2*b*d - 2*a*c, 2*c*d + 2*a*b, a^2 + d^2 - c^2 - b^2}};
QuaterionToRotationMatrix[{b_, c_, d_}] := QuaternionToRotationMatrix[
  {Sqrt[1.0 - (b^2 + c^2 + d^2)], b, c, d}];
Protect[QuaternionToRotationMatrix];

(* #BinaryStringFix *******************************************************************************)
BinaryStringFix[clist_] := StringJoin[TakeWhile[clist, ToCharacterCode[#] != {0}&]];
Protect[BinaryStringFix];

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

(* #NeuroticaPermanentDatum ***********************************************************************)
$NeuroticaPermanentDataPath = FileNameJoin[{$UserBaseDirectory, ".NeuroticaPermanentData.dat"}];
$NeuroticaPermanentData = If[FileExistsQ[$NeuroticaPermanentDataPath],
  Association @ Map[
    Function[#[[1]] -> StringJoin[Riffle[Rest[#], "\t"]]],
    Import[$NeuroticaPermanentDataPath, "Table"]],
  Association[]];
Protect[$NeuroticaPermamentDataPath, $NeuroticaPermamentData];

NeuroticaPermanentDataSave[] := Check[
  (If[Length[$NeuroticaPermanentData] == 0,
     If[FileExistsQ[$NeuroticaPermanentDataPath], DeleteFile[$NeuroticaPermanentDataPath]],
     Export[
       $NeuroticaPermanentDataPath,
       Map[(List @@ #)&, Normal[$NeuroticaPermanentData]],
       "Table"]];
   True),
  False];
Protect[NeuroticaPermanentDataSave];

NeuroticaPermanentDatum[name_String, value_String] := With[
  {cur = Lookup[$NeuroticaPermanentData, name]},
  Unprotect[$NeuroticaPermanentData];
  $NeuroticaPermanentData[name] = value;
  Protect[$NeuroticaPermanentData];
  NeuroticaPermanentDataSave[]];
NeuroticaPermanentDatum[name_String, None] := With[
  {cur = Lookup[$NeuroticaPermanentData, name]},
  Unprotect[$NeuroticaPermanentData];
  KeyDropFrom[$NeuroticaPermanentData, name];
  Protect[$NeuroticaPermanentData];
  NeuroticaPermanentDataSave[]];
NeuroticaPermanentDatum[name_String] := Replace[
  $NeuroticaPermanentData[name],
  Missing[___] -> None];
Protect[NeuroticaPermanentDatum];


End[];
EndPackage[];
