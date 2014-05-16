;    This file is part of XRDUA.
;
;    XRDUA is free software: you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation, either version 3 of the License, or
;    (at your option) any later version.
;
;    XRDUA is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;    GNU General Public License for more details.
;
;    You should have received a copy of the GNU General Public License
;    along with XRDUA.  If not, see <http://www.gnu.org/licenses/>.
;
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FUNCTION ProgramRootDir, OneUp=oneup, TwoUp=twoup, NoMark=nomark

   ; Return to caller on an error.
   On_Error, 2

   ; Get the current call stack.
   Help, Calls=callStack

   ; Get the name of the calling routine.
   callingRoutine = (StrSplit(StrCompress(callStack[1])," ", /Extract))[0]

   ; We don't know if the calling routine is a procedure or a function,
   ; and we don't have a way to get information without knowing this. So,
   ; we are going to try first to see if it is a procedure. If not, we
   ; will try it as a function. Unfortunately, if it is *not* a procedure,
   ; we will cause an error. We have to catch that and handle it silently.

   Catch, theError
   IF theError NE 0 THEN BEGIN
      Catch, /Cancel
      Message, /Reset
      thisRoutine = Routine_Info(callingRoutine, /Functions, /Source)
   ENDIF

   IF N_Elements(thisRoutine) EQ 0 THEN $
      thisRoutine = Routine_Info(callingRoutine, /Source)

   ; If there are no path separators, you are here.
   IF ( StrPos(thisRoutine.Path, Path_Sep()) ) EQ -1 THEN BEGIN
      CD, Current=thisDir
      sourcePath = FilePath(thisRoutine.Path, Root_Dir=thisDir)
   ENDIF ELSE BEGIN
      sourcePath = thisRoutine.Path
   ENDELSE

   ; Strip the root directory off the source path.
   root = StrMid(sourcePath, 0, StrPos(sourcePath, Path_Sep(), /Reverse_Search) + 1)

   ; If ONEUP is set, then climb up the root directory by one directory. This will
   ; be the *second* path separator, since the root directory has a path separator
   ; as its end.
   IF Keyword_Set(oneup) THEN BEGIN
      i = Where( Byte(root) EQ (Byte(Path_Sep()))[0], count)

      IF count GE 2 THEN BEGIN
         sourcePath = StrMid(root, 0, StrLen(root)-1)
         root = StrMid(sourcePath, 0, StrPos(sourcePath, Path_Sep(), /Reverse_Search) + 1)
      ENDIF
   ENDIF

   ; If TWOUP is set, then climb up the root directory by two directories. This will
   ; be the *third* path separator, since the root directory has a path separator
   ; as its end.
   IF Keyword_Set(twoup) THEN BEGIN
      i = Where( Byte(root) EQ (Byte(Path_Sep()))[0], count)

      IF count GE 3 THEN BEGIN
         sourcePath = StrMid(root, 0, StrLen(root)-1)
         root = StrMid(sourcePath, 0, StrPos(sourcePath, Path_Sep(), /Reverse_Search) + 1)
         sourcePath = StrMid(root, 0, StrLen(root)-1)
         root = StrMid(sourcePath, 0, StrPos(sourcePath, Path_Sep(), /Reverse_Search) + 1)
      ENDIF
   ENDIF

   ; Remove last path separation mark, if requested.
   IF Keyword_Set(nomark) THEN RETURN, StrMid(root, 0, StrLen(root)-1) ELSE RETURN, root
END;FUNCTION ProgramRootDir
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro CDProgramRootDir,_ref_extra=ex
; Set working dir to the directory of this file (or one up)
cd,ProgramRootDir(_extra=ex)
end;pro CDProgramRootDir
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
