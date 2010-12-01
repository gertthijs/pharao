/*******************************************************************************
options.cpp - Pharao
 
Copyright (C) 2005-2010 by Silicos NV
 
This file is part of the Open Babel project.
For more information, see <http://openbabel.sourceforge.net/>
 
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
*******************************************************************************/



#include "options.h"



Options::Options(): 
   cutOff(0.0),
   best(0),
   rankby(TANIMOTO),
   funcGroupVec(),
	noHybrid(false),
   epsilon(0.5),
   noAlign(false),
   withExclusion(false), 
   merge(false),
   noNormal(false),
   scoreOnly(false),
   isQuiet(false),
   version(false),
   help(false)
{
   molOutStream = NULL;
   molOutWriter = NULL;
   molOutFile = "";

   pharmOutStream = NULL;
   pharmOutWriter = NULL;
   pharmOutFile = "";

   scoreOutStream = NULL;
   scoreOutFile = "";

   dbInpStream = NULL;
   dbInpFile = "";
	dbInpType = UNKNOWN;

   refInpStream = NULL;
   refInpFile = "";
   refInpType = UNKNOWN;
}




Options::~Options(void)
{
   //
   // reference input
   if (!refInpFile.empty())
   {
      refInpFile = "";
   };
   if (refInpStream)
   {
      delete refInpStream;
      refInpStream = NULL;
   };


   //
   // Database input
   if (!dbInpFile.empty())
   {
      dbInpFile = "";
   };
   if (dbInpStream)
   {
      delete dbInpStream;
      dbInpStream = NULL;
   };


   //
   // Molecule output
   if (!molOutFile.empty())
   {
      molOutFile = "";
   };
   if (molOutStream)
   {
      delete molOutStream;
      molOutStream = NULL;
   };
   if (molOutWriter)
   {
      delete molOutWriter;
      molOutWriter = NULL;
   };
   
   
   //
   // Pharmacopore output
   if (!pharmOutFile.empty())
   {
      pharmOutFile = "";
   };
   if (pharmOutStream)
   {
      delete pharmOutStream;
      pharmOutStream = NULL;
   };
   if (pharmOutWriter)
   {
      delete pharmOutWriter;
      pharmOutWriter = NULL;
   };
   
   
   //
   // Score output
   if (!scoreOutFile.empty())
   {
      scoreOutFile = "";
   };
   if (scoreOutStream)
   {
      delete scoreOutStream;
      scoreOutStream = NULL;
   };


}
