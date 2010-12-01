/*******************************************************************************
solutionInfo.h - Pharao
 
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



#ifndef __SILICOS_PHARAO_SOLUTIONINFO_H__
#define __SILICOS_PHARAO_SOLUTIONINFO_H__



// General

// OpenBabel

// Pharao
#include "coordinate.h"
#include "siMath.h"



class SolutionInfo
{
   public:
   
		SiMath::Vector       rotor;
		double               volume;
		unsigned int         iterations;
		Coordinate           center1;
		Coordinate           center2;
		SiMath::Matrix       rotation1;
		SiMath::Matrix       rotation2;
		
		SolutionInfo(void);
      ~SolutionInfo(void);
};



#endif __SILICOS_PHARAO_SOLUTIONINFO_H__
