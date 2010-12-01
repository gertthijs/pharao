/*******************************************************************************
pharMerger.h - Pharao
 
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



#ifndef __SILICOS_PHARAO_PHARMERGER_H__
#define __SILICOS_PHARAO_PHARMERGER_H__



// General
#include <set>
#include <list>

// OpenBabel

// Pharao
#include "pharmacophore.h"
#include "utilities.h"



class PharMerger
{
   private:
        
      double  _deltaSigma;
      double  _threshold;
   
   public: 
  
      PharMerger();
  
      void merge(Pharmacophore& phar);
  };



#endif __SILICOS_PHARAO_PHARMERGER_H__
