//  
//       ,---.     ,--,    .---.     ,--,    ,---.    .-. .-. 
//       | .-'   .' .')   / .-. )  .' .'     | .-'    |  \| | 
//       | `-.   |  |(_)  | | |(_) |  |  __  | `-.    |   | | 
//       | .-'   \  \     | | | |  \  \ ( _) | .-'    | |\  | 
//       |  `--.  \  `-.  \ `-' /   \  `-) ) |  `--.  | | |)| 
//       /( __.'   \____\  )---'    )\____/  /( __.'  /(  (_) 
//      (__)              (_)      (__)     (__)     (__)     
//
//  This file is part of ECOGEN.
//
//  ECOGEN is the legal property of its developers, whose names 
//  are listed in the copyright file included with this source 
//  distribution.
//
//  ECOGEN is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published 
//  by the Free Software Foundation, either version 3 of the License, 
//  or (at your option) any later version.
//  
//  ECOGEN is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//  GNU General Public License for more details.
//  
//  You should have received a copy of the GNU General Public License
//  along with ECOGEN (file LICENSE).  
//  If not, see <http://www.gnu.org/licenses/>.

#include "Outils.h"

using namespace std;

Outils *BO;

//***********************************************************************

Outils::Outils() : ak(0), rhok(0), pk(0), akS(0), rhokS(0),eos(0)
{}

//***********************************************************************

Outils::Outils(const int &nombrePhases)
{
  m_nombrePhases = nombrePhases;
  ak = new double[nombrePhases];
  rhok = new double[nombrePhases];
  pk = new double[nombrePhases];
  akS = new double[nombrePhases];
  rhokS = new double[nombrePhases];
  eos = new Eos*[nombrePhases];
}

//***********************************************************************

Outils::~Outils()
{
  if (ak!=0) delete[] ak;
  if (rhok != 0) delete[] rhok;
  if (pk != 0) delete[] pk;
  if (akS != 0) delete[] akS;
  if (rhokS != 0) delete[] rhokS;
  if (eos != 0) delete[] eos;
}

//***********************************************************************

void Outils::majuscule(string &chaine)
{
  for (unsigned int c = 0; c < chaine.size(); c++){ chaine[c] = toupper(chaine[c]); }
}

//***********************************************************************

double Outils::pi()
{
  //return acos(-1.);
  return 3.14159;
}

//***********************************************************************

void Outils::alloueOutil(const int &nombrePhases)
{
  m_nombrePhases = nombrePhases;
  ak = new double[nombrePhases];
  rhok = new double[nombrePhases];
  pk = new double[nombrePhases];
  akS = new double[nombrePhases];
  rhokS = new double[nombrePhases];
  eos = new Eos*[nombrePhases];
}

//***********************************************************************