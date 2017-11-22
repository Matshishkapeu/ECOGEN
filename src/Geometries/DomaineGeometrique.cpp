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

#include "DomaineGeometrique.h"

using namespace std;

//******************************************************************

DomaineGeometrique::DomaineGeometrique(){}

//******************************************************************

DomaineGeometrique::DomaineGeometrique(string nom, vector<Phase*> vecPhases, Melange *melange, vector<Transport> vecTransports) : m_nom(nom)
{
  m_nombrePhases = vecPhases.size();
  m_nombreTransports = vecTransports.size();
  m_vecPhases = new Phase*[m_nombrePhases];
  for (int k = 0; k < m_nombrePhases; k++){
    vecPhases[k]->alloueEtCopiePhase(&m_vecPhases[k]);
  }
  melange->alloueEtCopieMelange(&m_melange);
  if (m_nombreTransports > 0) { m_vecTransports = new Transport[m_nombreTransports]; }
  for (int k = 0; k < m_nombreTransports; k++) {
    m_vecTransports[k].setValeur(vecTransports[k].getValeur());
  }
}

//******************************************************************

DomaineGeometrique::~DomaineGeometrique()
{
  for (int k = 0; k < m_nombrePhases; k++) {
    delete m_vecPhases[k];
  }
  delete[] m_vecPhases;
  delete m_melange;
  if (m_nombreTransports != 0) delete[] m_vecTransports;
}

//******************************************************************

void DomaineGeometrique::rempli(Cellule *cellule, const int &nombrePhases, const int &nombreTransports) const
{
  for (int k = 0; k < nombrePhases; k++) { cellule->copiePhase(k, m_vecPhases[k]); }
  cellule->copieMelange(m_melange);
  for (int k = 0; k < nombreTransports; k++) { cellule->setTransport(m_vecTransports[k].getValeur(), k); }
}