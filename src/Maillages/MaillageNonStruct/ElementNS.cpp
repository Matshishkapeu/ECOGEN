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

#include "ElementNS.h"

using namespace std;

//***********************************************************************

ElementNS::ElementNS(){}

//***********************************************************************

ElementNS::ElementNS(const int &typeGmsh, const int &nombreNoeuds, const int &nombreFaces, const int &typeVTK) :
m_typeGmsh(typeGmsh),
m_nombreNoeuds(nombreNoeuds),
m_nombreFaces(nombreFaces),
m_typeVTK(typeVTK),
m_isFantome(false),
m_isCommunicant(false)
{
  m_numNoeuds = new int[nombreNoeuds];
}

//***********************************************************************

ElementNS::~ElementNS()
{
  delete[] m_numNoeuds;
  delete[] m_autresCPU;
}

//***********************************************************************

void ElementNS::construitElement(const int *numNoeuds, const Coord *noeuds, const int numeroEntitePhysique, const int numeroEntiteGeometrique, int &indiceElement)
{
  m_indice = indiceElement;
  //Attribution des numero de noeud vis a vis du tableau de noeud global et calcul position du Centre de l'element
  m_position = 0.;
  for (int i = 0; i < m_nombreNoeuds; i++){
    m_numNoeuds[i] = numNoeuds[i];
    m_position += noeuds[i];
  }
  m_position /= static_cast<double>(m_nombreNoeuds);

  m_appartenancePhysique = numeroEntitePhysique;
  m_appartenanceGeometrique = numeroEntiteGeometrique;

  //Calculs des proprietes de l'element
  this->calculVolume(noeuds);
  this->calculLCFL(noeuds);
}

//***********************************************************************

void ElementNS::construitElementParallele(const Coord *noeuds)
{
  Coord *noeudslocal = new Coord[m_nombreNoeuds];
  for (int i = 0; i < m_nombreNoeuds; i++){ noeudslocal[i] = noeuds[m_numNoeuds[i]]; }

  //Attribution des numero de noeud vis a vis du tableau de noeud global et calcul position du Centre de l'element
  for (int i = 0; i < m_nombreNoeuds; i++)
  {
    m_position += noeudslocal[i];
  }
  m_position /= static_cast<double>(m_nombreNoeuds);

  //Calculs des proprietes de l'element
  this->calculVolume(noeudslocal);
  this->calculLCFL(noeudslocal);
}

//***********************************************************************

void ElementNS::setIndice(int &indice)
{
  m_indice = indice;
}

//***********************************************************************

void ElementNS::setAppartenancePhysique(int &appartenancePhysique)
{
  m_appartenancePhysique = appartenancePhysique;
}

//***********************************************************************

void ElementNS::setNumNoeud(int *numNoeuds)
{
  for (int i = 0; i < m_nombreNoeuds; i++){ m_numNoeuds[i] = numNoeuds[i]; }
}

//***********************************************************************

void ElementNS::setNumNoeud(int &noeud, int &numNoeud)
{
  m_numNoeuds[noeud] = numNoeud;
}

//***********************************************************************

void ElementNS::setIsFantome(bool isFantome)
{
  m_isFantome = isFantome;
}

//***********************************************************************

void ElementNS::setIsCommunicant(bool isCommunicant)
{
  m_isCommunicant = isCommunicant;
}

//***********************************************************************

void ElementNS::setAppartenanceCPU(const int *numCPU, const int &nombreCPU)
{
  m_CPU = numCPU[0] - 1;
  m_nombreautresCPU = nombreCPU - 1;
  m_autresCPU = new int[m_nombreautresCPU];
  for (int i = 1; i < nombreCPU; i++)
  {
    m_autresCPU[i - 1] = -numCPU[i] - 1;
  }
}

//***********************************************************************

void ElementNS::enleveCPUAutres(vector<int> &numCPU)
{
  //if (numCPU >= m_nombreautresCPU){ Erreurs::messageErreur("Probleme dans enleveCPUAutres"); }
  //Copie des anciens
  int *autresCPUTemp = new int[m_nombreautresCPU];
  for (int i = 0; i < m_nombreautresCPU; i++)
  {
    autresCPUTemp[i] = m_autresCPU[i];
  }

  //reperage des CPU a enlever
  bool *enleveCPU = new bool[m_nombreautresCPU];
  for (int i = 0; i < m_nombreautresCPU; i++)
  {
    enleveCPU[i] = false;
    for (unsigned int p = 0; p < numCPU.size(); p++)
    {
      if (i == numCPU[p]){ enleveCPU[i] = true; break; }
    }
  }

  //Construction nouveau
  delete[] m_autresCPU;
  m_autresCPU = new int[m_nombreautresCPU - numCPU.size()];
  int indiceNouveau(0);
  for (int i = 0; i < m_nombreautresCPU; i++)
  {
    if (!enleveCPU[i]){ m_autresCPU[indiceNouveau++] = autresCPUTemp[i]; }
  }
  m_nombreautresCPU = indiceNouveau;

  delete[] autresCPUTemp;
  delete[] enleveCPU;
}

//***********************************************************************

int ElementNS::getIndice() const
{
  return m_indice;
}

//***********************************************************************

int ElementNS::getNombreNoeuds() const
{
  return m_nombreNoeuds;
}

//***********************************************************************

int ElementNS::getNombreFaces() const
{
  return m_nombreFaces;
}

//***********************************************************************

int ElementNS::getTypeGmsh() const
{
  return m_typeGmsh;
}

//***********************************************************************

int ElementNS::getTypeVTK() const
{
  return m_typeVTK;
}

//***********************************************************************

int ElementNS::getNumNoeud(int &noeud) const
{
  return m_numNoeuds[noeud];
}

//***********************************************************************

int ElementNS::getAppartenancePhysique() const
{
  return m_appartenancePhysique;
}

//***********************************************************************

int ElementNS::getAppartenanceGeometrique() const
{
  return m_appartenanceGeometrique;
}

//***********************************************************************

int ElementNS::getCPU() const
{
  return m_CPU;
}

//***********************************************************************

int ElementNS::getNombreAutresCPU() const
{
  return m_nombreautresCPU;
}

//***********************************************************************

int ElementNS::getAutreCPU(const int &autreCPU) const
{
  if (sizeof(m_autresCPU) <= autreCPU)
  {
    Erreurs::messageErreur("probleme de dimension dans m_autresCPU");
  }
  return m_autresCPU[autreCPU];
}

//***********************************************************************

bool ElementNS::isFantome() const
{
  return m_isFantome;
}

//***********************************************************************

bool ElementNS::isCommunicant() const
{
  return m_isCommunicant;
}

//***********************************************************************

void ElementNS::afficheInfos() const
{
  cout << "-------------" << endl;
  //cout << "Infos Element" << endl;
  //for (int i = 0; i < m_nombreNoeuds; i++)
  //  cout << " " << m_numNoeuds[i];
  //cout << endl;
  cout << "centre : " << m_position.getX() << " " << m_position.getY() << " " << m_position.getZ() << endl;
  //cout << " volume : " << m_volume << endl;
  //cout << " lCfl : " << m_lCFL << endl;
  //cout << " CPU n : " << m_CPU << endl;
  //for (int i = 0; i < m_nombreautresCPU; i++)
  //  cout << " autre CPU : " << m_autresCPU[i] << endl;
}

//***********************************************************************