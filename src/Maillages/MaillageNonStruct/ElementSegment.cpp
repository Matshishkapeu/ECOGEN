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

#include "ElementSegment.h"

using namespace std;

const int ElementSegment::TYPEGMSH = 1;
const int ElementSegment::NOMBRENOEUDS = 2;
const int ElementSegment::NOMBREFACES = 2;
const int ElementSegment::TYPEVTK = 3;

//***********************************************************************

ElementSegment::ElementSegment() :
ElementNS(TYPEGMSH, NOMBRENOEUDS, NOMBREFACES, TYPEVTK)
{}

//***********************************************************************

ElementSegment::~ElementSegment(){}

//***********************************************************************

void ElementSegment::calculVolume(const Coord *noeuds)
{
   m_volume = (noeuds[1] - noeuds[0]).norme(); //longueur du segment
}

//***********************************************************************

void ElementSegment::calculLCFL(const Coord *noeuds)
{
  m_lCFL = (noeuds[1] - noeuds[0]).norme()/2.0; //demi longueur du segment
}

//***********************************************************************
// Nouvelle version beaucoup plus efficace avec recherche dans tableau temporaire
void ElementSegment::construitFaces(const Coord *noeuds, FaceNS **faces, int &iMax, int** facesTemp, int* sommeNoeudsTemp)
{
  //2 faces a traiter de type point
  int indiceFaceExiste(-1);
  int noeudAutre;
  for (int i = 0; i < NOMBREFACES; i++)
  {
    switch (i)
    {
      case 0: facesTemp[iMax][0] = m_numNoeuds[0]; noeudAutre = 1; break;
      case 1: facesTemp[iMax][0] = m_numNoeuds[1]; noeudAutre = 0; break;      
    }
    sommeNoeudsTemp[iMax] = facesTemp[iMax][0];
    //Existance face ?
    indiceFaceExiste = FaceNS::rechercheFace(facesTemp[iMax],sommeNoeudsTemp[iMax],facesTemp,sommeNoeudsTemp,1,iMax);
    //Creation face ou rattachement
    if (indiceFaceExiste==-1)
    {
      faces[iMax] = new FacePoint(facesTemp[iMax][0]); //pas besoin du tri ici
      faces[iMax]->construitFace(noeuds, m_numNoeuds[noeudAutre], this);
      iMax++;
    }
    else
    {
      faces[indiceFaceExiste]->ajouteElementVoisin(this);
    }
  }
}

//***********************************************************************
// Nouvelle version beaucoup plus efficace avec recherche dans tableau temporaire
void ElementSegment::construitFacesSimplifie(int &iMax, int** facesTemp, int* sommeNoeudsTemp)
{
  //2 faces a traiter de type point
  int indiceFaceExiste(-1);
  for (int i = 0; i < NOMBREFACES; i++)
  {
    switch (i)
    {
      case 0: facesTemp[iMax][0] = m_numNoeuds[0]; break;
      case 1: facesTemp[iMax][0] = m_numNoeuds[1]; break;
    }
    sommeNoeudsTemp[iMax] = facesTemp[iMax][0];
    //Existance face ?
    indiceFaceExiste = FaceNS::rechercheFace(facesTemp[iMax], sommeNoeudsTemp[iMax], facesTemp, sommeNoeudsTemp, 1, iMax);
    //Creation face ou rattachement
    if (indiceFaceExiste == -1)
    {
      iMax++;
    }
  }
}

//***********************************************************************

void ElementSegment::attributFaceLimite(const Coord *noeuds, FaceNS **faces, const int &indiceMaxFaces)
{
  int indiceFaceExiste(0);
  FaceSegment face(m_numNoeuds[0], m_numNoeuds[1]);
  if (face.faceExiste(faces, indiceMaxFaces, indiceFaceExiste))
  {
    faces[indiceFaceExiste]->ajouteElementVoisinLimite(this);
  }
  else
  {
    Erreurs::messageErreur("Probleme attribution des faces limites element Segment");
  }
}

//***********************************************************************

void ElementSegment::attributFaceCommunicante(const Coord *noeuds, FaceNS **faces, const int &indiceMaxFaces, const int &nombreNoeudsInternes)
{
  int indiceFaceExiste(0);
  //Verification face 1 :
  if (m_numNoeuds[0] < nombreNoeudsInternes)
  {
    FacePoint face(m_numNoeuds[0]);
    if (face.faceExiste(faces, indiceMaxFaces, indiceFaceExiste))
    {
      faces[indiceFaceExiste]->ajouteElementVoisinLimite(this);
      faces[indiceFaceExiste]->setEstComm(true);
    }
  }
  //Verification face 2 :
  if (m_numNoeuds[1] < nombreNoeudsInternes)
  {
    FacePoint face(m_numNoeuds[1]);
    if (face.faceExiste(faces, indiceMaxFaces, indiceFaceExiste))
    {
      faces[indiceFaceExiste]->ajouteElementVoisinLimite(this);
      faces[indiceFaceExiste]->setEstComm(true);
    }
  }
}

//***********************************************************************

int ElementSegment::compteFaceCommunicante(vector<int*> &facesTemp, vector<int> &sommeNoeudsTemp)
{
  //2 faces a traiter de type point
  int indiceFaceExiste(-1), nombreFacesCommunicante(0);
  int point;
  for (int i = 0; i < NOMBREFACES; i++)
  {
    switch (i)
    {
      case 0: point = m_numNoeuds[0]; break;
      case 1: point = m_numNoeuds[1]; break;
    }
    int iMax = sommeNoeudsTemp.size();
    //Recherche existance faces
    indiceFaceExiste = FaceNS::rechercheFace(&point, point, facesTemp, sommeNoeudsTemp, 1, iMax);
    if (indiceFaceExiste != -1)
    {
      nombreFacesCommunicante++;
    }
  }
  return nombreFacesCommunicante;
}

//***********************************************************************
//Nouvelle version plus efficace
int ElementSegment::compteFaceCommunicante(int &iMax, int **facesTemp, int *sommeNoeudsTemp)
{
  //2 faces a traiter de type point
  int indiceFaceExiste(-1), nombreFacesCommunicante(0);
  int point;
  for (int i = 0; i < NOMBREFACES; i++)
  {
    switch (i)
    {
      case 0: point = m_numNoeuds[0]; break;
      case 1: point = m_numNoeuds[1]; break;
    }
    //Recherche existance faces
    indiceFaceExiste = FaceNS::rechercheFace(&point, point, facesTemp, sommeNoeudsTemp, 1, iMax);
    if (indiceFaceExiste != -1)
    {
      nombreFacesCommunicante++;
    }
  }
  return nombreFacesCommunicante;
}

//***********************************************************************