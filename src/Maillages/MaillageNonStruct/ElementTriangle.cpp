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

#include "ElementTriangle.h"

using namespace std;

const int ElementTriangle::TYPEGMSH = 2;
const int ElementTriangle::NOMBRENOEUDS = 3;
const int ElementTriangle::NOMBREFACES = 3; /* ici il s'agit du nombre de segments*/
const int ElementTriangle::TYPEVTK = 5;

//***********************************************************************

ElementTriangle::ElementTriangle() :
ElementNS(TYPEGMSH, NOMBRENOEUDS, NOMBREFACES, TYPEVTK)
{}

//***********************************************************************

ElementTriangle::~ElementTriangle(){}

//***********************************************************************

void ElementTriangle::calculVolume(const Coord *noeuds)
{
  Coord v1(noeuds[1] - noeuds[0]), v2(noeuds[2] - noeuds[1]), v3(noeuds[0] - noeuds[2]);
  double a(v1.norme()), b(v2.norme()), c(v3.norme());
  double dp = 0.5*(a + b + c);
  m_volume = sqrt(dp*(dp - a)*(dp - b)*(dp - c)); //Aire du triangle
}

//***********************************************************************

void ElementTriangle::calculLCFL(const Coord *noeuds)
{
  Coord vec; m_lCFL = 1e10;
  vec = ((noeuds[0] + noeuds[1]) / 2.) - m_position;
  m_lCFL = min(m_lCFL, vec.norme());
  vec = ((noeuds[0] + noeuds[2]) / 2.) - m_position;
  m_lCFL = min(m_lCFL, vec.norme());
  vec = ((noeuds[1] + noeuds[2]) / 2.) - m_position;
  m_lCFL = min(m_lCFL, vec.norme());
}

//***********************************************************************
// Nouvelle version beaucoup plus efficace avec recherche dans tableau temporaire
void ElementTriangle::construitFaces(const Coord *noeuds, FaceNS **faces, int &iMax, int** facesTemp, int* sommeNoeudsTemp)
{
  //3 faces a traiter de type segment
  int indiceFaceExiste(-1);
  int noeudAutre;
  for (int i = 0; i < NOMBREFACES; i++)
  {
    switch (i)
    {
      case 0: facesTemp[iMax][0] = m_numNoeuds[0]; facesTemp[iMax][1] = m_numNoeuds[1]; noeudAutre = 2; break;
      case 1: facesTemp[iMax][0] = m_numNoeuds[1]; facesTemp[iMax][1] = m_numNoeuds[2]; noeudAutre = 0; break;      
      case 2: facesTemp[iMax][0] = m_numNoeuds[2]; facesTemp[iMax][1] = m_numNoeuds[0]; noeudAutre = 1; break;      
    }
    sommeNoeudsTemp[iMax] = facesTemp[iMax][0] + facesTemp[iMax][1];
    std::sort(facesTemp[iMax],facesTemp[iMax]+2);  //Tri des noeuds
    //Existance face ?
    indiceFaceExiste = FaceNS::rechercheFace(facesTemp[iMax],sommeNoeudsTemp[iMax],facesTemp,sommeNoeudsTemp,2,iMax);
    //Creation face ou rattachement
    if (indiceFaceExiste==-1)
    {
      faces[iMax] = new FaceSegment(facesTemp[iMax][0], facesTemp[iMax][1], 0); //pas besoin du tri ici
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
void ElementTriangle::construitFacesSimplifie(int &iMax, int** facesTemp, int* sommeNoeudsTemp)
{
  //3 faces a traiter de type segment
  int indiceFaceExiste(-1);
  for (int i = 0; i < NOMBREFACES; i++)
  {
    switch (i)
    {
      case 0: facesTemp[iMax][0] = m_numNoeuds[0]; facesTemp[iMax][1] = m_numNoeuds[1]; break;
      case 1: facesTemp[iMax][0] = m_numNoeuds[1]; facesTemp[iMax][1] = m_numNoeuds[2]; break;      
      case 2: facesTemp[iMax][0] = m_numNoeuds[2]; facesTemp[iMax][1] = m_numNoeuds[0]; break;      
    }
    sommeNoeudsTemp[iMax] = facesTemp[iMax][0] + facesTemp[iMax][1];
    std::sort(facesTemp[iMax],facesTemp[iMax]+2);  //Tri des noeuds
    //Existance face ?
    indiceFaceExiste = FaceNS::rechercheFace(facesTemp[iMax],sommeNoeudsTemp[iMax],facesTemp,sommeNoeudsTemp,2,iMax);
    //Creation face ou rattachement
    if (indiceFaceExiste==-1)
    {
      iMax++;
    }
  }
}

//***********************************************************************

void ElementTriangle::attributFaceLimite(const Coord *noeuds, FaceNS **faces, const int &indiceMaxFaces)
{
  int indiceFaceExiste(0);
  FaceTriangle face(m_numNoeuds[0], m_numNoeuds[1], m_numNoeuds[2]);
  if (face.faceExiste(faces, indiceMaxFaces, indiceFaceExiste))
  {
    faces[indiceFaceExiste]->ajouteElementVoisinLimite(this);
  }
  else
  {
    Erreurs::messageErreur("Probleme attribution des faces limites element Triangle");
  }
}

//***********************************************************************

void ElementTriangle::attributFaceCommunicante(const Coord *noeuds, FaceNS **faces, const int &indiceMaxFaces, const int &nombreNoeudsInternes)
{
  int indiceFaceExiste(0);
  //Verification face 1 :
  if (m_numNoeuds[0] < nombreNoeudsInternes && m_numNoeuds[1] < nombreNoeudsInternes)
  {
    FaceSegment face(m_numNoeuds[0], m_numNoeuds[1]);
    if (face.faceExiste(faces, indiceMaxFaces, indiceFaceExiste))
    {
      faces[indiceFaceExiste]->ajouteElementVoisinLimite(this);
      faces[indiceFaceExiste]->setEstComm(true);
    }
  }
  //Verification face 2 :
  if (m_numNoeuds[1] < nombreNoeudsInternes && m_numNoeuds[2] < nombreNoeudsInternes)
  {
    FaceSegment face(m_numNoeuds[1], m_numNoeuds[2]);
    if (face.faceExiste(faces, indiceMaxFaces, indiceFaceExiste))
    {
      faces[indiceFaceExiste]->ajouteElementVoisinLimite(this);
      faces[indiceFaceExiste]->setEstComm(true);
    }
  }
  //Verification face 3 :
  if (m_numNoeuds[2] < nombreNoeudsInternes && m_numNoeuds[0] < nombreNoeudsInternes)
  {
    FaceSegment face(m_numNoeuds[2], m_numNoeuds[0]);
    if (face.faceExiste(faces, indiceMaxFaces, indiceFaceExiste))
    {
      faces[indiceFaceExiste]->ajouteElementVoisinLimite(this);
      faces[indiceFaceExiste]->setEstComm(true);
    }
  }
}

//***********************************************************************

int ElementTriangle::compteFaceCommunicante(vector<int*> &facesTemp, vector<int> &sommeNoeudsTemp)
{
  //4 faces a traiter de type segment
  int indiceFaceExiste(-1), nombreFacesCommunicante(0);
  int face[2], sommeNoeuds;
  for (int i = 0; i < NOMBREFACES; i++)
  {
    switch (i)
    {
      case 0: face[0] = m_numNoeuds[0]; face[1] = m_numNoeuds[1]; break;
      case 1: face[0] = m_numNoeuds[1]; face[1] = m_numNoeuds[2]; break;
      case 2: face[0] = m_numNoeuds[2]; face[1] = m_numNoeuds[0]; break;     
    }
    int iMax = sommeNoeudsTemp.size();
    sommeNoeuds = face[0]+face[1];
    std::sort(face, face+2);
    //Recherche existance faces
    indiceFaceExiste = FaceNS::rechercheFace(face,sommeNoeuds,facesTemp,sommeNoeudsTemp,2,iMax);
    if (indiceFaceExiste!=-1)
    {
      nombreFacesCommunicante++;
    }
  }
  return nombreFacesCommunicante;
}

//***********************************************************************
//Nouvelle version plus efficace
int ElementTriangle::compteFaceCommunicante(int &iMax, int **facesTemp, int *sommeNoeudsTemp)
{
  //4 faces a traiter de type segment
  int indiceFaceExiste(-1), nombreFacesCommunicante(0);
  int face[2], sommeNoeuds;
  for (int i = 0; i < NOMBREFACES; i++)
  {
    switch (i)
    {
      case 0: face[0] = m_numNoeuds[0]; face[1] = m_numNoeuds[1]; break;
      case 1: face[0] = m_numNoeuds[1]; face[1] = m_numNoeuds[2]; break;
      case 2: face[0] = m_numNoeuds[2]; face[1] = m_numNoeuds[0]; break;       
    }
    sommeNoeuds = face[0]+face[1];
    std::sort(face, face+2);
    //Recherche existance faces
    indiceFaceExiste = FaceNS::rechercheFace(face,sommeNoeuds,facesTemp,sommeNoeudsTemp,2,iMax);
    if (indiceFaceExiste!=-1)
    {
      nombreFacesCommunicante++;
    }
  }
  return nombreFacesCommunicante;
}

//***********************************************************************