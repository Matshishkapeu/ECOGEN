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

#include "ElementQuadrangle.h"

using namespace std;

const int ElementQuadrangle::TYPEGMSH = 3;
const int ElementQuadrangle::NOMBRENOEUDS = 4;
const int ElementQuadrangle::NOMBREFACES = 4; /* ici il s'agit du nombre de segments*/
const int ElementQuadrangle::TYPEVTK = 9;

//***********************************************************************

ElementQuadrangle::ElementQuadrangle() :
ElementNS(TYPEGMSH, NOMBRENOEUDS, NOMBREFACES, TYPEVTK)
{}

//***********************************************************************

ElementQuadrangle::~ElementQuadrangle(){}

//***********************************************************************

void ElementQuadrangle::calculVolume(const Coord *noeuds)
{
  //une diagonale :
  Coord v0(noeuds[2] - noeuds[0]);
  double diagonale(v0.norme());
  //Les 4 cotes :
  Coord v1(noeuds[1] - noeuds[0]);
  Coord v2(noeuds[2] - noeuds[1]);
  Coord v3(noeuds[3] - noeuds[2]);
  Coord v4(noeuds[0] - noeuds[3]);
  double a(v1.norme()); double b(v2.norme()); double c(v3.norme()); double d(v4.norme());
  //Aire premier triangle
  double dp1 = 0.5*(a + b + diagonale);
  double surf1 = sqrt(dp1*(dp1 - a)*(dp1 - b)*(dp1 - diagonale));
  //Aire second triangle
  double dp2 = 0.5*(c + d + diagonale);
  double surf2 = sqrt(dp2*(dp2 - c)*(dp2 - d)*(dp2 - diagonale));
  //Air quadrangle
  m_volume = surf1 + surf2; //aire du quadrangle
}

//***********************************************************************

void ElementQuadrangle::calculLCFL(const Coord *noeuds)
{
  Coord vec; m_lCFL = 1e10;
  vec = ((noeuds[0] + noeuds[1]) / 2.) - m_position;
  m_lCFL = min(m_lCFL, vec.norme());
  vec = ((noeuds[1] + noeuds[2]) / 2.) - m_position;
  m_lCFL = min(m_lCFL, vec.norme());
  vec = ((noeuds[2] + noeuds[3]) / 2.) - m_position;
  m_lCFL = min(m_lCFL, vec.norme());
  vec = ((noeuds[3] + noeuds[0]) / 2.) - m_position;
  m_lCFL = min(m_lCFL, vec.norme());
}

//***********************************************************************
// Nouvelle version beaucoup plus efficace avec recherche dans tableau temporaire
void ElementQuadrangle::construitFaces(const Coord *noeuds, FaceNS **faces, int &iMax, int** facesTemp, int* sommeNoeudsTemp)
{
  //4 faces a traiter de type segment
  int indiceFaceExiste(-1);
  int noeudAutre;
  for (int i = 0; i < NOMBREFACES; i++)
  {
    switch (i)
    {
      case 0: facesTemp[iMax][0] = m_numNoeuds[0]; facesTemp[iMax][1] = m_numNoeuds[1]; noeudAutre = 2; break;
      case 1: facesTemp[iMax][0] = m_numNoeuds[1]; facesTemp[iMax][1] = m_numNoeuds[2]; noeudAutre = 3; break;      
      case 2: facesTemp[iMax][0] = m_numNoeuds[2]; facesTemp[iMax][1] = m_numNoeuds[3]; noeudAutre = 0; break;      
      case 3: facesTemp[iMax][0] = m_numNoeuds[3]; facesTemp[iMax][1] = m_numNoeuds[0]; noeudAutre = 1; break;      
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
void ElementQuadrangle::construitFacesSimplifie(int &iMax, int** facesTemp, int* sommeNoeudsTemp)
{
  //4 faces a traiter de type segment
  int indiceFaceExiste(-1);
  for (int i = 0; i < NOMBREFACES; i++)
  {
    switch (i)
    {
      case 0: facesTemp[iMax][0] = m_numNoeuds[0]; facesTemp[iMax][1] = m_numNoeuds[1]; break;
      case 1: facesTemp[iMax][0] = m_numNoeuds[1]; facesTemp[iMax][1] = m_numNoeuds[2]; break;      
      case 2: facesTemp[iMax][0] = m_numNoeuds[2]; facesTemp[iMax][1] = m_numNoeuds[3]; break;      
      case 3: facesTemp[iMax][0] = m_numNoeuds[3]; facesTemp[iMax][1] = m_numNoeuds[0]; break;      
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

void ElementQuadrangle::attributFaceLimite(const Coord *noeuds, FaceNS **faces, const int &indiceMaxFaces)
{
  int indiceFaceExiste(0);
  FaceQuadrangle face(m_numNoeuds[0], m_numNoeuds[1], m_numNoeuds[2], m_numNoeuds[3]);
  if (face.faceExiste(faces, indiceMaxFaces, indiceFaceExiste))
  {
    faces[indiceFaceExiste]->ajouteElementVoisinLimite(this);
  }
  else
  {
    Erreurs::messageErreur("Probleme attribution des faces limites element Quadrangle");
  }
}

//***********************************************************************

void ElementQuadrangle::attributFaceCommunicante(const Coord *noeuds, FaceNS **faces, const int &indiceMaxFaces, const int &nombreNoeudsInternes)
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
  if (m_numNoeuds[2] < nombreNoeudsInternes && m_numNoeuds[3] < nombreNoeudsInternes)
  {
    FaceSegment face(m_numNoeuds[2], m_numNoeuds[3]);
    if (face.faceExiste(faces, indiceMaxFaces, indiceFaceExiste))
    {
      faces[indiceFaceExiste]->ajouteElementVoisinLimite(this);
      faces[indiceFaceExiste]->setEstComm(true);
    }
  }
  //Verification face 4 :
  if (m_numNoeuds[3] < nombreNoeudsInternes && m_numNoeuds[0] < nombreNoeudsInternes)
  {
    FaceSegment face(m_numNoeuds[3], m_numNoeuds[0]);
    if (face.faceExiste(faces, indiceMaxFaces, indiceFaceExiste))
    {
      faces[indiceFaceExiste]->ajouteElementVoisinLimite(this);
      faces[indiceFaceExiste]->setEstComm(true);
    }
  }
}

//***********************************************************************

int ElementQuadrangle::compteFaceCommunicante(vector<int*> &facesTemp, vector<int> &sommeNoeudsTemp)
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
      case 2: face[0] = m_numNoeuds[2]; face[1] = m_numNoeuds[3]; break;     
      case 3: face[0] = m_numNoeuds[3]; face[1] = m_numNoeuds[0]; break;     
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
int ElementQuadrangle::compteFaceCommunicante(int &iMax, int **facesTemp, int *sommeNoeudsTemp)
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
      case 2: face[0] = m_numNoeuds[2]; face[1] = m_numNoeuds[3]; break;     
      case 3: face[0] = m_numNoeuds[3]; face[1] = m_numNoeuds[0]; break;     
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