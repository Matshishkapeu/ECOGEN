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

#include <iostream>
#include <algorithm>
#include <sstream>

#include "MaillageCartesien.h"

using namespace std;

//***********************************************************************

MaillageCartesien::MaillageCartesien(double lX, int nombreMaillesX,
  double lY, int nombreMaillesY,
  double lZ, int nombreMaillesZ) :
  m_lX(lX), m_nombreMaillesXGlobal(nombreMaillesX),
  m_lY(lY), m_nombreMaillesYGlobal(nombreMaillesY),
  m_lZ(lZ), m_nombreMaillesZGlobal(nombreMaillesZ)
{
  m_nombreCellulesCalcul = m_nombreMaillesXGlobal*m_nombreMaillesYGlobal*m_nombreMaillesZGlobal;
  m_dX = m_lX / static_cast<double>(m_nombreMaillesXGlobal);
  m_dY = m_lY / static_cast<double>(m_nombreMaillesYGlobal);
  m_dZ = m_lZ / static_cast<double>(m_nombreMaillesZGlobal);
  m_lCFL = 1.e10;
  m_geometrie = 0;
  if (nombreMaillesX != 1) { m_lCFL = min(m_lCFL, m_dX); m_geometrie += 1; }
  if (nombreMaillesY != 1) { m_lCFL = min(m_lCFL, m_dY); m_geometrie += 1; }
  if (nombreMaillesZ != 1) { m_lCFL = min(m_lCFL, m_dZ); m_geometrie += 1; }
  if (m_geometrie > 1) m_lCFL *= 0.6;
  m_volume = m_dX*m_dY*m_dZ;
  m_numNoeudDeb = 0;
  m_numNoeudFin = m_nombreMaillesXGlobal;
  m_type = REC;
}

//***********************************************************************

MaillageCartesien::~MaillageCartesien(){
  delete[] m_elements;
  delete[] m_faces;
  delete m_limXm;
  delete m_limXp;
  delete m_limYm;
  delete m_limYp;
  delete m_limZm;
  delete m_limZp;
}

//***********************************************************************

void MaillageCartesien::attributLimites(std::vector<CondLim*> &condLim)
{
  for (unsigned int i = 0; i < condLim.size(); i++) {
    switch (i) {
    case 0:
      m_limXm = condLim[i]; break;
    case 1:
      m_limXp = condLim[i]; break;
    case 2:
      m_limYm = condLim[i]; break;
    case 3:
      m_limYp = condLim[i]; break;
    case 4:
      m_limZm = condLim[i]; break;
    case 5:
      m_limZp = condLim[i]; break;
    default:
      Erreurs::messageErreur("Probleme de limites dans attributLimites"); break;
    }
  }

  for (int j = condLim.size(); j < 6; j++) {
    switch (j) {
    case 0:
      m_limXm = new CondLimAbs(j+1); break;
    case 1:
      m_limXp = new CondLimAbs(j+1); break;
    case 2:
      m_limYm = new CondLimAbs(j+1); break;
    case 3:
      m_limYp = new CondLimAbs(j+1); break;
    case 4:
      m_limZm = new CondLimAbs(j+1); break;
    case 5:
      m_limZp = new CondLimAbs(j+1); break;
    default:
      Erreurs::messageErreur("Probleme de limites dans attributLimites"); break;
    }
  }
}

//***********************************************************************

void MaillageCartesien::recupereIJK(const int &indice, int &i, int &j, int &k) const
{
  int reste;
  k = indice / (m_nombreMaillesX*m_nombreMaillesY);
  reste = indice % (m_nombreMaillesX*m_nombreMaillesY);
  j = reste / m_nombreMaillesX;
  i = reste % m_nombreMaillesX;
}

//***********************************************************************

void MaillageCartesien::construitIGlobal(const int &i, const int &j, const int &k, int &indice) const
{
  indice = 0;
  if (m_nombreMaillesX != 1) indice += i;
  if (m_nombreMaillesY != 1) indice += j*m_nombreMaillesX;
  if (m_nombreMaillesZ != 1) indice += k*m_nombreMaillesX*m_nombreMaillesY;
}

//***********************************************************************

int MaillageCartesien::initialiseGeometrie(Cellule ***cellules, BordDeMaille ***bord, bool pretraitementParallele, string ordreCalcul)
{
  if (Ncpu == 1)
  {
    this->initialiseGeometrieMonoCPU(cellules, bord, ordreCalcul);
  }
  else
  {
    this->initialiseGeometrieParallele(cellules, bord, ordreCalcul);
  }
  return m_geometrie;
}

//***********************************************************************

void MaillageCartesien::initialiseGeometrieMonoCPU(Cellule ***cellules, BordDeMaille ***bord, string ordreCalcul)
{
  int ix, iy, iz;
  int compteMaillesParallele(0);

  m_nombreMaillesX = m_nombreMaillesXGlobal;
  m_nombreMaillesY = m_nombreMaillesYGlobal;
  m_nombreMaillesZ = m_nombreMaillesZGlobal;
  m_nombreCellulesTotal = m_nombreCellulesCalcul;

  //Declaration du tableau de mailles et d'elements
  //-----------------------------------------------
  m_elements = new ElementCartesien[m_nombreCellulesCalcul];
  (*cellules) = new Cellule*[m_nombreCellulesCalcul];
  for (int i = 0; i < m_nombreCellulesCalcul; i++) { 
    if(ordreCalcul == "ORDRE1") { (*cellules)[i] = new Cellule; }
    else { (*cellules)[i] = new CelluleO2; }
    (*cellules)[i]->setElement(&m_elements[i], i); 
  }

  //Attribution des donnees geometriques cartesiennes aux mailles
  //-------------------------------------------------------------
  Coord tangente, normale, binormale;
  double surface(1.);
  double posX, posY, posZ;
  for (int i = 0; i < m_nombreCellulesCalcul; i++)
  {
    (*cellules)[i]->getElement()->setVolume(m_volume);
    (*cellules)[i]->getElement()->setLCFL(m_lCFL);
    this->recupereIJK(i, ix, iy, iz);
    //cout << i << " " << ix << " " << iy << " " << iz << endl;
    posX = (static_cast<double>(ix)+0.5)*m_dX;
    posY = (static_cast<double>(iy)+0.5)*m_dY;
    posZ = (static_cast<double>(iz)+0.5)*m_dZ;
    (*cellules)[i]->getElement()->setPos(posX, posY, posZ);
  }

  //Attribution des donnees geometriques cartesiennes aux bords de mailles
  //----------------------------------------------------------------------
  //Determination du nombre de faces selon calcul 1D/2D/3D
  //m_nombreFaces : nombre total de faces
  //m_nombreFacesLimites : nombre de faces sur les limites
  //nombreFacesInternes : nombre de faces communes a deux mailles
  m_nombreFacesTotal = 0;
  int m_nombreFacesLimites(0);
  if (m_nombreMaillesX != 1)
  {
    m_nombreFacesTotal += (m_nombreMaillesX + 1)*m_nombreMaillesY*m_nombreMaillesZ;
    m_nombreFacesLimites += 2 * m_nombreMaillesY*m_nombreMaillesZ;
  }
  if (m_nombreMaillesY != 1)
  {
    m_nombreFacesTotal += (m_nombreMaillesY + 1)*m_nombreMaillesX*m_nombreMaillesZ;
    m_nombreFacesLimites += 2 * m_nombreMaillesX*m_nombreMaillesZ;
  }
  if (m_nombreMaillesZ != 1)
  {
    m_nombreFacesTotal += (m_nombreMaillesZ + 1)*m_nombreMaillesX*m_nombreMaillesY;
    m_nombreFacesLimites += 2 * m_nombreMaillesX*m_nombreMaillesY;
  }
  int nombreFacesInternes(m_nombreFacesTotal - m_nombreFacesLimites);
  //Allocation globale
  (*bord) = new BordDeMaille*[m_nombreFacesTotal];
  m_faces = new FaceCartesien[m_nombreFacesTotal];

  //Initialisation des faces internes
  //*********************************
  int iMailleG, iMailleD, iFace(0), iTemp;
  Cellule *BGM=0, *BGP=0, *BDM=0, *BDP=0;
  //Faces selon X
  surface = m_dY*m_dZ;
  tangente.setXYZ(0., 1., 0.); normale.setXYZ(1., 0., 0.); binormale.setXYZ(0., 0., 1.);
  for (ix = 0; ix < m_nombreMaillesX - 1; ix++)
  {
    for (iy = 0; iy < m_nombreMaillesY; iy++)
    {
      for (iz = 0; iz < m_nombreMaillesZ; iz++)
      {
        if(ordreCalcul == "ORDRE1") { (*bord)[iFace] = new BordDeMaille; }
        else { (*bord)[iFace] = new BordDeMailleO2; }
        (*bord)[iFace]->setFace(&m_faces[iFace]);
        this->construitIGlobal(ix, iy, iz, iMailleG);
        this->construitIGlobal(ix + 1, iy, iz, iMailleD);
        (*bord)[iFace]->initialise((*cellules)[iMailleG], (*cellules)[iMailleD]);
        (*cellules)[iMailleG]->ajouteBord((*bord)[iFace]);
        (*cellules)[iMailleD]->ajouteBord((*bord)[iFace]);
        m_faces[iFace].initialiseAutres(surface, normale, tangente, binormale);
        posX = (static_cast<double>(ix) + 1.)*m_dX;
        posY = (static_cast<double>(iy) + 0.5)*m_dY;
        posZ = (static_cast<double>(iz) + 0.5)*m_dZ;
        m_faces[iFace].setPos(posX, posY, posZ);
        //Preparation des cellules pour ordre 2 multipentes
        BGM = 0; BGP = 0; BDM = 0; BDP = 0;
        if (ix != 0) { this->construitIGlobal(ix - 1, iy, iz, iTemp); BGM = (*cellules)[iTemp]; }
        if (ix < m_nombreMaillesX - 1) { this->construitIGlobal(ix + 1, iy, iz, iTemp); BGP = (*cellules)[iTemp]; }
        if (ix < m_nombreMaillesX - 2) { this->construitIGlobal(ix + 2, iy, iz, iTemp); BDM = (*cellules)[iTemp]; }
        this->construitIGlobal(ix, iy, iz, iTemp); BDP = (*cellules)[iTemp];
        //(*bord)[iFace]->setB(BG1M, BGM);
        //(*bord)[iFace]->setB(BG1P, BGP);
        //(*bord)[iFace]->setB(BD1M, BDM);
        //(*bord)[iFace]->setB(BD1P, BDP);
        //(*bord)[iFace]->setDistanceH(distanceHGM, m_dX);
        //(*bord)[iFace]->setDistanceH(distanceHGP, m_dX);
        //(*bord)[iFace]->setDistanceH(distanceHDM, m_dX);
        //(*bord)[iFace]->setDistanceH(distanceHDP, m_dX);
        iFace++;
      }
    }
  }
  //Faces selon Y
  surface = m_dX*m_dZ;
  tangente.setXYZ(-1., 0., 0.); normale.setXYZ(0., 1., 0.); binormale.setXYZ(0., 0., 1.);
  for (ix = 0; ix < m_nombreMaillesX; ix++)
  {
    for (iy = 0; iy < m_nombreMaillesY - 1; iy++)
    {
      for (iz = 0; iz < m_nombreMaillesZ; iz++)
      {
        if(ordreCalcul == "ORDRE1") { (*bord)[iFace] = new BordDeMaille; }
        else { (*bord)[iFace] = new BordDeMailleO2; }
        (*bord)[iFace]->setFace(&m_faces[iFace]);
        this->construitIGlobal(ix, iy, iz, iMailleG);
        this->construitIGlobal(ix, iy + 1, iz, iMailleD);
        (*bord)[iFace]->initialise((*cellules)[iMailleG], (*cellules)[iMailleD]);
        (*cellules)[iMailleG]->ajouteBord((*bord)[iFace]);
        (*cellules)[iMailleD]->ajouteBord((*bord)[iFace]);
        m_faces[iFace].initialiseAutres(surface, normale, tangente, binormale);
        posX = (static_cast<double>(ix) + 0.5)*m_dX;
        posY = (static_cast<double>(iy) + 1.)*m_dY;
        posZ = (static_cast<double>(iz) + 0.5)*m_dZ;
        m_faces[iFace].setPos(posX, posY, posZ);
        //Preparation des cellules pour ordre 2 multipentes
        BGM = 0; BGP = 0; BDM = 0; BDP = 0;
        if (iy != 0) { this->construitIGlobal(ix, iy - 1, iz, iTemp); BGM = (*cellules)[iTemp]; }
        if (iy < m_nombreMaillesY - 1) { this->construitIGlobal(ix, iy + 1, iz, iTemp); BGP = (*cellules)[iTemp]; }
        if (iy < m_nombreMaillesY - 2) { this->construitIGlobal(ix, iy + 2, iz, iTemp); BDM = (*cellules)[iTemp]; }
        this->construitIGlobal(ix, iy, iz, iTemp); BDP = (*cellules)[iTemp];
        //(*bord)[iFace]->setB(BG1M, BGM);
        //(*bord)[iFace]->setB(BG1P, BGP);
        //(*bord)[iFace]->setB(BD1M, BDM);
        //(*bord)[iFace]->setB(BD1P, BDP);
        //(*bord)[iFace]->setDistanceH(distanceHGM, m_dY);
        //(*bord)[iFace]->setDistanceH(distanceHGP, m_dY);
        //(*bord)[iFace]->setDistanceH(distanceHDM, m_dY);
        //(*bord)[iFace]->setDistanceH(distanceHDP, m_dY);
        iFace++;
      }
    }
  }
  //Faces selon Z
  surface = m_dX*m_dY;
  tangente.setXYZ(1., 0., 0.); normale.setXYZ(0., 0., 1.); binormale.setXYZ(0., 1., 0.);
  for (ix = 0; ix < m_nombreMaillesX; ix++)
  {
    for (iy = 0; iy < m_nombreMaillesY; iy++)
    {
      for (iz = 0; iz < m_nombreMaillesZ - 1; iz++)
      {
        if(ordreCalcul == "ORDRE1") { (*bord)[iFace] = new BordDeMaille; }
        else { (*bord)[iFace] = new BordDeMailleO2; }
        (*bord)[iFace]->setFace(&m_faces[iFace]);
        this->construitIGlobal(ix, iy, iz, iMailleG);
        this->construitIGlobal(ix, iy, iz + 1, iMailleD);
        (*bord)[iFace]->initialise((*cellules)[iMailleG], (*cellules)[iMailleD]);
        (*cellules)[iMailleG]->ajouteBord((*bord)[iFace]);
        (*cellules)[iMailleD]->ajouteBord((*bord)[iFace]);
        m_faces[iFace].initialiseAutres(surface, normale, tangente, binormale);
        posX = (static_cast<double>(ix) + 0.5)*m_dX;
        posY = (static_cast<double>(iy) + 0.5)*m_dY;
        posZ = (static_cast<double>(iz) + 1.)*m_dZ;
        m_faces[iFace].setPos(posX, posY, posZ);
        //Preparation des cellules pour ordre 2 multipentes
        BGM = 0; BGP = 0; BDM = 0; BDP = 0;
        if (iz != 0) { this->construitIGlobal(ix, iy, iz - 1, iTemp); BGM = (*cellules)[iTemp]; }
        if (iz < m_nombreMaillesZ - 1) { this->construitIGlobal(ix, iy, iz + 1, iTemp); BGP = (*cellules)[iTemp]; }
        if (iz < m_nombreMaillesZ - 2) { this->construitIGlobal(ix, iy, iz + 2, iTemp); BDM = (*cellules)[iTemp]; }
        this->construitIGlobal(ix, iy, iz, iTemp); BDP = (*cellules)[iTemp];
        //(*bord)[iFace]->setB(BG1M, BGM);
        //(*bord)[iFace]->setB(BG1P, BGP);
        //(*bord)[iFace]->setB(BD1M, BDM);
        //(*bord)[iFace]->setB(BD1P, BDP);
        //(*bord)[iFace]->setDistanceH(distanceHGM, m_dZ);
        //(*bord)[iFace]->setDistanceH(distanceHGP, m_dZ);
        //(*bord)[iFace]->setDistanceH(distanceHDM, m_dZ);
        //(*bord)[iFace]->setDistanceH(distanceHDP, m_dZ);
        iFace++;
      }
    }
  }

  //Initialisation des faces limites en X
  //---------------------------------------
  if (m_nombreMaillesX != 1)
  {
    surface = m_dY*m_dZ;
    //Limite X=0
    ix = 0;
    tangente.setXYZ(0., -1., 0.); normale.setXYZ(-1., 0., 0.); binormale.setXYZ(0., 0., 1.);
    for (iy = 0; iy < m_nombreMaillesY; iy++)
    {
      for (iz = 0; iz < m_nombreMaillesZ; iz++)
      {
        m_limXm->creeLimite(&(*bord)[iFace]);
        (*bord)[iFace]->setFace(&m_faces[iFace]);
        this->construitIGlobal(ix, iy, iz, iMailleG);
        iMailleD = iMailleG;
        (*bord)[iFace]->initialise((*cellules)[iMailleG], (*cellules)[iMailleD]);
        (*cellules)[iMailleG]->ajouteBord((*bord)[iFace]);
        m_faces[iFace].initialiseAutres(surface, normale, tangente, binormale); //FP//Q// Interet d'avoir deux fonctions initialise ?
        posX = 0.;
        posY = (static_cast<double>(iy) + 0.5)*m_dY;
        posZ = (static_cast<double>(iz) + 0.5)*m_dZ;
        m_faces[iFace].setPos(posX, posY, posZ);
        iFace++;
      }
    }
    //Limite X=lX
    ix = m_nombreMaillesX - 1;
    tangente.setXYZ(0., 1., 0.); normale.setXYZ(1., 0., 0.); binormale.setXYZ(0., 0., 1.);
    for (iy = 0; iy < m_nombreMaillesY; iy++)
    {
      for (iz = 0; iz < m_nombreMaillesZ; iz++)
      {
        m_limXp->creeLimite(&(*bord)[iFace]);
        (*bord)[iFace]->setFace(&m_faces[iFace]);
        this->construitIGlobal(ix, iy, iz, iMailleG);
        iMailleD = iMailleG;
        (*bord)[iFace]->initialise((*cellules)[iMailleG], (*cellules)[iMailleD]);
        (*cellules)[iMailleG]->ajouteBord((*bord)[iFace]);
        m_faces[iFace].initialiseAutres(surface, normale, tangente, binormale);
        posX = (static_cast<double>(ix) + 1.)*m_dX;
        posY = (static_cast<double>(iy) + 0.5)*m_dY;
        posZ = (static_cast<double>(iz) + 0.5)*m_dZ;
        m_faces[iFace].setPos(posX, posY, posZ);
        iFace++;
      }
    }
  }
  //Initialisation des faces limites en Y
  //---------------------------------------
  if (m_nombreMaillesY != 1)
  {
    surface = m_dX*m_dZ;
    //Limite Y=0
    iy = 0;
    tangente.setXYZ(1., 0., 0.); normale.setXYZ(0., -1., 0.); binormale.setXYZ(0., 0., 1.);
    for (ix = 0; ix < m_nombreMaillesX; ix++)
    {
      for (iz = 0; iz < m_nombreMaillesZ; iz++)
      {
        m_limYm->creeLimite(&(*bord)[iFace]);
        (*bord)[iFace]->setFace(&m_faces[iFace]);
        this->construitIGlobal(ix, iy, iz, iMailleG);
        iMailleD = iMailleG;
        (*bord)[iFace]->initialise((*cellules)[iMailleG], (*cellules)[iMailleD]);
        (*cellules)[iMailleG]->ajouteBord((*bord)[iFace]);
        m_faces[iFace].initialiseAutres(surface, normale, tangente, binormale);
        posX = (static_cast<double>(ix) + 0.5)*m_dX;
        posY = 0;
        posZ = (static_cast<double>(iz) + 0.5)*m_dZ;
        m_faces[iFace].setPos(posX, posY, posZ);
        iFace++;
      }
    }
    //Limite Y=lY
    iy = m_nombreMaillesY - 1;
    tangente.setXYZ(-1., 0., 0.); normale.setXYZ(0., 1., 0.); binormale.setXYZ(0., 0., 1.);
    for (ix = 0; ix < m_nombreMaillesX; ix++)
    {
      for (iz = 0; iz < m_nombreMaillesZ; iz++)
      {
        m_limYp->creeLimite(&(*bord)[iFace]);
        (*bord)[iFace]->setFace(&m_faces[iFace]);
        this->construitIGlobal(ix, iy, iz, iMailleG);
        iMailleD = iMailleG;
        (*bord)[iFace]->initialise((*cellules)[iMailleG], (*cellules)[iMailleD]);
        (*cellules)[iMailleG]->ajouteBord((*bord)[iFace]);
        m_faces[iFace].initialiseAutres(surface, normale, tangente, binormale);
        posX = (static_cast<double>(ix) + 0.5)*m_dX;
        posY = (static_cast<double>(iy) + 1.)*m_dY;
        posZ = (static_cast<double>(iz) + 0.5)*m_dZ;
        m_faces[iFace].setPos(posX, posY, posZ);
        iFace++;
      }
    }
  }
  //Initialisation des faces limites en Z
  //-------------------------------------
  if (m_nombreMaillesZ != 1)
  {
    surface = m_dX*m_dY;
    //Limite Z=0
    iz = 0;
    tangente.setXYZ(-1., 0., 0.); normale.setXYZ(0., 0., -1.); binormale.setXYZ(0., 1., 0.);
    for (ix = 0; ix < m_nombreMaillesX; ix++)
    {
      for (iy = 0; iy < m_nombreMaillesY; iy++)
      {
        m_limZm->creeLimite(&(*bord)[iFace]);
        (*bord)[iFace]->setFace(&m_faces[iFace]);
        this->construitIGlobal(ix, iy, iz, iMailleG);
        iMailleD = iMailleG;
        (*bord)[iFace]->initialise((*cellules)[iMailleG], (*cellules)[iMailleD]);
        (*cellules)[iMailleG]->ajouteBord((*bord)[iFace]);
        m_faces[iFace].initialiseAutres(surface, normale, tangente, binormale);
        posX = (static_cast<double>(ix) + 0.5)*m_dX;
        posY = (static_cast<double>(iy) + 0.5)*m_dY;
        posZ = 0.;
        m_faces[iFace].setPos(posX, posY, posZ);
        iFace++;
      }
    }
    //Limite Z=lZ
    iz = m_nombreMaillesZ - 1;
    tangente.setXYZ(1., 0., 0.); normale.setXYZ(0., 0., 1.); binormale.setXYZ(0., 1., 0.);
    for (ix = 0; ix < m_nombreMaillesX; ix++)
    {
      for (iy = 0; iy < m_nombreMaillesY; iy++)
      {
        m_limZp->creeLimite(&(*bord)[iFace]);
        (*bord)[iFace]->setFace(&m_faces[iFace]);
        this->construitIGlobal(ix, iy, iz, iMailleG);
        iMailleD = iMailleG;
        (*bord)[iFace]->initialise((*cellules)[iMailleG], (*cellules)[iMailleD]);
        (*cellules)[iMailleG]->ajouteBord((*bord)[iFace]);
        m_faces[iFace].initialiseAutres(surface, normale, tangente, binormale);
        posX = (static_cast<double>(ix) + 0.5)*m_dX;
        posY = (static_cast<double>(iy) + 0.5)*m_dY;
        posZ = (static_cast<double>(iz) + 1.)*m_dZ;
        m_faces[iFace].setPos(posX, posY, posZ);
        iFace++;
      }
    }
  }
}

//***********************************************************************

void MaillageCartesien::initialiseGeometrieParallele(Cellule ***cellules, BordDeMaille ***bord, string ordreCalcul)
{
  int ix, iy, iz;
  int compteMaillesParallele(0);

  m_nombreMaillesX = m_nombreMaillesXGlobal;
  m_nombreMaillesY = m_nombreMaillesYGlobal;
  m_nombreMaillesZ = m_nombreMaillesZGlobal;
  
  //Declaration du tableau de mailles
  //---------------------------------
  //Decoupage parallele du domaine geometrique
  this->decoupageParallele();
  m_elements = new ElementCartesien[m_nombreCellulesTotal];
  (*cellules) = new Cellule*[m_nombreCellulesTotal];
  for (int i = 0; i < m_nombreCellulesCalcul; i++) {
	  if (ordreCalcul == "ORDRE1") { (*cellules)[i] = new Cellule; }
	  else { (*cellules)[i] = new CelluleO2; }
    (*cellules)[i]->setElement(&m_elements[i], i);
  }
	for (int i = m_nombreCellulesCalcul; i < m_nombreCellulesTotal; i++) {
		if (ordreCalcul == "ORDRE1") { (*cellules)[i] = new Cellule; }
		else { (*cellules)[i] = new CelluleO2Ghost; }
		(*cellules)[i]->setElement(&m_elements[i], i);
	}
  
  //Attribution des donnees geometriques cartesiennes aux mailles
  //-------------------------------------------------------------
  Coord tangente, normale, binormale;
	Coord tangenteGhost, normaleGhost, binormaleGhost;
  double surface(1.);
  double posX, posY, posZ;
  for (int i = 0; i < m_nombreCellulesTotal; i++)
  {
    (*cellules)[i]->getElement()->setVolume(m_volume);
    (*cellules)[i]->getElement()->setLCFL(m_lCFL);
    this->recupereIJK(i, ix, iy, iz);
    //cout << i << " " << ix << " " << iy << " " << iz << endl;
    posX = (static_cast<double>(ix)+0.5)*m_dX + m_origine.getX();
    posY = (static_cast<double>(iy)+0.5)*m_dY + m_origine.getY();
    posZ = (static_cast<double>(iz)+0.5)*m_dZ + m_origine.getZ();
    (*cellules)[i]->getElement()->setPos(posX, posY, posZ);
  }

  //Attribution des donnees geometriques cartesiennes aux bords de mailles
  //----------------------------------------------------------------------
  //Determination du nombre de faces selon calcul 1D/2D/3D
  //m_nombreFaces : nombre total de faces
  //m_nombreFacesLimites : nombre de faces sur les limites
  //nombreFacesInternes : nombre de faces communes a deux mailles
  m_nombreFacesTotal = 0;
  int m_nombreFacesLimites(0);
  if (m_nombreMaillesX != 1)
  {
    m_nombreFacesTotal += (m_nombreMaillesX + 1)*m_nombreMaillesY*m_nombreMaillesZ;
    m_nombreFacesLimites += 2 * m_nombreMaillesY*m_nombreMaillesZ;
  }
  if (m_nombreMaillesY != 1)
  {
    m_nombreFacesTotal += (m_nombreMaillesY + 1)*m_nombreMaillesX*m_nombreMaillesZ;
    m_nombreFacesLimites += 2 * m_nombreMaillesX*m_nombreMaillesZ;
  }
  if (m_nombreMaillesZ != 1)
  {
    m_nombreFacesTotal += (m_nombreMaillesZ + 1)*m_nombreMaillesX*m_nombreMaillesY;
    m_nombreFacesLimites += 2 * m_nombreMaillesX*m_nombreMaillesY;
  }
  int nombreFacesInternes(m_nombreFacesTotal - m_nombreFacesLimites);
  //Allocation globale
  (*bord) = new BordDeMaille*[m_nombreFacesTotal];
  m_faces = new FaceCartesien[m_nombreFacesTotal];
  
  //Initialisation des faces internes
  //*********************************
  int iMailleG, iMailleD, iFace(0), iTemp;
  Cellule *BGM = 0, *BGP = 0, *BDM = 0, *BDP = 0;
  //Faces selon X
  surface = m_dY*m_dZ;
  tangente.setXYZ(0., 1., 0.); normale.setXYZ(1., 0., 0.); binormale.setXYZ(0., 0., 1.);
  for (ix = 0; ix < m_nombreMaillesX - 1; ix++)
  {
    for (iy = 0; iy < m_nombreMaillesY; iy++)
    {
      for (iz = 0; iz < m_nombreMaillesZ; iz++)
      {
        if(ordreCalcul == "ORDRE1") { (*bord)[iFace] = new BordDeMaille; }
        else { (*bord)[iFace] = new BordDeMailleO2; }
        (*bord)[iFace]->setFace(&m_faces[iFace]);
        this->construitIGlobal(ix, iy, iz, iMailleG);
        this->construitIGlobal(ix + 1, iy, iz, iMailleD);
        (*bord)[iFace]->initialise((*cellules)[iMailleG], (*cellules)[iMailleD]);
        (*cellules)[iMailleG]->ajouteBord((*bord)[iFace]);
        (*cellules)[iMailleD]->ajouteBord((*bord)[iFace]);
        m_faces[iFace].initialiseAutres(surface, normale, tangente, binormale);
        posX = (static_cast<double>(ix) + 1.)*m_dX + m_origine.getX();
        posY = (static_cast<double>(iy) + 0.5)*m_dY + m_origine.getY();
        posZ = (static_cast<double>(iz) + 0.5)*m_dZ + m_origine.getZ();
        m_faces[iFace].setPos(posX, posY, posZ);
        //Preparation des cellules pour ordre 2 multipentes
        BGM = 0; BGP = 0; BDM = 0; BDP = 0;
        if (ix != 0) { this->construitIGlobal(ix - 1, iy, iz, iTemp); BGM = (*cellules)[iTemp]; }
        if (ix < m_nombreMaillesX - 1) { this->construitIGlobal(ix + 1, iy, iz, iTemp); BGP = (*cellules)[iTemp]; }
        if (ix < m_nombreMaillesX - 2) { this->construitIGlobal(ix + 2, iy, iz, iTemp); BDM = (*cellules)[iTemp]; }
        this->construitIGlobal(ix, iy, iz, iTemp); BDP = (*cellules)[iTemp];
        //(*bord)[iFace]->setB(BG1M, BGM);
        //(*bord)[iFace]->setB(BG1P, BGP);
        //(*bord)[iFace]->setB(BD1M, BDM);
        //(*bord)[iFace]->setB(BD1P, BDP);
        //(*bord)[iFace]->setDistanceH(distanceHGM, m_dX);
        //(*bord)[iFace]->setDistanceH(distanceHGP, m_dX);
        //(*bord)[iFace]->setDistanceH(distanceHDM, m_dX);
        //(*bord)[iFace]->setDistanceH(distanceHDP, m_dX);
        iFace++;
      }
    }
  }
  //Faces selon Y
  surface = m_dX*m_dZ;
  tangente.setXYZ(-1., 0., 0.); normale.setXYZ(0., 1., 0.); binormale.setXYZ(0., 0., 1.);
  for (ix = 0; ix < m_nombreMaillesX; ix++)
  {
    for (iy = 0; iy < m_nombreMaillesY - 1; iy++)
    {
      for (iz = 0; iz < m_nombreMaillesZ; iz++)
      {
        if(ordreCalcul == "ORDRE1") { (*bord)[iFace] = new BordDeMaille; }
        else { (*bord)[iFace] = new BordDeMailleO2; }
        (*bord)[iFace]->setFace(&m_faces[iFace]);
        this->construitIGlobal(ix, iy, iz, iMailleG);
        this->construitIGlobal(ix, iy + 1, iz, iMailleD);
        (*bord)[iFace]->initialise((*cellules)[iMailleG], (*cellules)[iMailleD]);
        (*cellules)[iMailleG]->ajouteBord((*bord)[iFace]);
        (*cellules)[iMailleD]->ajouteBord((*bord)[iFace]);
        m_faces[iFace].initialiseAutres(surface, normale, tangente, binormale);
        posX = (static_cast<double>(ix) + 0.5)*m_dX + m_origine.getX();
        posY = (static_cast<double>(iy) + 1.)*m_dY + m_origine.getY();
        posZ = (static_cast<double>(iz) + 0.5)*m_dZ + m_origine.getZ();
        m_faces[iFace].setPos(posX, posY, posZ);
        //Preparation des cellules pour ordre 2 multipentes
        BGM = 0; BGP = 0; BDM = 0; BDP = 0;
        if (iy != 0) { this->construitIGlobal(ix, iy - 1, iz, iTemp); BGM = (*cellules)[iTemp]; }
        if (iy < m_nombreMaillesY - 1) { this->construitIGlobal(ix, iy + 1, iz, iTemp); BGP = (*cellules)[iTemp]; }
        if (iy < m_nombreMaillesY - 2) { this->construitIGlobal(ix, iy + 2, iz, iTemp); BDM = (*cellules)[iTemp]; }
        this->construitIGlobal(ix, iy, iz, iTemp); BDP = (*cellules)[iTemp];
        //(*bord)[iFace]->setB(BG1M, BGM);
        //(*bord)[iFace]->setB(BG1P, BGP);
        //(*bord)[iFace]->setB(BD1M, BDM);
        //(*bord)[iFace]->setB(BD1P, BDP);
        //(*bord)[iFace]->setDistanceH(distanceHGM, m_dY);
        //(*bord)[iFace]->setDistanceH(distanceHGP, m_dY);
        //(*bord)[iFace]->setDistanceH(distanceHDM, m_dY);
        //(*bord)[iFace]->setDistanceH(distanceHDP, m_dY);
        iFace++;
      }
    }
  }
  //Faces selon Z
  surface = m_dX*m_dY;
  tangente.setXYZ(1., 0., 0.); normale.setXYZ(0., 0., 1.); binormale.setXYZ(0., 1., 0.);
  for (ix = 0; ix < m_nombreMaillesX; ix++)
  {
    for (iy = 0; iy < m_nombreMaillesY; iy++)
    {
      for (iz = 0; iz < m_nombreMaillesZ - 1; iz++)
      {
        if(ordreCalcul == "ORDRE1") { (*bord)[iFace] = new BordDeMaille; }
        else { (*bord)[iFace] = new BordDeMailleO2; }
        (*bord)[iFace]->setFace(&m_faces[iFace]);
        this->construitIGlobal(ix, iy, iz, iMailleG);
        this->construitIGlobal(ix, iy, iz + 1, iMailleD);
        (*bord)[iFace]->initialise((*cellules)[iMailleG], (*cellules)[iMailleD]);
        (*cellules)[iMailleG]->ajouteBord((*bord)[iFace]);
        (*cellules)[iMailleD]->ajouteBord((*bord)[iFace]);
        m_faces[iFace].initialiseAutres(surface, normale, tangente, binormale);
        posX = (static_cast<double>(ix) + 0.5)*m_dX + m_origine.getX();
        posY = (static_cast<double>(iy) + 0.5)*m_dY + m_origine.getY();
        posZ = (static_cast<double>(iz) + 1.)*m_dZ + m_origine.getZ();
        m_faces[iFace].setPos(posX, posY, posZ);
        //Preparation des cellules pour ordre 2 multipentes
        BGM = 0; BGP = 0; BDM = 0; BDP = 0;
        if (iz != 0) { this->construitIGlobal(ix, iy, iz - 1, iTemp); BGM = (*cellules)[iTemp]; }
        if (iz < m_nombreMaillesZ - 1) { this->construitIGlobal(ix, iy, iz + 1, iTemp); BGP = (*cellules)[iTemp]; }
        if (iz < m_nombreMaillesZ - 2) { this->construitIGlobal(ix, iy, iz + 2, iTemp); BDM = (*cellules)[iTemp]; }
        this->construitIGlobal(ix, iy, iz, iTemp); BDP = (*cellules)[iTemp];
        //(*bord)[iFace]->setB(BG1M, BGM);
        //(*bord)[iFace]->setB(BG1P, BGP);
        //(*bord)[iFace]->setB(BD1M, BDM);
        //(*bord)[iFace]->setB(BD1P, BDP);
        //(*bord)[iFace]->setDistanceH(distanceHGM, m_dZ);
        //(*bord)[iFace]->setDistanceH(distanceHGP, m_dZ);
        //(*bord)[iFace]->setDistanceH(distanceHDM, m_dZ);
        //(*bord)[iFace]->setDistanceH(distanceHDP, m_dZ);
        iFace++;
      }
    }
  }
  //Calcul_Parallele.arreterCode();
  //Utile pour le parallele, ne change rien en mono-proc.
  compteMaillesParallele = m_nombreCellulesCalcul;

  //PROBLEME ADRESSE MULTI-D..... //FP//Q//?

  //Initialisation des faces limites en X
  //-------------------------------------
  if (m_nombreMaillesX != 1)
  {
    surface = m_dY*m_dZ;
    //Limite X=0
    ix = 0;
    tangente.setXYZ(0., -1., 0.); normale.setXYZ(-1., 0., 0.); binormale.setXYZ(0., 0., 1.);
    for (iy = 0; iy < m_nombreMaillesY; iy++)
    {
      for (iz = 0; iz < m_nombreMaillesZ; iz++)
      {
        //Traitement des limites communicantes
        if (m_decoupe == 'X' && rang > 0) {
					tangenteGhost.setXYZ(0., 1., 0.); normaleGhost.setXYZ(1., 0., 0.); binormaleGhost.setXYZ(0., 0., 1.);
					this->construitIGlobal(ix, iy, iz, iMailleD);
          iMailleG = compteMaillesParallele++;
          if (ordreCalcul == "ORDRE1") { (*bord)[iFace] = new BordDeMaille; }
          else { (*bord)[iFace] = new BordDeMailleO2; }
					(*cellules)[iMailleG]->getElement()->setPos((*cellules)[iMailleD]->getElement()->getPosition());
          (*cellules)[iMailleG]->getElement()->setPosX((*cellules)[iMailleD]->getElement()->getPosition().getX() - m_dX);
          (*cellules)[iMailleG]->ajouteBord((*bord)[iFace]);
					(*bord)[iFace]->setFace(&m_faces[iFace]);
					(*bord)[iFace]->initialise((*cellules)[iMailleG], (*cellules)[iMailleD]);
					(*cellules)[iMailleD]->ajouteBord((*bord)[iFace]);
					m_faces[iFace].initialiseAutres(surface, normaleGhost, tangenteGhost, binormaleGhost);
        } //On prend parmi les mailles fantomes dans l'ordre
        //Traitement des limites physiques
        else {
					this->construitIGlobal(ix, iy, iz, iMailleG);
          m_limXm->creeLimite(&(*bord)[iFace]);
          iMailleD = iMailleG;
					(*bord)[iFace]->setFace(&m_faces[iFace]);
					(*bord)[iFace]->initialise((*cellules)[iMailleG], (*cellules)[iMailleD]);
					(*cellules)[iMailleG]->ajouteBord((*bord)[iFace]);
					m_faces[iFace].initialiseAutres(surface, normale, tangente, binormale);
        }
        posX = m_origine.getX();
        posY = (static_cast<double>(iy) + 0.5)*m_dY + m_origine.getY();
        posZ = (static_cast<double>(iz) + 0.5)*m_dZ + m_origine.getZ();
        m_faces[iFace].setPos(posX, posY, posZ);
        iFace++;
      }
    }
    //Limite X=lX
    ix = m_nombreMaillesX - 1;
    tangente.setXYZ(0., 1., 0.); normale.setXYZ(1., 0., 0.); binormale.setXYZ(0., 0., 1.);
    for (iy = 0; iy < m_nombreMaillesY; iy++)
    {
      for (iz = 0; iz < m_nombreMaillesZ; iz++)
      {   
        this->construitIGlobal(ix, iy, iz, iMailleG);
        //Traitement des limites communicantes
        if (m_decoupe == 'X' && rang < Ncpu - 1) {
          iMailleD = compteMaillesParallele++;
          if (ordreCalcul == "ORDRE1") { (*bord)[iFace] = new BordDeMaille; }
          else { (*bord)[iFace] = new BordDeMailleO2; }
          (*cellules)[iMailleD]->getElement()->setPos((*cellules)[iMailleG]->getElement()->getPosition());
          (*cellules)[iMailleD]->getElement()->setPosX((*cellules)[iMailleG]->getElement()->getPosition().getX() + m_dX);
          (*cellules)[iMailleD]->ajouteBord((*bord)[iFace]);
        } //On prend parmi les mailles fantomes dans l'ordre
        //Traitement des limites physiques
        else {
          m_limXp->creeLimite(&(*bord)[iFace]);
          iMailleD = iMailleG;
        }
        (*bord)[iFace]->setFace(&m_faces[iFace]);
        (*bord)[iFace]->initialise((*cellules)[iMailleG], (*cellules)[iMailleD]);
        (*cellules)[iMailleG]->ajouteBord((*bord)[iFace]);
        m_faces[iFace].initialiseAutres(surface, normale, tangente, binormale);
        posX = (static_cast<double>(ix) + 1.)*m_dX + m_origine.getX();
        posY = (static_cast<double>(iy) + 0.5)*m_dY + m_origine.getY();
        posZ = (static_cast<double>(iz) + 0.5)*m_dZ + m_origine.getZ();
        m_faces[iFace].setPos(posX, posY, posZ);
        iFace++;
      }
    }
  }
  //Initialisation des faces limites en Y
  //-------------------------------------
  if (m_nombreMaillesY != 1)
  {
    surface = m_dX*m_dZ;
    //Limite Y=0
    iy = 0;
    tangente.setXYZ(1., 0., 0.); normale.setXYZ(0., -1., 0.); binormale.setXYZ(0., 0., 1.);
    for (ix = 0; ix < m_nombreMaillesX; ix++)
    {
      for (iz = 0; iz < m_nombreMaillesZ; iz++)
      {
        //Traitement des limites communicantes
        if (m_decoupe == 'Y' && rang > 0) {
					tangenteGhost.setXYZ(-1., 0., 0.); normaleGhost.setXYZ(0., 1., 0.); binormaleGhost.setXYZ(0., 0., 1.);
					this->construitIGlobal(ix, iy, iz, iMailleD);
          iMailleG = compteMaillesParallele++;
          if (ordreCalcul == "ORDRE1") { (*bord)[iFace] = new BordDeMaille; }
          else { (*bord)[iFace] = new BordDeMailleO2; }
          (*cellules)[iMailleG]->getElement()->setPos((*cellules)[iMailleD]->getElement()->getPosition());
          (*cellules)[iMailleG]->getElement()->setPosY((*cellules)[iMailleD]->getElement()->getPosition().getX() - m_dY);
          (*cellules)[iMailleG]->ajouteBord((*bord)[iFace]);
					(*bord)[iFace]->setFace(&m_faces[iFace]);
					(*bord)[iFace]->initialise((*cellules)[iMailleG], (*cellules)[iMailleD]);
					(*cellules)[iMailleD]->ajouteBord((*bord)[iFace]);
					m_faces[iFace].initialiseAutres(surface, normaleGhost, tangenteGhost, binormaleGhost);
        } //On prend parmi les mailles fantomes dans l'ordre
        //Traitement des limites physiques
        else {
					this->construitIGlobal(ix, iy, iz, iMailleG);
          m_limYm->creeLimite(&(*bord)[iFace]);
          iMailleD = iMailleG;
					(*bord)[iFace]->setFace(&m_faces[iFace]);
					(*bord)[iFace]->initialise((*cellules)[iMailleG], (*cellules)[iMailleD]);
					(*cellules)[iMailleG]->ajouteBord((*bord)[iFace]);
					m_faces[iFace].initialiseAutres(surface, normale, tangente, binormale);
        }
        posX = (static_cast<double>(ix) + 0.5)*m_dX + m_origine.getX();
        posY = m_origine.getY();
        posZ = (static_cast<double>(iz) + 0.5)*m_dZ + m_origine.getZ();
        m_faces[iFace].setPos(posX, posY, posZ);
        iFace++;
      }
    }
    //Limite Y=lY
    iy = m_nombreMaillesY - 1;
    tangente.setXYZ(-1., 0., 0.); normale.setXYZ(0., 1., 0.); binormale.setXYZ(0., 0., 1.);
    for (ix = 0; ix < m_nombreMaillesX; ix++)
    {
      for (iz = 0; iz < m_nombreMaillesZ; iz++)
      {
        this->construitIGlobal(ix, iy, iz, iMailleG);
        //Traitement des limites communicantes
        if (m_decoupe == 'Y' && rang < Ncpu - 1) {
          iMailleD = compteMaillesParallele++;
          if (ordreCalcul == "ORDRE1") { (*bord)[iFace] = new BordDeMaille; }
          else { (*bord)[iFace] = new BordDeMailleO2; }
          (*cellules)[iMailleD]->getElement()->setPos((*cellules)[iMailleG]->getElement()->getPosition());
          (*cellules)[iMailleD]->getElement()->setPosY((*cellules)[iMailleG]->getElement()->getPosition().getX() + m_dY);
          (*cellules)[iMailleD]->ajouteBord((*bord)[iFace]);
        } //On prend parmi les mailles fantomes dans l'ordre
        //Traitement des limites physiques
        else {
          m_limYp->creeLimite(&(*bord)[iFace]);
          iMailleD = iMailleG;
        }
        (*bord)[iFace]->setFace(&m_faces[iFace]);
        (*bord)[iFace]->initialise((*cellules)[iMailleG], (*cellules)[iMailleD]);
        (*cellules)[iMailleG]->ajouteBord((*bord)[iFace]);
        m_faces[iFace].initialiseAutres(surface, normale, tangente, binormale);
        posX = (static_cast<double>(ix) + 0.5)*m_dX + m_origine.getX();
        posY = (static_cast<double>(iy) + 1.)*m_dY + m_origine.getY();
        posZ = (static_cast<double>(iz) + 0.5)*m_dZ + m_origine.getZ();
        m_faces[iFace].setPos(posX, posY, posZ);
        iFace++;
      }
    }
  }
  //Initialisation des faces limites en Z //FP// Decoupage // selon Z a faire
  //-------------------------------------
  if (m_nombreMaillesZ != 1)
  {
    surface = m_dX*m_dY;
    //Limite Z=0
    iz = 0;
    tangente.setXYZ(-1., 0., 0.); normale.setXYZ(0., 0., -1.); binormale.setXYZ(0., 1., 0.);
    for (ix = 0; ix < m_nombreMaillesX; ix++)
    {
      for (iy = 0; iy < m_nombreMaillesY; iy++)
      {
        m_limZm->creeLimite(&(*bord)[iFace]);
        (*bord)[iFace]->setFace(&m_faces[iFace]);
        this->construitIGlobal(ix, iy, iz, iMailleG);
        iMailleD = iMailleG;
        (*bord)[iFace]->initialise((*cellules)[iMailleG], (*cellules)[iMailleD]);
        (*cellules)[iMailleG]->ajouteBord((*bord)[iFace]);
        (*cellules)[iMailleD]->ajouteBord((*bord)[iFace]);
        m_faces[iFace].initialiseAutres(surface, normale, tangente, binormale);
        posX = (static_cast<double>(ix) + 0.5)*m_dX + m_origine.getX();
        posY = (static_cast<double>(iy) + 0.5)*m_dY + m_origine.getY();
        posZ = m_origine.getZ();
        m_faces[iFace].setPos(posX, posY, posZ);
        iFace++;
      }
    }
    //Limite Z=lZ
    iz = m_nombreMaillesZ - 1;
    tangente.setXYZ(1., 0., 0.); normale.setXYZ(0., 0., 1.); binormale.setXYZ(0., 1., 0.);
    for (ix = 0; ix < m_nombreMaillesX; ix++)
    {
      for (iy = 0; iy < m_nombreMaillesY; iy++)
      {
        m_limZp->creeLimite(&(*bord)[iFace]);
        (*bord)[iFace]->setFace(&m_faces[iFace]);
        this->construitIGlobal(ix, iy, iz, iMailleG);
        iMailleD = iMailleG;
        (*bord)[iFace]->initialise((*cellules)[iMailleG], (*cellules)[iMailleD]);
        (*cellules)[iMailleG]->ajouteBord((*bord)[iFace]);
        (*cellules)[iMailleD]->ajouteBord((*bord)[iFace]);
        m_faces[iFace].initialiseAutres(surface, normale, tangente, binormale);
        posX = (static_cast<double>(ix) + 0.5)*m_dX + m_origine.getX();
        posY = (static_cast<double>(iy) + 0.5)*m_dY + m_origine.getY();
        posZ = (static_cast<double>(iz) + 1.)*m_dZ + m_origine.getZ();
        m_faces[iFace].setPos(posX, posY, posZ);
        iFace++;
      }
    }
  }
}

//***********************************************************************

double MaillageCartesien::getLCFL() const
{
  return m_lCFL;
}

//***********************************************************************

double MaillageCartesien::getVolume() const
{
  return m_volume;
}

//***********************************************************************

void MaillageCartesien::decoupageParallele()
{
  int ix, iy, iz;
  int maille_par_cpu, reste, iMaille, voisin;
  int nombreElements, compteMaillesParallele(0), countElements(0);
  int *elements_rec, *elements_env;

  //Selon X
  if (m_nombreMaillesXGlobal >= m_nombreMaillesYGlobal && m_nombreMaillesXGlobal >= m_nombreMaillesZGlobal)
  {
    m_decoupe = 'X'; //if (rang == 0){ cout << "decoupage parallele selon axe des X" << endl; }
    maille_par_cpu = m_nombreMaillesXGlobal / Ncpu;
    reste = m_nombreMaillesXGlobal % Ncpu;

    if (rang < reste){ ++maille_par_cpu; }
    m_nombreMaillesX = maille_par_cpu;
    m_nombreMaillesY = m_nombreMaillesYGlobal;
    m_nombreMaillesZ = m_nombreMaillesZGlobal;

    //Determination du premier et dernier noeud du CPU dans le domaine complet
    m_numNoeudDeb = 0; m_numNoeudFin = 0;
    for (int p = 0; p <= rang; p++) {
      m_numNoeudDeb = m_numNoeudFin;
      int nb_maille_cpu_p = m_nombreMaillesXGlobal / Ncpu;
      int reste = m_nombreMaillesXGlobal % Ncpu;
      if (p < reste) { ++nb_maille_cpu_p; }
      m_numNoeudFin += nb_maille_cpu_p;
    }

    double offset_x = static_cast<double>(rang)*static_cast<double>(maille_par_cpu);
    if (rang >= reste){ offset_x += static_cast<double>(reste); }
    m_origine.setX(offset_x*m_dX);
    //cout << m_origine.getX() << endl;

    //Nombre de mailles sur ce Cpu;
    m_nombreCellulesCalcul = m_nombreMaillesX*m_nombreMaillesY*m_nombreMaillesZ;

    //Determination du nombre d'element a envoyer/recevoir
    nombreElements = m_nombreMaillesY*m_nombreMaillesZ;

    //Nombre de mailles total en comptant les mailles necessaires au parallele;
    if (rang > 0 && rang < Ncpu - 1){ m_nombreCellulesTotal = m_nombreCellulesCalcul + 2 * nombreElements; }
    else{ m_nombreCellulesTotal = m_nombreCellulesCalcul + nombreElements; }

    ////***************Table de connectivite**********

    //Le numero des mailles paralleles est situe en dehors du nombre de maille "vraies"
    //donc commence a nombreMailles
    compteMaillesParallele = m_nombreCellulesCalcul;
    elements_rec = new int[nombreElements];
    elements_env = new int[nombreElements];
    countElements = 0;
    if (rang > 0)
    {
      voisin = rang - 1;
      Calcul_Parallele.setVoisin(voisin);
      ix = 0;
      for (iy = 0; iy < m_nombreMaillesY; iy++)
      {
        for (iz = 0; iz < m_nombreMaillesZ; iz++)
        {
          elements_rec[countElements] = compteMaillesParallele;
          this->construitIGlobal(ix, iy, iz, iMaille);
          elements_env[countElements] = iMaille;
          ++countElements;
          ++compteMaillesParallele;
        }
      }
      Calcul_Parallele.setElementsAEnvoyer(voisin, elements_env, nombreElements);
      Calcul_Parallele.setElementsARecevoir(voisin, elements_rec, nombreElements);
    }
    countElements = 0;
    if (rang < Ncpu - 1)
    {
      voisin = rang + 1;
      Calcul_Parallele.setVoisin(voisin);
      ix = m_nombreMaillesX - 1;
      for (iy = 0; iy < m_nombreMaillesY; iy++)
      {
        for (iz = 0; iz < m_nombreMaillesZ; iz++)
        {
          elements_rec[countElements] = compteMaillesParallele;
          this->construitIGlobal(ix, iy, iz, iMaille);
          elements_env[countElements] = iMaille;
          ++countElements;
          ++compteMaillesParallele;
        }
      }
      Calcul_Parallele.setElementsAEnvoyer(voisin, elements_env, nombreElements);
      Calcul_Parallele.setElementsARecevoir(voisin, elements_rec, nombreElements);
    }
  }
  else if (m_nombreMaillesYGlobal >= m_nombreMaillesXGlobal && m_nombreMaillesYGlobal >= m_nombreMaillesZGlobal)
  {
    m_decoupe = 'Y'; if (rang == 0){ cout << "decoupage parallele selon axe des Y" << endl; }
    maille_par_cpu = m_nombreMaillesYGlobal / Ncpu;
    reste = m_nombreMaillesYGlobal % Ncpu;

    if (rang < reste){ ++maille_par_cpu; }
    m_nombreMaillesY = maille_par_cpu;
    m_nombreMaillesX = m_nombreMaillesXGlobal;
    m_nombreMaillesZ = m_nombreMaillesZGlobal;

    double offset_y = static_cast<double>(rang)*static_cast<double>(maille_par_cpu);
    if (rang >= reste){ offset_y += static_cast<double>(reste); }
    m_origine.setY(offset_y*m_dY);

    //Definition des voisins et determination du nombre d'element a envoyer/recevoir
    nombreElements = m_nombreMaillesX*m_nombreMaillesZ;

    //Nombre de mailles totale en comptant les mailles necessaires au parallele;
    m_nombreCellulesCalcul = m_nombreMaillesX*m_nombreMaillesY*m_nombreMaillesZ;
    if (rang > 0 && rang < Ncpu - 1){ m_nombreCellulesTotal = m_nombreCellulesCalcul + 2 * nombreElements; }
    else{ m_nombreCellulesTotal = m_nombreCellulesCalcul + nombreElements; }

    ////***************Table de connectivite**********

    //Le numero des mailles paralleles est situe en dehors du nombre de maille "vraies"
    compteMaillesParallele = m_nombreMaillesX*m_nombreMaillesY*m_nombreMaillesZ;
    elements_rec = new int[nombreElements];
    elements_env = new int[nombreElements];
    countElements = 0;
    if (rang > 0)
    {
      voisin = rang - 1;
      Calcul_Parallele.setVoisin(voisin);
      iy = 0;
      for (ix = 0; ix < m_nombreMaillesX; ix++)
      {
        for (iz = 0; iz < m_nombreMaillesZ; iz++)
        {
          elements_rec[countElements] = compteMaillesParallele;
          this->construitIGlobal(ix, iy, iz, iMaille);
          elements_env[countElements] = iMaille;
          ++countElements;
          ++compteMaillesParallele;
        }
      }
      Calcul_Parallele.setElementsAEnvoyer(voisin, elements_env, nombreElements);
      Calcul_Parallele.setElementsARecevoir(voisin, elements_rec, nombreElements);
    }
    countElements = 0;
    if (rang < Ncpu - 1)
    {
      voisin = rang + 1;
      Calcul_Parallele.setVoisin(voisin);
      iy = m_nombreMaillesY - 1;
      for (ix = 0; ix < m_nombreMaillesX; ix++)
      {
        for (iz = 0; iz < m_nombreMaillesZ; iz++)
        {
          elements_rec[countElements] = compteMaillesParallele;
          this->construitIGlobal(ix, iy, iz, iMaille);
          elements_env[countElements] = iMaille;
          ++countElements;
          ++compteMaillesParallele;
        }
      }
      Calcul_Parallele.setElementsAEnvoyer(voisin, elements_env, nombreElements);
      Calcul_Parallele.setElementsARecevoir(voisin, elements_rec, nombreElements);
    }
  }
}

//***********************************************************************

string MaillageCartesien::quiSuisJe() const
{
  return "CARTESIEN";
}

//**************************************************************************
//******************************** ECRITURE ********************************
//**************************************************************************

//****************************************************************************

string MaillageCartesien::recupereChaineExtent(int rangLocal, bool global) const
{
  stringstream chaineExtent;

  //Determination du premier et dernier noeud du CPU dans le domaine complet
  int numNoeudDeb = 0; int numNoeudFin = 0;
  for (int p = 0; p <= rangLocal; p++) {
    numNoeudDeb = numNoeudFin;
    int nb_maille_cpu_p = m_nombreMaillesXGlobal / Ncpu;
    int reste = m_nombreMaillesXGlobal % Ncpu;
    if (p < reste) { ++nb_maille_cpu_p; }
    numNoeudFin += nb_maille_cpu_p;
  }
  //Gestion numNoeudDeb et Fin selon decoupage
  if (Ncpu == 1 || m_decoupe == 'X') {
    if(global) { chaineExtent << " 0 " << m_nombreMaillesXGlobal << " 0 " << m_nombreMaillesYGlobal << " 0 " << m_nombreMaillesZGlobal; }
    else{ chaineExtent << numNoeudDeb << " " << numNoeudFin << " 0 " << m_nombreMaillesY << " 0 " << m_nombreMaillesZ; }
  }
  else { Erreurs::messageErreur("decoupage en Y ou Z non gere dans ecritSolutionXML"); return 0; }
  return chaineExtent.str();
}

//****************************************************************************

void MaillageCartesien::recupereCoord(vector<Cellule *> *cellulesLvl, std::vector<double> &jeuDonnees, Axe axe) const
{
  int numCell;
  switch (axe) {
  case X:
    for (int i = 0; i < m_nombreMaillesX; i++)
    {
      construitIGlobal(i, 0, 0, numCell);
      jeuDonnees.push_back(cellulesLvl[0][numCell]->getElement()->getPosition().getX() - 0.5*m_dX);
    }
    jeuDonnees.push_back(cellulesLvl[0][numCell]->getElement()->getPosition().getX() + 0.5*m_dX);
    break;
  case Y:
    for (int j = 0; j < m_nombreMaillesY; j++)
    {
      construitIGlobal(0, j, 0, numCell);
      jeuDonnees.push_back(cellulesLvl[0][numCell]->getElement()->getPosition().getY() - 0.5*m_dY);
    }
    jeuDonnees.push_back(cellulesLvl[0][numCell]->getElement()->getPosition().getY() + 0.5*m_dY);
    break;
  case Z:
    for (int k = 0; k < m_nombreMaillesZ; k++)
    {
      construitIGlobal(0, 0, k, numCell);
      jeuDonnees.push_back(cellulesLvl[0][numCell]->getElement()->getPosition().getZ() - 0.5*m_dZ);
    }
    jeuDonnees.push_back(cellulesLvl[0][numCell]->getElement()->getPosition().getZ() + 0.5*m_dZ);
    break;
  }
}

//****************************************************************************
/*!
* \brief     Recuperation de donnees pour ecriture
* \details   Fonctions permettant de recuperer un esemble de donnees de phase ou de melange, 
*                  scalaire ou vectorielle
* \param    cellules         cellules contenant les donnees
* \param    var              numero de la variable demandee (>0 si scalaire, <0 si vectorielle)
* \param    phase            numero de la phase demandee (-1 si melange demande)
* \return   jeuDonnees       vecteur de double contenant toutes les donnees de la variables var de la phase phase a la suite
*/
void MaillageCartesien::recupereDonnees(vector<Cellule *> *cellulesLvl, std::vector<double> &jeuDonnees, const int var, int phase, int lvl) const
{
  int numCell;
  for (int k = 0; k < m_nombreMaillesZ; k++) {
    for (int j = 0; j < m_nombreMaillesY; j++) {
      for (int i = 0; i < m_nombreMaillesX; i++) {
        construitIGlobal(i, j, k, numCell);
        if (var > 0) { //On veut recuperer les donnees scalaires
          if (phase >= 0) { jeuDonnees.push_back(cellulesLvl[0][numCell]->getPhase(phase)->renvoieScalaire(var)); } //Donnees de phases
          else if (phase == -1) { jeuDonnees.push_back(cellulesLvl[0][numCell]->getMelange()->renvoieScalaire(var)); }               //Donnees de melange
          else if (phase == -2) { jeuDonnees.push_back(cellulesLvl[0][numCell]->getTransport(var - 1).getValeur()); }
          else if (phase == -3) { jeuDonnees.push_back(cellulesLvl[0][numCell]->getXi()); }
          else if (phase == -4) { jeuDonnees.push_back(cellulesLvl[0][numCell]->getGradient()); }
          else { Erreurs::messageErreur("MaillageCartesien::recupereDonnees : numero de phase inconnu : ", phase); }
        }
        else { //On veut recuperer les donnees vectorielles
          if (phase >= 0) { //Donnees de phases
            jeuDonnees.push_back(cellulesLvl[0][numCell]->getPhase(phase)->renvoieVecteur(-var)->getX());
            jeuDonnees.push_back(cellulesLvl[0][numCell]->getPhase(phase)->renvoieVecteur(-var)->getY());
            jeuDonnees.push_back(cellulesLvl[0][numCell]->getPhase(phase)->renvoieVecteur(-var)->getZ());
          }
          else if (phase == -1) {  //Donnees de melange
            jeuDonnees.push_back(cellulesLvl[0][numCell]->getMelange()->renvoieVecteur(-var)->getX());
            jeuDonnees.push_back(cellulesLvl[0][numCell]->getMelange()->renvoieVecteur(-var)->getY());
            jeuDonnees.push_back(cellulesLvl[0][numCell]->getMelange()->renvoieVecteur(-var)->getZ());
          }
          else { Erreurs::messageErreur("MaillageCartesien::recupereDonnees : numero de phase inconnu : ", phase); }
        } //Fin vecteur
      } // Fin X
    } //Fin Y
  } //Fin Z
}

//****************************************************************************
//*****************************Accesseurs*************************************
//****************************************************************************

double MaillageCartesien::getdX() const
{
  return m_dX;
}

//***********************************************************************

double MaillageCartesien::getdY() const
{
  return m_dY;
}

//***********************************************************************

double MaillageCartesien::getdZ() const
{
  return m_dZ;
}

//***********************************************************************