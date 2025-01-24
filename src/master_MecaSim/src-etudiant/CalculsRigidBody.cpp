/*
 * CalculsRigidBody.cpp :
 * Copyright (C) 2016 Florence Zara, LIRIS
 *               florence.zara@liris.univ-lyon1.fr
 *               http://liris.cnrs.fr/florence.zara/
 *
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 2.1 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/** \file CalculsMSS.cpp
 Programme calculant pour chaque particule i d un MSS son etat au pas de temps suivant
 (methode d 'Euler semi-implicite) : principales fonctions de calculs.
 \brief Fonctions de calculs de la methode semi-implicite sur un systeme masses-ressorts.
 */

#include <stdio.h>
#include <math.h>
#include <vector>
#include <iostream>
#include <numeric>

#include "vec.h"
#include "ObjetSimule.h"
#include "ObjetSimuleRigidBody.h"
#include "Viewer.h"

using namespace std;




/*
 * Calcul de la masse de l objet rigide
 */
void ObjetSimuleRigidBody::CalculMasse()
{
    _Mass = std::accumulate(M.begin(), M.end(), 0.0);

    _BaryCentre = Vector();

    for (int i = 0; i < _Nb_Sommets; ++i)
    {
        _BaryCentre = _BaryCentre + M[i] * P[i];
    }

    _BaryCentre = _BaryCentre / _Mass;
}


/*
 * Calcul du tenseur d inertie de l objet rigide - - partie constante Ibody
 * et remplissage tableau roi (position particules dans repere objet)
 */
void ObjetSimuleRigidBody::CalculIBody()
{
    _ROi.resize(_Nb_Sommets);

    for (unsigned int i = 0; i < _Nb_Sommets; i++)
    {
        // Calculer la position de la particule i dans l'objet
        _ROi[i] = P[i] - _Position;
    }

    Matrix t;

    for (unsigned int i = 0; i < _Nb_Sommets; i++)
    {
        float xi = _ROi[i].x;
        float yi = _ROi[i].y;
        float zi = _ROi[i].z;

        // Tenseur d'inertie
        t.m_Values[0] = xi * xi;
        t.m_Values[1] = xi * yi;
        t.m_Values[2] = xi * zi;
        t.m_Values[3] = yi * xi;
        t.m_Values[4] = yi * yi;
        t.m_Values[5] = yi * zi;
        t.m_Values[6] = zi * xi;
        t.m_Values[7] = zi * yi;
        t.m_Values[8] = zi * zi;

        _Ibody += ((xi * xi + yi * yi + zi * zi) * Matrix::UnitMatrix() - t);
        _IbodyInv = _Ibody.InverseConst();
    }

}


/*
 * Calcul de l etat de l objet rigide.
 */
void ObjetSimuleRigidBody::CalculStateX()
{
    _Position = _Position - _BaryCentre;

    // Calcul de l'inverse du tenseur d inertie
    _InertieTenseurInv = _Rotation * _IbodyInv * _Rotation.TransposeConst();
    
    // Calcul de la vitesse angulaire
    _VitesseAngulaire = _InertieTenseurInv * _MomentCinetique;
}



/*
 * Calcul de la derivee de l etat : d/dt X(t).
 */
void ObjetSimuleRigidBody::CalculDeriveeStateX(Vector gravite)
{
    // v(t) = P(t) / mass
    _Vitesse = _QuantiteMouvement / _Mass;

    // R'(t)
    _RotationDerivee = StarMatrix(_VitesseAngulaire) * _Rotation;
    
    // P'(t) = F(t)
    // gravite = normalize(Vector(1.f, 1.f, 2.f) - _Position) * 9.81f;
    _Force = _Mass * gravite * 0.001f;
    
    // L'(t) = moment total de la force (torque)
    _Torque = Vector();
    for (size_t i = 0; i < _Ri.size(); ++i)
    {
        _Torque = _Torque + cross(_ROi[i] - _Position, _Force); 
    }
}


/**
 * Schema integration pour obbtenir X(t+dt) en fct de X(t) et d/dt X(t)
 */
void ObjetSimuleRigidBody::Solve(float visco)
{
    // x(t+dt) - position
    _Position = _Position + _Vitesse * _delta_t;
    
    // R(t+dt) - rotation
    _Rotation = _Rotation + _RotationDerivee * _delta_t;
    
    // L(t+dt) - moment cinetique
    _MomentCinetique = _MomentCinetique + _Torque * _delta_t;
    
    // P(t+dt) - quantite mouvement    
    _QuantiteMouvement = _QuantiteMouvement + _Force * _delta_t;

    _VitesseAngulaire = _InertieTenseurInv * _MomentCinetique;

    // ri(t)
    for (size_t i = 0; i < P.size(); ++i)
    {
        _Ri[i] = _Rotation * _ROi[i] + _Position;
    }

}//void



/**
 * Gestion des collisions avec le sol.
 */
void ObjetSimuleRigidBody::Collision()
{
    for (size_t i = 0; i < _Ri.size(); ++i)
    {
        if (_Ri[i].y < -3.f)
        {
            _Position.y += 6.f;
            return;
        }
    }
}// void

