/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "hronTurekInlet.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::hronTurekInlet::
hronTurekInlet
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    maxValue_(0),
    n_(1, 0, 0),
    y_(0, 1, 0),
    rampTime_(0)
{}


Foam::hronTurekInlet::
hronTurekInlet
(
    const hronTurekInlet& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    maxValue_(ptf.maxValue_),
    n_(ptf.n_),
    y_(ptf.y_),
    rampTime_(ptf.rampTime_)
{}


Foam::hronTurekInlet::
hronTurekInlet
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict),
    maxValue_(readScalar(dict.lookup("maxValue"))),
    n_(dict.lookup("n")),
    y_(dict.lookup("y")),
    rampTime_(readScalar(dict.lookup("rampTime")))
{}


Foam::hronTurekInlet::
hronTurekInlet
(
    const hronTurekInlet& pivpvf
)
:
    fixedValueFvPatchVectorField(pivpvf),
    maxValue_(0),
    n_(1, 0, 0),
    y_(0, 1, 0),
    rampTime_(0)
{}


Foam::hronTurekInlet::
hronTurekInlet
(
    const hronTurekInlet& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF),
    maxValue_(pivpvf.maxValue_),
    n_(pivpvf.n_),
    y_(pivpvf.y_),
    rampTime_(pivpvf.rampTime_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::hronTurekInlet::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    scalar curMaxValue = maxValue_;

    if(this->db().time().value() < rampTime_)
    {
        scalar t = this->db().time().value();
        curMaxValue *= 0.5*(1.0-cos(M_PI*t/rampTime_));
    }

    // get range and orientation
    boundBox bb(patch().patch().localPoints(), true);

    vector ctr = 0.5*(bb.max() + bb.min());
 
    const vectorField& c = patch().Cf();

    // Calculate local 1-D coordinate for the parabolic profile
    scalarField coord = 2*((c-ctr) & y_)/((bb.max() - bb.min()) & y_);

    vectorField::operator=(n_*curMaxValue*(1.0 - sqr(coord)));

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::hronTurekInlet::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("maxValue") <<  maxValue_ << token::END_STATEMENT << nl;
    os.writeKeyword("n") <<  n_ << token::END_STATEMENT << nl;
    os.writeKeyword("y") <<  y_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::hronTurekInlet::operator=
(
    const fvPatchField<vector>& pvf
)
{
    fvPatchField<vector>::operator=(patch().nf()*(patch().nf() & pvf));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        hronTurekInlet
    );
}

// ************************************************************************* //
