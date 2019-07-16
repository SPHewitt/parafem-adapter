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

Class
    Foam::parabolicInletVelocityFvPatchVectorField

Description
    This velocity inlet boundary condition is applied to patches where the
    Velocity is specified.  The inflow velocity is obtained from the flux with
    a direction normal to the patch faces.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            parabolicInletVelocity;
        phi             phi;
        rho             rho;
        value           uniform 0;
    }
    \endverbatim

Note
    If reverse flow is possible or expected use
    the parabolicInletOutletVelocityFvPatchVectorField condition instead.

See also
    Foam::fixedValueFvPatchField
    Foam::parabolicInletOutletVelocityFvPatchVectorField

SourceFiles
    parabolicInletVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef parabolicInletVelocityFvPatchVectorField_H
#define parabolicInletVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class parabolicInletVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class parabolicInletVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Peak velocity magnitude 
        scalar maxValue_;

        //- Flow direction 
        vector n_;

        //- Direction of the y-coordinate
        vector y_;

public:

    //- Runtime type information
    TypeName("parabolicInletVelocity");


    // Constructors

        //- Construct from patch and internal field
        parabolicInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        parabolicInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given parabolicInletVelocityFvPatchVectorField
        //  onto a new patch
        parabolicInletVelocityFvPatchVectorField
        (
            const parabolicInletVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        parabolicInletVelocityFvPatchVectorField
        (
            const parabolicInletVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new parabolicInletVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        parabolicInletVelocityFvPatchVectorField
        (
            const parabolicInletVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new parabolicInletVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Attributes

            //- Return true: this patch field is altered by assignment
            virtual bool assignable() const
            {
                return true;
            }


        // Access

            //- Return maxValue 
            const scalar& maxValue() const
            {
                return maxValue_;
            }

            //- Return flow direction 
            vector& n()
            {
                return n_;
            }

            //- Return y direction 
            vector& y()
            {
                return y_;
            }


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;


    // Member operators

        virtual void operator=(const fvPatchField<vector>& pvf);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
