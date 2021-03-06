#ComputeAxisymmetricRZSmallStrain
!syntax description /Materials/ComputeAxisymmetricRZSmallStrain


## Description
!include docs/content/documentation/modules/tensor_mechanics/common/supplementalAxisymmetricRZStrain.md

Only minimal changes are required to adapt a plane strain problem to an axisymmetric problem; thus,
`ComputeAxisymmetricRZISmallStrain` inherits from `Compute2DSmallStrain`.  The finite strain axisymmetric code overwrites the method used to calculate the strain component $\nabla u_{22}$, before calculating the total strain measure with the small strain assumption, as described in the [Strains](tensor_mechanics/Strains.md) page.

!syntax parameters /Materials/ComputeAxisymmetricRZSmallStrain

!syntax inputs /Materials/ComputeAxisymmetricRZSmallStrain

!syntax children /Materials/ComputeAxisymmetricRZSmallStrain
