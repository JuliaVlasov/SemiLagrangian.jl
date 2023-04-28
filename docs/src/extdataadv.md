# ExtDataAdv

The `AbstractExtDataAdv` interface is an interface which allows to
define the treatments which will make it possible to obtain the
values to be applied to the advections.

##  Functions that needs to be implemented

# Methods to define

- `initcoef!(parext::AbstractExtDataAdv, self::AdvectionData)` : this method called at the beginning of each advection to initialize parext data. The `self.parext` mutable structure is the only data that initcoef! can modify otherwise it leads to unpredictable behaviour.
- `getalpha(parext::AbstractExtDataAdv, self::AdvectionData, ind)` : return the alpha number that is used for interpolation.
- `getperm(parext::AbstractExtDataAdv, advd::AdvectionData)` : get the permutation of the dimension as a function of the current state, the dimension where advection occurs must be first, the dimensions used to compute alpha must be at the end.
