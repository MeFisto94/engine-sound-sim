# Rust Engine Sound Simulation Library
This project is an approach to procedural audio / physical simulation of the sound of internal combustion engines.
It is primarily targeted at real-time applications and hence does not intend to use rigid body physics, as compared
to other implementations.

It is currently in a very early stage, not much more than a scratch-book of formulae, but I got lost in the 
thermodynamics and as such I thought sharing it in this stage may help get more eyes on the maths. 
_We can always refactor later™_.

Currently, the code mostly consists of the test inside lib.rs that tries to create a csv file for plotting and a wave 
file. The wave file currently sounds broken, which is somewhat related due to normalizing the pressure samples into 
the [0, 1] domain. Still, I am not satisfied with the pressure over time, since it falls too steep (if you plot p-V, the
pressure is mostly constant for some time in real engines). This may be related to the shape of the function calculating
the pressure coming from ignition (it's completely simplified as in separating compression and ignition pressure and
adding them later), which is a linear curve instead of something that is used in real world scenarios
Vibe Curve [\[1\]](https://www.researchgate.net/figure/Vibe-curves-and-significant-burning-points-for-different-fuel-mixtures-at-constant-spark_fig1_314193040)
[\[2\]](https://doc.simulationx.com/4.0/1033/Content/Libraries/PowerTransmission/CombustionEngines/ExcitationModels/Vibe.htm)

### Open Questions
1. With the Exhaust Pressure Wave being "filtered" by the exhaust valve, is it that only the exhaust stroke is relevant
for the sound? (Unless we'd need the maximum combustion pressure for the exhaust stroke pressure modelling)
   1. The exhaust stroke pressure curve seems easier than the actual combustion/work stroke ones and thus maybe
   my worries about the curve form may be completely irrelevant and it's just a simple expansion thermodynamics wise.
