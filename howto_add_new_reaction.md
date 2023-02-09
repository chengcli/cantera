# Steps to add freezing reaction

1. create a header file under kinetics, FreezingRate.h
1.1. derive class FreezingData from ReactionData
1.2. derive class FreezingRate from ReactionRate

2.2. add the header file to ReactionFactory.cpp
```
  #include "cantera/kinetics/FreezingRate.h"
```
2.2. register the reaction in ReactionRateFactory.cpp
```
  reg("freezing", []()) {
  }
```

