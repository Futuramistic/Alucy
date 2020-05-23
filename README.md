
# Alucy
Mesh simplification and Loop subdivision project for 3D modelling and animatiom (CS3242) module at NUS.

## Mesh simplification
Mesh simplification using iterative decimation with edge contraction (Criterion: curvature and length of each edge)
<div style="display:inline-block">
  <img style="vertical-align:top" width="300" height="300" src="https://github.com/Futuramistic/Alucy/blob/master/Images/rabbit.PNG">
  <img style="vertical-align:top" width="300" height="300" src="https://github.com/Futuramistic/Alucy/blob/master/Images/rabbit2500.PNG">
  <img style="vertical-align:top" width="300" height="300" src="https://github.com/Futuramistic/Alucy/blob/master/Images/rabbit1000.PNG">
</div>
<i>Figures: Rabbit with 5000 faces / Rabbit with 2500 faces / Rabbit with 1000 faces respectively.</i>

### Issues encountered with implementation:
* Relatively slow
* Link condition not implemented <b>(!)</b>

## Loop subdivision
Implementation of Loop Subdivision algorithm
<div style="display:inline-block">
  <img style="vertical-align:top" width="200" height="300" src="https://github.com/Futuramistic/Alucy/blob/master/Images/skull15kV.PNG">
  <img style="vertical-align:top" width="200" height="300" src="https://github.com/Futuramistic/Alucy/blob/master/Images/skull55kV.PNG">
  <img style="vertical-align:top" width="200" height="300" src="https://github.com/Futuramistic/Alucy/blob/master/Images/skull220kV.PNG">
</div>

### Issues encountered with implementation:
* Problems with sharp/crease vertices - not recognized by algorithm
