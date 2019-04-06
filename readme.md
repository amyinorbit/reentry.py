# Python Reentry Simulations

This is the python codes I use to simulate spacecraft reentries. Requires Python 3 & Matplotlib.

## Method

The program solves the newtonian equations of motion and computes the trajectory of the reentering spacecraft using the forces applied to it (gravity, aerodynamic drag & lift) using a very basic Euler integrator. You can read more about that on [the blog post](https://amyparent.com/post/crew-dragon-reentry/).

The atmospheres are modelled using the scale height method. It's simplistic and could benefit from better models (NASA Standard Atmosphere, maybe even Mars GRAM and co.) but the gains are likely to be dwarfed by errors induced by the simple assumptions made everywhere (constant lift-drag ratio, constant drag coefficient…)
