# Python Reentry Simulations

This is the python codes I use to simulate spacecraft reentries. Requires Python 3 & Matplotlib.

## Usage

The script is very basic at the moment. Run data are defined in a `.json` file. At the moment, only `Earth` and `Mars` are valid planets, but more could be added fairly easily in the script.

````
{
    "craft": {
        "name": "MSL",
        "ballistic_coef": 146,
        "lift_drag": 0.24
    },
    "planet": "Mars",
    "sim": {
        "max_it": 1000000,
        "delta_t": 0.1,
        "entry_interface": 135e3,
        "fpa": -15.5,
        "velocity": 5.6e3,
        "stop_alt": 0
    }
}
````

Once you have a run file, just call the script on it and it'll generate your plots

````bash
$ ./reentry.py msl.json
$ ls

    MSL-dtg.png
    MSL-load_alt.png
    MSL-load_time.png
    MSL-traj.png
    MSL-vel.png
````

## Method

The program solves the newtonian equations of motion and computes the trajectory of the reentering spacecraft using the forces applied to it (gravity, aerodynamic drag & lift) using a very basic Euler integrator. You can read more about that on [the blog post](https://amyparent.com/post/crew-dragon-reentry/).

The atmospheres are modelled using the scale height method. It's simplistic and could benefit from better models (NASA Standard Atmosphere, maybe even Mars GRAM and co.) but the gains are likely to be dwarfed by errors induced by the simple assumptions made everywhere (constant lift-drag ratio, constant drag coefficient…)
