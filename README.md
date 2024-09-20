# Small asteroid fragmentation model
In recent years, nine small near-Earth asteroids (NEAs), with a diameter of the order of one meter1, have been discovered a few hours before hitting our planet, from 2008 TC3 on 7 October 2008, (Farnocchia et al., 2017; Shaddad37 et al., 2010), to 2024 BX1 on 21 January 2024 (Spurný et al., 2024), ending with 2024 RW1 dropped in early September 2024. Asteroids of such small dimensions can be discovered thanks to the increased detection efficiency of NEA surveys over the last 20 years. The discovery of small asteroids before they hit the Earth offers a good opportunity to understand the physical origin of these bodies compared to the more imprecise satellite observations (Peña-Asensio et al., 2022). Small NEAs of this size enter the atmosphere at hypervelocity and generally experience ablation, fragmentation and airburst, generating brilliant fireballs. The major fragments surviving the airburst can enter the dark flight phase and generate a fall event. In this case, several meteorites may be found on the so-called “strewn field” on the ground, allowing the chemical-physical characterisation of the parent body. Finally, knowledge of the orbit allows us to reconstruct, within certain limits, the dynamic history of the NEA.\

When a small asteroid falls into the terrestrial atmosphere, the path can be qualitatively divided into three main phases. The first is when the asteroid enters the upper atmosphere with a speed of the order of 10-20 km/s; the drag starts, and a shock wave is formed in the front of the body, the air is compressed and heated, and mass ablation begins: this is the fireball phase. The second phase starts when the pressure of the shock wave is higher than the body’s strength, and fragmentation occurs. After fragmentation, the body’s mass spreads over a greater area, the quantity of intercepted atmosphere increases, and thus, braking and ablation increase: the body loses kinetic energy very rapidly, i.e. explosively, and there is a little airburst. Afterwards, the fragments decouple, developing their shock wave as independent fireballs.\

Usually, the meteoroids fragmentation model assumes that the process starts when the aerodynamic pressure in front of the body is equal or superior to mechanical strength 𝑆 of the body. A priori, we do not know the masses of the fragments in the strewn field, i.e. the meteorites’ masses. So, in our model, we assume the following final masses with spherical shape: 1, 0.3, 0.2, 0.1, 0.05, 0.02, 0.005 and 0.001 kg. In practice, we use some “sample” particles to show the final strewn field’s possible range and position. With good approximation, these masses for the fragments correspond to the typical masses of meteorites on the ground. In all cases, the sum of the fragment’s mass is inferior to the meteoroid mass after the fragmentation phase, so we assume that the great part of the asteroid’s mass was lost in ablation and debris cloud. Note that after the airburst, there is a residual ablation phase, which tends to decrease the mass of the fragments, which we consider in the computation. So, the original fragment mass after the fragmentation is larger than our assumed meteorite’s mass that “samples” the strewn field. There is no guarantee that any meteorites corresponding to these assumed masses will be produced, even if they are the mass of typical meteorites.\

For these fragments, we will assume that the speed is equal to that of the fragmentation phase, with the same inclination and direction as the original body. With this last assumption, we simplified the model because we neglect the possible lateral expansion speed imprinted during fragmentation. This speed component is usually very small, so pieces from a fragmented fireball, in most cases, continue along their original trajectory.\

A good atmospheric profile is one of the critical points, along with the value of strength, for computing a realistic strewn field. The essential quantities it must report are the air density at different heights (or the pressure and temperature from which the density can be calculated) and the speed and direction of the wind. Solving the classical drag and ablation equations numerically with a Runge-Kutta 4th/5th order solver and using the World Geodetic System 84 (WGS 84) provides the position and velocity of the meteoroid as it falls towards the ground, and the intersection of the trajectory with the ground gives the impact point.\

This is a Matlab software for describing the fall and the possible strewn field of a small asteroid using a simple fragmentation model. There are input data for 2023 CX1 and 2024 BX1 with initial conditions in the Settings.txt file. There are also the atmospheric profile from the starting height of the dark flight to the ground. As outputs there are .txt files with numerical data about height, speed, latitude, longitude, mass and time of the main body and fragments paths. There is a file with the strewn field coordinates and some plots of the results.
