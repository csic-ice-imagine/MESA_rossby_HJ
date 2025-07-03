This repository contains the necessary MESA files to obtain the results shown in Elias-López et al 2025 b "Rossby number regime, convection suppression, and dynamo-generated magnetism in inflated hot Jupiters".

MESA version 24.08.1 was used. Additional Python scripts for analysis can be given upon request.

The files in planet_create_and_core/ can create a planet and add a core with a relaxation stage. This can be used to generate many initial HJs of different radii and masses.

The files in planet_evolve_heat/ can evolve an already created planet with irradiation and possible internal heating using the Thorngen & Fortney 2018 Bayesian formula. The script run.py can run a series of models for different irradiation and thus different heating efficiencies.
