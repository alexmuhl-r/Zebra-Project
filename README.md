# Zebra Stripe Salience Project

This repository includes code and data from our zebra stripe salience project. You can read about this project in our preprint available here: https://www.biorxiv.org/content/10.1101/2021.04.16.440148v1.

## Human and Computational Model Salience Estimates

The first part of this project involves human and computational model estimates of salience in response to images of zebra. These images included four stripe patterns from Grevy's, mountain, plains zebra with shadow stripes and plains zebra without shadow stripes. These images were adjusted to simulate motion, distance and predator vision (colour and acuity). Human salience estimates were gathered online using the Gorilla Experiment Builder: https://app.gorilla.sc/openmaterials/207836. Model estimates were computed using Matlab (details in the preprint).

## Lion-Zebra Pursuit Simulations

The second part of this project involves simulations of predator-prey pursuit, included in a Matlab script.

The prey (zebra) runs along the x axis from point (0,0), the predator (lion) starts at given distance in y and at x=0. Starting distance and predator max speed are altered sequentially to model the difference in chase time when the predator aims at the front of the zebra versus the rear of the zebra (details in the preprint). Many parameters can be adjusted including: Ppredator max speed, prey max speed, predator acceleration, prey acceleration, predator speed increase increment, predator starting distance, predator distance increase increment, max chase time.
