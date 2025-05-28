# R-project_Tocilizumab
This is a cost-utility model to evaluate the cost effectiveness of Tocilizumab for people with Giant Cell Arthritis (GCA). 
The model is based on the HTA submitted by Roche in 2017. 

## Interventions
Intervention arm: Tocilizumab for 2 years + 26 weeks tapering steroids
Comparison arm: Placebo + 52 weeks tapering steroids
Cumulative GC dose based on GiACTA trial and real-world data: 1) prior to the first flare; 2) during secondary remission; 3) during relapse/flare

## Overview of Model Structure
A semi-Markov model using a weekly cycle length and life-time (30 years) time horizon:
- On remission and on steroid
- On remission and off steroid
- On relapse/flare
- On remission and on maintenance steroids (escape)
- GCA-related complications
- Steroid-related AEs
- Dead

## Data source
GiACTA trial data to estimate the impact of tocilizumab on disease control (time in remission and number of flares)
Real world data to estimate the relationship between cumulative prednisone dose and risk of steroid related adverse events (fractures and diabetes)
Previous literature for the risk of GCA-related complications (vision loss and stroke) among patients at relapse/flare state

## Survival analyses based on GiACTA
Transition probabilities from the initial remission state to the first flare/relapse event were based on individually fitted parametric models using patient-level data. 
Transition probabilities from the subsequent remission state to flare were based on a separate poisson regression. 

## Model assumptions
Only focus on newly diagnosed patients who starts with remission on steroids
30-year time horizon and weekly cycle length

## Our contributions
Update the real-world evidence on steroid adverse events



