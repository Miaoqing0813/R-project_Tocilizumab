# R-project_Tocilizumab
- This is a cost-utility model to evaluate the cost effectiveness of Tocilizumab for people with Giant Cell Arthritis (GCA). 
- The model is based on the HTA model developed by Roche in 2017. 

## Interventions
- Intervention arm: Tocilizumab for 2 years + 26 weeks tapering steroids (prednisolone rather than prednisone)
- Comparison arm: Placebo + 52 weeks tapering steroids (typically as long as 18-24 months in routine clinical practice)
- Cumulative GC dose based on GiACTA trial and real-world data: 1) prior to the first flare; 2) during secondary remission; 3) during relapse/flare

## Population
- Overall ITT population in the GiACTA trial
- Subgroups: 1) newly diagnosed; 2) relapsed

## Overview of Model Structure
A semi-Markov model using a monthly cycle and life-time (30 years) time horizon:
- On initial remission (on and off steroids) 
- On relapse/flare
- On subsequent remission 
- Dead

- GCA-related complications (visual loss, stroke) are modeled as events which only impact the health state values assigned to a proportion of patients at the relapse/flare states.  
- Steroid-related AEs (fractures, diabetes, heart failure, MI, infection) are modeled as events which only impact the health state values assigned to a proportion of patients at the remission/flare states.  

## Transition probabilities
- Transition probabilities from the initial remission state to the first flare/relapse event are time-dependent and based on individually fitted parametric models using patient-level data from the GiACTA trial. 
- Transition probabilities from the subsequent remission state to flare were based on a separate poisson regression using the subgroup of patients following an initial flare from the GiACTA trial. 
- The probabilities of GCA-related complications and GC-related adverse events are based on previous literature: GCA complications are only linked to relapse/flare state; GC-related AEs can be experienced by a patient in any of the remission or flare state. 
- The model assumes no excess mortality risk relative to the general population other than that arising due to stroke. 

## Utilities
- Utilities for the remission and relapse/flare were sourced from a mixed effect regression model based on EQ5D data from GiACTA trial.
- Data was combined across the separate arms given the lack of significant difference by treatment arm. 
- The relapse/flare utility was applied for a 4-week duration (28 days). 
- Utility decrement for GCA and GC-related complications were sourced from the external literature. 

## Costs
- The treatment costs of tocilizumab and GC were based on published prices only. 
- No additional administration costs were assumed for tocilizumab.
- The cost of GC treatment was based on published prices for prednisone. 
- The costs of GCA related complications and GC related AEs were derived from the published literature. 

## Our contributions
- Following York ERG's comment, we extended the weekly cycle to monthly so that the cycle length is consistent with health state values.
- Add more steroid-related adverse events (ie. MI, heart failure, infection, stroke)
- Vary transition probabilities by current and cumulative steroid dosages based on the latest real-world evidence, while the HTA model only used real-world evidence based on cumulative steroid dosages and applied the same probability to all states (did not link the individual model states with steroid dosages).   
- Apply excessive mortality for GCA patients based on the latest real-world evidence, rather than accounting for stroke-related mortality only.
- Following the York ERG comment, 'On remission+on steroid' and 'On remission+off steroid' are combined into a single remission state. 

## Questions to discuss
- The model starts with the initial remission state and did not account for period prior to remission, although ERG did not consider this leads to any significant bias as the evidence does not suggest that this period is longer with tocilizumab. 
- Costs assigned to the flare/relapse state were considered to represent additional costs that would be incurred on top of the background management costs applied to the remission states. However, the ERG notes that uncertainty remains regarding whether it is appropriate to include these background costs in addition to the 3-month event cost assigned to the relapse/flare state. 
- The company used 52-week prednisone tapering regimen while the ERG pointed that clinicians often used a longer tapering regimen in routine clinical practice (typically 18-24 months), which could lead to a greater cumulative steroid burden. Shall we explore additional scenario? 
